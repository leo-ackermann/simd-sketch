use std::path::PathBuf;

use clap::Parser;
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use log::info;
use packed_seq::{AsciiSeqVec, SeqVec};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use regex::bytes::RegexBuilder;
use simd_sketch::SketchParams;

/// Compute the sketch distance between two fasta files.
#[derive(clap::Parser)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

/// TODO: Support for writing sketches to disk.
#[derive(clap::Subcommand)]
enum Command {
    /// Compute the distance between two sequences.
    Dist {
        #[command(flatten)]
        params: SketchParams,
        /// First input fasta file.
        path_a: PathBuf,
        /// Second input fasta file.
        path_b: PathBuf,
    },
    /// Takes paths to fasta files, and outputs a Phylip distance matrix to stdout.
    Triangle {
        #[command(flatten)]
        params: SketchParams,
        /// Paths to (directories of) (gzipped) fasta files.
        paths: Vec<PathBuf>,
        /// Write phylip distance matrix here, or default to stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
}

fn main() {
    env_logger::init();

    let args = Args::parse();

    let (params, paths) = match &args.command {
        Command::Dist {
            params,
            path_a,
            path_b,
        } => (params, vec![path_a.clone(), path_b.clone()]),
        Command::Triangle { params, paths, .. } => (params, collect_paths(&paths)),
    };

    let q = paths.len();

    let sketcher = params.build();

    let style = indicatif::ProgressStyle::with_template(
        "{msg:.bold} [{elapsed_precise:.cyan}] {bar} {pos}/{len} ({percent:>3}%)",
    )
    .unwrap()
    .progress_chars("##-");

    let start = std::time::Instant::now();

    let sketches: Vec<_> = paths
        .par_iter()
        .progress_with_style(style.clone())
        .with_message("Sketching")
        .with_finish(indicatif::ProgressFinish::AndLeave)
        .map(|path| {
            let mut seqs = vec![];
            let mut reader = needletail::parse_fastx_file(&path).unwrap();
            while let Some(r) = reader.next() {
                let seq: &[u8] = &r.unwrap().seq().into_owned();
                let re = RegexBuilder::new(r"[N]+")
                    .case_insensitive(true)
                    .unicode(false)
                    .build()
                    .unwrap();
                let subseqs_without_n: Vec<&[u8]> = re.split(seq).collect();

                for swn in subseqs_without_n {
                    seqs.push(AsciiSeqVec::from_ascii(swn));
                }
            }
            let slices = seqs.iter().map(|s| s.as_slice()).collect_vec();
            sketcher.sketch_seqs(&slices)
        })
        .collect();
    let t_sketch = start.elapsed();

    info!(
        "Sketching {q} seqs took {t_sketch:?} ({:?} avg)",
        t_sketch / q as u32
    );

    let num_pairs = q * (q - 1) / 2;
    let mut pairs = Vec::with_capacity(num_pairs);
    for i in 0..q {
        for j in 0..i {
            pairs.push((i, j));
        }
    }
    let start = std::time::Instant::now();
    let dists: Vec<_> = pairs
        .into_par_iter()
        .progress_with_style(style.clone())
        .with_message("Distances")
        .with_finish(indicatif::ProgressFinish::AndLeave)
        .map(|(i, j)| sketches[i].mash_distance(&sketches[j]))
        .collect();
    let t_dist = start.elapsed();

    let cnt = q * (q - 1) / 2;
    info!(
        "Computing {cnt} dists took {t_dist:?} ({:?} avg)",
        t_dist / cnt.max(1) as u32
    );

    match &args.command {
        Command::Dist { .. } => {
            println!("Distance: {:.4}", dists[0]);
            return;
        }
        Command::Triangle { output, .. } => {
            use std::io::Write;

            // Output Phylip triangle format.
            let mut out = Vec::new();
            writeln!(out, "{q}").unwrap();
            let mut d = dists.iter();
            for i in 0..q {
                write!(out, "{}", paths[i].to_str().unwrap()).unwrap();
                for _ in 0..i {
                    write!(out, "\t{:.7}", d.next().unwrap()).unwrap();
                }
                writeln!(out).unwrap();
            }

            match output {
                Some(output) => std::fs::write(output, out).unwrap(),
                None => println!("{}", str::from_utf8(&out).unwrap()),
            }
        }
    }
}

fn collect_paths(paths: &Vec<PathBuf>) -> Vec<PathBuf> {
    let mut res = vec![];
    for path in paths {
        if path.is_dir() {
            res.extend(path.read_dir().unwrap().map(|entry| entry.unwrap().path()));
        } else {
            res.push(path.clone());
        }
    }
    res.sort();
    res
}
