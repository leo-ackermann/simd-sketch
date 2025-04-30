use std::path::PathBuf;

use clap::Parser;
use packed_seq::{AsciiSeqVec, SeqVec};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use simd_sketch::SketchParams;
use tracing::info;

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
        #[clap(long)]
        output: Option<PathBuf>,
    },
}

fn main() {
    init_trace();

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

    let start = std::time::Instant::now();

    let sketches: Vec<_> = paths
        .par_iter()
        .map(|path| {
            let mut seq = AsciiSeqVec::default();
            let mut reader = needletail::parse_fastx_file(&path).unwrap();
            while let Some(r) = reader.next() {
                // FIXME: Skip adjacent k-mers across fasta records.
                seq.push_ascii(&r.unwrap().seq());
            }
            sketcher.sketch(seq.as_slice())
        })
        .collect();
    let t_sketch = start.elapsed();
    info!(
        "Sketching {q} seqs took {t_sketch:?} ({:?} avg)",
        t_sketch / q as u32
    );

    let mut pairs = Vec::with_capacity(q * (q - 1) / 2);
    for i in 0..q {
        for j in 0..i {
            pairs.push((i, j));
        }
    }
    let start = std::time::Instant::now();
    let dists: Vec<_> = pairs
        .into_par_iter()
        .map(|(i, j)| sketches[i].mash_dist(&sketches[j]))
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

fn init_trace() {
    use tracing::level_filters::LevelFilter;
    use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};

    tracing_subscriber::registry()
        .with(tracing_subscriber::fmt::layer().with_writer(std::io::stderr))
        .with(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(LevelFilter::TRACE.into())
                .from_env_lossy(),
        )
        .init();
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
