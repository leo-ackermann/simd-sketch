use std::path::PathBuf;

use clap::Parser;
use itertools::Itertools;
use packed_seq::{AsciiSeqVec, SeqVec};
use std::io::Write;
use tracing::{info, trace};

#[derive(clap::Parser, Debug, Clone)]
struct Args {
    paths: Vec<PathBuf>,
    #[clap(long)]
    bucket: bool,

    /// k-mer length
    #[clap(short, default_value_t = 31)]
    k: usize,

    /// Sketch size
    #[clap(short, default_value_t = 10000)]
    s: usize,
    /// Store bottom-b bits of each element. Must be multiple of 8.
    #[clap(short, default_value_t = 32)]
    b: usize,

    #[clap(long)]
    stats: Option<PathBuf>,
}

fn main() {
    init_trace();

    let args = Args::parse();
    let paths = collect_paths(&args.paths);
    let q = paths.len();

    let k = args.k;
    let s = args.s;
    let b = args.b;

    let mut sketcher = simd_sketch::Sketcher::new_rc(k, s, b);
    sketcher.filter_empty = true;

    let mut bottom_sketches = vec![];
    let mut bucket_sketches = vec![];
    let start = std::time::Instant::now();

    for path in paths {
        trace!("Sketching {path:?}");
        let mut seq = AsciiSeqVec::default();
        let mut reader = needletail::parse_fastx_file(path).unwrap();
        let start = std::time::Instant::now();
        while let Some(r) = reader.next() {
            // let record = r
            //     .unwrap()
            //     .seq();
            // .iter()
            // .filter_map(|&b| if b == b'N' { None } else { Some(b) })
            // .collect::<Vec<_>>();
            // seq.push_ascii(&record);
            seq.push_ascii(&r.unwrap().seq());
            // FIXME: Skip adjacent k-mers.
        }
        trace!("Reading & filtering took {:?}", start.elapsed());
        let start = std::time::Instant::now();
        if args.bucket {
            bucket_sketches.push(sketcher.sketch(seq.as_slice()));
        } else {
            bottom_sketches.push(sketcher.bottom_sketch(seq.as_slice()));
        };
        trace!("sketching itself took {:?}", start.elapsed());
    }
    let t_sketch = start.elapsed();
    info!(
        "Sketching {q} seqs took {t_sketch:?} ({:?} avg)",
        t_sketch / q as u32
    );

    let start = std::time::Instant::now();
    let dists = if args.bucket {
        bucket_sketches
            .iter()
            .tuple_combinations()
            .map(|(s1, s2)| s1.similarity(s2))
            .collect_vec()
    } else {
        bottom_sketches
            .iter()
            .tuple_combinations()
            .map(|(s1, s2)| s1.similarity(s2))
            .collect_vec()
    };
    let t_dist = start.elapsed();
    let cnt = q * (q - 1) / 2;
    info!(
        "Computing {cnt} dists took {t_dist:?} ({:?} avg)",
        t_dist / cnt.max(1) as u32
    );
    info!(
        "Params {:?}",
        Args {
            paths: vec![],
            ..args.clone()
        }
    );

    if let Some(stats) = &args.stats {
        let mut writer = std::fs::File::options()
            .create(true)
            .append(true)
            .write(true)
            .open(stats)
            .unwrap();
        writeln!(
            writer,
            "SimdSketch {} {q} {k} {s} {b} {} {}",
            if args.bucket { "bucket" } else { "bottom" },
            t_sketch.as_secs_f32(),
            t_dist.as_secs_f32()
        )
        .unwrap();
    }

    for dist in dists {
        println!("{dist}");
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
