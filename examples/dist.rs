use std::path::PathBuf;

use clap::Parser;
use itertools::Itertools;
use packed_seq::{AsciiSeqVec, SeqVec};
use simd_sketch::SketchParams;
use std::io::Write;
use tracing::{info, trace};

#[derive(clap::Parser, Debug, Clone)]
struct Args {
    #[clap(flatten)]
    params: SketchParams,

    paths: Vec<PathBuf>,

    #[clap(long)]
    stats: Option<PathBuf>,
}

fn main() {
    init_trace();

    let args = Args::parse();
    let paths = collect_paths(&args.paths);
    let q = paths.len();

    let k = args.params.k;
    let s = args.params.s;
    let b = args.params.b;

    let sketcher = SketchParams {
        alg: simd_sketch::SketchAlg::Bucket,
        rc: true,
        k,
        s,
        b,
        filter_empty: true,
    }
    .build();

    let mut sketches = vec![];
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
        sketches.push(sketcher.sketch(seq.as_slice()));
        trace!("sketching itself took {:?}", start.elapsed());
    }
    let t_sketch = start.elapsed();
    info!(
        "Sketching {q} seqs took {t_sketch:?} ({:?} avg)",
        t_sketch / q as u32
    );

    let start = std::time::Instant::now();
    let dists = sketches
        .iter()
        .tuple_combinations()
        .map(|(s1, s2)| s1.similarity(s2))
        .collect_vec();
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
            "SimdSketch {:?} {q} {k} {s} {b} {} {}",
            args.params.alg,
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
