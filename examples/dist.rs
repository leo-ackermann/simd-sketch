use std::path::PathBuf;

use clap::Parser;
use itertools::Itertools;
use packed_seq::{PackedSeqVec, SeqVec};
use tracing::{info, trace};

#[derive(clap::Parser)]
struct Args {
    paths: Vec<PathBuf>,
    #[clap(short, long)]
    bin: bool,

    #[clap(short, default_value_t = 31)]
    k: usize,
    #[clap(short, default_value_t = 10000)]
    h: usize,
}

fn main() {
    init_trace();

    let args = Args::parse();
    let s = args.paths.len();

    let k = args.k;
    let h = args.h;

    let masher = simd_mash::Masher::new(k, h);

    let mut bottom_mashes = vec![];
    let mut bin_mashes = vec![];
    let start = std::time::Instant::now();
    for path in args.paths {
        trace!("Sketching {path:?}");
        let mut seq = PackedSeqVec::default();
        let mut reader = needletail::parse_fastx_file(path).unwrap();
        while let Some(r) = reader.next() {
            let record = r
                .unwrap()
                .seq()
                .iter()
                .filter_map(|&b| if b == b'N' { None } else { Some(b) })
                .collect::<Vec<_>>();
            seq.push_ascii(&record);
        }
        let start = std::time::Instant::now();
        if args.bin {
            bin_mashes.push(masher.bin_mash(seq.as_slice()));
        } else {
            bottom_mashes.push(masher.bottom_mash(seq.as_slice()));
        };
        eprintln!("sketching took {:?}", start.elapsed());
    }
    let t = start.elapsed();
    info!("Sketching {s} seqs took {t:?} ({:?} avg)", t / s as u32);

    let start = std::time::Instant::now();
    let dists = if args.bin {
        bin_mashes
            .iter()
            .tuple_combinations()
            .map(|(s1, s2)| s1.similarity(s2))
            .collect_vec()
    } else {
        bottom_mashes
            .iter()
            .tuple_combinations()
            .map(|(s1, s2)| s1.similarity(s2))
            .collect_vec()
    };
    let t = start.elapsed();
    let cnt = s * (s - 1) / 2;
    info!(
        "Computing {cnt} dists took {t:?} ({:?} avg)",
        t / cnt.max(1) as u32
    );
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
