use std::path::PathBuf;

use clap::Parser;
use packed_seq::{PackedSeqVec, SeqVec};

#[derive(clap::Parser)]
struct Args {
    p1: PathBuf,
    p2: PathBuf,
    #[clap(short, long)]
    bin: bool,
}

fn main() {
    init_trace();

    let args = Args::parse();

    let k = 31;
    let h = 10_000;

    let sketches = [args.p1, args.p2].map(|path| {
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
        let mash = if args.bin {
            simd_mash::bin_mash::<false, _>(seq.as_slice(), k, h)
        } else {
            simd_mash::mash::<false, _>(seq.as_slice(), k, h)
        };

        let elapsed = start.elapsed();
        tracing::info!("{h} hashes in {elapsed:?}");
        mash
    });

    let overlap = if args.bin {
        simd_mash::bin_intersection(&sketches[0], &sketches[1])
    } else {
        simd_mash::set_intersection_size(&sketches[0], &sketches[1])
    };
    tracing::info!("overlap = {}", overlap);
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
