use std::hint::black_box;

use packed_seq::{PackedSeqVec, SeqVec};

fn main() {
    init_trace();
    // Compute rolling 32bit hash

    let k = 31;
    let h = 10_000;
    let n = 2_000_000;

    let seq = PackedSeqVec::random(n);

    let start = std::time::Instant::now();
    let masher = simd_mash::Masher::new(k, h);
    let mash = masher.bottom_mash(seq.as_slice());
    let elapsed = start.elapsed();

    tracing::info!("{h} hashes in {elapsed:?}");
    black_box(mash);
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
