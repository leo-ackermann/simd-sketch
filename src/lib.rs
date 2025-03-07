//! # Compute a 32-bit bottom-h MASH.
//!
//! ```
//! use packed_seq::SeqVec;
//!
//! // Generate a random sequence of 2M characters.
//! let n = 2_000_000;
//! let seq = packed_seq::PackedSeqVec::random(n);
//!
//! // Bottom h=10000 sketch of k=31-mers.
//! let k = 31;
//! let h = 10_000;
//!
//! // Return 32-bit hashes.
//! // `false`: use a forward ntHash. `true`: use canonical ntHash instead.
//! let mash: Vec<u32> = simd_mash::mash::<false, _>(seq.as_slice(), k, h);
//!
//! assert_eq!(mash.len(), h);
//! // E.g.:
//! // mash == [332, 904, 1929, 2091, 3058, 5121, 8768, 11738, 13923, 17572, ...]
//! ```
//!
//! This crate uses `packed-seq` for iterating characters of a sequence in parallel using SIMD instructions.
mod intrinsics;

use packed_seq::{u32x8, Seq};
use simd_minimizers::private::nthash::NtHasher;
use tracing::info;

/// Return the `h` smallest `u32` k-mer hashes.
///
/// Set `RC` to `true` for using canonical ntHash.
pub fn mash<'s, const RC: bool, S: Seq<'s>>(seq: S, k: usize, h: usize) -> Vec<u32> {
    // Iterate all kmers and compute 32bit nthashes.
    let n = seq.len();
    let kmers = n - k + 1;
    assert!(
        kmers >= h,
        "Sequence of length n={n} has {kmers} kmers, which is less than h={h}."
    );

    // 32-bit random hashes.
    let mut shift = 0;
    let mut out = vec![];
    loop {
        let target = ((u32::MAX as usize / n * h) << shift).min(u32::MAX as usize) as u32;

        // Should be fine. In practice 10% overlap is probably good enough.
        let bound = target + target / 2;
        tracing::trace!(
            "h {h} target {target} bound {bound}, expected count {}",
            h + h / 2
        );
        let simd_bound = u32x8::splat(bound);

        let (hashes_head, _hashes_tail) =
            simd_minimizers::private::nthash::nthash_seq_simd::<RC, S, NtHasher>(seq, k, 1);

        out.clear();
        out.resize(2 * h, 0);
        let mut write_idx = 0;
        for (i, hashes) in hashes_head.enumerate() {
            let mask = hashes.cmp_lt(simd_bound);
            assert!(
                write_idx < out.len() - 8,
                "i {i} n {n} Buffer len {} is not sufficient for write_idx {write_idx}",
                out.len()
            );
            unsafe { intrinsics::append_from_mask(hashes, mask, &mut out, &mut write_idx) };
        }

        tracing::trace!("Got {}, needed {h}, expected {}", write_idx, h + h / 2);

        if write_idx >= h {
            out.resize(write_idx, 0);
            let start = std::time::Instant::now();
            out.sort_unstable();
            tracing::trace!("Sorting took {:?}", start.elapsed());
            out.resize(h, 0);
            break out;
        }
        shift += 1;
        info!("Doing another iteration of mash because found only {write_idx} values with hash < {bound}, while h={h} are needed and {} were expected.", k+k/2);
    }
}
