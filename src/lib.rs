//! # Compute a 32-bit bottom-h MASH.
//!
//! This library works by partitioning the input sequence into 8 chunks,
//! and processing those in paralle using SIMD. This is based on the `packed-seq` and `simd-minimizers` crates.
//!
//! The largest hash should be around `target = u32::MAX * h / n`.
//! To ensure a branch-free algorithm, we first collect all hashes up to `bound = 1.5 * target`.
//! Then we sort the collected hashes, and return the bottom `h`.
//!
//! ### Example
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
        let bound = target.saturating_add(target / 2);
        tracing::trace!(
            "h {h} target {target} bound {bound}, expected count {}",
            h + h / 2
        );
        let write_idx = collect_up_to_bound::<RC, S>(seq, k, h, bound, &mut out);

        tracing::trace!("Got {}, needed {h}, expected {}", out.len(), h + h / 2);

        if out.len() >= h {
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

/// For each remainder modulo h, return the smallest seen value.
///
/// Set `RC` to `true` for using canonical ntHash.
pub fn bin_mash<'s, const RC: bool, S: Seq<'s>>(seq: S, k: usize, h: usize) -> Vec<u32> {
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
    'l: loop {
        let target = ((u32::MAX as usize / n * h) << shift).min(u32::MAX as usize) as u32;

        // Should be fine. In practice 10% overlap is probably good enough.
        let bound = target * h.ilog2() * 2;
        tracing::trace!(
            "h {h} target {target} bound {bound}, expected count {}",
            h + h / 2
        );
        let write_idx = collect_up_to_bound::<RC, S>(seq, k, h, bound, &mut out);

        tracing::trace!("Got {}, needed {h}, expected {}", out.len(), h + h / 2);

        if out.len() >= h {
            let m = FM32::new(h as u32);
            let start = std::time::Instant::now();
            let mut buckets = vec![u32::MAX; h];
            for &hash in &out {
                let bucket = m.reduce(hash);
                buckets[bucket] = buckets[bucket].min(hash);
            }
            let mut empty = 0;
            for &x in &buckets {
                if x == u32::MAX {
                    empty += 1;
                }
            }
            if empty > 0 {
                info!("Doing another iteration of mash because found {empty} empty buckets, while h={h} are needed and {} were expected.", k+k/2);
                shift += 1;
                continue 'l;
            }
            tracing::trace!("Sorting took {:?}", start.elapsed());
            break buckets;
        }
        shift += 1;
        info!("Doing another iteration of mash because found only {write_idx} values with hash < {bound}, while h={h} are needed and {} were expected.", k+k/2);
    }
}

fn collect_up_to_bound<'s, const RC: bool, S: Seq<'s>>(
    seq: S,
    k: usize,
    h: usize,
    bound: u32,
    out: &mut Vec<u32>,
) -> usize {
    let simd_bound = u32x8::splat(bound);

    let (hashes_head, hashes_tail) =
        simd_minimizers::private::nthash::nthash_seq_simd::<RC, S, NtHasher>(seq, k, 1);

    out.clear();
    out.resize(2 * h, 0);
    let mut write_idx = 0;
    for hashes in hashes_head {
        let mask = hashes.cmp_lt(simd_bound);
        if write_idx + 8 >= out.len() {
            out.resize(write_idx * 3 / 2 + 8, 0);
        }
        unsafe { intrinsics::append_from_mask(hashes, mask, out, &mut write_idx) };
    }

    out.resize(write_idx, 0);

    for hash in hashes_tail {
        if hash <= bound {
            out.push(hash);
        }
    }
    write_idx
}

/// FastMod32, using the low 32 bits of the hash.
/// Taken from https://github.com/lemire/fastmod/blob/master/include/fastmod.h
#[derive(Copy, Clone, Debug)]
pub struct FM32 {
    d: u64,
    m: u64,
}
impl FM32 {
    fn new(d: u32) -> Self {
        Self {
            d: d as u64,
            m: u64::MAX / d as u64 + 1,
        }
    }
    fn reduce(self, h: u32) -> usize {
        let lowbits = self.m * (h as u64);
        ((lowbits as u128 * self.d as u128) >> 64) as usize
    }
}

/// The number of common elements between two sorted lists.
pub fn set_intersection_size(a: &[u32], b: &[u32]) -> usize {
    assert_eq!(a.len(), b.len());
    let mut count = 0;
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        count += (a[i] == b[j]) as usize;
        let di = (a[i] <= b[j]) as usize;
        let dj = (a[i] >= b[j]) as usize;
        i += di;
        j += dj;
    }
    count
}

/// The number of common elements between two sorted lists.
pub fn bin_intersection(a: &[u32], b: &[u32]) -> usize {
    assert_eq!(a.len(), b.len());
    std::iter::zip(a, b).map(|(a, b)| (a == b) as usize).sum()
}

#[cfg(test)]
#[test]
fn test_overlap() {
    let count = set_intersection_size(&[1, 2, 3, 5, 9, 10], &[1, 3, 5, 7, 9, 11]);
    assert_eq!(count, 4);
}

#[cfg(test)]
#[test]
fn test() {
    use packed_seq::SeqVec;

    let k = 31;
    for n in 31..100 {
        let h = n - k + 1;
        let seq = packed_seq::PackedSeqVec::random(n);
        let mash = crate::mash::<false, _>(seq.as_slice(), k, h);
        assert_eq!(mash.len(), h);
        assert!(mash.is_sorted());

        let h = h.min(10);
        let seq = packed_seq::PackedSeqVec::random(n);
        let mash = crate::mash::<false, _>(seq.as_slice(), k, h);
        assert_eq!(mash.len(), h);
        assert!(mash.is_sorted());
    }
}

#[cfg(test)]
#[test]
fn rc() {
    use packed_seq::SeqVec;

    for k in (0..10).map(|_| rand::random_range(1..100)) {
        for n in (0..10).map(|_| rand::random_range(k..1000)) {
            for h in (0..10).map(|_| rand::random_range(0..n - k + 1)) {
                let seq = packed_seq::AsciiSeqVec::random(n);
                let mash = crate::mash::<true, _>(seq.as_slice(), k, h);
                assert_eq!(mash.len(), h);
                assert!(mash.is_sorted());

                let seq_rc = packed_seq::AsciiSeqVec::from_ascii(
                    &seq.seq
                        .iter()
                        .rev()
                        .map(|c| packed_seq::complement_char(*c))
                        .collect::<Vec<_>>(),
                );

                let mash_rc = crate::mash::<true, _>(seq_rc.as_slice(), k, h);
                assert_eq!(mash, mash_rc);
            }
        }
    }
}
