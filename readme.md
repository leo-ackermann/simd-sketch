# simd-mash

[![crates.io](https://img.shields.io/crates/v/simd-mash.svg)](https://crates.io/crates/simd-mash)
[![docs.rs](https://img.shields.io/docsrs/simd-mash.svg)](https://docs.rs/simd-mash)

A SIMD-accelerated library to compute the bottom $h$ sketch, or _Mash_ of a sequence.

Takes 4 seconds for a 3Gbp human genome. This library returns 32-bit `u32`
hashes. This means that currently it may not be very suitable for sequences that are
too close to 1Gbp in length, since the bottom hash values will be relatively dense.

The largest hash should be around `target = u32::MAX * h / n`.
To ensure a branch-free algorithm, we first collect all hashes up to `bound = 1.5 * target`.
Then we sort the collected hashes, and return the bottom `h`.

The underlying streaming and hashing algorithms are described in the following [preprint](https://doi.org/10.1101/2025.01.27.634998):

- SimdMinimizers: Computing random minimizers, fast.
  Ragnar Groot Koerkamp, Igor Martayan
  bioRxiv 2025.01.27 [doi.org/10.1101/2025.01.27.634998](https://doi.org/10.1101/2025.01.27.634998)


## Usage example
Full documentation can be found on [docs.rs](https://docs.rs/simd-mash).

```rust
use packed_seq::SeqVec;

// Generate a random sequence of 2M characters.
let n = 2_000_000;
let seq = packed_seq::PackedSeqVec::random(n);

// Bottom h=10000 sketch of k=31-mers.
let k = 31;
let h = 10_000;

// Return 32-bit hashes.
// `false`: use a forward ntHash. `true`: use canonical ntHash instead.
let mash: Vec<u32> = simd_mash::mash::<false, _>(seq.as_slice(), k, h);

assert_eq!(mash.len(), h);
// E.g.:
// mash == [332, 904, 1929, 2091, 3058, 5121, 8768, 11738, 13923, 17572, ...]
```
