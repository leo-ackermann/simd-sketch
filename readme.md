# SimdSketch

[![crates.io](https://img.shields.io/crates/v/simd-sketch.svg)](https://crates.io/crates/simd-sketch)
[![docs.rs](https://img.shields.io/docsrs/simd-sketch.svg)](https://docs.rs/simd-sketch)

A SIMD-accelerated library to compute two types of sketches:
- Classic bottom-$h$ sketch, containing the $h$ smallest distinct k-mer hashes.
- Bucket sketch, which partitions the hashes into $h$ parts and returns the smallest
  hash in each part.


Sketching takes 2 seconds for a 3Gbp human genome. This library returns 32-bit `u32`
hashes. This means that currently it may not be very suitable for sequences that are
too close to 1Gbp in length, since the bottom hash values will be relatively dense.

Good performance is mostly achieved by using a branch-free implementation: all
hashes are computed using 8 parallel streams, and appended to a vector when they
are sufficiently small to likely be part of the sketch. That small subset
is sorted and deduplicated. If not enough hashes remain, the initial threshold
is increased.

The underlying streaming and hashing algorithms are described in the following [preprint](https://doi.org/10.1101/2025.01.27.634998):

- SimdMinimizers: Computing random minimizers, fast.
  Ragnar Groot Koerkamp, Igor Martayan
  bioRxiv 2025.01.27 [doi.org/10.1101/2025.01.27.634998](https://doi.org/10.1101/2025.01.27.634998)


## Usage
Please see [lib.rs](src/lib.rs) and [docs.rs](https://docs.rs/simd-sketch) for
full docs.

```rust
use packed_seq::SeqVec;

// Bottom h=10000 sketch of k=31-mers.
let k = 31;
let h = 10_000;

// Use `new_rc` for a canonical version instead.
let sketch = simd_sketch::Sketcher::new(k, h);

// Generate two random sequences of 2M characters.
let n = 2_000_000;
let seq1 = packed_seq::PackedSeqVec::random(n);
let seq2 = packed_seq::PackedSeqVec::random(n);

let sketch1: simd_sketch::BottomSketch = sketcher.bottom_sketch(seq1.as_slice());
let sketch2: simd_sketch::BottomSketch = sketcher.bottom_sketch(seq2.as_slice());

// Value between 0 and 1, estimating the fraction of shared k-mers.
let similarity = sketch1.similarity(&sketch2);

// Bucket sketch variant

let sketch1: simd_sketch::BinSketch = sketcher.bucket_sketch(seq1.as_slice());
let sketch2: simd_sketch::BinSketch = sketcher.bucket_sketch(seq2.as_slice());

// Value between 0 and 1, estimating the fraction of shared k-mers.
let similarity = sketch1.similarity(&sketch2);
```
