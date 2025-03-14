# SimdSketch

[![crates.io](https://img.shields.io/crates/v/simd-sketch.svg)](https://crates.io/crates/simd-sketch)
[![docs.rs](https://img.shields.io/docsrs/simd-sketch.svg)](https://docs.rs/simd-sketch)

A SIMD-accelerated library to compute two types of sketches:
- Classic bottom $s$ sketch, containing the $s$ smallest distinct k-mer hashes.
- Bucket sketch, which partitions the hashes into $s$ parts and returns the smallest
  hash in each part. (Introduced as *one permutation hashing* in Li, Owen, Zhang 2012.)

See the corresponding [blog post](https://curiouscoding.nl/posts/simd-sketch/)
for background and evaluation.

Sketching takes 2 seconds for a 3Gbp human genome. This library returns 32-bit `u32`
hashes. This means that currently it may not be very suitable for sequences that are
too close to 1Gbp in length, since the bottom hash values will be relatively dense.

**Algorithm.**
For the bottom $s$ sketch, we first collect all ``sufficiently small'' hashes
into a vector. Then, that vector is sorted and deduplicated, and the smallest
$s$ values are returned. This ensures that the runtime is $O(n + s \log s)$ when
the number of duplicate k-mers is limited.

For the bucket sketch, the classic method is to partition hashes linearly, e.g.,
for $s=2$ into the bottom and top half. Then, a single value is kept per part,
and each hash is compared against the rolling minimum of its bucket.

Instead, here we make buckets by the remainder modulo $s$. This way, we can
again pre-filter for ``sufficiently small'' values, and then only scan those for
the minimum.

In both variants, we double the ``smallness'' threshold until either $s$
distinct values are found or all $s$ buckets have a value in them.


**Implementation notes.**
Good performance is mostly achieved by using a branch-free implementation: all
hashes are computed using 8 parallel streams using SIMD, and appended to a vector when they
are sufficiently small to likely be part of the sketch.

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

let sketch1: simd_sketch::BinSketch = sketcher.sketch(seq1.as_slice());
let sketch2: simd_sketch::BinSketch = sketcher.sketch(seq2.as_slice());

// Value between 0 and 1, estimating the fraction of shared k-mers.
let similarity = sketch1.similarity(&sketch2);
```

**TODO:** If you would like a binary instead of a library, please create an
issue and propose an API.
