use packed_seq::{PackedSeqVec, SeqVec};
use simd_sketch::BitSketch;

fn main() {
    let k = 31;
    let n = 4_000_000;

    let s = 32768;
    let b = 8;

    let mut counts = [0; 256];
    let mut bottom = 0.0;
    let mut bucket = 0.0;
    for _ in 0..10 {
        let seq1 = PackedSeqVec::random(n);
        let seq2 = PackedSeqVec::random(n);

        let sketcher = simd_sketch::Sketcher::new_rc(k, s, b);
        let s1 = sketcher.bottom_sketch(seq1.as_slice());
        let s2 = sketcher.bottom_sketch(seq2.as_slice());
        bottom += s1.similarity(&s2);
        let s1 = sketcher.sketch(seq1.as_slice());
        let s2 = sketcher.sketch(seq2.as_slice());
        bucket += s1.similarity(&s2);

        if let BitSketch::B8(x) = s1.buckets {
            for &v in x.iter() {
                counts[v as usize] += 1;
            }
        }
    }
    println!("AVG Bottom: {}", bottom / 10.0);
    println!("AVG Bucket: {}", bucket / 10.0);
    // if b == 8 {
    //     for (i, v) in counts.iter().enumerate() {
    //         println!("{i:>3} {v:>3}");
    //     }
    // }
}
