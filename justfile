build:
    cargo build -r --example dist

bench: build

simd_bot s: build
    ./target/release/examples/dist input -s {{s}}                --stats py/stats > output/simd_bot_s{{s}}.dist
simd_bin s b: build
    ./target/release/examples/dist input -s {{s}} -b {{b}} --bin --stats py/stats > output/simd_bin_s{{s}}_b{{b}}.dist

simd_bot_all: (simd_bot "128") (simd_bot "1024") (simd_bot "8192") (simd_bot"65536")
simd_bin_all_b s: (simd_bin s "1") (simd_bin s "8") (simd_bin s "16") (simd_bin s "32")
simd_bin_all: (simd_bin_all_b "128") (simd_bin_all_b "1024") (simd_bin_all_b "8192")

bindash_bot s:
    time -f %U ./bindash sketch input/*fa --minhashtype=0 --kmerlen=31 --sketchsize={{s}} --outfname=output/tmp 2>&1 | tee >(cat 1>&2) | tail -1 > output/tmp_sketch_time
    time -f %U ./bindash dist output/tmp --outfname=output/bindash_bot_s{{s}}.dist 2>&1 | tee >(cat 1>&2) | tail -1 > output/tmp_dist_time
    echo BinDash bottom 1000 31 {{s}} 64 `cat output/tmp_sketch_time` `cat output/tmp_dist_time` >> py/stats
bindash_bin s b:
    time -f %U ./bindash sketch input/*fa --minhashtype=2 --kmerlen=31 --sketchsize={{s}} --bbits={{b}} --outfname=output/tmp 2>&1 | tee >(cat 1>&2) | tail -1 > output/tmp_sketch_time
    time -f %U ./bindash dist output/tmp --outfname=output/bindash_bin_s{{s}}_b{{b}}.dist 2>&1 | tee >(cat 1>&2) | tail -1 > output/tmp_dist_time
    echo BinDash bin 1000 31 {{s}} {{b}} `cat output/tmp_sketch_time` `cat output/tmp_dist_time` >> py/stats


bindash_bot_all: (bindash_bot "128") (bindash_bot "1024") (bindash_bot "8192")
bindash_bin_all_b s: (bindash_bin s "1") (bindash_bin s "8") (bindash_bin s "16") (bindash_bin s "32")
bindash_bin_all: (bindash_bin_all_b "128") (bindash_bin_all_b "1024") (bindash_bin_all_b "8192")

bindashrs s:
    ./bindash-rs -t 6 -k 31 --sketch_size {{s}} -q input_files -r input_files --stats py/stats -o output/bindashrs_s{{s}}.dist


bindashrs_all: (bindashrs "128") (bindashrs "1024") (bindashrs "8192")
