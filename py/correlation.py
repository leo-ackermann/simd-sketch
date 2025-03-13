#!/usr/bin/env python3

import matplotlib.pyplot as plt
from pathlib import Path
import sys
import re
import random
import numpy as np
import matplotlib.patches as mpatches

plt.close()
paths = sys.argv[1:]
n = 50000

dir = Path("../output/")


def key(s, _nsre=re.compile(r"(\d+)")):
    return [
        int(text) if text.isdigit() else text.lower() for text in _nsre.split(str(s))
    ]


def correlation(a, b):
    return np.corrcoef(a, b)[0, 1]


baseline = dir / "simd_bottom_s65536.dist"
groups = [
    (
        sorted(
            list(x for x in dir.glob("simd_bottom_*.dist") if x != baseline), key=key
        ),
        "SimdSketch bottom",
    ),
    (sorted(list(dir.glob("bindash_bottom_s*.dist")), key=key), "BinDash bottom"),
    (sorted(list(dir.glob("bindashrs_*.dist")), key=key), "BinDash-rs bucket"),
    # (
    #     sorted(list(dir.glob("bindash_bottom_fixed*.dist")), key=key),
    #     "BinDash bottom (fixed)",
    # ),
    (
        sorted(list(dir.glob("simd_bucket_*b32.dist")), key=key),
        "SimdSketch bucket b=32",
    ),
    (sorted(list(dir.glob("simd_bucket_*b8.dist")), key=key), "SimdSketch bucket b=8"),
    (sorted(list(dir.glob("simd_bucket_*b1.dist")), key=key), "SimdSketch bucket b=1"),
    (
        sorted(list(dir.glob("bindash_bucket_*b32.dist")), key=key),
        "BinDash bucket b=32",
    ),
    (sorted(list(dir.glob("bindash_bucket_*b8.dist")), key=key), "BinDash bucket b=8"),
    (sorted(list(dir.glob("bindash_bucket_*b1.dist")), key=key), "BinDash bucket b=1"),
]


def read(p):
    print("Reading", p)
    return [float(x) for x in Path(p).read_text().splitlines()]


d0 = read(baseline)

# Sample n random lines from each file
indices = random.sample(range(len(d0)), n)

# one subfigure for each group
for i, (group, title) in enumerate(groups):
    print(*group)
    names = [Path(p).stem for p in group]
    dists = [read(p) for p in group]
    plt.subplot(3, 3, i + 1)
    for name, d in zip(names, dists):
        print(f"plotting len")
        # extract value of s from name, _s\d+_
        c = correlation(d0, d)
        plt.scatter(
            [d0[idx] for idx in indices],
            [d[idx] for idx in indices],
            label=f"{c:.4f}",
            alpha=0.4,
            s=2,
        )
    # plt.xlabel(baseline.stem)
    leg = plt.legend()
    for lh in leg.legend_handles:
        lh.set_alpha(1)
        lh.set_sizes([50] * 4)
    plt.title(title)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xticks([0, 1])
    plt.yticks([0, 1])
    plt.plot([0, 1], [0, 1], color="black", linestyle="-", lw=0.5)

# Add legend under the plot mapping s to each colour.
# s=128: blue
# s=1024: orange
# s=8192: green
# s=65536: red

# Build manual legend
plt.subplot(4, 3, 12)
plt.axis("off")
# Red circle for legend

handles = [
    mpatches.Patch(color="blue", label="s = 128"),
    mpatches.Patch(color="orange", label="s = 1024"),
    mpatches.Patch(color="green", label="s = 8192"),
    mpatches.Patch(color="red", label="s = 16384"),
    mpatches.Patch(color="purple", label="s = 32768"),
]

plt.figlegend(
    handles=handles,
    loc="lower center",
    ncol=5,
    labelspacing=0.0,
    bbox_to_anchor=(0.5, 0.05),
)

# Set the figure size
plt.gcf().set_size_inches(15, 10)

plt.tight_layout()
plt.savefig("plots/correlation.png", dpi=300, bbox_inches="tight")
