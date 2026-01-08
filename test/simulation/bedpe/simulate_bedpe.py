#!/usr/bin/env python3
"""
Simulate BEDPE contact pairs (e.g., Hi-C read-pairs).

Output columns (BEDPE-ish):
chrom1 start1 end1 chrom2 start2 end2 name score

Notes
- Provide a chrom.sizes file: two columns -> chrom <TAB> size
- Cis distances follow a truncated power-law to mimic distance decay.
- Coordinates are binned by --bin-size (e.g. 5000 for 5 kb).
"""

from __future__ import annotations
import argparse
import gzip
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


def read_chrom_sizes(path: str) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, size = line.split()[:2]
            sizes[chrom] = int(size)
    if not sizes:
        raise ValueError("chrom.sizes appears empty.")
    return sizes


def choose_weighted(items: List[str], weights: np.ndarray, rng: np.random.Generator) -> str:
    idx = rng.choice(len(items), p=weights)
    return items[idx]


def truncated_power_law(
    n: int,
    xmin: float,
    xmax: float,
    alpha: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample from p(x) ~ x^(-alpha) on [xmin, xmax], alpha > 1.
    Inverse-CDF for truncated power law.
    """
    if alpha <= 1:
        raise ValueError("alpha must be > 1 for a proper power law.")
    u = rng.random(n)
    a1 = 1.0 - alpha
    # CDF(x) = (x^(1-alpha) - xmin^(1-alpha)) / (xmax^(1-alpha) - xmin^(1-alpha))
    # Invert:
    x = (u * (xmax**a1 - xmin**a1) + xmin**a1) ** (1.0 / a1)
    return x


def bin_interval(pos: int, bin_size: int) -> Tuple[int, int]:
    start = (pos // bin_size) * bin_size
    end = start + bin_size
    return start, end


def simulate_bedpe(
    chrom_sizes: Dict[str, int],
    n_pairs: int,
    cis_fraction: float,
    bin_size: int,
    cis_min_dist: int,
    cis_max_dist: int,
    cis_alpha: float,
    seed: Optional[int],
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    chroms = list(chrom_sizes.keys())
    lengths = np.array([chrom_sizes[c] for c in chroms], dtype=float)
    weights = lengths / lengths.sum()

    # decide cis vs trans
    is_cis = rng.random(n_pairs) < cis_fraction

    out = {
        "chrom1": [],
        "start1": [],
        "end1": [],
        "chrom2": [],
        "start2": [],
        "end2": [],
        "name": [],
        "score": [],
    }

    # pre-sample cis distances
    n_cis = int(is_cis.sum())
    if n_cis > 0:
        dists = truncated_power_law(
            n=n_cis,
            xmin=float(cis_min_dist),
            xmax=float(cis_max_dist),
            alpha=float(cis_alpha),
            rng=rng,
        ).astype(int)
    else:
        dists = np.array([], dtype=int)

    cis_i = 0
    for i in range(n_pairs):
        if is_cis[i]:
            chrom = choose_weighted(chroms, weights, rng)
            L = chrom_sizes[chrom]

            # choose first position uniformly
            p1 = int(rng.integers(0, max(L, 1)))
            dist = int(dists[cis_i])
            cis_i += 1

            # random direction
            if rng.random() < 0.5:
                p2 = p1 + dist
            else:
                p2 = p1 - dist

            # if out of bounds, resample a bit (fast + simple)
            # (for realistic sims you might reflect; resampling is fine here)
            tries = 0
            while (p2 < 0 or p2 >= L) and tries < 10:
                p1 = int(rng.integers(0, max(L, 1)))
                if rng.random() < 0.5:
                    p2 = p1 + dist
                else:
                    p2 = p1 - dist
                tries += 1
            # clamp as last resort
            p2 = min(max(p2, 0), L - 1)

            s1, e1 = bin_interval(p1, bin_size)
            s2, e2 = bin_interval(p2, bin_size)
            # ensure end doesn't exceed chrom length (optional strictness)
            e1 = min(e1, L)
            e2 = min(e2, L)

            out["chrom1"].append(chrom)
            out["start1"].append(s1)
            out["end1"].append(e1)
            out["chrom2"].append(chrom)
            out["start2"].append(s2)
            out["end2"].append(e2)
        else:
            # trans: pick two *different* chromosomes weighted by size
            c1 = choose_weighted(chroms, weights, rng)
            c2 = choose_weighted(chroms, weights, rng)
            tries = 0
            while c2 == c1 and tries < 10:
                c2 = choose_weighted(chroms, weights, rng)
                tries += 1
            if c2 == c1:
                # in case only 1 chrom in file
                c2 = c1

            L1 = chrom_sizes[c1]
            L2 = chrom_sizes[c2]
            p1 = int(rng.integers(0, max(L1, 1)))
            p2 = int(rng.integers(0, max(L2, 1)))

            s1, e1 = bin_interval(p1, bin_size)
            s2, e2 = bin_interval(p2, bin_size)
            e1 = min(e1, L1)
            e2 = min(e2, L2)

            out["chrom1"].append(c1)
            out["start1"].append(s1)
            out["end1"].append(e1)
            out["chrom2"].append(c2)
            out["start2"].append(s2)
            out["end2"].append(e2)

        out["name"].append(f"readpair_{i}")
        out["score"].append(int(rng.integers(1, 1000)))  # placeholder "score"

    df = pd.DataFrame(out)

    # BEDPE convention often expects start<=end, and for each end start<end
    # Our binning guarantees start<=end; ensure non-negative.
    for col in ["start1", "end1", "start2", "end2"]:
        df[col] = df[col].clip(lower=0).astype(int)

    return df


def write_tsv(df: pd.DataFrame, path: str) -> None:
    if path.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as f:
            df.to_csv(f, sep="\t", index=False, header=False)
    else:
        df.to_csv(path, sep="\t", index=False, header=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom-sizes", required=True, help="Path to chrom.sizes (chrom<TAB>length).")
    ap.add_argument("--n", type=int, default=1_000_000, help="Number of pairs.")
    ap.add_argument("--cis-fraction", type=float, default=0.85, help="Fraction of cis contacts.")
    ap.add_argument("--bin-size", type=int, default=5000, help="Bin size (bp).")
    ap.add_argument("--cis-min-dist", type=int, default=1000, help="Min cis distance (bp).")
    ap.add_argument("--cis-max-dist", type=int, default=2_000_000, help="Max cis distance (bp).")
    ap.add_argument("--cis-alpha", type=float, default=1.4, help="Power-law exponent (>1). Larger => more short-range.")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed.")
    ap.add_argument("-o", "--out", default="sim.bedpe.gz", help="Output BEDPE (.tsv or .gz).")
    args = ap.parse_args()

    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    df = simulate_bedpe(
        chrom_sizes=chrom_sizes,
        n_pairs=args.n,
        cis_fraction=args.cis_fraction,
        bin_size=args.bin_size,
        cis_min_dist=args.cis_min_dist,
        cis_max_dist=args.cis_max_dist,
        cis_alpha=args.cis_alpha,
        seed=args.seed,
    )
    write_tsv(df, args.out)


if __name__ == "__main__":
    main()
