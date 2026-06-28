"""Lead variant clustering by genomic distance.
"""

from __future__ import annotations

import numpy as np


def cluster_ids_from_positions(
    chrom: np.ndarray,
    pos: np.ndarray,
    window_bp: int,
) -> np.ndarray:
    """Assign cluster IDs from chromosome blocks and distance gaps.

Parameters
----------
chrom : numpy.ndarray
    Chromosome labels in row order.
pos : numpy.ndarray
    Positions in row order.
window_bp : int
    Maximum within-cluster distance in base pairs.
Returns
-------
numpy.ndarray
    Integer cluster ID per row (0-based cumulative breaks).

Notes
-----
    Implementation note: GWASLab lead clustering. Pandas groupby orchestration
    lives in ``util.util_in_get_sig._collect_leads_generic``.
"""
    chrom = np.asarray(chrom)
    pos_arr = np.asarray(pos)
    if pos_arr.size and np.issubdtype(pos_arr.dtype, np.floating):
        if np.isnan(pos_arr).any():
            raise ValueError(
                "cluster_ids_from_positions: position array contains NaN; "
                "filter invalid coordinates before clustering"
            )
    pos = np.asarray(pos_arr, dtype=np.int64)
    n = pos.size
    if n == 0:
        return np.array([], dtype=np.int64)

    new_chr = np.empty(n, dtype=bool)
    new_chr[0] = True
    new_chr[1:] = chrom[1:] != chrom[:-1]

    pos_diff = np.zeros(n, dtype=np.int64)
    pos_diff[1:] = pos[1:] - pos[:-1]
    pos_diff[new_chr] = 0

    breaks = new_chr | (pos_diff > window_bp)
    return np.cumsum(breaks) - 1
