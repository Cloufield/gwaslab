"""Sliding-window signal density on sorted genomic positions.
"""

from __future__ import annotations

from typing import TypedDict

import numpy as np


class DensitySummary(TypedDict):
    mean: float
    median: float
    std: float
    max: float
    max_index: int


def density_all_variants(positions: np.ndarray, window_bp: int) -> np.ndarray:
    """Count variants within a symmetric window for each sorted position.

Parameters
----------
positions : numpy.ndarray
    Genomic positions sorted ascending on one chromosome.
window_bp : int
    Half-window size in base pairs.
Returns
-------
numpy.ndarray
    Variant counts within ``±window_bp``, excluding self.

Notes
-----
    Implementation note: GWASLab sliding-window density using ``searchsorted``.
    Orchestration lives in ``util.util_in_get_density``.
"""
    positions = np.asarray(positions, dtype=np.int64)
    n = positions.size
    if n == 0:
        return np.array([], dtype=np.int32)
    right_idx = np.searchsorted(positions, positions + window_bp, side="right")
    left_idx = np.searchsorted(positions, positions - window_bp, side="left")
    return (right_idx - left_idx - 1).astype(np.int32)


def density_from_signals(
    all_positions: np.ndarray,
    signal_positions: np.ndarray,
    window_bp: int,
) -> np.ndarray:
    """Increment density counts for variants near each signal position.

Parameters
----------
all_positions : numpy.ndarray
    All variant positions on one chromosome (sorted).
signal_positions : numpy.ndarray
    Significant variant positions (sorted not required).
window_bp : int
    Half-window size in base pairs.
Returns
-------
numpy.ndarray
    Count of signal-centered windows covering each variant.

Notes
-----
    Implementation note: GWASLab conditional density for ``sig_sumstats`` mode.
"""
    all_positions = np.asarray(all_positions, dtype=np.int64)
    signal_positions = np.asarray(signal_positions, dtype=np.int64)
    counts = np.zeros(all_positions.shape[0], dtype=np.int32)
    for sig_pos in signal_positions:
        left_idx = np.searchsorted(all_positions, sig_pos - window_bp, side="left")
        right_idx = np.searchsorted(all_positions, sig_pos + window_bp, side="right")
        counts[left_idx:right_idx] += 1
    return counts


def density_summary(counts: np.ndarray) -> DensitySummary:
    """Summarize density count distribution.

Parameters
----------
counts : numpy.ndarray
    Density values (may contain NaN).
Returns
-------
DensitySummary
    Mean, median, std, max, and argmax index.
"""
    valid = np.asarray(counts, dtype=float)
    valid = valid[np.isfinite(valid)]
    if valid.size == 0:
        return DensitySummary(mean=float("nan"), median=float("nan"), std=float("nan"), max=float("nan"), max_index=-1)
    max_index = int(np.nanargmax(valid))
    std = float(valid.std()) if valid.size > 1 else 0.0
    return DensitySummary(
        mean=float(np.nanmean(valid)),
        median=float(np.nanmedian(valid)),
        std=std,
        max=float(np.nanmax(valid)),
        max_index=max_index,
    )
