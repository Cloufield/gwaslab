"""Effective sample size formulas."""

from __future__ import annotations

import numpy as np


def metal_effective_sample_size(
    n_case: np.ndarray,
    n_control: np.ndarray,
) -> np.ndarray:
    """
    Compute METAL effective sample size for case-control studies.

    Parameters
    ----------
    n_case : numpy.ndarray
        Case counts per variant or study row.
    n_control : numpy.ndarray
        Control counts per variant or study row.

    Returns
    -------
    numpy.ndarray
        Effective sample sizes ``N_EFF``.

    Notes
    -----
    Orchestration (column attach, logging) lives in ``util.util_in_estimate_ess``.

    References
    ----------
    Willer, C. J., Li, Y., & Abecasis, G. R. (2010). METAL: fast and efficient
    meta-analysis of genomewide association scans. Bioinformatics, 26(17),
    2190-2191.
    """
    n_case = np.asarray(n_case, dtype=np.float64)
    n_control = np.asarray(n_control, dtype=np.float64)
    return 4.0 / (1.0 / n_case + 1.0 / n_control)
