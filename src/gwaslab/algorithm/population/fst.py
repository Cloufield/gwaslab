"""Population genetics statistics."""

from __future__ import annotations


def fst_from_allele_frequencies(p_1: float, p_2: float) -> float:
    """
    Compute Wright's F_ST between two populations from allele frequencies.

    Parameters
    ----------
    p_1 : float
        Allele frequency in population 1.
    p_2 : float
        Allele frequency in population 2.

    Returns
    -------
    float
        F_ST estimate in [0, 1] for valid inputs with ``ht > 0``.

    Notes
    -----
    Orchestration lives in ``util.util_ex_infer_ancestry``.

    References
    ----------
    See fixation-index tutorials such as BIOS1140 F_ST notes for the
    heterozygosity formulation used here.
    """
    q_1 = 1.0 - p_1
    q_2 = 1.0 - p_2
    p_t = (p_1 + p_2) / 2.0
    q_t = 1.0 - p_t
    hs_1 = 2.0 * p_1 * q_1
    hs_2 = 2.0 * p_2 * q_2
    hs = (hs_1 + hs_2) / 2.0
    ht = 2.0 * p_t * q_t
    if ht == 0:
        return float("nan")
    return (ht - hs) / ht
