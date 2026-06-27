"""Allele string operations."""

from __future__ import annotations

_COMPLEMENT = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})


def reverse_complement(allele: str) -> str:
    """
    Return the reverse-complement of a nucleotide allele string.

    Parameters
    ----------
    allele : str
        Allele sequence using A/C/G/T (case preserved through translation).

    Returns
    -------
    str
        Reverse-complemented allele.

    Notes
    -----
    Orchestration in ``qc.qc_fix_sumstats`` and ``hm.hm_harmonize_sumstats``
    delegates here.
    """
    return allele[::-1].translate(_COMPLEMENT)
