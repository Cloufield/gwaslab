"""Palindromic SNP detection."""

from __future__ import annotations

import numpy as np


def is_palindromic_pair(ea: str, nea: str) -> bool:
    """
    Test whether an allele pair is strand-ambiguous (A/T or C/G).

    Parameters
    ----------
    ea : str
        Effect allele.
    nea : str
        Non-effect allele.

    Returns
    -------
    bool
        True if the pair is A/T, T/A, C/G, or G/C.
    """
    pair = {ea.upper(), nea.upper()}
    return pair == {"A", "T"} or pair == {"C", "G"}


def is_palindromic_alleles(
    effect_allele: np.ndarray,
    other_allele: np.ndarray,
) -> np.ndarray:
    """
    Vectorized palindromic test for allele columns.

    Parameters
    ----------
    effect_allele : numpy.ndarray
        Effect alleles (string-like).
    other_allele : numpy.ndarray
        Non-effect alleles (string-like).

    Returns
    -------
    numpy.ndarray
        Boolean array indicating palindromic pairs.
    """
    ea = np.asarray(effect_allele, dtype=str)
    nea = np.asarray(other_allele, dtype=str)
    ea_u = np.char.upper(ea)
    nea_u = np.char.upper(nea)
    gc = (ea_u == "G") & (nea_u == "C")
    cg = (ea_u == "C") & (nea_u == "G")
    at = (ea_u == "A") & (nea_u == "T")
    ta = (ea_u == "T") & (nea_u == "A")
    return gc | cg | at | ta
