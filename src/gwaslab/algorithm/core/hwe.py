"""Exact Hardy-Weinberg equilibrium test.
"""

from __future__ import annotations

import numpy as np


def snphwe(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
    """Compute exact Hardy-Weinberg equilibrium P-value for a SNP.

    Uses the Wigginton et al. (2005) exact test via a combinatorial
    recurrence (Python port of Jeremy McRae's C++ implementation).

Parameters
----------
obs_hets : int
    Observed heterozygote count.
obs_hom1 : int
    Observed homozygote count for the first allele.
obs_hom2 : int
    Observed homozygote count for the second allele.
Returns
-------
float
    Two-sided HWE P-value capped at 1.0.

References
----------
    Wigginton, J. E., Cutler, D. J., & Abecasis, G. R. (2005). A note on exact
    tests of Hardy-Weinberg equilibrium. American Journal of Human Genetics,
    76(5), 887-893.
"""
    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
        raise ValueError("snphwe: negative allele count")

    obs_homr = min(obs_hom1, obs_hom2)
    obs_homc = max(obs_hom1, obs_hom2)

    rare = 2 * obs_homr + obs_hets
    genotypes = obs_hets + obs_homc + obs_homr

    if genotypes == 0:
        raise ValueError("snphwe: zero genotypes")

    probs = np.zeros(rare + 1, dtype=float)

    mid = rare * (2 * genotypes - rare) // (2 * genotypes)
    if mid % 2 != rare % 2:
        mid += 1

    probs[mid] = 1.0
    sum_p = 1.0
    curr_homr = (rare - mid) // 2
    curr_homc = genotypes - mid - curr_homr

    for curr_hets in range(mid, 1, -2):
        probs[curr_hets - 2] = (
            probs[curr_hets]
            * curr_hets
            * (curr_hets - 1.0)
            / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        )
        sum_p += probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1

    curr_homr = (rare - mid) // 2
    curr_homc = genotypes - mid - curr_homr

    for curr_hets in range(mid, rare - 1, 2):
        probs[curr_hets + 2] = (
            probs[curr_hets]
            * 4.0
            * curr_homr
            * curr_homc
            / ((curr_hets + 2.0) * (curr_hets + 1.0))
        )
        sum_p += probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1

    target = probs[obs_hets]
    p_hwe = 0.0
    for p in probs:
        if p <= target:
            p_hwe += p / sum_p

    return min(p_hwe, 1.0)
