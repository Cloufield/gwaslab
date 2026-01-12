"""
MR-MEGA (Meta-Regression for Genome-wide Association Studies) Extension Module

This module provides the MR-MEGA algorithm for meta-regression of GWAS summary statistics.
MR-MEGA uses a Multidimensional Scaling (MDS) based approach with principal components
derived from EAF correlation distances between cohorts.

Example usage:
    from gwaslab.extension.mrmega import meta_regress_mrmega
    
    # Create SumstatsMulti object
    multi = gl.SumstatsMulti(sumstats_list)
    
    # Perform MR-MEGA meta-regression (MDS-based)
    result = meta_regress_mrmega(multi, num_pcs=2, use_genomic_control=True)
"""

from .mrmega import meta_regress_mrmega, _check_mrmega_available

__all__ = [
    "meta_regress_mrmega",
    "_check_mrmega_available",
]
