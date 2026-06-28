"""Pure GWAS calculation helpers for GWASLab.

Subpackages: ``core``, ``heterogeneity``, ``finemapping``, ``allele``,
``leads``, ``density``. See ``algorithm/README.md`` and ``algorithm/STYLE.md``.
"""

from gwaslab.algorithm.allele.complement import reverse_complement
from gwaslab.algorithm.core.genomic_control import (
    lambda_gc_from_chisq,
    lambda_gc_from_mlog10p,
    lambda_gc_from_p,
    lambda_gc_from_z,
)
from gwaslab.algorithm.core.hwe import snphwe
from gwaslab.algorithm.core.power import (
    genotype_disease_probabilities,
    infer_beta_range,
    infer_maf_range,
    or_to_rr,
    power_binary,
    power_contour_at_target,
    power_grid_binary,
    power_grid_quantitative,
    power_quantitative,
    resolve_genotype_rr,
)
from gwaslab.algorithm.density.signal import (
    density_all_variants,
    density_from_signals,
    density_summary,
)
from gwaslab.algorithm.finemapping.abf import log_abf_from_summary, pip_from_log_abf
from gwaslab.algorithm.finemapping.winners_curse import winners_curse_correct
from gwaslab.algorithm.heterogeneity.jackknife import jackknife_correlation_se

__all__ = [
    "reverse_complement",
    "lambda_gc_from_p",
    "lambda_gc_from_mlog10p",
    "lambda_gc_from_z",
    "lambda_gc_from_chisq",
    "snphwe",
    "or_to_rr",
    "resolve_genotype_rr",
    "genotype_disease_probabilities",
    "power_quantitative",
    "power_binary",
    "power_grid_quantitative",
    "power_grid_binary",
    "power_contour_at_target",
    "infer_maf_range",
    "infer_beta_range",
    "density_all_variants",
    "density_from_signals",
    "density_summary",
    "log_abf_from_summary",
    "pip_from_log_abf",
    "winners_curse_correct",
    "jackknife_correlation_se",
]
