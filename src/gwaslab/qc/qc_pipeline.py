"""Unified QC pipeline shared by basic_check() and harmonize(basic_check=True).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, Optional

import gc

from gwaslab.io.io_process_kwargs import remove_overlapping_kwargs
from gwaslab.qc.qc_fix_sumstats import (
    _fix_allele,
    _fix_chr,
    _fix_ID,
    _fix_pos,
    _parallelize_normalize_allele,
    _remove_dup,
    _sort_column,
    _sort_coordinate,
)
from gwaslab.qc.qc_sanity_check import _check_data_consistency, _sanity_check_stats

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats


@dataclass
class QcRunReport:
    """Per-step execution record from _run_qc_core.
"""

    steps: Dict[str, bool] = field(default_factory=dict)
    parameters: Dict[str, Any] = field(default_factory=dict)

    def mark(self, step: str, ran: bool) -> None:
        self.steps[step] = bool(ran)


def _run_qc_core(
    sumstats_obj: "Sumstats",
    *,
    remove: bool = False,
    remove_dup: bool = False,
    normalize: bool = True,
    threads: int = 1,
    include_consistency: bool = True,
    verbose: bool = True,
    fix_id_kwargs: Optional[Dict[str, Any]] = None,
    fix_chr_kwargs: Optional[Dict[str, Any]] = None,
    fix_pos_kwargs: Optional[Dict[str, Any]] = None,
    fix_allele_kwargs: Optional[Dict[str, Any]] = None,
    sanity_check_stats_kwargs: Optional[Dict[str, Any]] = None,
    consistency_check_kwargs: Optional[Dict[str, Any]] = None,
    normalize_allele_kwargs: Optional[Dict[str, Any]] = None,
    remove_dup_kwargs: Optional[Dict[str, Any]] = None,
) -> QcRunReport:
    """Canonical QC step sequence used by basic_check and harmonize.

    Order: fix_id → fix_chr → mapper refresh → fix_pos → fix_allele → sanity
    → [consistency] → [normalize] → [remove_dup] → sort_coord → sort_column.
"""
    report = QcRunReport(
        parameters={
            "remove": remove,
            "remove_dup": remove_dup,
            "normalize": normalize,
            "threads": threads,
            "include_consistency": include_consistency,
            "fix_id_kwargs": fix_id_kwargs or {},
            "fix_chr_kwargs": fix_chr_kwargs or {},
            "fix_pos_kwargs": fix_pos_kwargs or {},
            "fix_allele_kwargs": fix_allele_kwargs or {},
            "sanity_check_stats_kwargs": sanity_check_stats_kwargs or {},
            "consistency_check_kwargs": consistency_check_kwargs or {},
            "normalize_allele_kwargs": normalize_allele_kwargs or {},
            "remove_dup_kwargs": remove_dup_kwargs or {},
        }
    )
    log = sumstats_obj.log

    fix_id_kwargs = remove_overlapping_kwargs(fix_id_kwargs or {}, {"log", "remove", "verbose"})
    sumstats_obj.data = _fix_ID(sumstats_obj, log=log, verbose=verbose, **fix_id_kwargs)
    report.mark("id", True)

    fix_chr_kwargs = remove_overlapping_kwargs(fix_chr_kwargs or {}, {"log", "remove", "verbose"})
    sumstats_obj.data = _fix_chr(sumstats_obj, log=log, remove=remove, verbose=verbose, **fix_chr_kwargs)
    report.mark("chr", True)
    sumstats_obj._update_mapper_from_data()

    fix_pos_kwargs = remove_overlapping_kwargs(fix_pos_kwargs or {}, {"log", "remove", "verbose"})
    sumstats_obj.data = _fix_pos(sumstats_obj, log=log, remove=remove, verbose=verbose, **fix_pos_kwargs)
    report.mark("pos", True)

    fix_allele_kwargs = remove_overlapping_kwargs(fix_allele_kwargs or {}, {"log", "remove", "verbose"})
    sumstats_obj.data = _fix_allele(sumstats_obj, log=log, remove=remove, verbose=verbose, **fix_allele_kwargs)
    report.mark("allele", True)

    sanity_check_stats_kwargs = remove_overlapping_kwargs(
        sanity_check_stats_kwargs or {}, {"log", "remove", "verbose"}
    )
    sumstats_obj.data = _sanity_check_stats(
        sumstats_obj, log=log, verbose=verbose, **sanity_check_stats_kwargs
    )
    report.mark("sanity", True)

    if include_consistency:
        consistency_check_kwargs = remove_overlapping_kwargs(
            consistency_check_kwargs or {}, {"log", "remove", "verbose"}
        )
        _check_data_consistency(sumstats_obj, log=log, verbose=verbose, **consistency_check_kwargs)
        report.mark("consistency", True)
    else:
        report.mark("consistency", False)

    if normalize:
        normalize_allele_kwargs = remove_overlapping_kwargs(
            normalize_allele_kwargs or {}, {"log", "remove", "verbose", "threads"}
        )
        sumstats_obj.data = _parallelize_normalize_allele(
            sumstats_obj, threads=threads, verbose=verbose, log=log, **normalize_allele_kwargs
        )
        report.mark("normalize", True)
    else:
        report.mark("normalize", False)

    if remove_dup:
        remove_dup_kwargs = remove_overlapping_kwargs(remove_dup_kwargs or {}, {"log", "remove", "verbose"})
        sumstats_obj.data = _remove_dup(sumstats_obj, log=log, verbose=verbose, **remove_dup_kwargs)
        report.mark("remove_dup", True)
    else:
        report.mark("remove_dup", False)

    sumstats_obj.data = _sort_coordinate(sumstats_obj, verbose=verbose, log=log)
    report.mark("sort_coord", True)

    sumstats_obj.data = _sort_column(sumstats_obj, verbose=verbose, log=log)
    report.mark("sort_column", True)

    gc.collect()
    return report
