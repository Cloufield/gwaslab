from typing import TYPE_CHECKING, Union, Optional
import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_ancestry_ref import resolve_ancestry_af
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.algorithm.population.fst import fst_from_allele_frequencies as calculate_fst

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

SUPER_POPS = ["EUR", "EAS", "AMR", "SAS", "AFR"]
SUB_POPS = [
    "GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV",
    "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB",
    "JPT", "LWK", "ASW", "MXL", "TSI", "GIH",
]
ALL_POPS = SUB_POPS + SUPER_POPS


def _mean_fst_section(mean_fst: pd.Series, pops: list[str]) -> pd.Series:
    cols = [f"FST_{p}" for p in pops]
    return mean_fst.loc[cols].sort_values()


def _write_fst_section(
    log: Log,
    mean_fst: pd.Series,
    pops: list[str],
    title: str,
    verbose: bool,
) -> None:
    log.write(f"  -{title}", verbose=verbose)
    for col, value in _mean_fst_section(mean_fst, pops).items():
        log.write(f"  -{col} : {value}", verbose=verbose)


@with_logging(
        start_to_msg="infer ancestry based on Fst",
        finished_msg="inferring ancestry",
        start_cols=["CHR","POS","EA","NEA","EAF"],
        start_function="infer_ancestry",
        must_kwargs=["build"]
)
def _infer_ancestry(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    ancestry_af: Optional[str] = None,
    build: Optional[str] = None,
    _core: bool = False,
    log: Log = Log(),
    verbose: bool = True
) -> str:
    """Infer ancestry based on Fst values from effective allele frequencies.

    A high Fst value indicates that populations are genetically distinct. This function
    compares the effective allele frequencies from the sumstats with those from 1kg data
    to determine the closest ancestry. Inconsistency may suggest mislabeling of EAF.

Parameters
----------
sumstats_or_dataframe : Sumstats or pd.DataFrame
    Sumstats object or DataFrame to process.
ancestry_af : str, optional
    Path to allele frequency file, or keywords ``1kg_hm3_hg19_eaf`` / ``1kg_hm3_hg38_eaf``.
    If None, uses downloaded full PAN when available, otherwise the builtin core panel.
build : str, optional
    Genome build version. Options are "19" or "38". Required when ``ancestry_af`` is None.
_core : bool, optional
    Internal flag. If True, force use of the builtin core EAF panel and skip downloaded
    reference lookup.
verbose : bool, optional
    If True, write log messages. Default is True.
Returns
-------
str
    The closest ancestry determined by the minimum average Fst value, derived from
    the header name of the corresponding column.

Notes
-----
    This function internally uses `calculate_fst` to compute Fst values for each variant.
"""
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data

    ancestry_af = resolve_ancestry_af(
        ancestry_af, build, _core=_core, log=log, verbose=verbose
    )

    ref_af = pd.read_csv(ancestry_af, sep="\t")

    data_af = pd.merge(sumstats[["CHR","POS","EA","NEA","EAF"]] ,ref_af,on=["CHR","POS"],how="inner")

    log.write(f"  -Estimating Fst using {len(data_af)} variants...", verbose=verbose)

    is_filp = data_af["EA"] == data_af["ALT"]
    data_af.loc[is_filp, ["EA","NEA"]] = data_af.loc[is_filp, ["NEA","EA"]]
    data_af.loc[is_filp, "EAF"] = 1 - data_af.loc[is_filp, "EAF"]

    headers = [f"FST_{p}" for p in ALL_POPS]
    for pop in ALL_POPS:
        data_af[f"FST_{pop}"] = data_af.apply(
            lambda x, p=pop: calculate_fst(x["EAF"], x[p]), axis=1
        )

    mean_fst = data_af[headers].mean()
    _write_fst_section(
        log, mean_fst, SUPER_POPS, "Superpopulation (mean Fst):", verbose
    )
    _write_fst_section(
        log, mean_fst, SUB_POPS, "Population (mean Fst):", verbose
    )

    closest_super = _mean_fst_section(mean_fst, SUPER_POPS).idxmin()
    closest_pop = _mean_fst_section(mean_fst, SUB_POPS).idxmin()
    log.write(
        f"  -Closest superpopulation: {closest_super.split('_')[1]}",
        verbose=verbose,
    )
    log.write(
        f"  -Closest population: {closest_pop.split('_')[1]}",
        verbose=verbose,
    )

    closest_ancestry = mean_fst.idxmin()
    log.write("Finished inferring ancestry.", verbose=verbose)
    return closest_ancestry.split("_")[1]
