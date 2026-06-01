from typing import TYPE_CHECKING, Union, Optional
import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_ancestry_ref import resolve_ancestry_af
from gwaslab.qc.qc_decorator import with_logging

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

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

    ancestry_af = resolve_ancestry_af(ancestry_af, build, log=log, verbose=verbose)

    ref_af = pd.read_csv(ancestry_af, sep="\t")

    data_af = pd.merge(sumstats[["CHR","POS","EA","NEA","EAF"]] ,ref_af,on=["CHR","POS"],how="inner")

    log.write(f"  -Estimating Fst using {len(data_af)} variants...", verbose=verbose)

    is_filp = data_af["EA"] == data_af["ALT"]
    data_af.loc[is_filp, ["EA","NEA"]] = data_af.loc[is_filp, ["NEA","EA"]]
    data_af.loc[is_filp, "EAF"] = 1 - data_af.loc[is_filp, "EAF"]

    headers = []
    pop_cols = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX',
        'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB', 'GWD', 'ESN', 'BEB', 'MSL',
        'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK', 'ASW', 'MXL', 'TSI',
        'GIH', 'EUR', 'EAS', 'AMR', 'SAS', 'AFR']
    for i in pop_cols:
        headers.append(f"FST_{i}")
        data_af[f"FST_{i}"] = data_af.apply(lambda x: calculate_fst(x["EAF"], x[i]), axis=1)

    mean_fst = data_af[headers].mean().sort_values()
    for i, value in mean_fst.items():
        log.write( f"  -{i} : {value}", verbose=verbose)

    closest_ancestry = mean_fst.idxmin()

    log.write(f"  -Closest Ancestry: {closest_ancestry.split('_')[1]}", verbose=verbose)
    log.write("Finished inferring ancestry.", verbose=verbose)
    return closest_ancestry.split("_")[1]

def calculate_fst(p_1: float, p_2: float) -> float:
    # https://bios1140.github.io/understanding-fst-the-fixation-index.html
    # calculate q1 and q2
    q_1 = 1 - p_1
    q_2 = 1 - p_2

    # calculate total allele frequency
    p_t = (p_1 + p_2)/2
    q_t = 1 - p_t

    # calculate expected heterozygosity
    # first calculate expected heterozygosity for the two populations
    # pop1
    hs_1 = 2*p_1*q_1
    # pop2
    hs_2 = 2*p_2*q_2
    # then take the mean of this
    hs = (hs_1 + hs_2)/2

    # next calculate expected heterozygosity for the metapopulations
    ht = 2*p_t*q_t

    # calculate fst
    fst = (ht - hs)/ht

    # return output
    return fst
