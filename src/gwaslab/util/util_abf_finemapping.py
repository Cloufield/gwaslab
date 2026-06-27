from typing import TYPE_CHECKING, Union, Optional, Tuple, Any
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_filter_value import _get_flanking_by_chrpos
from gwaslab.util.util_in_filter_value import _get_flanking_by_id

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

from gwaslab.algorithm.finemapping.abf import log_abf_from_summary, pip_from_log_abf


def calc_abf(
    insumstats: pd.DataFrame,
    w: float = 0.2,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs: Any
) -> pd.DataFrame:
    log.write("Start to calculate approximate Bayesian factor for {} variants".format(len(insumstats)), verbose=verbose)
    log.write(" - Priors for the standard deviation W of the effect size parameter β : {} ".format(w), verbose=verbose)
    insumstats = insumstats.copy()
    insumstats.loc[:, "log_ABF"] = log_abf_from_summary(
        insumstats["BETA"].to_numpy(),
        insumstats["SE"].to_numpy(),
        w=w,
    )
    return insumstats


def calc_PIP(
    insumstats: pd.DataFrame,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs: Any
) -> pd.DataFrame:
    insumstats = insumstats.copy()
    log.write("Start to calculate PIP for {} variants".format(len(insumstats)), verbose=verbose)
    insumstats.loc[:, "PIP"] = pip_from_log_abf(insumstats["log_ABF"].to_numpy())
    insumstats.loc[:, "log_PIP"] = np.log(insumstats["PIP"].clip(lower=np.finfo(float).tiny))
    return insumstats

def _abf_finemapping(
    insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    region: Optional[Tuple[int, int, int]] = None,
    chrpos: Optional[Tuple[int, int]] = None,
    snpid: Optional[str] = None,
    log: Log = Log(),
    **kwargs: Any
) -> pd.DataFrame:
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        insumstats = insumstats_or_dataframe.data.copy()

    if region is not None:
        region_data = insumstats[(insumstats["CHR"] == region[0]) & (insumstats["POS"] >= region[1]) & (insumstats["POS"] <= region[2])]
    elif chrpos is not None:
        region_data = _get_flanking_by_chrpos(insumstats, chrpos=chrpos,**kwargs)
    elif snpid is not None:
        region_data = _get_flanking_by_id(insumstats, snpid=snpid,**kwargs)

    region_data = calc_abf(region_data,log=log,**kwargs)
    region_data = calc_PIP(region_data,log=log,**kwargs)
    return region_data

def _make_cs(
    insumstats: pd.DataFrame,
    threshold: float = 0.95,
    log: Log = Log(),
    verbose: bool = True
) -> pd.DataFrame:
    insumstats = insumstats.sort_values(by="PIP",ascending=False)
    pip_sum = 0
    cs = pd.DataFrame()
    for index, row in insumstats.iterrows():
        cs = pd.concat([cs,pd.DataFrame(row).T])
        pip_sum += row["PIP"]
        if pip_sum > threshold:
            break
    log.write("Finished constructing a {}% credible set with {} variant(s)".format(str(threshold * 100),str(len(cs))),verbose=verbose)
    return cs
