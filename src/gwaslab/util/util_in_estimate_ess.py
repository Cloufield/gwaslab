from typing import TYPE_CHECKING, Union, Optional
import numpy as np
import pandas as pd
from scipy.stats import norm
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

@with_logging(
        start_to_msg="estimate effective sample size (N_EFF)",
        finished_msg="estimating effective sample size (N_EFF)"
)
def _get_ess(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    method: Union[str, float] = "metal",
    log: Log = Log(),
    verbose: bool = True
) -> pd.DataFrame:
    """
    Estimate effective sample size (N_EFF) for GWAS summary statistics. Summary statistics DataFrame containing N_CASE and N_CONTROL columns.
    
    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    method : str or float, optional
        Method for ESS calculation:
        - "metal": Uses formula from Willer et al. (2010)
        - float: Directly uses the provided value
    
    Returns
    -------
    pandas.DataFrame
        Modified sumstats DataFrame with N_EFF column added
    
    References
    ----------
    Willer, C. J., Li, Y., & Abecasis, G. R. (2010). 
    METAL: fast and efficient meta-analysis of genomewide association scans. 
    Bioinformatics, 26(17), 2190-2191.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    if type(method) is str:
        if method =="metal":
            log.write(" - Method: {} ".format(method), verbose=verbose)
            log.write(" - Referencec: {} ".format("Willer, C. J., Li, Y., & Abecasis, G. R. (2010)"), verbose=verbose)
            log.write(" - Equation: {} ".format(" N_EFF = 4 * N_CASE * N_CONTROL / (N_CASE + N_CONTROL)"), verbose=verbose)
            # Willer, C. J., Li, Y., & Abecasis, G. R. (2010). METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics, 26(17), 2190-2191.
            sumstats["N_EFF"] =  4 / (1/sumstats["N_CASE"] + 1/sumstats["N_CONTROL"])
    else:
        sumstats["N_EFF"] =  method
    return sumstats
