from typing import TYPE_CHECKING, Union, Optional
import pandas as pd
import numpy as np
import scipy as sp
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

#20220312
def _lambda_GC(insumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
                include_chrXYMT: bool = True, 
                x: Union[int, str] = 23,
                y: Union[int, str] = 24, 
                mt: Union[int, str] = 25, 
                mode: Optional[str] = None,
                level: float = 0.5,
                verbose: bool = True,
                log: Log = Log()) -> float:
    """Calculate the Genomic Inflation Factor (LambdaGC) for genomic control in GWAS.
    
Parameters
----------
insumstats_or_dataframe : Sumstats or pd.DataFrame
    Sumstats object or DataFrame to process. Can be a full DataFrame or a subset with CHR and mode columns.
include_chrXYMT : bool, optional
    If False, exclude sex chromosomes (X, Y) and mitochondrial (MT) from calculation
    x, y, mt : int or str, optional
    Identifiers for sex and mitochondrial chromosomes (default: 23, 24, 25)
mode : {'P', 'MLOG10P', 'Z', 'CHISQ'} or None, optional
    Input data type to use for calculation. If None, will auto-detect based on available columns:
    - 'P': p-values (default if available)
    - 'MLOG10P': -log10(p-values)
    - 'Z': Z-scores
    - 'CHISQ': Chi-squared statistics
level : float, optional default=0.5
    Quantile level for calculation, default value is 0.5 which is median
verbose : bool, optional
    If True, write progress messages to log
Returns
-------
float
    Genomic inflation factor (LambdaGC), calculated as the ratio of observed to expected median chi-squared statistics
    
References
----------
    Devlin, B., & Roeder, K. (1999). Genomic control for association studies. 
    Biometrics, 55(4), 964-975.
"""
    # Extract DataFrame if Sumstats object
    if isinstance(insumstats_or_dataframe, pd.DataFrame):
        insumstats = insumstats_or_dataframe
    else:
        # Assume it's a Sumstats object
        try:
            insumstats = insumstats_or_dataframe.data
        except AttributeError:
            # Not a Sumstats object, pass through as-is
            insumstats = insumstats_or_dataframe

    # Auto-detect mode if not provided
    if mode is None:
        if "P" in insumstats.columns:
            mode = "P"
        elif "Z" in insumstats.columns:
            mode = "Z"
        elif "CHISQ" in insumstats.columns:
            mode = "CHISQ"
        elif "MLOG10P" in insumstats.columns:
            mode = "MLOG10P"
        else:
            log.write("  -No available columns (P, Z, CHISQ, or MLOG10P) for calculation.", verbose=verbose)
            return np.nan

    mode=mode.upper()
    sumstats=insumstats.loc[:,["CHR",mode]]
    
    if include_chrXYMT is False:
        log.write(" -Excluding chrX, chrY, chrMT from lambda GC calculation.", verbose=verbose)
        xymt= [x,y,mt,"chrx","chry","chrmt","chrX","chrY","chrMT","chrM","M","x","y","mt","X","Y","MT"]
        sumstats = sumstats.loc[~sumstats["CHR"].isin(xymt),:]

    indata = sumstats[mode].values
    if len(indata) == 0:
        log.write("  -No available variants to use for calculation.", verbose=verbose)
        return np.nan

    from gwaslab.algorithm.core.genomic_control import (
        lambda_gc_from_chisq,
        lambda_gc_from_mlog10p,
        lambda_gc_from_p,
        lambda_gc_from_z,
    )

    if mode in ("p", "P"):
        lambdagc = lambda_gc_from_p(indata, quantile=level)
        log.write(" -Lambda GC (P mode) at " + str(1 - level) + " is", " ", "{:.5f}".format(lambdagc), verbose=verbose)
    elif mode in ("mlog10p", "MLOG10P"):
        lambdagc = lambda_gc_from_mlog10p(indata, quantile=level)
        log.write(" -Lambda GC (MLOG10P mode) at " + str(1 - level) + " is", " ", "{:.5f}".format(lambdagc), verbose=verbose)
    elif mode in ("z", "Z"):
        lambdagc = lambda_gc_from_z(indata, quantile=level)
        if verbose:
            log.write(" -Lambda GC (Z mode) at " + str(1 - level) + " is", " ", "{:.5f}".format(lambdagc), verbose=verbose)
    elif mode in ("chi2", "CHISQ"):
        lambdagc = lambda_gc_from_chisq(indata, quantile=level)
        log.write(" -Lambda GC (CHISQ mode) at " + str(1 - level) + " is", " ", "{:.5f}".format(lambdagc), verbose=verbose)
    else:
        return np.nan
    return lambdagc
