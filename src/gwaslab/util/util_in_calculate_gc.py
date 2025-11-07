import pandas as pd
import numpy as np
import scipy as sp
from gwaslab.g_Log import Log
#20220312
def lambdaGC(insumstats,include_chrXYMT=True, x=23 ,y=24, mt=25, mode="P",level=0.5,verbose=True,log=Log()):
    """
    Calculate the Genomic Inflation Factor (LambdaGC) for genomic control in GWAS.
    
    Parameters
    ----------
    insumstats : pandas.DataFrame
        Input summary statistics DataFrame containing 'CHR' and a column with test statistics
    include_chrXYMT : bool, optional
        If False, exclude sex chromosomes (X, Y) and mitochondrial (MT) from calculation
    x, y, mt : int or str, optional
        Identifiers for sex and mitochondrial chromosomes (default: 23, 24, 25)
    mode : {'P', 'MLOG10P', 'Z', 'CHISQ'}, optional
        Input data type to use for calculation:
        - 'P': p-values (default)
        - 'MLOG10P': -log10(p-values)
        - 'Z': Z-scores
        - 'CHISQ': Chi-squared statistics
    level : float, optional
        Quantile level for calculation (default: 0.5, median)
    verbose : bool, optional
        If True, write progress messages to log
    log : gwaslab.g_Log.Log object, optional
        Logging object for output messages
    
    Returns
    -------
    float
        Genomic inflation factor (LambdaGC), calculated as the ratio of observed to expected median chi-squared statistics
    
    References
    ----------
    Devlin, B., & Roeder, K. (1999). Genomic control for association studies. 
    Biometrics, 55(4), 964-975.
    """

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
    if mode=="p" or mode=="P":
        observedMedianChi2 = sp.stats.chi2.isf(np.nanmedian(indata),1)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        log.write(" -Lambda GC (P mode) at "+ str(1 - level)+ " is"," ","{:.5f}".format(lambdagc), verbose=verbose)
    elif mode=="mlog10p" or mode=="MLOG10P":
        observedMedianChi2 = sp.stats.chi2.isf( np.nanmedian(np.power(10,-indata)) ,1)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        log.write(" -Lambda GC (MLOG10P mode) at "+ str(1- level)+ " is"," ","{:.5f}".format(lambdagc), verbose=verbose)
    elif mode=="z" or mode=="Z":
        observedMedianChi2 = np.median((indata)**2)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:log.write(" -Lambda GC (Z mode) at "+ str(1- level)+ " is"," ","{:.5f}".format(lambdagc), verbose=verbose)
    elif mode=="chi2" or mode=="CHISQ":
        observedMedianChi2 = np.median(indata)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        log.write(" -Lambda GC (CHISQ mode) at "+ str(1- level)+ " is"," ","{:.5f}".format(lambdagc), verbose=verbose)
    else:
        return np.nan
    return lambdagc
