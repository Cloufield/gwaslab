import pandas as pd
import numpy as np
import scipy as sp
from gwaslab.Log import Log
#20220312
def lambdaGC(insumstats,include_chrXYMT=True, x=23 ,y=24, mt=25, mode="P",level=0.5,verbose=True,log=Log()):
    '''
    Calculate Genomic inflation factor for genomic control (LambdaGC)
    '''

    mode=mode.upper()
    sumstats=insumstats.loc[:,["CHR",mode]]
    
    if include_chrXYMT is False:
        if verbose: log.write(" -Excluding chrX, chrY, chrMT from lambda GC calculation.")
        xymt= [x,y,mt,"chrx","chry","chrmt","chrX","chrY","chrMT","chrM","M","x","y","mt","X","Y","MT"]
        sumstats = sumstats.loc[~sumstats["CHR"].isin(xymt),:]

    indata = sumstats[mode].values
    if len(indata) == 0: 
        if verbose: log.write("  -No available variants to use for calculation.")
        return np.nan   
    if mode=="p" or mode=="P":
        observedMedianChi2 = sp.stats.chi2.isf(np.nanmedian(indata),1)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose: log.write(" -Lambda GC (P mode) at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="mlog10p" or mode=="MLOG10P":
        observedMedianChi2 = sp.stats.chi2.isf( np.nanmedian(np.power(10,-indata)) ,1)
        expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose: log.write(" -Lambda GC (MLOG10P mode) at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="z" or mode=="Z":
        observedMedianChi2 = np.median((indata)**2)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:log.write(" -Lambda GC (Z mode) at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="chi2" or mode=="CHISQ":
        observedMedianChi2 = np.median(indata)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:log.write(" -Lambda GC (CHISQ mode) at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    else:
        return np.nan
    return lambdagc
