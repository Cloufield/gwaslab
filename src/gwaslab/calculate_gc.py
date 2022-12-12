import pandas as pd
import numpy as np
import scipy as sp

#20220312
def lambdaGC(insumstats,include_chrXYMT=[True,True,True],x=23 ,y=24, mt=25, mode="P",level=0.5,verbose=True):
    '''
    Calculate Genomic inflation factor for genomic control (LambdaGC)
    '''

    mode=mode.upper()
    sumstats=insumstats.loc[:,["CHR",mode]].dropna()
    
    if include_chrXYMT[0]!=True:
        sumstats = sumstats.loc[~sumstats["CHR"]==x,:]
    if include_chrXYMT[1]!=True:
        sumstats = sumstats.loc[~sumstats["CHR"]==y,:]
    if include_chrXYMT[2]!=True:
        sumstats = sumstats.loc[~sumstats["CHR"]==mt,:]

    indata = sumstats[mode]
    if len(indata) == 0: 
        if verbose: print("No available variants to use for calculation.")
        return np.nan   
    if mode=="p" or mode=="P":
        observedMedianChi2 = sp.stats.chi2.isf(np.median(indata),1)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose: print("(P mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="mlog10p" or mode=="MLOG10P":
        observedMedianChi2 = sp.stats.chi2.isf( np.median(np.power(10,-indata)) ,1)
        expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose: print("(MLOG10P mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="z" or mode=="Z":
        observedMedianChi2 = np.median((indata)**2)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:print("(Z mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="chi2" or mode=="CHISQ":
        observedMedianChi2 = np.median(indata)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:print("(CHISQ mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    else:
        return np.nan
    return lambdagc
