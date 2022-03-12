import pandas as pd
import numpy as np
import scipy as sp

#20220312
def gc(insumstats,mode="p",level=0.5,verbose=True):
    indata=insumstats.dropna()
    if mode=="p" or mode=="P":
        observedMedianChi2 = sp.stats.chi2.isf(np.median(indata),1)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose: print("(p mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="z" or mode=="Z":
        observedMedianChi2 = np.median((indata)**2)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:print("(z mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="chi2" or mode=="chi2":
        observedMedianChi2 = np.median(indata)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        if verbose:print("(chi2 mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    else:
        return None
    return lambdagc

    