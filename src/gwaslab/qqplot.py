import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp

def qqplot(insumstats,p,scaled=False,title="QQ plot",figsize=(8,8)):
    sumstats=insumstats.loc[~insumstats[p].isna(),p].astype(np.float64)
    print("Processing "+ str(len(sumstats)) +" variants...")
    observedMedianChi2 = sp.stats.chi2.isf(np.median(sumstats),1)
    expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
    lambdagc=observedMedianChi2/expectedMedianChi2
    if scaled:
        print("P values are already -log10 scaled.")
        p_toplot = sumstats
    else:
        print("P values are being scaled to -log10(P)...")
        p_toplot = -np.log10(sumstats)    
    
    minit=1/len(sumstats)
    observed = p_toplot.sort_values(ascending=False)
    expected = -np.log10(np.linspace(minit,1,len(observed)))
    
    fig, axs = plt.subplots(1, 1,figsize=figsize)
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
    plt.scatter(expected,observed,s=3)
    
    plt.plot([0,-np.log10(minit)],[0,-np.log10(minit)],linestyle="--",color="grey")
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title(title)
    plt.text(0.05, 0.95,"Lambda GC = "+"{:.4f}".format(lambdagc),horizontalalignment='left',verticalalignment='top',transform=axs.transAxes)
    
def maf_qqplot():
    return 0

def gc(insumstats,mode="p",level=0.5):
    indata=insumstats.dropna()
    if mode=="p" or mode=="P":
        observedMedianChi2 = sp.stats.chi2.isf(np.median(indata),1)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        print("(p mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="z" or mode=="Z":
        observedMedianChi2 = np.median((indata)**2)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        print("(z mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    elif mode=="chi2" or mode=="chi2":
        observedMedianChi2 = np.median(indata)
        expectedMedianChi2 = sp.stats.chi2.ppf(level,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        print("(chi2 mode) Lambda GC at "+ str(level)+ " is"," ","{:.5f}".format(lambdagc))
    else:
        return None
    return lambdagc

    