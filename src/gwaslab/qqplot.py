import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
def qqplot(insumstats,
           p,
           title="QQ plot",
           figsize=(8,8),
           fontsize=15,
           color="#000042",
           s=10,
           gc=True,
           save=None,
           saveargs={"dpi":300,"facecolor":"white"},
           sig_line_color="grey",
           verbose=True,
           scaled=False,
           rtspine=True
          ):
    sumstats=insumstats.loc[~insumstats[p].isna(),p].astype(np.float64)
    if verbose: print(" - Processing "+ str(len(sumstats)) +" non-missing variants...")
    
    if scaled:
        if verbose:print(" - P values are already -log10 scaled.")
        p_toplot = sumstats
    else:
        if verbose:print(" - P values are being scaled to -log10(P)...")
        p_toplot = -np.log10(sumstats)    

    # sort x,y for qq plot
    minit=1/len(p_toplot)
    observed = p_toplot.sort_values(ascending=False)
    expected = -np.log10(np.linspace(minit,1,len(observed)))
    
    fig, ax2 = plt.subplots(1, 1,figsize=figsize)
    ax2.scatter(expected,observed,s=s,color=color)
    ax2.plot([0,-np.log10(minit)],[0,-np.log10(minit)],linestyle="--",color=sig_line_color)
    ax2.set_xlabel("Expected -log10(p)",fontsize=fontsize)
    ax2.set_ylabel("Observed -log10(p)",fontsize=fontsize)
    if not rtspine:
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(True)
    
    if gc:
        observedMedianChi2 = sp.stats.chi2.isf( np.median(np.power(10,-p_toplot)) ,1)
        expectedMedianChi2 = sp.stats.chi2.ppf(0.5,1)
        lambdagc=observedMedianChi2/expectedMedianChi2
        ax2.text(0.05, 0.95,"Lambda GC = "+"{:.4f}".format(lambdagc),
                 horizontalalignment='left',
                 verticalalignment='top',
                 transform=ax2.transAxes,
                 fontsize=fontsize)

    if verbose: print(" - Created QQ plot successfully!")
    if title:
        plt.title(title,fontsize=fontsize,pad=10)
    if save:
        if save==True:
            plt.savefig("./qqplot.png",bbox_inches="tight",**saveargs)
            if verbose: print(" - Saved to "+ "./qqplot.png" + " successfully!" )
        else:    
            plt.savefig(save,bbox_inches="tight",**saveargs)
            if verbose: print(" - Saved to "+ save + " successfully!" )
    return ax2
    
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

    