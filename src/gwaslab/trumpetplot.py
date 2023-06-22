import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
from gwaslab.Log import Log
from gwaslab.calculate_power import get_beta
from gwaslab.calculate_power import get_beta_binary
from gwaslab.fill import filldata
from matplotlib.collections import LineCollection
import matplotlib.colors as mc
import matplotlib

def plottrumpet(mysumstats,
                mode="q",
                prevalence=None,
                scase=None,
                scontrol=None, 
                raweaf="EAF",
                eaf="MAF",
                p="P",
                sig_level=5e-8,
                p_level=5e-8,
                beta="BETA",
                n="N",
                eaf_range=None,
                beta_range=None, 
                ts=None,
                verbose=True,
                n_matrix=1000,
                cmap="cool",
                markercolor="#349beb",
                fontsize=12,
                fontfamily="Arial",
                sizes=None,
                save=False,
                saveargs=None,
                log=Log()):
    
    if sizes is None:
        sizes = (20,80)

    matplotlib.rc('font', family=fontfamily) 
    if verbose: log.write("Start to create trumpet plot...")
    if (beta not in mysumstats.columns) or (raweaf not in mysumstats.columns):
        if verbose:
            log.write(" -No EAF or BETA columns. Skipping...")
        return None
    if mode=="b":
        if scase is None or scontrol is None:
            if verbose:
                log.write(" -No scase or scontrol. Skipping...")
            return None
        if prevalence is None:
                prevalence= scase / (scase + scontrol)

    if p in mysumstats.columns:
        sumstats = mysumstats.loc[mysumstats[p]< p_level, [beta,raweaf,n, p]].copy()
        if verbose: log.write("Excluding variants with P values > {}".format(p_level))
    else:
        sumstats = mysumstats[[beta,raweaf,n]].copy()

    if verbose: log.write("Plotting {} variants...".format(len(sumstats)))

    if eaf not in sumstats.columns:
        sumstats = filldata(sumstats,to_fill=["MAF"])
    
    if n == "N":
        n = sumstats["N"].median() 
    elif n == "max":
        n = sumstats["N"].max() 
    elif n == "min":
        n = sumstats["N"].min() 
    elif n == "median":
        n = sumstats["N"].median() 
    elif n == "mean":
        n = sumstats["N"].mean() 

    if verbose:
        if mode=="q":
            log.write("N for power calculation: {}".format(n))

    if eaf_range is None:
        eaf_range=(min(sumstats[eaf].min(),0.0001),0.5)
    if beta_range is None:
        if sumstats[beta].max()>3:
            beta_range=(0.0001,sumstats[beta].max())
        else:
            beta_range=(0.0001,3)

    if ts is None:
        ts=[0.3,0.5,0.8]

    cmap_to_use = plt.cm.get_cmap(cmap)
    if cmap_to_use.N >100:
        rgba = cmap_to_use(ts)
    else:
        rgba = cmap_to_use(range(len(ts)))
    
    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))
    output_hex_colors

    fig, ax = plt.subplots(figsize=(10,10))
    if mode=="q":
        for i,t in enumerate(ts):
            xpower = get_beta(mode="q",          
                            eaf_range=eaf_range,
                            beta_range=beta_range, 
                            n=n,
                            t=t,
                            sig_level=sig_level,
                            n_matrix=n_matrix)
            xpower2 = xpower.copy()
            xpower2[1] = -xpower2[1] 
            lines = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i],zorder=0)
            ax.add_collection(lines)
            #ax.plot(,label=t)
            #ax.plot(xpower[0],-xpower[1],label=t)
    else:
        for i,t in enumerate(ts):
            xpower = get_beta_binary(        
                            eaf_range=eaf_range,
                            beta_range=beta_range, 
                            prevalence=prevalence,
                            scase=scase, 
                            scontrol=scontrol, 
                            t=t,
                            sig_level=sig_level,
                            n_matrix=n_matrix)
            xpower2 = xpower.copy()
            xpower2[1] = -xpower2[1] 
            lines = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i])
            ax.add_collection(lines)
    
    sumstats["ABS_BETA"] = sumstats[beta].abs()
    
    sns.scatterplot(data=sumstats,
                    x=eaf,
                    y=beta,
                    size="ABS_BETA", 
                    ax=ax, 
                    sizes=sizes,
                    color=markercolor,
                    legend=False, 
                    alpha=0.6)
    
    ax.set_xscale('log')
    ax.set_xticks([0.001,0.01,0.05,0.1,0.2,0.5],[0.001,0.01,0.05,0.1,0.2,0.5],fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    leg = ax.legend(title="Power",fontsize =fontsize,title_fontsize=fontsize)
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    ax.axhline(y=0,color="grey",linestyle="dashed")
    ax.set_xlim(min(sumstats[eaf].min()/2,0.001/2),0.5)
    ax.set_ylabel("Effect size",fontsize=fontsize)
    ax.set_xlabel("Minor allele frequency",fontsize=fontsize)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    
    if save:
        if verbose: log.write("Saving plot:")
        if save==True:
            fig.savefig("./trumpet_plot.png",bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ "./trumpet_plot.png" + " successfully!" )
        else:
            fig.savefig(save,bbox_inches="tight",**saveargs)
            log.write(" -Saved to "+ save + " successfully!" )

    if verbose: log.write("Finished creating trumpet plot!")
    return fig