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
from adjustText import adjust_text

def plottrumpet(mysumstats,
                snpid="SNPID",
                mode="q",
                anno=None,
                prevalence=None,
                scase=None,
                scontrol=None, 
                raweaf="EAF",
                eaf="MAF",
                p="P",
                sig_level=5e-8,
                p_level=5e-8,
                beta="BETA",
                anno_y = 1,
                anno_x = 0.01,
                n="N",
                eaf_range=None,
                beta_range=None, 
                ts=None,
                verbose=True,
                n_matrix=1000,
                xscale="log",
                yscale_factor=1,
                cmap="cool",
                markercolor="#349beb",
                fontsize=12,
                fontfamily="Arial",
                sizes=None,
                save=False,
                saveargs=None,
                ann_args=None,
                log=Log()):
    
    if sizes is None:
        sizes = (20,80)
    if ann_args is None:
        anno_args={"fontsize":10}
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
    
    cols_to_use = [beta,raweaf,n, p]
    if anno is not None:
        cols_to_use.append(anno)
    
    if p in mysumstats.columns:
        sumstats = mysumstats.loc[mysumstats[p]< p_level,cols_to_use ].copy()
        if verbose: log.write("Excluding variants with P values > {}".format(p_level))
    else:
        cols_to_use.remove(p)
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
            xpower2[1] = xpower2[1] * yscale_factor
            xpower[1] = xpower[1] * yscale_factor
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
            xpower2[1] = xpower2[1] * yscale_factor
            xpower[1] = xpower[1] * yscale_factor
            lines = LineCollection([xpower2,xpower], label=t,color=output_hex_colors[i])
            ax.add_collection(lines)
    
    sumstats["ABS_BETA"] = sumstats[beta].abs()
    

    sumstats[beta] = sumstats[beta]*yscale_factor

    sns.scatterplot(data=sumstats,
                    x=eaf,
                    y=beta,
                    size="ABS_BETA", 
                    ax=ax, 
                    sizes=sizes,
                    color=markercolor,
                    legend=False, 
                    alpha=0.6)
    


    if xscale== "log":
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

    if anno is not None:
        if anno in sumstats.columns:
            variants_toanno = sumstats.dropna(axis=0)
            variants_toanno = variants_toanno.loc[ variants_toanno[beta].abs() > anno_y,:]
            variants_toanno = variants_toanno.loc[ variants_toanno[eaf] < anno_x,:]
            texts_u=[]
            texts_d=[]
            for index, row in variants_toanno.iterrows():
                if row[beta] >0 :
                    texts_u.append(ax.annotate(row[anno], xy=(row[eaf], row[beta]),xytext=(row[eaf], row[beta]*1.1),arrowprops=dict(arrowstyle="-|>"),ha="left",va="bottom",fontsize=anno_args["fontsize"]))
                    #texts_u.append(plt.text(row[eaf], row[beta], row[anno],ha="right",va="bottom"))
                else:
                    texts_d.append(ax.annotate(row[anno], xy=(row[eaf], row[beta]),xytext=(row[eaf], row[beta]*1.1),arrowprops=dict(arrowstyle="-|>"),ha="left",va="top",fontsize=anno_args["fontsize"]))
                    #texts_d.append(plt.text(row[eaf], row[beta], row[anno],ha="left",va="top"))

            adjust_text(texts_u + texts_d, autoalign =False,
                        precision =0.01,lim=1000, 
                        expand_text=(0.5,0.5), 
                        expand_points=(0.5,0.5),
                        force_objects=(0.5,0.5), 
                        #arrowprops=dict(arrowstyle='-|>', color='grey'),
                        ax=ax)
            #adjust_text(texts_d, arrowprops=dict(arrowstyle='-|>', color='grey'),ax=ax)  
    
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

####################################################################
