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
from gwaslab.getsig import annogene
from gwaslab.annotateplot import annotate_single

def plottrumpet(mysumstats,
                snpid="SNPID",
                mode="q",
                chrom="CHR",
                pos="POS",
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
                fontsize=15,
                font_family="Arial",
                sizes=None,
                save=False,
                saveargs=None,
                ann_args=None,
                build="99",
                anno_set=None,
                anno_alias=None,
                anno_d=None,
                anno_args=None,
                anno_source = "ensembl",
                log=Log()):
    
    matplotlib.rc('font', family=font_family)
    if sizes is None:
        sizes = (20,80)
    if ann_args is None:
        anno_args={"fontsize":12}
    if anno_set is None:
        anno_set=list()
    if anno_alias is None:
        anno_alias=dict()
    if anno_d is None:
        anno_d=dict()
    if anno_args is None:
        anno_args=dict()


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
    
    cols_to_use = [snpid, beta,raweaf,n, p]
    if anno is not None:
        if anno != "GENENAME":
            cols_to_use.append(anno)
        else:
            cols_to_use.append(pos)
            cols_to_use.append(chrom)
            
    
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
    
    sumstats["i"] = sumstats[eaf] * 5000

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
        if anno in sumstats.columns or anno=="GENENAME" :
            variants_toanno = sumstats.dropna(axis=0)
            variants_toanno = variants_toanno.loc[ variants_toanno[beta].abs() > anno_y,:]
            variants_toanno = variants_toanno.loc[ variants_toanno[eaf] < anno_x,:]
            if (variants_toanno.empty is not True) and anno=="GENENAME":
                variants_toanno = annogene(variants_toanno,
                                    id=snpid,
                                    chrom=chrom,
                                    pos=pos,
                                    log=log,
                                    build=build,
                                    source=anno_source,
                                    verbose=verbose).rename(columns={"GENE":"GENENAME"})

            texts_u=[]
            texts_d=[]
            variants_toanno["scaled_P"] = variants_toanno[beta]
            if len(variants_toanno)>0:
                offsety = 0.2 * max(variants_toanno[beta].abs().max(),1.5)
                offsetx = 0.01
                for index, row in variants_toanno.iterrows():

                    if row[beta] >0 :
                        texts_u.append(ax.annotate(row[anno], xy=(row[eaf], row[beta]),xytext=(row[eaf]+offsetx, row[beta]+offsety),arrowprops=dict(arrowstyle="-|>"),ha="left",va="bottom",fontsize=anno_args["fontsize"]))
                #        #texts_u.append(plt.text(row[eaf], row[beta], row[anno],ha="right",va="bottom"))
                    else:
                        texts_d.append(ax.annotate(row[anno], xy=(row[eaf], row[beta]),xytext=(row[eaf]+offsetx, row[beta]-offsety),arrowprops=dict(arrowstyle="-|>"),ha="left",va="top",fontsize=anno_args["fontsize"]))
                #        #texts_d.append(plt.text(row[eaf], row[beta], row[anno],ha="left",va="top"))

                adjust_text(texts_u + texts_d, 
                            autoalign =True,
                            precision =0.001,
                            lim=1000, 
                            expand_text=(0.5,0.5), 
                            expand_points=(0.5,0.5),
                            force_objects=(0.1,0.1), 
                            ax=ax)
    
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
