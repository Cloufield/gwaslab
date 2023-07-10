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
from gwaslab.textreposition import adjust_text_position
from gwaslab.figuresave import save_figure

def plottrumpet(mysumstats,
                snpid="SNPID",
                mode="q",
                chrom="CHR",
                pos="POS",
                n="N",
                p="P",
                eaf="MAF",
                raweaf="EAF",
                beta="BETA",
                ts=None,
                anno=None,
                prevalence=None,
                scase=None,
                scontrol=None, 
                sig_level=5e-8,
                p_level=5e-8,
                anno_y = 1,
                anno_x = 0.01,                
                eaf_range=None,
                beta_range=None, 
                n_matrix=1000,
                xscale="log",
                yscale_factor=1,
                cmap="Reds",
                ylim=None,
                markercolor="#597FBD",
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
                anno_style="expand",
                anno_source = "ensembl",
                anno_max_iter=100,
                arm_scale=1,
                repel_force=0.05,
                ylabel="Effect size",
                xlabel="Minor allele frequency",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                sort="beta",
                verbose=True,
                log=Log()):
    
    #Checking columns#################################################################################################################
    matplotlib.rc('font', family=font_family)
    if sizes is None:
        sizes = (20,80)
    if ann_args is None:
        anno_args={"fontsize":12,"fontstyle":"italic"}
    if anno_set is None:
        anno_set=list()
    if anno_alias is None:
        anno_alias=dict()
    if anno_d is None:
        anno_d=dict()
    if xticks is None:
        if xscale== "log":
            xticks = [0.001,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks
        else:
            xticks = [0,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks            

    #Checking columns#################################################################################################################
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
                prevalence= scase/(scase + scontrol)
                log.write(" -Prevalence is not given. Estimating based on scase and scontrol :{}...".format(prevalence))

    #loading columns #################################################################################################################
    cols_to_use = [snpid, beta,raweaf,n, p]
    if anno is not None:
        if anno != "GENENAME":
            if anno!=True:
                log.write(" -Loading column {} for annotation...".format(anno))
                cols_to_use.append(anno)
        else:
            cols_to_use.append(pos)
            cols_to_use.append(chrom)
            
    #filter by p #################################################################################################################
    if p in mysumstats.columns:
        sumstats = mysumstats.loc[mysumstats[p]< p_level,cols_to_use ].copy()
        if verbose: log.write("Excluding variants with P values > {}".format(p_level))
    else:
        cols_to_use.remove(p)
        sumstats = mysumstats[[beta,raweaf,n]].copy()
    if verbose: log.write("Plotting {} variants...".format(len(sumstats)))
    
    #add maf column #################################################################################################################
    if eaf not in sumstats.columns:
        sumstats = filldata(sumstats,to_fill=["MAF"])
    
    #configure n #################################################################################################################
    if mode=="q":
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
            log.write("N for power calculation: {}".format(n))

    #configure beta and maf range ###################################################################################################
    if eaf_range is None:
        eaf_min_power = np.floor( -np.log10(sumstats[eaf].min())) + 1
        eaf_range=(min(np.power(10.0,-eaf_min_power),np.power(10.0,-4)),0.5)
    if beta_range is None:
        if sumstats[beta].max()>3:
            beta_range=(0.0001,sumstats[beta].max())
        else:
            beta_range=(0.0001,3)
    
    #configure power threshold###################################################################################################
    if ts is None:
        ts=[0.3,0.5,0.8]
    
    #configure colormap##########################################################################################################
    cmap_to_use = plt.cm.get_cmap(cmap)
    if cmap_to_use.N >100:
        rgba = cmap_to_use(ts)
    else:
        rgba = cmap_to_use(range(len(ts)))
    
    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))
    output_hex_colors

    ##################################################################################################
    fig, ax = plt.subplots(figsize=(10,10))
    
    ##creating power line############################################################################################
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
    ###################################################################################################
    # get abs  and convert using scaling factor
    sumstats[beta] = sumstats[beta]*yscale_factor
    sumstats["ABS_BETA"] = sumstats[beta].abs()

    ##################################################################################################
    
    sns.scatterplot(data=sumstats,
                    x=eaf,
                    y=beta,
                    size="ABS_BETA", 
                    ax=ax, 
                    sizes=sizes,
                    color=markercolor,
                    legend=False, 
                    edgecolor="black",
                    alpha=0.8)
    
    ##################################################################################################

    ax.tick_params(axis='y', labelsize=fontsize)
    leg = ax.legend(title="Power",fontsize =fontsize,title_fontsize=fontsize)

    for line in leg.get_lines():
        line.set_linewidth(5.0)
    ax.axhline(y=0,color="grey",linestyle="dashed")
    
    if xscale== "log":
        ax.set_xscale('log')
        rotation=0
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(min(sumstats[eaf].min()/2,0.001/2),0.52)
    else:
        rotation=90    
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(-0.02,0.52)
    
    if ylim is not None:
        ax.set_ylim(ylim)
    if yticks is not None:
        ax.set_yticks(yticks, yticklabels)

    ax.set_ylabel(ylabel,fontsize=fontsize)
    ax.set_xlabel(xlabel,fontsize=fontsize)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)

    ############  Annotation ##################################################################################################
    if anno is not None:
        if anno in sumstats.columns or anno=="GENENAME" :
            variants_toanno = sumstats.copy()
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

            if len(variants_toanno)>0:
                
                maxy = max(variants_toanno[beta].abs().max(),1.5)
                variants_toanno["ADJUSTED_i"] = np.nan 
                y_span = 0.5
                
                if sort == "beta" : 
                    variants_toanno = variants_toanno.sort_values(by=beta, key= np.abs, ascending = False)
                else:
                    variants_toanno = variants_toanno.sort_values(by=eaf, key= np.abs, ascending = True)
                
                if anno_style == "expand":
                    if len(variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"])>1:
                        variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"] = adjust_text_position(variants_toanno.loc[variants_toanno[beta]>0,eaf].values.copy(), 
                                                                                y_span, 
                                                                                repel_force=repel_force,
                                                                                max_iter=anno_max_iter,
                                                                                log=log,
                                                                                amode=xscale,
                                                                                verbose=verbose)

                    if len(variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"])>1:
                        variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"] = adjust_text_position(variants_toanno.loc[variants_toanno[beta]<0,eaf].values.copy(), 
                                                                y_span, 
                                                                repel_force=repel_force,
                                                                max_iter=anno_max_iter,
                                                                log=log,
                                                                amode=xscale,
                                                                verbose=verbose)

                last_pos = min(variants_toanno[eaf])/2
                for index, row in variants_toanno.iterrows():
                    
                    armB_length_in_point = ax.transData.transform((0,1.1*maxy))[1]-ax.transData.transform((0, abs(row[beta])))[1]
                    armB_length_in_point = armB_length_in_point*arm_scale

                    if anno_style == "right" :
                        #right style
                        if row[eaf]>last_pos*(repel_force+1):
                            last_pos=row[eaf]
                        else:
                            last_pos*= (repel_force+1)
                    elif anno_style == "expand" :
                        last_pos = row["ADJUSTED_i"]

                    if anno_style == "right"  or anno_style == "expand":
                        if row[beta] >0 :
                            texts_u.append(ax.annotate(row[anno], xy=(row[eaf], row[beta]),xytext=(last_pos , 1.2*maxy),
                                                    arrowprops=dict(relpos=(0,0),arrowstyle="-|>",facecolor='black',connectionstyle="arc,angleA=-90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)),rotation=90,
                                                    ha="left",va="bottom",**anno_args))
                        else:
                            texts_d.append(ax.annotate(row[anno], xy=(row[eaf], row[beta]),xytext=(last_pos , -1.2*maxy),
                                                    arrowprops=dict(relpos=(0,1),arrowstyle="-|>",facecolor='black',connectionstyle="arc,angleA=90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)),rotation=90,
                                                    ha="left",va="top",**anno_args))
                    
                    if anno_style=="tight":
                        texts_d.append(ax.text(row[eaf], row[beta], row[anno]))
                        adjust_text(texts_d, 
                                    autoalign =True,
                                    precision =0.001,
                                    lim=1000, 
                                    expand_text=(0.5,0.5), 
                                    expand_points=(0.5,0.5),
                                    force_objects=(0.5,0.5), 
                                    ax=ax)
    ############  Annotation ##################################################################################################
    if mode=="q":
        save_figure(fig, save, keyword="trumpet_q",saveargs=saveargs, log=log, verbose=verbose)
    elif mode=="b":
        save_figure(fig, save, keyword="trumpet_b",saveargs=saveargs, log=log, verbose=verbose)

    if verbose: log.write("Finished creating trumpet plot!")
    return fig

####################################################################
