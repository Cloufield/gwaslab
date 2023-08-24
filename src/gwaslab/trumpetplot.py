import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
from gwaslab.Log import Log
from gwaslab.calculate_power import get_beta
from gwaslab.calculate_power import get_power
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
                maf="MAF",
                eaf="EAF",
                beta="BETA",
                ts=None,
                anno=None,
                prevalence=None,
                or_to_rr=False,
                ncase=None,
                ncontrol=None, 
                sig_level=5e-8,
                p_level=5e-8,
                anno_y = 1,
                anno_x = 0.01,                
                maf_range=None,
                beta_range=None, 
                n_matrix=1000,
                xscale="log",
                yscale_factor=1,
                cmap="cool",
                ylim=None,
                xlim=None,
                markercolor="#597FBD",
                fontsize=15,
                font_family="Arial",
                size= "ABS_BETA",
                sizes=None,
                save=False,
                saveargs=None,
                figargs=None,
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
    if anno_args is None:
        anno_args={"fontsize":12,"fontstyle":"italic"}
    if anno_set is None:
        anno_set=list()
    if anno_alias is None:
        anno_alias=dict()
    if anno_d is None:
        anno_d=dict()
    if ts is None:
        ts = [0.2,0.4,0.6,0.8]
    if xticks is None:
        if xscale== "log":
            xticks = [0.001,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks
        else:
            xticks = [0,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks            
    if figargs is None:
        figargs={"figsize":(10,15)}
    #Checking columns#################################################################################################################
    if verbose: log.write("Start to create trumpet plot...")
    
    #parameter check##################################################################################################################
    if (beta not in mysumstats.columns) or (eaf not in mysumstats.columns):
        if verbose:
            log.write(" -No EAF or BETA columns. Skipping...")
        return None
    if mode=="b":
        if ncase is None or ncontrol is None:
            if verbose:
                log.write(" -No scase or scontrol. Skipping...")
            return None
        if prevalence is None:
                prevalence= ncase/(ncase + ncontrol)
                log.write(" -Prevalence is not given. Estimating based on scase and scontrol :{}...".format(prevalence))
    
    #print settings##################################################################################################################
    if verbose:
        log.write(" -Settings:")
        log.write("  -Mode: {}".format(mode))
        if mode == "q" :
            log.write("  -N: {}".format(n))
        if mode == "b" :
            log.write("  -N_CASE: {}".format(ncase))
            log.write("  -N_CONTROL: {}".format(ncontrol))
            log.write("  -PREVALENCE: {}".format(prevalence))
        log.write("  -BETA: {}".format(beta))
        log.write("  -Significance level: {}".format(sig_level))
        log.write("  -Power thresholds: {}".format(ts))
        log.write("  -Power line smoothness: {}".format(n_matrix))
    
    #loading columns #################################################################################################################
    cols_to_use = [snpid, beta, eaf, n, p]
    
    if anno is not None:
        if anno != "GENENAME":
            if anno!=True:
                log.write(" -Loading column {} for annotation...".format(anno))
                if anno not in cols_to_use:
                    cols_to_use.append(anno)
        else:
            cols_to_use.append(pos)
            cols_to_use.append(chrom)

    if size != "ABS_BETA":
        if size not in cols_to_use:
            cols_to_use.append(size)
    #filter by p #################################################################################################################
    if p in mysumstats.columns:
        sumstats = mysumstats.loc[mysumstats[p]< p_level,cols_to_use ].copy()
        if verbose: log.write(" -Excluding variants with P values > {}".format(p_level))
    else:
        cols_to_use.remove(p)
        sumstats = mysumstats[[beta,eaf,n]].copy()
    if verbose: log.write(" -Plotting {} variants...".format(len(sumstats)))
    
    #add maf column #################################################################################################################
    if maf not in sumstats.columns:
        sumstats = filldata(sumstats,to_fill=["MAF"],verbose=False)
        is_filpped = (sumstats["MAF"] < sumstats[eaf]) & (sumstats[eaf] > 0.5)& (sumstats["MAF"] < 0.5)
        if verbose: log.write(" -Flipping {} variants...".format(sum(is_filpped)))
        sumstats.loc[is_filpped, beta] = -sumstats.loc[is_filpped, beta]
    
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
            log.write(" -N for power calculation: {}".format(n))

    #configure beta and maf range ###################################################################################################
    if maf_range is None:
        maf_min_power = np.floor( -np.log10(sumstats[maf].min())) + 1
        maf_range=(min(np.power(10.0,-maf_min_power),np.power(10.0,-4)),0.5)
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
    fig, ax = plt.subplots(**figargs)
    
    ##creating power line############################################################################################
    if mode=="q":
        for i,t in enumerate(ts):
            xpower = get_beta(mode="q",          
                            eaf_range=maf_range,
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
                            eaf_range=maf_range,
                            beta_range=beta_range, 
                            prevalence=prevalence,
                            or_to_rr = or_to_rr,
                            ncase=ncase, 
                            ncontrol=ncontrol, 
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
    
    dots = sns.scatterplot(data=sumstats,
                    x=maf,
                    y=beta,
                    size=size, 
                    ax=ax, 
                    sizes=sizes,
                    color=markercolor,
                    legend=True, 
                    edgecolor="black",
                    alpha=0.8)
    
    #second_legend = ax.legend(title="Power", loc="upper right",fontsize =fontsize,title_fontsize=fontsize)
    
    h,l = ax.get_legend_handles_labels()
    if len(ts)>0:
        l1 = ax.legend(h[:int(len(ts))],l[:int(len(ts))], title="Power", loc="upper right",fontsize =fontsize,title_fontsize=fontsize)
        for line in l1.get_lines():
            line.set_linewidth(5.0)
    l2 = ax.legend(h[int(len(ts)):],l[int(len(ts)):], title=size, loc="lower right",fontsize =fontsize,title_fontsize=fontsize)
    if len(ts)>0:
        ax.add_artist(l1)
    #first_legend = ax.legend(handles=dots, loc="lower right" ,title=size,fontsize =fontsize,title_fontsize=fontsize)
    #ax.add_artist(first_legend)
    ##################################################################################################

    ax.tick_params(axis='y', labelsize=fontsize)

    ax.axhline(y=0,color="grey",linestyle="dashed")
    
    if xscale== "log":
        ax.set_xscale('log')
        rotation=0
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(min(sumstats[maf].min()/2,0.001/2),0.52)
    else:
        rotation=90    
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(-0.02,0.52)
        
    if xlim is not None:
        ax.set_xlim(xlim)
    
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
            variants_toanno = variants_toanno.loc[ variants_toanno[maf] < anno_x,:]
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
                    variants_toanno = variants_toanno.sort_values(by=maf, key= np.abs, ascending = True)
                
                if anno_style == "expand":
                    if len(variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"])>1:
                        variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"] = adjust_text_position(variants_toanno.loc[variants_toanno[beta]>0,maf].values.copy(), 
                                                                                y_span, 
                                                                                repel_force=repel_force,
                                                                                max_iter=anno_max_iter,
                                                                                log=log,
                                                                                amode=xscale,
                                                                                verbose=verbose)

                    if len(variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"])>1:
                        variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"] = adjust_text_position(variants_toanno.loc[variants_toanno[beta]<0,maf].values.copy(), 
                                                                y_span, 
                                                                repel_force=repel_force,
                                                                max_iter=anno_max_iter,
                                                                log=log,
                                                                amode=xscale,
                                                                verbose=verbose)

                
                for variants_toanno_half in [variants_toanno.loc[variants_toanno[beta]<0,:], variants_toanno.loc[variants_toanno[beta]>0,:]]:
                    last_pos = min(variants_toanno_half[maf])/2
                    for index, row in variants_toanno_half.iterrows():
                        
                        armB_length_in_point = ax.transData.transform((0,1.1*maxy))[1]-ax.transData.transform((0, abs(row[beta])))[1]
                        armB_length_in_point = armB_length_in_point*arm_scale

                        if anno_style == "right" :
                            #right style
                            if row[maf]>last_pos*(repel_force+1):
                                last_pos=row[maf]
                            else:
                                last_pos*= (repel_force+1)
                        elif anno_style == "expand" :
                            last_pos = row["ADJUSTED_i"]

                        if anno_style == "right"  or anno_style == "expand":
                            if row[beta] >0 :
                                texts_u.append(ax.annotate(row[anno], xy=(row[maf], row[beta]),xytext=(last_pos , 1.2*maxy),
                                                        arrowprops=dict(relpos=(0,0),arrowstyle="-|>",facecolor='black',connectionstyle="arc,angleA=-90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)),rotation=90,
                                                        ha="left",va="bottom",**anno_args))
                            else:
                                texts_d.append(ax.annotate(row[anno], xy=(row[maf], row[beta]),xytext=(last_pos , -1.2*maxy),
                                                        arrowprops=dict(relpos=(0,1),arrowstyle="-|>",facecolor='black',connectionstyle="arc,angleA=90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)),rotation=90,
                                                        ha="left",va="top",**anno_args))
                        
                        if anno_style=="tight":
                            texts_d.append(ax.text(row[maf], row[beta], row[anno]))
                if anno_style=="tight":
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


def plot_power( ns=1000,
                mode="q",
                ts=None,
                prevalences=0.1,
                or_to_rr=False,
                ncases=5000,
                ncontrols=5000, 
                sig_levels=5e-8,       
                maf_range=None,
                beta_range=None, 
                n_matrix=1000,
                xscale="log",
                yscale_factor=1,
                cmap="cool",
                ylim=None,
                fontsize=15,
                font_family="Arial",
                sizes=None,
                save=False,
                saveargs=None,
                ylabel="Effect size",
                xlabel="Minor allele frequency",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                verbose=True,
                log=Log()):
    
    #Checking columns#################################################################################################################
    matplotlib.rc('font', family=font_family)
    if sizes is None:
        sizes = (20,80)
    if ts is None:
        ts=[0.2,0.4,0.6,0.8]
    if xticks is None:
        if xscale== "log":
            xticks = [0.001,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks
        else:
            xticks = [0,0.01,0.05,0.1,0.2,0.5]
            xticklabels = xticks            

    #Checking columns#################################################################################################################
    if verbose: log.write("Start to create trumpet plot...")
    
    if mode=="b":
        if ncases is None or ncontrols is None:
            if verbose:
                log.write(" -No scase or scontrol. Skipping...")
            return None
        
    #configure beta and maf range ###################################################################################################
    if maf_range is None:
        maf_range=(np.power(10.0,-3),0.5)
    if beta_range is None:
        beta_range=(0.0001,3)
    
    #configure power threshold###################################################################################################

    if type(ns) is list:
        var_to_change = ns
        legend_title = "N"
    if type(ts) is list:
        var_to_change = ts
        legend_title = "Power"
    if type(ncases) is list:
        var_to_change = ncases
        legend_title = "Number of cases"    
    if type(ncontrols) is list:
        var_to_change = ncontrols
        legend_title = "Number of controls"  
    if type(sig_levels) is list:
        var_to_change = sig_levels
        legend_title = "Significance level" 
    if type(prevalences) is list:
        var_to_change = prevalences
        legend_title = "Prevalence" 
    
    #Print settings#############################################################################
    
    #configure colormap##########################################################################################################
    cmap_to_use = plt.cm.get_cmap(cmap)
    if cmap_to_use.N >100:
        max_value = max(var_to_change)
        min_value = min(var_to_change)
        norm = lambda x: x/max_value
        if legend_title == "Significance level":
            norm = lambda x: np.log10(x)/np.log10(min_value)
        rgba = cmap_to_use(list(map(norm, var_to_change)))
    else:
        rgba = cmap_to_use(range(len(var_to_change)))
    
    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))
    output_hex_colors

    ##################################################################################################
    fig, ax = plt.subplots(figsize=(10,10))
    
    ##creating power line############################################################################################
    if mode=="q":
        for i,value in enumerate(var_to_change):
            
            n = ns
            t = ts
            sig_level = sig_levels
            
            if legend_title == "N":
                n = value
            elif legend_title == "Power":
                t = value
            elif legend_title == "Significance level":
                sig_level = value
                    
            xpower = get_beta(mode="q",          
                            eaf_range=maf_range,
                            beta_range=beta_range, 
                            n=n,
                            t=t,
                            sig_level=sig_level,
                            n_matrix=n_matrix)
            xpower2 = xpower.copy()
            xpower2[1] = -xpower2[1] 
            xpower2[1] = xpower2[1] * yscale_factor
            xpower[1] = xpower[1] * yscale_factor
            lines = LineCollection([xpower2,xpower], label=value,color=output_hex_colors[i],zorder=0)
            ax.add_collection(lines)
    else:
        for i,value in enumerate(var_to_change):
            
            ncase = ncases
            ncontrol = ncontrols
            t = ts
            prevalence = prevalences
            sig_level = sig_levels
            
            if legend_title == "Prevalence":
                prevalence = value
            elif legend_title == "Power":
                t = value
            elif legend_title == "Significance level":
                sig_level = value
            elif legend_title == "Number of cases":
                ncase = value
            elif legend_title == "Number of controls":
                ncontrol = value

            xpower = get_beta_binary(        
                            eaf_range=maf_range,
                            beta_range=beta_range, 
                            prevalence=prevalence,
                            or_to_rr = or_to_rr,
                            ncase=ncase, 
                            ncontrol=ncontrol, 
                            t=t,
                            sig_level=sig_level,
                            n_matrix=n_matrix)
            xpower2 = xpower.copy()
            xpower2[1] = -xpower2[1] 
            xpower2[1] = xpower2[1] * yscale_factor
            xpower[1] = xpower[1] * yscale_factor
            lines = LineCollection([xpower2,xpower], label=value,color=output_hex_colors[i])
            ax.add_collection(lines)
    ###################################################################################################

    ax.tick_params(axis='y', labelsize=fontsize)
    leg = ax.legend(title=legend_title,fontsize =fontsize,title_fontsize=fontsize)

    for line in leg.get_lines():
        line.set_linewidth(5.0)
    ax.axhline(y=0,color="grey",linestyle="dashed")
    
    if xscale== "log":
        ax.set_xscale('log')
        rotation=0
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(maf_range[0]/2,0.52)
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
 
    if mode=="q":
        save_figure(fig, save, keyword="power_q",saveargs=saveargs, log=log, verbose=verbose)
    elif mode=="b":
        save_figure(fig, save, keyword="power_b",saveargs=saveargs, log=log, verbose=verbose)

    if verbose: log.write("Finished creating trumpet plot!")
    return fig


def plot_power_x(
                x=None, 
                ts=None,
                mode="q",
                ns=10000,
                mafs=0.1,
                betas=0.1,
                prevalences=0.1,
                or_to_rr=False,
                ncases=5000,
                ncontrols=5000, 
                sig_levels=5e-8,       
                maf_range=None,
                beta_range=None, 
                n_range=None,
                prevalence_range=None,
                n_matrix=1000,
                xscale="log",
                yscale_factor=1,
                cmap="cool",
                ylim=None,
                fontsize=15,
                font_family="Arial",
                save=False,
                saveargs=None,
                ylabel="Power",
                xlabel="N",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                verbose=True,
                log=Log()):
    
    #Checking columns#################################################################################################################
    if verbose: log.write("Start to create power plot...")
    matplotlib.rc('font', family=font_family)
    if verbose:
        log.write(" -Settings:")
        log.write("  -Mode: {}".format(mode))
        if mode == "q" :
            log.write("  -X axis: {}".format(x))
            if x!="N":
                log.write("  -N: {}".format(ns))
            if x!="MAF":
                log.write("  -MAF: {}".format(mafs))
        if mode == "b" :
            log.write("  -X axis: {}".format(x))
            if x!="N_CASE":
                log.write("  -N_CASE: {}".format(ncases))
            if x!="N_CASE":
                log.write("  -N_CONTROL: {}".format(ncontrols))
            if x!="PREVALENCE":
                log.write("  -PREVALENCE: {}".format(prevalences))
        if x!="BETA":
            log.write("  -BETA: {}".format(betas))

        log.write(" -Significance level: {}".format(sig_levels))

    if x is None:
        if mode=="b":
            x = "N_CASE"
        else:
            x = "N"

    if ylim is None:
        ylim=(0,1)

    if ts is None:
        ts = [0.8]
    #Checking columns#################################################################################################################
    
    
    if mode=="b":
        if ncases is None or ncontrols is None:
            if verbose:
                log.write(" -No scase or scontrol. Skipping...")
            return None

    #configure beta and maf range ###################################################################################################
    if n_range is None:
        n_range = np.linspace(1,50000,n_matrix)
    else:
        n_range = np.linspace(n_range[0],n_range[1],n_matrix)

    if beta_range is None:
        beta_range = np.linspace(0.0001,3,n_matrix)
    else:
        beta_range = np.linspace(beta_range[0],beta_range[1],n_matrix)
    
    if maf_range is None:
        if mode=="q":
            maf_range = np.linspace(0.0001,0.5,n_matrix)
        else:
            maf_range = np.linspace(0.0001,0.9999,n_matrix)
    else:
        maf_range = np.linspace(maf_range[0],maf_range[1],n_matrix)
    
    if prevalence_range is None:
        prevalence_range = np.linspace(0.01,0.99,n_matrix)
    else:
        prevalence_range = np.linspace(prevalence_range[0],prevalence_range[1],n_matrix)
    #configure power threshold###################################################################################################
    
    if type(ns) is list:
        var_to_change = ns
        legend_title = "N"
    if type(betas) is list:
        var_to_change = betas
        legend_title = "BETA"
    if type(mafs) is list:
        var_to_change = mafs
        legend_title = "MAF"
    if type(ncases) is list:
        var_to_change = ncases
        legend_title = "Number of cases"    
    if type(ncontrols) is list:
        var_to_change = ncontrols
        legend_title = "Number of controls"  
    if type(sig_levels) is list:
        var_to_change = sig_levels
        legend_title = "Significance level" 
    if type(prevalences) is list:
        var_to_change = prevalences
        legend_title = "Prevalence" 
    
    #Print settings#############################################################################
    
    #configure colormap##########################################################################################################
    cmap_to_use = plt.cm.get_cmap(cmap)
    if cmap_to_use.N >100:
        max_value = max(var_to_change)
        min_value = min(var_to_change)
        norm = lambda x: x/max_value
        if legend_title == "Significance level":
            norm = lambda x: np.log10(x)/np.log10(min_value)
        rgba = cmap_to_use(list(map(norm, var_to_change)))
    else:
        rgba = cmap_to_use(range(len(var_to_change)))
    
    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))
    output_hex_colors

    ##################################################################################################
    fig, ax = plt.subplots(figsize=(10,10))
    
    ##creating power line############################################################################################
    if mode=="q":
        for i,value in enumerate(var_to_change):
            
            n = ns
            beta = betas
            maf = mafs
            sig_level = sig_levels
            
            if legend_title == "BETA":
                beta = value
            elif legend_title == "MAF":
                maf = value
            elif legend_title == "N":
                n = value
            elif legend_title == "Significance level":
                sig_level = value


            if x == "N":
                x_values = n_range
                n = x_values
            elif x=="MAF":
                x_values = maf_range
                maf = x_values
            elif x=="BETA":
                x_values = beta_range
                beta = x_values

            xpower = get_power(mode="q",          
                              eaf=maf,
                              beta=beta, 
                              n=n,
                              sig_level=sig_level,
                              verbose=False)
            
            ax.plot(x_values, xpower, label=value,color=output_hex_colors[i],zorder=0)

    else:
        for i,value in enumerate(var_to_change):
            
            beta = betas
            maf = mafs
            ncase = ncases
            ncontrol = ncontrols
            prevalence = prevalences
            sig_level = sig_levels
            
            if legend_title == "BETA":
                beta = value
            elif legend_title == "MAF":
                maf = value
            elif legend_title == "Prevalence":
                prevalence = value
            elif legend_title == "Significance level":
                sig_level = value
            elif legend_title == "Number of cases":
                ncase = value
            elif legend_title == "Number of controls":
                ncontrol = value
            
            if x == "N_CASE":
                x_values = n_range
                ncase = x_values
            elif x== "N_CONTROL":
                x_values = n_range
                ncontrol = x_values
            elif x=="MAF":
                x_values = maf_range
                maf = x_values
            elif x=="BETA":
                x_values = beta_range
                beta = x_values
            elif x=="PREVALENCE":
                x_values = prevalence_range
                prevalence = x_values

            xpower = get_power(mode="b",          
                              beta = beta, 
                              daf = maf,
                              ncase=ncase,
                              ncontrol=ncontrol,
                              prevalence=prevalence,
                              sig_level=sig_level,
                              or_to_rr = or_to_rr,
                              verbose=False)
            
            ax.plot(x_values, xpower, label=value,color=output_hex_colors[i],zorder=0)
    ###################################################################################################

    ax.tick_params(axis='y', labelsize=fontsize)
    leg = ax.legend(title=legend_title,fontsize =fontsize,title_fontsize=fontsize)

    for line in leg.get_lines():
        line.set_linewidth(5.0)
    ax.axhline(y=0,color="grey",linestyle="dashed")
    
    for i in ts:
        ax.axhline(y=i,color="grey",linestyle="dashed")
    

    #if xscale== "log":
    #    ax.set_xscale('log')
    #    rotation=0
    #    ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
    #    ax.set_xlim(0,max(ns))
    #else:
    #    rotation=90    
    #    ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
    #    ax.set_xlim(0,max(ns))
    
    if ylim is not None:
        ax.set_ylim(ylim)
    if yticks is not None:
        ax.set_yticks(yticks, yticklabels)

    ax.set_ylabel(ylabel,fontsize=fontsize)
    
    if x == "N_CASE":
        xlabel = "Number of cases"
    elif x== "N_CONTROL":
        xlabel = "Number of controls"
    elif x== "N":
        xlabel = "Number of samples"
    elif x=="MAF":
        xlabel = "Minor allele frequency"
    elif x=="BETA":
        xlabel = "Risk allele effect size"
    elif x=="PREVALENCE":
        xlabel = "Prevalence"
    ax.set_xlabel(xlabel,fontsize=fontsize)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)

    ############  Annotation ##################################################################################################
 
    if mode=="q":
        save_figure(fig, save, keyword="power_xq",saveargs=saveargs, log=log, verbose=verbose)
    elif mode=="b":
        save_figure(fig, save, keyword="power_xb",saveargs=saveargs, log=log, verbose=verbose)

    if verbose: log.write("Finished creating power plot!")
    return fig
