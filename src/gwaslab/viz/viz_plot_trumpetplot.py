import matplotlib
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
from adjustText import adjust_text
from matplotlib.collections import LineCollection
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_calculate_power import get_beta
from gwaslab.util.util_in_calculate_power import get_beta_binary
from gwaslab.util.util_in_calculate_power import get_power
from gwaslab.util.util_in_fill_data import _fill_data
from gwaslab.util.util_in_get_sig import _anno_gene
from gwaslab.viz.viz_aux_annotate_plot import annotate_single
from gwaslab.viz.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_plot_mqqplot import _process_highlight

def _plot_trumpet(mysumstats,
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
                hue=None,
                highlight = None,
                highlight_chrpos = False,
                highlight_color="#CB132D",
                highlight_windowkb = 500,
                highlight_anno_kwargs = None,
                highlight_lim = None,
                highlight_lim_mode = "absolute",
                pinpoint= None,
                pinpoint_color ="red",
                scatter_kwargs=None,
                fontsize=15,
                font_family="Arial",
                size= None,
                sizes=None,
                save=False,
                save_kwargs=None,
                fig_kwargs=None,
                build="99",
                anno_set=None,
                anno_alias=None,
                anno_d=None,
                anno_kwargs=None,
                anno_style="expand",
                anno_source = "ensembl",
                anno_max_iter=100,
                arm_scale=1,
                repel_force=0.01,
                ylabel="Effect size",
                xlabel="Minor allele frequency",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                sort="beta",
                verbose=True,
                title=None,
                title_fontsize=15,
                log=Log()):
    """
    Create trumpet plot for genetic association results.
    
    This function generates a trumpet plot to visualize effect sizes and minor allele frequencies
    with power curves and optional highlighting of significant variants.
    
    Parameters
    ----------
    mode : str, required
        Analysis mode: "q" (quantitative trait) or "b" (binary trait). Default is "q"
    n : str, optional
        Column name for sample size for quantitative trait analysis. Used in "q" mode. Default is "N"
    ts : list, optional
        Power thresholds to plot. Default is [0.2,0.4,0.6,0.8]
    anno : str or None, optional
        Column name for variant annotation. Default is None
    prevalence : float, optional
        Prevalence for binary trait analysis. Used in "b" mode. Default is None
    or_to_rr : bool, optional
        Convert odds ratio to risk ratio. Used in "b" mode. Default is False
    ncase : int, optional
        Number of cases for binary trait analysis. Used in "b" mode. Default is None
    ncontrol : int, optional
        Number of controls for binary trait analysis.Used in "b" mode. Default is None
    sig_level : float, optional
        Significance threshold. Default is 5e-8
    p_level : float, optional
        P-value threshold for variant inclusion. Default is 5e-8
    anno_y : float, optional
        Annotation y-axis threshold. Default is 1
    anno_x : float, optional
        Annotation x-axis threshold. Default is 0.01
    maf_range : tuple, optional
        MAF range for power calculation. If None, auto-detected. Default is None
    beta_range : tuple, optional
        Beta range for power calculation. If None, auto-detected. Default is None
    n_matrix : int, optional
        Power curve smoothness parameter. Default is 1000
    xscale : str, optional
        X-axis scale: "log" or "linear". Default is "log"
    yscale_factor : float, optional
        Effect size scaling factor. Default is 1
    cmap : str, optional
        Colormap for power curves. Default is "cool"
    ylim : tuple, optional
        Y-axis limits. For example, [-3, 3]. If None, auto-detected. Default is None.
    xlim : tuple, optional
        X-axis limits. Default is None
    markercolor : str, optional
        Color for scatter points. Default is "#597FBD"
    hue : str, optional
        Column name for color grouping. Only used when you want to color variants differently. Default is None
    highlight : list, optional
        Variants to highlight. Default is None
    pinpoint : list, optional
        Variants to pinpoint. Default is None
    pinpoint_color : str, optional
        Color for pinpointed variants. Default is "red"
    scatter_kwargs : dict, optional
        Additional seaborn scatter plot arguments. Excluding x, y, size, ax,sizes, size_norm, legend, edgecolor, alpha, zorder. Default is None
    fontsize : int, optional
        Font size. Default is 15
    font_family : str, optional
        Font family. Default is "Arial"
    size : str, optional
        Column name for point size. Default is None
    sizes : tuple, optional
        Size range for points. Default is None
    save : bool, optional
        Save the plot. Default is False
    save_kwargs : dict, optional
        Arguments for saving the plot. Default is None
    fig_kwargs : dict, optional
        Figure creation arguments passed to plt.subplots. Default is None
    build : str, optional
        Reference genome build. Default is "99"
    anno_set : list, optional
        Annotation set. Default is None
    anno_alias : dict, optional
        Annotation alias mapping. Default is None
    anno_d : dict, optional
        Annotation dictionary. Default is None
    anno_kwargs : dict, optional
        Annotation arguments. Default is None
    anno_style : str, optional
        Annotation style. Default is "expand"
    anno_source : str, optional
        Annotation source. Default is "ensembl"
    sort : str, optional
        Sorting method for annotations. Default is "beta"
    verbose : bool, optional
        Verbose output. Default is True
    
    Returns
    -------
    matplotlib.figure.Figure
        The generated trumpet plot figure
    
    Less used parameters
    -------
    ylabel : str, optional
        Y-axis label. Default is "Effect size"
    xlabel : str, optional
        X-axis label. Default is "Minor allele frequency"
    xticks : list, optional
        X-axis tick positions. Default is None
    xticklabels : list, optional
        X-axis tick labels. Default is None
    yticks : list, optional
        Y-axis tick positions. Default is None
    yticklabels : list, optional
        Y-axis tick labels. Default is None
    anno_max_iter : int, optional
        Maximum annotation iterations. Default is 100
    arm_scale : float, optional
        Annotation arm scaling. Default is 1
    repel_force : float, optional
        Repulsion force for annotations. Default is 0.01
    highlight_chrpos : bool, optional
        Use chromosome position for highlighting. Default is False
    highlight_windowkb : int, optional
        Window size for highlighting. Default is 500
    highlight_anno_kwargs : dict, optional
        Annotation arguments for highlights. Default is None
    highlight_lim : tuple, optional
        Highlight limits. Default is None
    highlight_lim_mode : str, optional
        Highlight limit mode. Default is "absolute"
    """
    
    # Extract dataframe if Sumstats object is passed
    if hasattr(mysumstats, 'data') and not isinstance(mysumstats, pd.DataFrame):
        mysumstats = mysumstats.data
    
    #Checking columns#################################################################################################################
    matplotlib.rc('font', family=font_family)
    if sizes is None:
        sizes = (20,80)
    if anno_kwargs is None:
        anno_kwargs={"fontsize":12,"fontstyle":"italic"}
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
    style = set_plot_style(
        plot="plot_trumpet",
        mode=mode,
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs if save_kwargs is not None else save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", style.get("fig_kwargs", {}))
    save_kwargs = style.get("save_kwargs", style.get("save_kwargs", {}))
    scatter_kwargs = style.get("scatter_kwargs", style.get("scatter_kwargs", {}))
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    if "figsize" not in fig_kwargs:
        fig_kwargs = {**{"figsize":(10,8)}, **fig_kwargs}
    if scatter_kwargs is None:
        scatter_kwargs ={}
    if hue is not None:
        scatter_kwargs["hue"]=hue
    if markercolor is not None:
        scatter_kwargs["color"]=markercolor
    if highlight is None:
        highlight = list()
    if pinpoint is None:
        pinpoint = list()  
    #Checking columns#################################################################################################################
    log.write("Start to create trumpet plot...", verbose=verbose)
    
    #parameter check##################################################################################################################
    if (beta not in mysumstats.columns) or (eaf not in mysumstats.columns):
        log.write(" -No EAF or BETA columns. Skipping...", verbose=verbose)
        return None
    if mode=="b":
        if ncase is None or ncontrol is None:
            log.write(" -No scase or scontrol. Skipping...", verbose=verbose)
            return None
        if prevalence is None:
                prevalence= ncase/(ncase + ncontrol)
                log.write(" -Prevalence is not given. Estimating based on scase and scontrol :{}...".format(prevalence), verbose=verbose)
    
    #print settings##################################################################################################################

    log.write(" -Settings:", verbose=verbose)
    log.write("  -Mode: {}".format(mode), verbose=verbose)
    if mode == "q" :
        log.write("  -N: {}".format(n), verbose=verbose)
    if mode == "b" :
        log.write("  -N_CASE: {}".format(ncase), verbose=verbose)
        log.write("  -N_CONTROL: {}".format(ncontrol), verbose=verbose)
        log.write("  -PREVALENCE: {}".format(prevalence), verbose=verbose)
    log.write("  -BETA: {}".format(beta), verbose=verbose)
    log.write("  -Significance level: {}".format(sig_level), verbose=verbose)
    log.write("  -Power thresholds: {}".format(ts), verbose=verbose)
    log.write("  -Power line smoothness: {}".format(n_matrix), verbose=verbose)
    
    #loading columns #################################################################################################################
    cols_to_use = [snpid, beta, eaf, n, p]
    
    if len(highlight)>0: 
        cols_to_use.append(pos)
        cols_to_use.append(chrom)
    
    if anno is not None:
        if anno != "GENENAME":
            if anno!=True:
                log.write(" -Loading column {} for annotation...".format(anno), verbose=verbose)
                if anno not in cols_to_use:
                    if anno!=False:
                        cols_to_use.append(anno)
        else:
            cols_to_use.append(pos) if pos not in cols_to_use else cols_to_use
            cols_to_use.append(chrom) if chrom not in cols_to_use else cols_to_use

    if size != "ABS_BETA":
        if size is not None and size not in cols_to_use:
            cols_to_use.append(size)
    if "hue" in scatter_kwargs.keys():
        if scatter_kwargs["hue"] not in cols_to_use:
            cols_to_use.append(scatter_kwargs["hue"]) 
    #filter by p #################################################################################################################
    if p in mysumstats.columns:
        sumstats = mysumstats.loc[mysumstats[p]< p_level,cols_to_use ].copy()
        log.write(" -Excluding variants with P values > {}".format(p_level), verbose=verbose)
    else:
        cols_to_use.remove(p)
        sumstats = mysumstats[[beta,eaf,n]].copy()
    log.write(" -Plotting {} variants...".format(len(sumstats)), verbose=verbose)
    
    #add maf column #################################################################################################################
    if maf not in sumstats.columns:
        sumstats = _fill_data(sumstats,to_fill=["MAF"],verbose=False)
        is_filpped = (sumstats["MAF"] < sumstats[eaf]) & (sumstats[eaf] > 0.5)& (sumstats["MAF"] < 0.5)
        log.write(" -Flipping {} variants...".format(sum(is_filpped)), verbose=verbose)
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
        else:
            n = n
        log.write(" -N for power calculation: {}".format(n), verbose=verbose)

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
    cmap_to_use = matplotlib.colormaps.get_cmap(cmap)
    if cmap_to_use.N >100:
        rgba = cmap_to_use(ts)
    else:
        rgba = cmap_to_use(range(len(ts)))
    
    output_hex_colors=[]
    for i in range(len(rgba)):
        output_hex_colors.append(mc.to_hex(rgba[i]))
    output_hex_colors

    if len(highlight)>0:
        sumstats["HUE"] = pd.NA
        sumstats["HUE"] = sumstats["HUE"].astype("Int64")
        sumstats = _process_highlight(sumstats=sumstats, 
                                                    highlight=highlight, 
                                                    highlight_chrpos=highlight_chrpos, 
                                                    highlight_windowkb=highlight_windowkb, 
                                                    highlight_lim=highlight_lim,
                                                    highlight_lim_mode=highlight_lim_mode,
                                                    snpid=snpid, 
                                                    chrom=chrom, 
                                                    pos=pos)
    ##################################################################################################

    fig, ax = plt.subplots(**fig_kwargs)
    
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
    if size is None:
        size_norm = None
    else:
        size_norm = (sumstats[size].min(), sumstats[size].max())
    ## if highlight  ##################################################################################################
    explicit = {"zorder","alpha","edgecolor","legend","size_norm","data","x","y","size","sizes"}
    scatter_kwargs = {k: v for k, v in scatter_kwargs.items() if k not in explicit}
    log.write(" -Creating scatter plot...", verbose=verbose)
    dots = sns.scatterplot(data=sumstats,
                    x=maf,
                    y=beta,
                    size=size, 
                    ax=ax, 
                    sizes=sizes,
                    size_norm=size_norm,
                    legend=True, 
                    edgecolor="black",
                    alpha=0.8,
                    zorder=2,
                    **scatter_kwargs)
    log.write(" -Finished screating scatter plot...", verbose=verbose)
    if len(highlight) >0:
        
        legend = None
        style=None
        linewidth=0
        edgecolor="black"

        if pd.api.types.is_list_like(highlight[0]) and highlight_chrpos==False:
            for i, highlight_set in enumerate(highlight):
                scatter_kwargs["color"]=highlight_color[i%len(highlight_color)]
                log.write(" -Highlighting set {} target loci...".format(i+1),verbose=verbose)
                sns.scatterplot(data=sumstats.loc[sumstats["HUE"]==i], x=maf,
                    y=beta,
                    legend=legend,
                    style=style,
                    size=size,
                    sizes=sizes,
                    size_norm=size_norm,
                    linewidth=linewidth,
                    zorder=3+i,
                    ax=ax,edgecolor=edgecolor,**scatter_kwargs)  

        else:
            log.write(" -Highlighting target loci...",verbose=verbose)
            scatter_kwargs["color"]=highlight_color
            sns.scatterplot(data=sumstats.loc[sumstats["HUE"]==0], x=maf,
                    y=beta,
                legend=legend,
                size=size,
                sizes=sizes,
                size_norm=size_norm,
                zorder=3,
                ax=ax,
                edgecolor="black",
                **scatter_kwargs)  
    ####################################################################################################################
    if len(pinpoint)>0:
        legend = None
        style=None
        linewidth=0
        edgecolor="black"
        if pd.api.types.is_list_like(pinpoint[0]):

            for i, pinpoint_set in enumerate(pinpoint):
                scatter_kwargs["color"]=pinpoint_color[i%len(pinpoint_color)]
                if sum(sumstats[snpid].isin(pinpoint_set))>0:
                    to_pinpoint = sumstats.loc[sumstats[snpid].isin(pinpoint_set),:]
                    log.write(" -Pinpointing set {} target variants...".format(i+1),verbose=verbose)
                    sns.scatterplot(data=to_pinpoint, 
                    x=maf,
                    y=beta,
                    legend=legend,
                    size=size,
                    sizes=sizes,
                    size_norm=size_norm,
                    zorder=3,
                    ax=ax,
                    edgecolor="black",
                    **scatter_kwargs)  
                    #ax.scatter(to_pinpoint[maf],to_pinpoint[beta],color=pinpoint_color[i%len(pinpoint_color)],zorder=100,s=to_pinpoint[size])
                else:
                    log.write(" -Target variants to pinpoint were not found. Skip pinpointing process...",verbose=verbose)
        else:
            scatter_kwargs["color"]=pinpoint_color
            if sum(sumstats[snpid].isin(pinpoint))>0:
                to_pinpoint = sumstats.loc[sumstats[snpid].isin(pinpoint),:]
                log.write(" -Pinpointing target variants...",verbose=verbose)
                sns.scatterplot(data=to_pinpoint, x=maf,
                    y=beta,
                    legend=legend,
                    size=size,
                    sizes=sizes,
                    size_norm=size_norm,
                    zorder=3,
                    ax=ax,
                    edgecolor="black",
                    **scatter_kwargs)  
                #ax.scatter(to_pinpoint[maf],to_pinpoint[beta],color=pinpoint_color[i%len(pinpoint_color)],zorder=100,s=to_pinpoint[size])
            else:
                log.write(" -Target variants to pinpoint were not found. Skip pinpointing process...",verbose=verbose)
    
    ####################################################################################################################
    
    #second_legend = ax.legend(title="Power", loc="upper right",fontsize =fontsize,title_fontsize=fontsize)
    log.write(" -Creating legends...")
    # curve, size,  hue
    h,l = ax.get_legend_handles_labels()
    
    if len(ts)>0:
        # power curves
        l1 = ax.legend(h[:int(len(ts))],l[:int(len(ts))], title="Power", loc="upper right",fontsize =fontsize,title_fontsize=fontsize)
        for line in l1.get_lines():
            line.set_linewidth(5.0)
    ## hue
    if hue is not None or size is not None:
    #    # sizes
        if hue is None and size is not None:
            hue_size_legend_title = size
        elif hue is not None and size is None:
            hue_size_legend_title = hue
        else:
            hue_size_legend_title = None
        l2 = ax.legend(h[int(len(ts)):],l[int(len(ts)):], title=hue_size_legend_title, loc="lower right",fontsize =fontsize,title_fontsize=fontsize)
    
    ## hue
    if len(ts)>0:
        ax.add_artist(l1)
    if hue is not None or size is not None:
        ax.add_artist(l2)
    #first_legend = ax.legend(handles=dots, loc="lower right" ,title=size,fontsize =fontsize,title_fontsize=fontsize)
    #ax.add_artist(first_legend)
    ##################################################################################################

    ax.tick_params(axis='y', labelsize=fontsize)

    ax.axhline(y=0,color="grey",linestyle="dashed")
    
    if xscale== "log":
        ax.set_xscale('log')
        rotation=0
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(sumstats[maf].min()/2,0.52)
    else:
        rotation=90    
        ax.set_xticks(xticks,xticklabels,fontsize=fontsize,rotation=rotation)
        ax.set_xlim(-0.01,0.52)
        
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
                variants_toanno = _anno_gene(variants_toanno,
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
                maxy = variants_toanno[beta].abs().max()
                #maxy = max(variants_toanno[beta].abs().max(),1.5)
                variants_toanno["ADJUSTED_i"] = np.nan 
                y_span = 0.5
                
                if sort == "beta" : 
                    variants_toanno = variants_toanno.sort_values(by=beta, key= np.abs, ascending = False)
                else:
                    variants_toanno = variants_toanno.sort_values(by=maf, key= np.abs, ascending = True)
                
                if anno_style == "expand":

                    min_factor=None
                    
                    if len(variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"])>1:
                        variants_toanno.loc[variants_toanno[beta]>0, "ADJUSTED_i"] = adjust_text_position(variants_toanno.loc[variants_toanno[beta]>0,maf].values.copy(), 
                                                                                y_span, 
                                                                                repel_force=repel_force,
                                                                                max_iter=anno_max_iter,
                                                                                log=log,
                                                                                amode=xscale,
                                                                                verbose=verbose,min_factor=min_factor)

                    if len(variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"])>1:
                        variants_toanno.loc[variants_toanno[beta]<0, "ADJUSTED_i"] = adjust_text_position(variants_toanno.loc[variants_toanno[beta]<0,maf].values.copy(), 
                                                                y_span, 
                                                                repel_force=repel_force,
                                                                max_iter=anno_max_iter,
                                                                log=log,
                                                                amode=xscale,
                                                                verbose=verbose,min_factor=min_factor)

                
                for variants_toanno_half in [variants_toanno.loc[variants_toanno[beta]<0,:], variants_toanno.loc[variants_toanno[beta]>0,:]]:
                    if len(variants_toanno_half)<1:
                        continue
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
                                                        ha="left",va="bottom",**anno_kwargs))
                            else:
                                texts_d.append(ax.annotate(row[anno], xy=(row[maf], row[beta]),xytext=(last_pos , -1.2*maxy),
                                                        arrowprops=dict(relpos=(0,1),arrowstyle="-|>",facecolor='black',connectionstyle="arc,angleA=90,armA={},angleB=0,armB=0,rad=0".format(armB_length_in_point)),rotation=90,
                                                        ha="left",va="top",**anno_kwargs))
                        
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
    if title:
        ax.set_title(title,fontsize=title_fontsize,family=font_family)
    
    if mode=="q":
        save_figure(fig, save, keyword="trumpet_q",save_kwargs=save_kwargs, log=log, verbose=verbose)
    elif mode=="b":
        save_figure(fig, save, keyword="trumpet_b",save_kwargs=save_kwargs, log=log, verbose=verbose)

    log.write("Finished creating trumpet plot!", verbose=verbose)
    return fig

####################################################################


def plot_power(ns=1000,
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
                save_kwargs=None,
                ylabel="Effect size",
                xlabel="Minor allele frequency",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                verbose=True,
                log=Log()):
    """
    Generate power curves for genetic association studies.
    
    This function creates power curves visualizing the relationship between
    effect sizes and minor allele frequencies under various parameter settings.
    
    Parameters
    ----------
    ns : int or list, optional
        Sample size(s) for quantitative trait analysis. Default is 1000
    mode : str, optional
        Analysis mode: "q" (quantitative) or "b" (binary). Default is "q"
    ts : list or None, optional
        Power thresholds to plot. Default is [0.2, 0.4, 0.6, 0.8]
    prevalences : float or list, optional
        Prevalence for binary trait analysis. Default is 0.1
    or_to_rr : bool, optional
        Convert odds ratio to risk ratio. Default is False
    ncases : int or list, optional
        Number of cases for binary trait analysis. Default is 5000
    ncontrols : int or list, optional
        Number of controls for binary trait analysis. Default is 5000
    sig_levels : float or list, optional
        Significance threshold(s). Default is 5e-8
    maf_range : tuple or None, optional
        MAF range for power calculation. Default is (0.001, 0.5)
    beta_range : tuple or None, optional
        Beta range for power calculation. Default is (0.0001, 3)
    n_matrix : int, optional
        Power curve smoothness parameter. Default is 1000
    xscale : str, optional
        X-axis scale: "log" or "linear". Default is "log"
    yscale_factor : float, optional
        Effect size scaling factor. Default is 1
    cmap : str, optional
        Colormap for power curves. Default is "cool"
    ylim : tuple or None, optional
        Y-axis limits. Default is None
    fontsize : int, optional
        Font size for labels. Default is 15
    font_family : str, optional
        Font family for text. Default is "Arial"
    sizes : tuple or None, optional
        Size range for points. Default is (20, 80)
    save : bool, optional
        Save the plot. Default is False
    save_kwargs : dict or None, optional
        Arguments for saving the plot. Default is None
    ylabel : str, optional
        Y-axis label. Default is "Effect size"
    xlabel : str, optional
        X-axis label. Default is "Minor allele frequency"
    xticks : list or None, optional
        X-axis tick positions. Default is None
    xticklabels : list or None, optional
        X-axis tick labels. Default is None
    yticks : list or None, optional
        Y-axis tick positions. Default is None
    yticklabels : list or None, optional
        Y-axis tick labels. Default is None
    verbose : bool, optional
        Verbose output. Default is True
    log : gwaslab.g_Log.Log, optional
        Logger object. Default is new Log()
    
    Returns
    -------
    matplotlib.figure.Figure
        The generated power curve plot figure
    """
    
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
    log.write("Start to create trumpet plot...", verbose=verbose)
    
    if mode=="b":
        if ncases is None or ncontrols is None:
            log.write(" -No scase or scontrol. Skipping...", verbose=verbose)
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
    cmap_to_use = matplotlib.colormaps.get_cmap(cmap)
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
        ax.set_xlim(-0.01,0.52)
    
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
        save_figure(fig, save, keyword="power_q",save_kwargs=save_kwargs, log=log, verbose=verbose)
    elif mode=="b":
        save_figure(fig, save, keyword="power_b",save_kwargs=save_kwargs, log=log, verbose=verbose)

    log.write("Finished creating trumpet plot!", verbose=verbose)
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
                save_kwargs=None,
                ylabel="Power",
                xlabel="N",
                xticks = None,
                xticklabels = None,
                yticks = None,
                yticklabels=None,
                verbose=True,
                log=Log()):
    
    #Checking columns#################################################################################################################
    """
    Create power cure with user-specified x axis
    q mode:
        N
        MAF
    b mode:
        N_CASE
        N_CONTROL
        PREVALENCE
    BETA
    """

    log.write("Start to create power plot...", verbose=verbose)
    matplotlib.rc('font', family=font_family)

    log.write(" -Settings:", verbose=verbose)
    log.write("  -Mode: {}".format(mode), verbose=verbose)
    
    if mode == "q" :
        log.write("  -X axis: {}".format(x), verbose=verbose)
        if x!="N":
            log.write("  -N: {}".format(ns), verbose=verbose)
        if x!="MAF":
            log.write("  -MAF: {}".format(mafs), verbose=verbose)

    if mode == "b" :
        log.write("  -X axis: {}".format(x), verbose=verbose)
        if x!="N_CASE":
            log.write("  -N_CASE: {}".format(ncases), verbose=verbose)
        if x!="N_CASE":
            log.write("  -N_CONTROL: {}".format(ncontrols), verbose=verbose)
        if x!="PREVALENCE":
            log.write("  -PREVALENCE: {}".format(prevalences), verbose=verbose)
    if x!="BETA":
        log.write("  -BETA: {}".format(betas), verbose=verbose)

    log.write(" -Significance level: {}".format(sig_levels), verbose=verbose)

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
            log.write(" -No scase or scontrol. Skipping...", verbose=verbose)
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
    # for hue
    # which is the variable
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
    cmap_to_use = matplotlib.colormaps.get_cmap(cmap)
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
            # iterate through variables
            n = ns
            beta = betas
            maf = mafs
            sig_level = sig_levels
            
            # update the variable
            if legend_title == "BETA":
                beta = value
            elif legend_title == "MAF":
                maf = value
            elif legend_title == "N":
                n = value
            elif legend_title == "Significance level":
                sig_level = value

            # update X
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
        save_figure(fig, save, keyword="power_xq",save_kwargs=save_kwargs, log=log, verbose=verbose)
    elif mode=="b":
        save_figure(fig, save, keyword="power_xb",save_kwargs=save_kwargs, log=log, verbose=verbose)

    log.write("Finished creating power plot!", verbose=verbose)
    return fig
