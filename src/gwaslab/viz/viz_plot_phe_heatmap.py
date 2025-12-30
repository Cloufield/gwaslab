from typing import TYPE_CHECKING, Optional, Dict, Any, Union, Tuple, List
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as sp
import copy
from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version

from gwaslab.viz.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz.viz_aux_quickfix import _get_largenumber
from gwaslab.viz.viz_aux_quickfix import _quick_fix_p_value
from gwaslab.viz.viz_aux_quickfix import _quick_fix_pos
from gwaslab.viz.viz_aux_quickfix import _quick_fix_chr
from gwaslab.viz.viz_aux_quickfix import _quick_fix_eaf
from gwaslab.viz.viz_aux_quickfix import _quick_fix_mlog10p
from gwaslab.viz.viz_aux_quickfix import _dropna_in_cols
from gwaslab.viz.viz_plot_mqqplot import _process_p_value
from gwaslab.viz.viz_plot_mqqplot import _configure_fig_save_kwargs
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes

def _gwheatmap(
    insumstats: Union['Sumstats', pd.DataFrame],
    chrom: str = "CHR",
    pos: str = "POS",
    ref_chrom: str = "REF_CHR",
    ref_pos: str = "REF_START",
    p: str = "P",
    scaled: bool = False,
    sizes: Tuple[int, int] = (10, 50),
    alpha: float = 0.5,
    mlog10p: str = "MLOG10P",
    snpid: str = "SNPID",
    eaf: Optional[str] = None,
    group: str = "CIS/TRANS",
    ea: str = "EA",
    nea: str = "NEA",
    colors: Optional[List[str]] = None,
    check: bool = True,
    chr_dict: Optional[Dict[str, int]] = None,
    xchrpad: float = 0,
    ychrpad: float = 0,
    use_rank: bool = False,
    xtick_chr_dict: Optional[Dict[int, str]] = None,
    ytick_chr_dict: Optional[Dict[int, str]] = None,
    fontsize: float = 10,
    add_b: bool = False,
    log: Log = Log(),
    fig_kwargs: Optional[Dict[str, Any]] = None,
    scatter_kwargs: Optional[Dict[str, Any]] = None,
    height_ratios: Optional[List[float]] = None,
    hspace: float = 0.1,
    font_family: str = "Arial",
    cis_windowsizekb: int = 100,
    verbose: bool = True,
    save: Union[bool, str] = True,
    save_kwargs: Optional[Dict[str, Any]] = None,
    grid_linewidth: float = 1,
    grid_linecolor: str = "grey",
    **mqq_kwargs: Any
) -> Tuple['Figure', Log]:
    # Extract dataframe if Sumstats object is passed
    if hasattr(insumstats, 'data') and not isinstance(insumstats, pd.DataFrame):
        insumstats = insumstats.data
    
    log.write("Start to create genome-wide scatter plot...{}:".format(_get_version()),verbose=verbose)
    if height_ratios is None:
        height_ratios = [1, 2]
    if xtick_chr_dict is None:         
        xtick_chr_dict = get_number_to_chr()
    if ytick_chr_dict is None:         
        ytick_chr_dict = get_number_to_chr()
    if chr_dict is None:          
        chr_dict = get_chr_to_number()
    if colors is None:
        colors=["#CB132D","#597FBD"]
    style = set_plot_style(
        plot="plot_gwheatmap",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        colors=colors,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    scatter_kwargs = style.get("scatter_kwargs", {})
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    colors = style["colors"]

    fig_kwargs, scatter_kwargs, qq_scatter_kwargs, save_kwargs = _configure_fig_save_kwargs(save=save, 
                                                                                    fig_kwargs = fig_kwargs, 
                                                                                    scatter_kwargs = scatter_kwargs, 
                                                                                    qq_scatter_kwargs = dict(), 
                                                                                    save_kwargs = save_kwargs)

    sumstats = insumstats.copy()
    
    # Data QC and format
    if check ==True:
        sumstats[pos] = _quick_fix_pos(sumstats[pos], verbose=verbose, log=log)
        sumstats[chrom] = _quick_fix_chr(sumstats[chrom], chr_dict=chr_dict, verbose=verbose, log=log)
        sumstats[ref_pos] = _quick_fix_pos(sumstats[ref_pos], verbose=verbose, log=log)
        sumstats[ref_chrom] = _quick_fix_chr(sumstats[ref_chrom], chr_dict=chr_dict, verbose=verbose, log=log)
        sumstats = _dropna_in_cols(sumstats, [pos, chrom, ref_pos, ref_chrom], log=log, verbose=verbose)
    
    # dropna
    sumstats = sumstats.sort_values(by=group)

    if scaled is True:
        sumstats["raw_P"] = pd.to_numeric(sumstats[mlog10p], errors='coerce')
    else:
        sumstats["raw_P"] = sumstats[p].astype("float64")

    sumstats =  _process_p_value(sumstats=sumstats, 
                                mode="m",
                                p=p, 
                                mlog10p=mlog10p, 
                                scaled=scaled, 
                                log=log, 
                                verbose=verbose )
    
    
    
    if add_b ==False:
        fig, ax1 = plt.subplots(**fig_kwargs)
    else:
        fig, (ax2, ax1) = plt.subplots( nrows=2 ,sharex=True, gridspec_kw={'height_ratios': height_ratios }, **fig_kwargs)
        plt.subplots_adjust(hspace=hspace)   

    ## assign i for variants
    sumstats, chrom_df_x = _quick_assign_i_with_rank(sumstats, 
                                                chrpad=xchrpad, 
                                                use_rank=use_rank, 
                                                chrom=chrom,
                                                pos=pos,
                                                verbose=verbose)
    chrom_df_b = chrom_df_x
    sumstats = sumstats.rename(columns={"i":"i_x"})
    add_x_unique = list(sumstats["_ADD"].unique()) 
    
    ## determine grouping methods for Y
    ## assign i for Y group
    sumstats, chrom_df_y = _quick_assign_i_with_rank(sumstats, 
                              chrpad=ychrpad, 
                              use_rank=use_rank, 
                              chrom=ref_chrom,
                              pos=ref_pos,
                              verbose=verbose)
    
    sumstats = sumstats.rename(columns={"i":"i_y"})
    add_y_unique = list(sumstats["_ADD"].unique()) 
    
    if add_b == True:
        sumstats["i"] = sumstats["i_x"]
        fig,log = _mqqplot(sumstats,
                        chrom=chrom,
                        pos=pos,
                        p=p,
                        mlog10p=mlog10p,
                        snpid=snpid,
                        scaled=scaled,
                        log=log, 
                        mode="b",
                        figax=(fig,ax2),
                        _chrom_df_for_i = chrom_df_b,
                        _invert=False,
                        _if_quick_qc=False,
                        **mqq_kwargs
                        )
    ## 
    #min_xy = min(min(sumstats["i_x"]),min(sumstats["i_y"]))
    #max_xy = max(max(sumstats["i_x"]),max(sumstats["i_y"]))
    
    ## determine color 

    ## determine dot size

    ## plot  
    legend = True
    style=None
    linewidth=0
    edgecolor="black"

    palette = sns.color_palette(colors,n_colors=sumstats[group].nunique()) 
    
    #for index,g in enumerate(sumstats[group].unique()):
    #    
    #    palette = sns.color_palette("dark:{}".format(colors[index]), as_cmap=True)
    #    
    #    plot = sns.scatterplot(data=sumstats.loc[sumstats[group]==g,:], x='i_x', y='i_y',
    #            hue="scaled_P",
    #            palette=palette,
    #            size="scaled_P",
    #            alpha=alpha,
    #            sizes=sizes,
    #            legend=legend,
    #            style=style,
    #            linewidth=linewidth,
    #            edgecolor = edgecolor,
    #            zorder=2,
    #            ax=ax1)   
    
    plot = sns.scatterplot(data=sumstats, x='i_x', y='i_y',
            hue=group,
            palette=palette,
            size="scaled_P",
            alpha=alpha,
            sizes=sizes,
            legend=legend,
            style=style,
            linewidth=linewidth,
            edgecolor = edgecolor,
            zorder=2,
            ax=ax1, **scatter_kwargs)   
    
    handles, labels = ax1.get_legend_handles_labels()
    new_labels = []
    ncol = len(labels)
    for i in labels:
        if i==group:
            new_labels.append("Group")
        elif i=="scaled_P":
            new_labels.append("$-log_{10}(P)$")
        else:
            new_labels.append(i)
    
    ax1.legend(labels = new_labels,  handles=handles, loc="lower center", bbox_to_anchor=(.45, -0.17), 
                    ncol=ncol, scatterpoints=2, title=None, frameon=False)

    ## add vertical line
    for i in add_x_unique:
        ax1.axvline(x = i+0.5, linewidth = grid_linewidth,color=grid_linecolor,zorder=1000 )
    for i in add_y_unique:
        ax1.axhline(y = i+0.5,  linewidth = grid_linewidth,color=grid_linecolor,zorder=1000 )

    
    ## add X tick label 
    ax1 = _process_xtick(ax1, chrom_df_x, xtick_chr_dict, fontsize, font_family, log=log,verbose=True)
    ## add Y tick label 
    ax1 = _process_ytick(ax1, chrom_df_y, ytick_chr_dict, fontsize, font_family, log=log,verbose=True)
    
    ## set x y lim
    ax1.set_ylim([0.5,sumstats["i_y"].max()+1])
    ax1.set_xlim([0.5,sumstats["i_x"].max()+1])

    ## set x y label

    xlabel = "pQTL position"
    ax1.set_xlabel(xlabel,fontsize=fontsize,family=font_family)
    ylabel = "location of the gene encoding the target protein"
    ax1.set_ylabel(ylabel,fontsize=fontsize,family=font_family)
    
    save_figure(fig = fig, save = save, keyword="gwheatmap",  save_kwargs=save_kwargs, log = log, verbose=verbose)

    return fig, log

################################################################################################################
def _process_xtick(
    ax1: 'Axes',
    chrom_df: pd.Series,
    xtick_chr_dict: Dict[int, str],
    fontsize: float,
    font_family: str,
    log: Log = Log(),
    verbose: bool = True
) -> 'Axes':
    log.write(" -Processing X ticks...",verbose=False)
    ax1.set_xticks(chrom_df.astype("float64"))
    ax1.set_xticklabels(chrom_df.index.astype("Int64").map(xtick_chr_dict),fontsize=fontsize,family=font_family)    
    return ax1

def _process_ytick(
    ax1: 'Axes',
    chrom_df: pd.Series,
    ytick_chr_dict: Dict[int, str],
    fontsize: float,
    font_family: str,
    log: Log = Log(),
    verbose: bool = True
) -> 'Axes':
    log.write(" -Processing Y ticks...",verbose=verbose)
    ax1.set_yticks(chrom_df.astype("float64"))
    ax1.set_yticklabels(chrom_df.index.astype("Int64").map(ytick_chr_dict),fontsize=fontsize,family=font_family)    
    return ax1
