from typing import TYPE_CHECKING, Optional, List, Tuple, Dict, Any, Union
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as sp
from math import ceil
from gwaslab.info.g_Log import Log

from gwaslab.viz.viz_aux_quickfix import _set_yticklabels
from gwaslab.util.util_in_calculate_gc import _lambda_GC

if TYPE_CHECKING:
    from matplotlib.axes import Axes

# qq plot module for mqqplot
from gwaslab.viz.viz_aux_save_figure import safefig
@safefig
def _plot_qq(
    sumstats: pd.DataFrame,
    p_toplot_raw: pd.DataFrame,
    ax2: 'Axes',
    maxticker: int,
    marker_size: Tuple[int, int],
    gc: bool,
    cut: Optional[float],
    cutfactor: Optional[float],
    cut_log: bool,
    skip: float,
    maxy: float,
    ystep: float,
    colors: List[str],
    qq_line_color: str,
    stratified: bool,
    eaf_raw: pd.Series,
    maf_bins: List[Tuple[float, float]],
    maf_bin_colors: List[str],
    fontsize: float,
    font_family: str,
    qtitle: Optional[str],
    qtitle_pad: float,
    title_fontsize: float,
    include_chrXYMT: bool,
    cut_line_color: str,
    linewidth: float,
    ytick3: Optional[List[float]],
    ylabels: Optional[List[str]],
    xlabels: Optional[List[float]],
    xlim: Optional[Tuple[float, float]],
    ylabels_converted: Optional[List[str]],
    qq_scatter_kwargs: Dict[str, Any],
    expected_min_mlog10p: float,
    verbose: bool = True,
    log: Log = Log()
) -> 'Axes':
    """
    Generate a QQ plot to visualize the distribution of p-values from GWAS results.
    
    Parameters
    ----------
    gc : bool
        Whether to calculate and display genomic inflation factor.
    stratified : bool, default=False
        Whether to stratify by MAF (Minor Allele Frequency).
    maf_bins : list of tuples, optional
        Bins for MAF stratification [[lower1, upper1], ...].
    maf_bin_colors : list, optional
        Colors for each MAF bin.
    qtitle : str, optional
        Title for the QQ plot.
    qtitle_pad : float, optional
        Padding for the title.
    include_chrXYMT : bool, optional
        Whether to include chrX, chrY, and MT in GC calculation.
    expected_min_mlog10p : float, optional
        Expected minimum -log10(p) value for theoretical distribution.
    fig_kwargs : dict, optional
        Figure arguments for subplots. Default is None. 
    colors : list, default=['#597FBD','#74BAD3']
        Color palette for plot. 
    marker_size : tuple, default=(5,20)
        Size range for markers. marker_size[1] will be used for qq plot.
    use_rank : bool, default=False
        Whether to use rank for plotting.
    verbose : bool, default=True
        Whether to show progress. 
    build : str, optional
        Genomic build version. Default is None. 
    dpi : int, default=200
        Dots per inch for figure resolution. 
    save : str, optional
        File path to save plot. Default is None. 
    save_kwargs : dict, default={"dpi":600,"transparent":True}
        Arguments for saving the plot. 
    Returns
    -------
    matplotlib.axes.Axes
        Modified axes object with the QQ plot.
    """
            
    # QQ plot #########################################################################################################
    # ax2 qqplot
    log.write("Start to create QQ plot with "+str(len(sumstats))+" variants:",verbose=verbose )
    
    # plotting qq plots using processed data after cut and skip
    
    # select -log10 scaled p to plot
    p_toplot = sumstats["scaled_P"]
    
    # min p value for uniform distribution 
    minit=1/len(p_toplot)
    
    
    upper_bound_p = np.power(10.0, -expected_min_mlog10p)

    explicit = {"color","marker_size"}
    qq_scatter_kwargs = {k: v for k, v in qq_scatter_kwargs.items() if k not in explicit}
    if stratified is False:
        log.write(" -Plotting all variants...",verbose=verbose)
        # sort x,y for qq plot
        # high to low
        observed = p_toplot.sort_values(ascending=False)
        
        # uniform distribution using raw number -> -log10 -> observed number (omit variants with low -log10p)
        #expected = -np.log10(np.linspace(minit,1,len(p_toplot_raw)))[:len(observed)]
    
        expected_all = -np.log10(np.linspace(minit,upper_bound_p,len(p_toplot_raw)))[:len(observed)]

        log.write(" -Expected range of P: (0,{})".format(upper_bound_p),verbose=verbose)
        #p_toplot = sumstats["scaled_P"]
        ax2.scatter(expected_all,observed,s=marker_size[1],color=colors[0],**qq_scatter_kwargs)

    else:
        # stratified qq plot
        log.write(" -Plotting variants stratified by MAF...",verbose=verbose)
        observed = p_toplot.sort_values(ascending=False)
        expected_all = -np.log10(np.linspace(minit,upper_bound_p,len(p_toplot_raw)))[:len(observed)]

        for i,(lower, upper) in enumerate(maf_bins):
            # extract data for a maf_bin
            
            databin = sumstats.loc[(sumstats["MAF"]>lower) &( sumstats["MAF"]<=upper),["MAF","scaled_P"]]
            # raw data : variants with maf(eaf_raw) in maf_bin
            
            databin_raw = eaf_raw[(eaf_raw>lower) & (eaf_raw<=upper)]
            # sort x,y for qq plot
            # high to low
            
            observed = databin["scaled_P"].sort_values(ascending=False)
            # uniform distribution using raw number -> -log10 -> observed number (omit variants with low -log10p)
            
            expected = -np.log10(np.linspace(minit,upper_bound_p,max(len(databin_raw),len(databin))))[:len(observed)]

            label ="( "+str(lower)+","+str(upper) +" ]"
            ax2.scatter(expected,observed,s=marker_size[1],color=maf_bin_colors[i],label=label,**qq_scatter_kwargs)
            ax2_legend= ax2.legend(loc="best",fontsize=fontsize,markerscale=1.5,frameon=False, handletextpad=0.4, borderaxespad=0.4)
            plt.setp(ax2_legend.texts, family=font_family)

    qq_x = max(skip, expected_min_mlog10p)
    ax2.plot([qq_x,-np.log10(minit)],[qq_x,-np.log10(minit)],linestyle="--",color=qq_line_color)
    ax2.set_xlabel("Expected $\mathregular{-log_{10}(P)}$",fontsize=fontsize,family=font_family)
    ax2.set_ylabel("Observed $\mathregular{-log_{10}(P)}$",fontsize=fontsize,family=font_family)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(True)
    
    # calculate genomic inflation factor and add annotation
    if gc == True:
        # gc was calculated using raw data (before cut and skip)
        p_toplot_raw = p_toplot_raw.rename(columns={"scaled_P":"MLOG10P"})
        
        level = 0.5

        if expected_min_mlog10p!=0:
            level = 1 -  np.power(10.0,-np.nanmedian(expected_all))
            log.write(" -Level for calculating lambda GC : {}".format(1 - level),verbose=verbose)

        if not include_chrXYMT : log.write(" -Excluding chrX,Y, MT from calculation of lambda GC.",verbose=verbose)
        lambdagc = _lambda_GC(p_toplot_raw, 
                            mode="MLOG10P", 
                            level=level, 
                            include_chrXYMT=include_chrXYMT,
                            log=log,
                            verbose=verbose)
        
        # annotate lambda gc to qq plot
        ax2.text(0.10, 1.03,"$\mathregular{\\lambda_{GC}}$ = "+"{:.4f}".format(lambdagc),
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    fontsize=fontsize,family=font_family)
    
    ax2 = _set_yticklabels(cut=cut,
                        cutfactor=cutfactor,
                        cut_log=cut_log,
                        ax1=ax2,
                        skip=skip,
                        maxy=maxy,
                        maxticker=maxticker,
                        ystep=ystep,
                        cut_line_color=cut_line_color,
                        fontsize=fontsize,
                        sc_linewidth = linewidth,
                        font_family=font_family,
                        ylabels=ylabels,
                        ytick3=ytick3,
                        ylabels_converted=ylabels_converted,
                        log=log,
                        verbose=verbose
                        )
    if xlim is not None:
        ax2.set_xlim(xlim)
        
    if xlabels is not None:
        ax2.set_xticks(xlabels)
        ax2.set_xticklabels(xlabels,fontsize=fontsize,family=font_family)

    ax2.tick_params(axis='both', which='both', labelsize=fontsize,labelfontfamily=font_family)
    #
    if qtitle:
        pad=(ax2.transData.transform((skip, qtitle_pad*maxy))[1]-ax2.transData.transform((skip, maxy)))[1]
        ax2.set_title(qtitle,fontsize=title_fontsize,pad=pad,family=font_family)

    log.write("Finished creating QQ plot successfully!",verbose=verbose)
    
    # Creating QQ plot Finished #############################################################################################
    return ax2
