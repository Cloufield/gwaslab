import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as sp
from math import ceil
from gwaslab.viz_aux_quickfix import _set_yticklabels
from gwaslab.g_Log import Log
from gwaslab.util_in_calculate_gc import lambdaGC
# qq plot module for mqqplot
def _plot_qq(
    sumstats,
    p_toplot_raw,
    ax2,
    maxticker,
    marker_size,
    gc,
    cut,
    cutfactor,
    cut_log,
    skip,
    maxy,
    ystep,
    colors,
    qq_line_color,
    stratified,
    eaf_raw,
    maf_bins,
    maf_bin_colors,
    fontsize,
    font_family,
    qtitle,
    title_fontsize,
    include_chrXYMT,
    cut_line_color,
    linewidth,
    ytick3,
    ylabels,
    ylabels_converted,
    qq_scatter_args,
    expected_min_mlog10p,
    verbose=True,
    log=Log()
):

    # QQ plot #########################################################################################################
    # ax2 qqplot
    log.write("Start to create QQ plot with "+str(len(sumstats))+" variants:",verbose=verbose )
    
    # plotting qq plots using processed data after cut and skip
    
    # select -log10 scaled p to plot
    p_toplot = sumstats["scaled_P"]
    
    # min p value for uniform distribution 
    minit=1/len(p_toplot)
    
    
    upper_bound_p = np.power(10.0, -expected_min_mlog10p)

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
        ax2.scatter(expected_all,observed,s=marker_size[1],color=colors[0],**qq_scatter_args)

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

            label ="("+str(lower)+","+str(upper) +"]"
            ax2.scatter(expected,observed,s=marker_size[1],color=maf_bin_colors[i],label=label,**qq_scatter_args)
            ax2_legend= ax2.legend(loc="best",fontsize=fontsize,markerscale=3,frameon=False)
            plt.setp(ax2_legend.texts, family=font_family)

    qq_x = max(skip, expected_min_mlog10p)
    ax2.plot([qq_x,-np.log10(minit)],[qq_x,-np.log10(minit)],linestyle="--",color=qq_line_color)
    ax2.set_xlabel("Expected $-log_{10}(P)$",fontsize=fontsize,family=font_family)
    ax2.set_ylabel("Observed $-log_{10}(P)$",fontsize=fontsize,family=font_family)
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
        lambdagc = lambdaGC(p_toplot_raw, 
                            mode="MLOG10P", 
                            level=level, 
                            include_chrXYMT=include_chrXYMT,
                            log=log,
                            verbose=verbose)
        
        # annotate lambda gc to qq plot
        ax2.text(0.10, 1.03,"$\\lambda_{GC}$ = "+"{:.4f}".format(lambdagc),
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

    #if cut == 0: 
    #    ax2.set_ylim(skip, ceil(maxy*1.2) )
    #    
    #if cut:
    #    qcutline = ax2.axhline(y=cut, linewidth = linewidth, linestyle="--",color=cut_line_color,zorder=1)
    #    step=2
    #    if ((maxticker-cut)/cutfactor + cut) > cut:
    #        if ystep == 0:
    #            if (cut - skip ) // step > 10:
    #                step = (cut - skip ) // 10
    #        else:
    #            step = ystep
#
    #        ax2.set_yticks([x for x in range(skip,cut-step,step)]+[cut]+[(maxticker-cut)/cutfactor + cut])
    #        ax2.set_yticklabels([x for x in range(skip,cut-step,step)]+[cut]+[maxticker],fontsize=fontsize,family=font_family)
    #    else:
    #        if ystep == 0:
    #            if (cut - skip ) // step > 10:
    #                step = (cut - skip ) // 10
    #        else:
    #            step = ystep
#
    #        ax2.set_yticks([x for x in range(skip,cut-step,step)]+[cut])
    #        ax2.set_yticklabels([x for x in range(skip,cut-step,step)]+[cut],fontsize=fontsize,family=font_family)
    #    ax2.set_ylim(bottom = skip)

    ax2.tick_params(axis='both', which='both', labelsize=fontsize)
    #
    if qtitle:
        ax2.set_title(qtitle,fontsize=title_fontsize,pad=10,family=font_family)

    log.write("Finished creating QQ plot successfully!",verbose=verbose)
    
    # Creating QQ plot Finished #############################################################################################
    return ax2


