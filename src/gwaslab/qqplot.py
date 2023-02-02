import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy as sp
from gwaslab.Log import Log
from gwaslab.calculate_gc import lambdaGC

# qq plot module for mqqplot
def _plot_qq(
    sumstats,
    p_toplot_raw,
    ax2,
    maxticker,
    gc,
    cut,
    cutfactor,
    skip,
    maxy,
    colors,
    sig_line_color,
    stratified,
    eaf_raw,
    maf_bins,
    maf_bin_colors,
    fontsize,
    qtitle,
    title_fontsize,
    include_chrXYMT,
    cut_line_color,
    verbose=True,
    log=Log()
):

    # QQ plot #########################################################################################################
    # ax2 qqplot
    if verbose:log.write("Start to create QQ plot with "+str(len(sumstats))+" variants:")
    
    # plotting qq plots using processed data after cut and skip
    
    # select -log10 scaled p to plot
    p_toplot = sumstats["scaled_P"]
    
    # min p value for uniform distribution 
    minit=1/len(p_toplot)
    
    if stratified is False:
        # sort x,y for qq plot
        # high to low
        observed = p_toplot.sort_values(ascending=False)
        
        # uniform distribution using raw number -> -log10 -> observed number (omit variants with low -log10p)
        expected = -np.log10(np.linspace(minit,1,len(p_toplot_raw)))[:len(observed)]
        
        #p_toplot = sumstats["scaled_P"]
        ax2.scatter(expected,observed,s=8,color=colors[0])
    else:
        # stratified qq plot
        for i,(lower, upper) in enumerate(maf_bins):
            # extract data for a maf_bin
            databin = sumstats.loc[(sumstats["MAF"]>lower) &( sumstats["MAF"]<=upper),["MAF","scaled_P"]]
            # raw data : varaints with maf(eaf_raw) in maf_bin
            databin_raw = eaf_raw[(eaf_raw>lower) & (eaf_raw<=upper)]
            # sort x,y for qq plot
            # high to low
            observed = databin["scaled_P"].sort_values(ascending=False)
            # uniform distribution using raw number -> -log10 -> observed number (omit variants with low -log10p)
            expected = -np.log10(np.linspace(minit,1,max(len(databin_raw),len(databin))))[:len(observed)]
            label ="("+str(lower)+","+str(upper) +"]"
            ax2.scatter(expected,observed,s=8,color=maf_bin_colors[i],label=label)
            ax2_legend= ax2.legend(loc="best",fontsize=fontsize,markerscale=3,frameon=False)
            plt.setp(ax2_legend.texts, family="sans-serif")
    
    ax2.plot([skip,-np.log10(minit)],[skip,-np.log10(minit)],linestyle="--",color=sig_line_color)
    ax2.set_xlabel("Expected $-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
    ax2.set_ylabel("Observed $-log_{10}(P)$",fontsize=fontsize,family="sans-serif")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(True)
    
    # calculate genomic inflation factor and add annotation
    if gc == True:
        # gc was calculated using raw data (before cut and skip)
        p_toplot_raw = p_toplot_raw.rename(columns={"scaled_P":"MLOG10P"})
        if verbose and not include_chrXYMT : log.write(" -Excluding chrX,Y, MT from calculation of lambda GC.")
        lambdagc = lambdaGC(p_toplot_raw, mode="MLOG10P", include_chrXYMT=include_chrXYMT,log=log,verbose=False)
        if verbose: log.write(" -Calculating lambda GC:",lambdagc)
        
        # annotate lambda gc to qq plot
        ax2.text(0.10, 1.03,"$\\lambda_{GC}$ = "+"{:.4f}".format(lambdagc),
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    fontsize=fontsize,family="sans-serif")
    
    # set set_yticks and set_yticklabels 
    if cut == 0: 
        ax2.set_ylim(skip,maxy*1.2)
    if cut:
        qcutline=ax2.axhline(y=cut, linewidth = 2,linestyle="--",color=cut_line_color,zorder=1)
        if ((maxticker-cut)/cutfactor + cut) > cut:
            ax2.set_yticks([x for x in range(skip,cut+1,2)]+[(maxticker-cut)/cutfactor + cut])
            ax2.set_yticklabels([x for x in range(skip,cut+1,2)]+[maxticker],fontsize=fontsize,family="sans-serif")
        else:
            ax2.set_yticks([x for x in range(skip,cut+1,2)])
            ax2.set_yticklabels([x for x in range(skip,cut+1,2)],fontsize=fontsize,family="sans-serif")
    
    ax2.tick_params(axis='both', which='both', labelsize=fontsize)
    #
    if qtitle:
        ax2.set_title(qtitle,fontsize=title_fontsize,pad=10,family="sans-serif")

    if verbose: log.write("Finished creating QQ plot successfully!")
    
    # Creating QQ plot Finished #############################################################################################
    return ax2