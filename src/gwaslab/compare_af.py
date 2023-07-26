import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns
from gwaslab.Log import Log
from gwaslab.figuresave import save_figure

################################################################################################################################
def plotdaf(sumstats,
             eaf="EAF",
             daf="DAF",
             threshold=0.16,
             xlabel="Alternative Allele Frequency in Reference Population (RAF)",
             ylabel="Effect Allele Frequency in Sumstats (EAF)",
             is_reg=True,
             r2=True,
             is_45_helper_line=True,
             is_threshold=True,
             helper_line_args=None,
             threshold_line_args=None,
             reg_line_args=None,
             plt_args=None,
             scatter_args=None,
             scatter_args_outlier =None,
             histplot_args=None,
             font_args=None,
             r2_args=None,
             legend1=True,
             legend2=True,
             save=False,
             save_args=None,
             saveargs=None,
             verbose=True,
             log=Log()
           ):
    
    if font_args is None:
        font_args={'family':'sans','fontname':'Arial','fontsize':8}
    if scatter_args is None:
        scatter_args={"s":1}
    if scatter_args_outlier is None:
        scatter_args_outlier={"s":3,"c":"red"}
    if plt_args is None:
        plt_args={"figsize":(8,4),"dpi":300}
    if histplot_args is None:
        histplot_args={"log_scale":(False,True)}
    if reg_line_args is None:
        reg_line_args={"color":'#cccccc', "linestyle":'--'}
    if threshold_line_args is None:
        threshold_line_args={"color":'#cccccc', "linestyle":'dotted'}
    if helper_line_args is None:
        helper_line_args={"color":'black', "linestyle":'-',"lw":1}
    if r2_args is None:
        r2_args = {"va":"bottom","ha":"right"}
    if saveargs is None:
        if save_args is None:
            saveargs = save_args = {}
        else:
            saveargs = save_args

    if verbose: log.write("Start to plot Reference frequency vs Effect allele frequency plot...")
    if not ((eaf in sumstats.columns) and (daf in sumstats.columns)):
        raise ValueError("EAF and/or DAF columns were not detected.")
    
    if "SNPID" in sumstats.columns:
        snpid = "SNPID"
    else:
        snpid = "rsID"
    
    alleles =[]
    if "EA" in sumstats.columns:
        alleles.append("EA")
    if "NEA" in sumstats.columns:
        alleles.append("NEA")
    

    sumstats = sumstats.loc[(~sumstats[eaf].isna())&(~sumstats[daf].isna()),[snpid,eaf,daf]+alleles].copy()
    sumstats.loc[:,daf] = sumstats.loc[:,daf].astype("float")
    sumstats.loc[:,eaf] = sumstats.loc[:,eaf].astype("float")
    if verbose: log.write(" -Plotting valriants:" + str(len(sumstats)))
    
    sumstats.loc[:,"RAF"]=sumstats[eaf] - sumstats[daf]
    sns.set_style("ticks")
    fig, (ax1, ax2) = plt.subplots(1, 2,**plt_args)
    ax1.scatter(sumstats["RAF"],sumstats[eaf],label="Non-outlier", **scatter_args)
    
    if is_threshold is True:
        is_outliers = sumstats[daf].abs() > threshold 
        if sum(is_outliers)>0:
            ax1.scatter(sumstats.loc[is_outliers, "RAF"],sumstats.loc[is_outliers, eaf],label="Outlier", **scatter_args_outlier)
    
    if legend1 ==True:
        ax1.legend()
    
    if is_reg is True:
        if verbose: log.write(" -Plotting regression line...")
        reg = ss.linregress(sumstats["RAF"],sumstats[eaf])
        if verbose:log.write(" -Beta = ", reg[0])
        if verbose:log.write(" -Intercept = ", reg[1])
        if verbose:log.write(" -R2 = ", reg[2])
        ax1.axline(xy1=(0,reg[1]),slope=reg[0],zorder=1,**reg_line_args)
        if r2 is True:
            ax1.text(0.98,0.02, "$R^2 = {:.3f}$".format(reg[2]), transform=ax1.transAxes, **r2_args)

    if is_threshold is True:
        if verbose: log.write(" -Threshold : " + str(threshold))
        num = sum(np.abs(sumstats[daf])>threshold )
        if verbose: log.write(" -Variants with relatively large DAF : ",num )
        if verbose: log.write(" -Percentage for variants with relatively large DAF : ",num/len(sumstats) )
        ax1.axline(xy1=(0,threshold),slope=1,zorder=1,**threshold_line_args)
        ax1.axline(xy1=(threshold,0),slope=1,zorder=1,**threshold_line_args)
    
    xl,xh=ax1.get_xlim()
    yl,yh=ax1.get_ylim()
    
    if is_45_helper_line is True:
        ax1.axline([0,0], [1,1],zorder=1, **helper_line_args)
    
    ax1.set_xlabel(xlabel,**font_args)
    ax1.set_ylabel(ylabel,**font_args)
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    

    sumstats.loc[:,"ID"] = sumstats.index
    
    to_plot = pd.melt(sumstats,id_vars=['ID'], value_vars=['EAF',"RAF"], var_name='Types', value_name='Allele Frequency')
    
    sns.histplot(data=to_plot, x="Allele Frequency", hue="Types", fill=True, ax=ax2, legend = legend2 ,**histplot_args)
    ax2.set_xlabel("Allele Frequency",**font_args)


    plt.tight_layout()
    save_figure(fig, save, keyword="afc",saveargs=saveargs, log=log, verbose=verbose)

    #if save:
    #    if verbose: log.write("Saving plot:")
    #    if save==True:
    #        fig.savefig("./allele_frequency_comparison.png",bbox_inches="tight",**save_args)
    #        log.write(" -Saved to "+ "./allele_frequency_comparison.png" + " successfully!" )
    #    else:
    #        fig.savefig(save,bbox_inches="tight",**save_args)
    #        log.write(" -Saved to "+ save + " successfully!" )
    sumstats = sumstats.drop(columns="ID")
    return fig, sumstats[is_outliers].copy()
    
