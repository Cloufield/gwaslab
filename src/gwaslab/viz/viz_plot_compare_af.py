import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

################################################################################################################################
def plotdaf(sumstats,
             eaf="EAF",
             daf="DAF",
             raf="RAF",
             threshold=0.16,
             xlabel="Alternative Allele Frequency in Reference Population (RAF)",
             ylabel="Effect Allele Frequency in Sumstats (EAF)",
             is_reg=True,
             r2=True,
             is_45_helper_line=True,
             is_threshold=True,
             helper_line_kwargs=None,
             threshold_line_kwargs=None,
             reg_line_kwargs=None,
            fig_kwargs=None,
             scatter_kwargs=None,
             scatter_kwargs_outlier =None,
             histplot_kwargs=None,
             font_kwargs=None,
             r2_kwargs=None,
             legend1=True,
             legend2=True,
             save=False,
             save_kwargs=None,
             verbose=True,
             log=Log()
           ):
    
    # Extract dataframe if Sumstats object is passed
    if hasattr(sumstats, 'data') and not isinstance(sumstats, pd.DataFrame):
        sumstats = sumstats.data
    
    if font_kwargs is None:
        font_kwargs={'family':'sans','fontname':'Arial','fontsize':8}
    if scatter_kwargs is None:
        scatter_kwargs={"s":1}
    if scatter_kwargs_outlier is None:
        scatter_kwargs_outlier={"s":3,"c":"red"}
    if histplot_kwargs is None:
        histplot_kwargs={"log_scale":(False,False)}
    if reg_line_kwargs is None:
        reg_line_kwargs={"color":'#cccccc', "linestyle":'--'}
    if threshold_line_kwargs is None:
        threshold_line_kwargs={"color":'#cccccc', "linestyle":'dotted'}
    if helper_line_kwargs is None:
        helper_line_kwargs={"color":'black', "linestyle":'-',"lw":1}
    if r2_kwargs is None:
        r2_kwargs = {"va":"bottom","ha":"right"}
    style = set_plot_style(
        plot="plot_daf",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs,
        scatter_kwargs_outlier=scatter_kwargs_outlier,
        histplot_kwargs=histplot_kwargs,
        helper_line_kwargs=helper_line_kwargs,
        threshold_line_kwargs=threshold_line_kwargs,
        reg_line_kwargs=reg_line_kwargs,
        r2_kwargs=r2_kwargs,
        fontsize=font_kwargs.get('fontsize') if font_kwargs else None,
        fontfamily=font_kwargs.get('fontname') if font_kwargs else None,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    scatter_kwargs = style.get("scatter_kwargs", {})
    scatter_kwargs_outlier = style.get("scatter_kwargs_outlier", {})
    histplot_kwargs = style.get("histplot_kwargs", {})
    helper_line_kwargs = style["helper_line_kwargs"]
    threshold_line_kwargs = style["threshold_line_kwargs"]
    reg_line_kwargs = style["reg_line_kwargs"]
    r2_kwargs = style["r2_kwargs"]
    if font_kwargs is None:
        font_kwargs = {'family': 'sans', 'fontname': style["font_family"], 'fontsize': style["fontsize"]}

    log.write("Start to plot allele frequency comparison plot...", verbose=verbose)
    
    if not ((eaf in sumstats.columns) and ((daf in sumstats.columns)) or (raf in sumstats.columns)):
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
    
    if daf not in sumstats.columns:
        sumstats[daf] = sumstats[eaf] - sumstats[raf]

    sumstats = sumstats.loc[(~sumstats[eaf].isna())&(~sumstats[daf].isna()),[snpid,eaf,daf]+alleles].copy()
    sumstats[daf] = sumstats[daf].astype("float")
    sumstats[eaf] = sumstats[eaf].astype("float")
    log.write(" -Plotting valriants:" + str(len(sumstats)), verbose=verbose)
    if raf not in sumstats.columns:
        sumstats[raf] = sumstats[eaf] - sumstats[daf]
    sns.set_style("ticks")
    fig, [ax1, ax2] = plt.subplots(1, 2, **fig_kwargs)
    ax1.scatter(sumstats[raf],sumstats[eaf],label="Non-outlier", **scatter_kwargs)
    
    if is_threshold is True:
        is_outliers = sumstats[daf].abs() > threshold 
        if sum(is_outliers)>0:
            ax1.scatter(sumstats.loc[is_outliers, raf],sumstats.loc[is_outliers, eaf],label="Outlier", **scatter_kwargs_outlier)
    
    if legend1 ==True:
        ax1.legend()
    
    if is_reg is True:
        log.write(" -Plotting regression line...", verbose=verbose)
        reg = ss.linregress(sumstats[raf],sumstats[eaf])
        log.write(" -Beta = ", reg[0], verbose=verbose)
        log.write(" -Intercept = ", reg[1], verbose=verbose)
        log.write(" -R2 = ", reg[2], verbose=verbose)
        ax1.axline(xy1=(0,reg[1]),slope=reg[0],zorder=1,**reg_line_kwargs)
        if r2 is True:
            ax1.text(0.98,0.02, "$R^2 = {:.3f}$".format(reg[2]), transform=ax1.transAxes, **r2_kwargs)

    if is_threshold is True:
        log.write(" -Threshold : " + str(threshold), verbose=verbose)
        num = sum(np.abs(sumstats[daf])>threshold )
        log.write(" -Variants with relatively large DAF : ",num , verbose=verbose)
        log.write(" -Percentage for variants with relatively large DAF : ",num/len(sumstats) , verbose=verbose)
        ax1.axline(xy1=(0,threshold),slope=1,zorder=1,**threshold_line_kwargs)
        ax1.axline(xy1=(threshold,0),slope=1,zorder=1,**threshold_line_kwargs)
    
    xl,xh=ax1.get_xlim()
    yl,yh=ax1.get_ylim()
    
    if is_45_helper_line is True:
        ax1.axline([0,0], [1,1],zorder=1, **helper_line_kwargs)
    
    ax1.set_xlabel(xlabel,**font_kwargs)
    ax1.set_ylabel(ylabel,**font_kwargs)
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    

    sumstats["ID"] = sumstats.index
    
    to_plot = pd.melt(sumstats,id_vars=['ID'], value_vars=[eaf,raf], var_name='Types', value_name='Allele Frequency').dropna()

    sns.histplot(data=to_plot, x="Allele Frequency", 
                 hue="Types", fill=True, 
                 ax=ax2, legend = legend2,
                 **histplot_kwargs)
 
    ax2.set_xlabel("Allele Frequency",**font_kwargs)

    plt.tight_layout()
    save_figure(fig, save, keyword="afc", save_kwargs=save_kwargs, log=log, verbose=verbose)
    sumstats = sumstats.drop(columns="ID")

    return fig, sumstats[is_outliers].copy()
    
