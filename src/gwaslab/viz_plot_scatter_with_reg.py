import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as ss
import seaborn as sns
import gc
import math
import scipy.stats as ss
from matplotlib.patches import Rectangle
from adjustText import adjust_text
from gwaslab.viz_aux_save_figure import save_figure
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_get_sig import annogene
from gwaslab.g_Log import Log
from gwaslab.util_in_correct_winnerscurse import wc_correct
from gwaslab.util_in_correct_winnerscurse import wc_correct_test
from gwaslab.g_Sumstats import Sumstats
from gwaslab.io_process_args import _merge_and_sync_dic
from gwaslab.io_process_args import _extract_kwargs

def scatter(df, 
            x, 
            y,
            mode="0",
            reg_box=None,
            is_reg=True,
            fdr=False,
            allele_match=False,
            r_se=False,
            is_45_helper_line=False,
            plt_args=None,
            xylabel_prefix="Per-allele effect size in ",
            helper_line_args=None,
            font_args=None,
            fontargs=None,
            build="19",
            r_or_r2="r",
            err_kwargs=None,
            legend_args=None,
            log = Log(),
            save=False,
            reg_xmin=None,
            verbose=True,
            save_args=None,
            scatter_kwargs=None,
            font_kwargs=None,
            plt_kwargs=None,
            null_beta=0,
            engine="plt",
            **kwargs):
    
    if save_args is None:
        save_args = {"dpi":300,"facecolor":"white"}
    if reg_box is None:
        reg_box = dict(boxstyle='round', facecolor='white', alpha=1,edgecolor="None")
    if err_kwargs is None:
        err_kwargs={"ecolor":"#cccccc","elinewidth":1}
    if font_kwargs is None:
        font_kwargs={'fontsize':12,'family':'sans','fontname':'Arial'}
    if helper_line_args is None:
        helper_line_args={"color":'black', "linestyle":'-',"lw":1}
    if plt_kwargs is None:
        plt_kwargs={"figsize":(8,8),"dpi":300}
    if scatter_kwargs is None:
        scatter_kwargs={"s":20}
    if reg_xmin is None:
        reg_xmin = df[x].min()
    
    save_kwargs =      _extract_kwargs("save", save_args, locals())
    err_kwargs =       _extract_kwargs("err", err_kwargs, locals())
    plt_kwargs =       _extract_kwargs("plt", plt_kwargs,  locals())
    scatter_kwargs =   _extract_kwargs("scatter", scatter_kwargs, locals())
    font_kwargs =      _extract_kwargs("font",font_kwargs, locals())

    log.write("Start to create scatter plot...", verbose=verbose)
    fig,ax = plt.subplots(**plt_kwargs) 

   # plot x=0,y=0, and a 45 degree line
    xl,xh=ax.get_xlim()
    yl,yh=ax.get_ylim()

    #ax.axhline(y=0, zorder=1,**helper_line_args)
    #ax.axvline(x=0, zorder=1,**helper_line_args)
    
    #for spine in ['top', 'right']:
    #    ax.spines[spine].set_visible(False)
    
    log.write(" -Creating scatter plot : {} - {}...".format(x, y), verbose=verbose)
    if engine=="plt":
        ax.scatter(df[x],df[y],**scatter_kwargs)
    elif engine=="sns":
        sns.scatterplot(data=df,x=x,y=y,ax=ax,**scatter_kwargs)
    ###regression line##############################################################################################################################
    ax, reg = confire_regression_line(x, y,
                                 is_reg,
                                 reg_box, 
                                 df, 
                                 ax, 
                                 mode,
                                 xl,
                                 yl,
                                 xh,
                                 yh, 
                                 null_beta, 
                                 r_se, 
                                is_45_helper_line,
                                helper_line_args, 
                                font_kwargs,
                                log, 
                                verbose, reg_xmin)
    
    save_figure(fig = fig, save = save, keyword="scatter", save_args=save_args, log = log, verbose=verbose)

    return fig, ax, reg


def confire_regression_line(x, y, is_reg, reg_box, df,  ax, mode,xl,yl,xh,yh, null_beta, r_se, 
                            is_45_helper_line,helper_line_args, font_kwargs,
                            log, verbose, reg_xmin):
    # if N <3
    if len(df)<3: 
        is_reg=False
    
    if is_reg is True:
        # reg 
        # slope, intercept, r, p, slope_se, intercept_se
        if mode=="0":
            reg = ss.linregress(df[x],df[y])
            # estimate se for r
            if r_se==True:
                log.write(" -Estimating SE for rsq using Jackknife method.", verbose=verbose)
                r_se_jackknife = jackknife_r(df,x,y,log,verbose)
                r_se_jackknife_string = " ({:.2f})".format(r_se_jackknife)
            else:
                r_se_jackknife_string= ""
        else:
            reg = ss.linregress(df[x],df[y])
            r_se_jackknife_string= ""

        #### calculate p values based on selected value , default = 0 
        create_reg_log(reg, log, verbose)

        reg_string = create_reg_string(reg, 
                    r_se_jackknife_string)
        
        ax.text(0.99,0.01, reg_string, va="bottom",ha="right",transform=ax.transAxes,bbox=reg_box,**font_kwargs)
        
        ax = create_helper_line(ax, reg[0], is_45_helper_line, helper_line_args, reg_xmin=reg_xmin)
        ax = create_reg_line(ax, reg, reg_xmin=reg_xmin)

    return ax, reg

#############################################################################################################################################################################
def create_reg_log(reg,log, verbose):
    #t_score = (reg[0]-null_beta) / reg[4]
    #degree = len(df.dropna())-2
    p =  reg[3]
    #ss.t.sf(abs(t_score), df=degree)*2
    log.write(" -Beta = ", reg[0], verbose=verbose)
    log.write(" -Beta_se = ", reg[4], verbose=verbose)
    log.write(" -H0 beta =  0",", default p = ", "{:.2e}".format(reg[3]), verbose=verbose)
    log.write(" -Peason correlation coefficient =  ", "{:.2f}".format(reg[2]), verbose=verbose)
    log.write(" -r2 =  ", "{:.2f}".format(reg[2]**2), verbose=verbose)

def create_helper_line(ax, 
                       slope,
                       is_45_helper_line, 
                       helper_line_args,
                       reg_xmin=0):
    
    if is_45_helper_line is True:
        xl,xh=ax.get_xlim()
        yl,yh=ax.get_ylim()
        if slope >0:
            ax.axline([min(xl,yl),min(xl,yl)], [max(xh, yh),max(xh, yh)],zorder=1,**helper_line_args)
        else:
            ax.axline([min(xl,yl),-min(xl,yl)], [max(xh, yh),-max(xh, yh)],zorder=1,**helper_line_args)

    return ax

def create_reg_line(ax, reg, reg_xmin=0):
    xy1 = (reg_xmin,reg[0]*reg_xmin+reg[1])
    ax.axline(xy1=xy1,slope=reg[0],color="#cccccc",linestyle='--',zorder=1)
    return ax

def create_reg_string(reg,  
                      r_se_jackknife_string):
    p = reg[2]
    try:
        p12=str("{:.2e}".format(p)).split("e")[0]
        pe =str(int("{:.2e}".format(p).split("e")[1]))
    except:
        p12="0"
        pe="0"

    p_text="$p = " + p12 + " \\times  10^{"+pe+"}$"
    p_latex= f'{p_text}'

    reg_string = "$y =$ "+"{:.2f}".format(reg[1]) +" $+$ "+ "{:.2f}".format(reg[0])+" $x$, "+ p_latex + ", $r =$" +"{:.2f}".format(reg[2])+r_se_jackknife_string
    
    return reg_string

def jackknife_r(df,x,y,log,verbose):
    """Jackknife estimation of se for rsq
    """

    # dropna
    df_nona = df.loc[:,[x,y]].dropna()
    # non-empty entries
    n=len(df)
    # assign row number
    df_nona["_NROW"] = range(n)
    # a list to store r2
    r_list=[]
    # estimate r
    for i in range(n):
        # exclude 1 record
        records_to_use = df_nona["_NROW"]!=i
        # estimate r
        reg_jackknife = ss.linregress(df_nona.loc[records_to_use, x],df_nona.loc[records_to_use,y])
        # add r_i to list
        r_list.append(reg_jackknife[2])

    # convert list to array
    rs = np.array(r_list)
    # https://en.wikipedia.org/wiki/Jackknife_resampling
    r_se = np.sqrt( (n-1)/n * np.sum((rs - np.mean(rs))**2) )
    log.write(" -R se (jackknife) = {:.2e}".format(r_se), verbose=verbose)
    return r_se