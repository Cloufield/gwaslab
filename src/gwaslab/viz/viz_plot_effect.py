import pandas as pd
import matplotlib.pyplot as plt
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.info.g_Log import Log
from gwaslab.io.io_process_kwargs import _extract_kwargs
import seaborn as sns
from gwaslab.viz.viz_aux_style_options import set_plot_style

def _plot_effect(to_plot, 
                 y=None, 
                 y_sort=None, 
                 group=None, 
                 x="BETA", 
                 se="SE", 
                 eaf="EAF", 
                 snpr2="SNPR2", 
                 ylabel="Variant",
                 eaf_panel=True, 
                 snpvar_panel=True, 
                 rename_dic=None, 
                 err_kwargs=None,
                 font_kwargs=None,
                 save=None,
                 title=None,
                 save_kwargs=None,
                 eaf_kwargs=None,
                 snpr2_kwargs=None,
                 fig_kwargs=None,
                 scatter_kwargs=None,
                 effect_label=None,
                 eaf_label=None,
                 snpr2_label=None,
                 xlim_eaf=None,
                 xlim_snpr2 = None,
                 log=Log(),
                 verbose=True,
                 legend_mode=1,
                 ncol=2,
                 gap=1,
                 fontsize=12,
                 font_family="Arial",
                 size=None,
                 hue=None,
                 style_col=None,
                 sort_kwargs=None,
                 **args):
    """Plot effect sizes with optional panels for EAF and SNPR2.

    Parameters
    ----------
    to_plot : pandas.DataFrame
        DataFrame containing the data to plot.
    y : list of str, optional
        Columns to use for y-axis labeling. Default is None.
    y_sort : list of str, optional
        Columns to sort by for y-axis order. Default is None.
    group : list of str, optional
        Columns to group by. Default is None.
    x : str, default="BETA"
        Column name for x-axis (effect size).
    se : str, default="SE"
        Column name for standard error.
    eaf : str, default="EAF"
        Column name for effect allele frequency.
    snpr2 : str, default="SNPR2"
        Column name for SNP R-squared.
    ylabel : str, default="Variant"
        Label for y-axis.
    eaf_panel : bool, default=True
        Whether to include EAF panel.
    snpvar_panel : bool, default=True
        Whether to include SNP variance panel.
    rename_dic : dict, optional
        Dictionary for renaming labels. Default is None.
    err_kwargs : dict, optional
        Arguments for error bars. Default is None.
    font_kwargs : dict, optional
        Font-related arguments. Default is None.
    save : str, optional
        File path to save the figure. Default is None.
    title : str, optional
        Title of the plot. Default is None.
    save_kwargs : dict, optional
        Additional arguments for saving. Default is None.
    eaf_kwargs : dict, optional
        Arguments for EAF panel. Default is None.
    snpr2_kwargs : dict, optional
        Arguments for SNPR2 panel. Default is None.
    fig_kwargs : dict, optional
        Figure size and DPI settings. Default is None.
    scatter_kwargs : dict, optional
        Scatter plot arguments. Default is None.
    effect_label : str, optional
        Custom label for effect size axis. Default is None.
    eaf_label : str, optional
        Custom label for EAF axis. Default is None.
    snpr2_label : str, optional
        Custom label for SNPR2 axis. Default is None.
    xlim_eaf : tuple, optional
        X-axis limits for EAF panel. Default is None.
    xlim_snpr2 : tuple, optional
        X-axis limits for SNPR2 panel. Default is None.
    log : gwaslab.g_Log.Log, default=Log()
        Logging object.
    verbose : bool, default=True
        Whether to print progress messages.
    legend_mode : int, default=1
        Legend positioning mode.
    ncol : int, default=2
        Number of legend columns.
    gap : int, default=1
        Gap between groups.
    fontsize : int, default=12
        Font size for labels.
    font_family : str, default="Arial"
        Font family for text.
    size : str, optional
        Column name for marker size. Default is None.
    hue : str, optional
        Column name for color encoding. Default is None.
    style_col : str, optional
        Column name for marker style. Default is None.
    sort_kwargs : dict, optional
        Additional sorting arguments. Default is None.
    **args : dict
        Additional arguments passed to scatterplot.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure object.
    """
    # Extract dataframe if Sumstats object is passed
    if hasattr(to_plot, 'data') and not isinstance(to_plot, pd.DataFrame):
        to_plot = to_plot.data
    
    style = set_plot_style(
        plot="plot_effect",
        fig_kwargs=fig_kwargs if fig_kwargs is not None else fig_kwargs,
        save_kwargs=save_kwargs if save_kwargs is not None else save_kwargs,
        save=save,
        scatter_kwargs=scatter_kwargs if scatter_kwargs is not None else scatter_kwargs,
        err_kwargs=err_kwargs,
        eaf_kwargs=eaf_kwargs,
        snpr2_kwargs=snpr2_kwargs,
        font_kwargs=font_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", style.get("fig_kwargs", {}))
    save_kwargs = style.get("save_kwargs", style.get("save_kwargs", {}))
    scatter_kwargs = style.get("scatter_kwargs", style.get("scatter_kwargs", {}))
    err_kwargs = style["err_kwargs"]
    eaf_kwargs = style["eaf_kwargs"]
    snpr2_kwargs = style["snpr2_kwargs"]
    font_kwargs = style["font_kwargs"]
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    if sort_kwargs is None:
        sort_kwargs={}

    legend_titles=[]
    if hue is not None:
        args["hue"] = hue
        legend_titles.append(hue)

    if size is not None:
        args["size"] = size
        legend_titles.append(size)

    if style_col is not None:
        args["style_col"] = style_col
        legend_titles.append(style_col)

        
    save_kwargs =      _extract_kwargs("save", save_kwargs, locals())
    err_kwargs =       _extract_kwargs("err", err_kwargs, locals())
    scatter_kwargs =   _extract_kwargs("scatter", scatter_kwargs, locals())
    font_kwargs =      _extract_kwargs("font", font_kwargs, locals())
    
    def concat_cols(cols):
        string = "-".join(map(str,cols))
        return string
    
    if len(to_plot)>100:
        return "Too many variants to plot"
    
    if y is list:
        y_name = "-".join(y)
        to_plot[y_name] = to_plot[y].apply(lambda x: concat_cols(x), axis=1)
    else:
        y_name = y
    

    # sort y 
    if y_sort is None:
        y_sort = [i  for i in ["CHR","POS","STUDY"] if i in to_plot.columns]
    
    #to_plot = to_plot.sort_values(by=y_sort)

    if group is None:
        group = ["CHR","POS"] + y_sort
    
    sort_columns= group + y_sort

    to_plot = to_plot.sort_values(by=sort_columns,**sort_kwargs)
    
    # calculate cum sum
    cum_sizes = to_plot.groupby(group).size()
    cum_sizes = cum_sizes +  gap 
    cum_sizes = cum_sizes.cumsum() 
    
    # create index for y axis
    to_plot['_GROUP_CUMSUM'] = to_plot.set_index(group).index.map(cum_sizes)
    to_plot['_VAR_INDEX'] = to_plot.groupby(group).cumcount()
    to_plot["_VAR_INDEX"]=  to_plot['_GROUP_CUMSUM'] - to_plot['_VAR_INDEX'] 

    y="_VAR_INDEX"

    if rename_dic is None:
        rename_dic = {
            "BETA":"Per-allele effect size",
            "STUDY":"Study"
                      }
    ncols=1
    if eaf_panel and eaf in to_plot.columns:
        ncols+=1
    if snpvar_panel and snpr2 in to_plot.columns:
        ncols+=1

    if ncols==1:
        fig,ax1 = plt.subplots(ncols=ncols, **fig_kwargs)
    elif ncols==2:
        if eaf_panel==True and eaf in to_plot.columns:
            fig,axes = plt.subplots(ncols=ncols,sharey=True, **fig_kwargs)
            ax1=axes[0]
            ax2=axes[1]
        else:
            fig,axes = plt.subplots(ncols=ncols,sharey=True, **fig_kwargs)
            ax1=axes[0]
            ax3=axes[1]
    else:
        fig,axes = plt.subplots(ncols=ncols,sharey=True, **fig_kwargs)
        ax1=axes[0]
        ax2=axes[1]
        ax3=axes[2]

    # 上方已通过 set_plot_style 合并，无需再次调用

    sns.scatterplot(data=to_plot, x=x, y=y, ax=ax1, zorder=100, **{**args, **scatter_kwargs})

    ax1.errorbar(y=to_plot[y], x=to_plot[x], xerr=to_plot[se], 
                  **err_kwargs)
    
    ax1.axvline(x=0,linestyle="dashed",c="grey")
    ax1.set_yticks(to_plot[y], labels = to_plot[y_name], fontsize=fontsize, family=font_family)
    ax1.set_ylabel(ylabel, fontsize=fontsize, family=font_family) 
    ax1.set_xlabel(x, fontsize=fontsize, family=font_family) 

    if title is not None:
        ax1.set_title(title,fontsize=fontsize, family=font_family)

    if eaf_panel==True and eaf in to_plot.columns:
        ax2.barh(y=to_plot[y], width=to_plot[eaf], zorder=100, **eaf_kwargs)
        ax2.set_xlabel(eaf, fontsize=fontsize, family=font_family)
        if xlim_eaf is not None:
            ax3.set_xlim(xlim_eaf)

    if snpvar_panel==True and snpr2 in to_plot.columns:
        ax3.barh(y=to_plot[y], width=to_plot[snpr2], zorder=100,**snpr2_kwargs)
        ax3.set_xlabel(snpr2, fontsize=fontsize, family=font_family)
        if xlim_snpr2 is not None:
            ax3.set_xlim(xlim_snpr2)
    #try:
    if legend_mode==1:
        #if ncols==1:
        if len(legend_titles)>0: 
            sns.move_legend(
                ax1, "upper left",
                bbox_to_anchor=(1, 1), title=None, frameon=False, bbox_transform = axes[-1].transAxes, 
                title_fontproperties={"size":fontsize,"family":font_family},
                prop={"size":fontsize,"family":font_family}
                )
            #else:
##
            #    sns.move_legend(
            #        ax1, "lower left",
            #        bbox_to_anchor=(0, ncols), title=None, frameon=False,
            #    )
        #elif legend_mode==2:
        #    sns.move_legend(
        #        ax1, "lower center",
        #        bbox_to_anchor=(0, 1), ncol=ncol, title=None, frameon=False,
        #        )
    #except:
    #    pass


    #handles, labels = ax1.get_legend_handles_labels()
    #if len(labels)>0:
    #    #new_labels = []
    #    #ncol = len(labels)
    #    max_col=0
    #    new_labels=[]
    #    new_labels_i = []
    #    previous_i = 0
    #    max_string_len=0
    #    for i in range(len(labels)):
    #        if len(labels[i]) > max_string_len:
    #            max_string_len = len(labels[i])
    #        if labels[i] in legend_titles:
    #            new_labels_i.append(i)
    #            col_number = i - previous_i
    #            if col_number > max_col:
    #                max_col = col_number
    #            previous_i = i
    #    for i in labels:
    #        new_labels.append(str(i).ljust(max_string_len))
    #    print(new_labels)
    #    new_labels_i.append(len(labels))
#
    #    legend_rows = []
    #    #new_labels_i[index+1] - i
    #    for index, i in enumerate(new_labels_i):
    #        if index<len(new_labels_i)-1:
    #            legend_row = ax1.legend(labels = new_labels[i:new_labels_i[index+1]],  
    #                                    handles= handles[i:new_labels_i[index+1]],
    #                                    loc="lower left", 
    #                                    bbox_to_anchor=(-0.2, 1.02 + 0.05*index), 
    #                                    ncol=max_col, 
    #                                    scatterpoints=1, 
    #                                    title=None, 
    #                                    borderpad=0,
    #                                    handletextpad=0.1,
    #                                    handlelength=0.7,
    #                                    borderaxespad =0,
    #                                    alignment = "left",
    #                                    fontsize=8,
    #                                    frameon=False)
    #            legend_rows.append(legend_row)
    #    for legend_row in legend_rows[:-1]:
    #        ax1.add_artist(legend_row)



    ax1.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family) 
    
    if effect_label is not None:
        ax1.set_xlabel(effect_label, fontsize=fontsize, family=font_family)
        ax1.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family)
    if eaf_label is not None and eaf in to_plot.columns:
        ax2.set_xlabel(eaf_label, fontsize=fontsize, family=font_family)
        ax2.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family)
    if snpr2_label is not None and snpr2 in to_plot.columns:
        ax3.set_xlabel(snpr2_label, fontsize=fontsize, family=font_family)
        ax3.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family)
    save_figure(fig, save, keyword="forest",save_kwargs=save_kwargs, log=log, verbose=verbose)   
    
    return fig
