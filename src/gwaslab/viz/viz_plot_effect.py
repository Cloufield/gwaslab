import pandas as pd
import matplotlib.pyplot as plt
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.g_Log import Log
from gwaslab.io.io_process_args import _extract_kwargs

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
                 err_args=None,
                 font_args=None,
                 save=None,
                 title=None,
                 save_args=None,
                 eaf_args=None,
                 snpr2_args=None,
                 fig_args=None,
                 scatter_args=None,
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
                 style=None,
                 sort_args=None,
                 **args):

    if err_args is None:
        err_args={"ecolor":"#cccccc",
                  "linewidth":0,
                  "zorder":90,
                  "elinewidth":1}
    if eaf_args is None:
        eaf_args={"color":"#74BAD3"}
    if snpr2_args is None:
        snpr2_args={"color":"#74BAD3"}
    if font_args is None:
        font_args={'fontsize':12,'family':'sans','fontname':'Arial'}
    if fig_args is None:
        fig_args={"figsize":(8,8),"dpi":300}
    if scatter_args is None:
        scatter_args={"s":20}
    if sort_args is None:
        sort_args={}

    legend_titles=[]
    if hue is not None:
        args["hue"] = hue
        legend_titles.append(hue)

    if size is not None:
        args["size"] = size
        legend_titles.append(size)

    if style is not None:
        args["style"] = style
        legend_titles.append(style)

        
    save_kwargs =      _extract_kwargs("save", save_args, locals())
    err_kwargs =       _extract_kwargs("err", err_args, locals())
    scatter_kwargs =   _extract_kwargs("scatter", scatter_args, locals())
    font_kwargs =      _extract_kwargs("font",font_args, locals())
    
    def concat_cols(cols):
        string = "-".join(map(str,cols))
        return string
    
    y_name = "-".join(y)
    
    to_plot[y_name] = to_plot[y].apply(lambda x: concat_cols(x), axis=1)
    

    # sort y 
    if y_sort is None:
        y_sort = ["CHR","POS","STUDY"]
    
    #to_plot = to_plot.sort_values(by=y_sort)

    if group is None:
        group = ["CHR","POS"] + y_sort
    
    sort_columns= group + y_sort

    to_plot = to_plot.sort_values(by=sort_columns,**sort_args)
    
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
    if eaf_panel:
        ncols+=1
    if snpvar_panel:
        ncols+=1

    if ncols==1:
        fig,ax1 = plt.subplots(ncols=ncols, **fig_args)
    elif ncols==2:
        if eaf_panel==True:
            fig,axes = plt.subplots(ncols=ncols, dpi=400,sharey=True, **fig_args)
            ax1=axes[0]
            ax2=axes[1]
        else:
            fig,axes = plt.subplots(ncols=ncols, dpi=400,sharey=True, **fig_args)
            ax1=axes[0]
            ax3=axes[1]
    else:
        fig,axes = plt.subplots(ncols=ncols, dpi=400,sharey=True, **fig_args)
        ax1=axes[0]
        ax2=axes[1]
        ax3=axes[2]

    sns.scatterplot(data=to_plot, x=x, y=y, ax=ax1, zorder=100, **args)

    ax1.errorbar(y=to_plot[y], x=to_plot[x], xerr=to_plot[se], 
                  **err_kwargs)
    
    ax1.axvline(x=0,linestyle="dashed",c="grey")
    ax1.set_yticks(to_plot[y], labels = to_plot[y_name], fontsize=fontsize, family=font_family)
    ax1.set_ylabel(ylabel, fontsize=fontsize, family=font_family) 
    ax1.set_xlabel(x, fontsize=fontsize, family=font_family) 

    if title is not None:
        ax1.set_title(title,fontsize=fontsize, family=font_family)

    if eaf_panel==True:
        ax2.barh(y=to_plot[y], width=to_plot[eaf], zorder=100, **eaf_args)
        ax2.set_xlabel(eaf, fontsize=fontsize, family=font_family)
        if xlim_eaf is not None:
            ax3.set_xlim(xlim_eaf)

    if snpvar_panel==True:
        ax3.barh(y=to_plot[y], width=to_plot[snpr2], zorder=100,**snpr2_args)
        ax3.set_xlabel(snpr2, fontsize=fontsize, family=font_family)
        if xlim_snpr2 is not None:
            ax3.set_xlim(xlim_snpr2)
    #try:
    if legend_mode==1:
        #if ncols==1:
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
    if eaf_label is not None:
        ax2.set_xlabel(eaf_label, fontsize=fontsize, family=font_family)
        ax2.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family)
    if snpr2_label is not None:
        ax3.set_xlabel(snpr2_label, fontsize=fontsize, family=font_family)
        ax3.tick_params(axis='x', 
                        labelsize=fontsize,
                        labelfontfamily=font_family)
    save_figure(fig, save, keyword="forest",save_args=save_kwargs, log=log, verbose=verbose)   
    
    return fig