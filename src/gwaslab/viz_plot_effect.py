import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
import copy
from math import ceil
from shutil import which
from pyensembl import EnsemblRelease
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import MaxNLocator
import gc as garbage_collect
from adjustText import adjust_text
from gwaslab.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz_aux_annotate_plot import annotate_single
from gwaslab.viz_plot_qqplot import _plot_qq
from gwaslab.hm_harmonize_sumstats import auto_check_vcf_chr_dict
from gwaslab.viz_plot_regional2 import _plot_regional
from gwaslab.viz_plot_regional2 import process_vcf
from gwaslab.viz_plot_regional2 import _get_lead_id
from gwaslab.viz_aux_quickfix import _get_largenumber
from gwaslab.viz_aux_quickfix import _quick_fix_p_value
from gwaslab.viz_aux_quickfix import _quick_fix_pos
from gwaslab.viz_aux_quickfix import _quick_fix_chr
from gwaslab.viz_aux_quickfix import _quick_fix_eaf
from gwaslab.viz_aux_quickfix import _quick_fix_mlog10p
from gwaslab.viz_aux_quickfix import _quick_add_tchrpos
from gwaslab.viz_aux_quickfix import _quick_merge_sumstats
from gwaslab.viz_aux_quickfix import _quick_assign_i
from gwaslab.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz_aux_quickfix import _quick_extract_snp_in_region
from gwaslab.viz_aux_quickfix import _quick_assign_highlight_hue_pair
from gwaslab.viz_aux_quickfix import _quick_assign_marker_relative_size
from gwaslab.viz_aux_quickfix import _cut
from gwaslab.viz_aux_quickfix import _set_yticklabels
from gwaslab.viz_aux_quickfix import _jagged_y
from gwaslab.viz_aux_save_figure import save_figure
from gwaslab.g_Log import Log
from gwaslab.util_in_calculate_gc import lambdaGC
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_get_sig import annogene
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_recombination_rate
from gwaslab.bd_common_data import get_gtf
from gwaslab.util_in_filter_value import _filter_region
from gwaslab.g_version import _get_version
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import to_hex
from gwaslab.io_process_args import _extract_kwargs

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
                 log=Log(),
                 verbose=True,
                 legend_mode=1,
                 ncol=2,
                 size=None,
                 hue=None,
                 style=None,
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
    
    to_plot = to_plot.sort_values(by=y_sort)

    if group is None:
        group = ["CHR","POS"]
        to_plot = to_plot.sort_values(by=group)
    # Assign group IDs based on the sorted 'score'
        
    to_plot['_VAR_GROUP'] = to_plot.groupby(group).ngroup() + 1

    to_plot["_VAR_INDEX"] = range(len(to_plot))
    to_plot["_VAR_INDEX"]= to_plot["_VAR_INDEX"] + to_plot['_VAR_GROUP']

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
            fig,axes = plt.subplots(ncols=ncols, dpi=400,sharey=True)
            ax1=axes[0]
            ax2=axes[1]
        else:
            fig,axes = plt.subplots(ncols=ncols, dpi=400,sharey=True)
            ax1=axes[0]
            ax3=axes[1]
    else:
        fig,axes = plt.subplots(ncols=ncols, dpi=400,sharey=True)
        ax1=axes[0]
        ax2=axes[1]
        ax3=axes[2]

    sns.scatterplot(data=to_plot, x=x, y=y, ax=ax1, zorder=100, **args)

    ax1.errorbar(y=to_plot[y], x=to_plot[x], xerr=to_plot[se], 
                  **err_kwargs)
    
    ax1.axvline(x=0,linestyle="dashed",c="grey")
    ax1.set_yticks(to_plot[y], labels = to_plot[y_name])
    ax1.set_ylabel(ylabel) 

    if title is not None:
        ax1.set_title(title)

    if eaf_panel==True:
        ax2.barh(y=to_plot[y], width=to_plot[eaf], zorder=100, **eaf_args)
        ax2.set_xlabel(eaf)

    if snpvar_panel==True:
        ax3.barh(y=to_plot[y], width=to_plot[snpr2], zorder=100,**snpr2_args)
        ax3.set_xlabel(snpr2)
    
    #try:
    if legend_mode==1:
        #if ncols==1:
        sns.move_legend(
            ax1, "upper left",
            bbox_to_anchor=(1, 1), title=None, frameon=False, bbox_transform = axes[-1].transAxes
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



    if effect_label is not None:
        ax1.set_xlabel(effect_label)
    if eaf_label is not None:
        ax2.set_xlabel(eaf_label)
    if snpr2_label is not None:
        ax3.set_xlabel(snpr2_label)

    save_figure(fig, save, keyword="forest",save_args=save_kwargs, log=log, verbose=verbose)   
    
    return fig