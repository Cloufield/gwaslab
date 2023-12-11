import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import scipy as sp
import gwaslab as gl
from pyensembl import EnsemblRelease
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gc as garbage_collect
from adjustText import adjust_text
from gwaslab.g_Log import Log
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_get_sig import annogene
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_recombination_rate
from gwaslab.bd_common_data import get_gtf
from gwaslab.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz_aux_quickfix import _quick_fix
from gwaslab.viz_aux_quickfix import _get_largenumber
from gwaslab.viz_aux_quickfix import _quick_add_tchrpos
from gwaslab.viz_aux_quickfix import _quick_merge_sumstats
from gwaslab.viz_aux_quickfix import _quick_assign_i
from gwaslab.viz_aux_quickfix import _quick_assign_i_with_rank
from gwaslab.viz_aux_quickfix import _quick_extract_snp_in_region
from gwaslab.viz_aux_quickfix import _quick_assign_highlight_hue_pair
from gwaslab.viz_aux_quickfix import _quick_assign_marker_relative_size
from gwaslab.viz_aux_annotate_plot import annotate_pair
from gwaslab.io_to_pickle import load_pickle
from gwaslab.io_to_pickle import load_data_from_pickle
from gwaslab.g_Sumstats import Sumstats
from gwaslab.viz_aux_save_figure import save_figure
from gwaslab.viz_plot_mqqplot import mqqplot

def plot_stacked_mqq(objects,
                        vcfs,
                        mode="r",
                        region=None,
                        region_titles= None,
                        gtf=None,
                        gene_track_height=0.5,
                        fig_args=None,
                        region_hspace=0.05,
                        subplot_height=4,
                        region_lead_grid_line=None,
                        build="99", 
                        save=None,
                        save_args=None,
                        verbose=True,
                        **mqq_args
                        ):
    

    if fig_args is None:
        fig_args = {"dpi":200}
    if region_lead_grid_line is None:
        region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}
    # create figure and axes
    if mode=="r":
        if len(vcfs)==1:
            vcfs = vcfs *len(objects)
        n_plot = len(objects)
        n_plot_plus_gene_track = n_plot + 1

        fig_args["figsize"] = [16,subplot_height*n_plot_plus_gene_track]
        fig, axes = plt.subplots(n_plot_plus_gene_track, 1, sharex=True, 
                             gridspec_kw={'height_ratios': [1 for i in range(n_plot_plus_gene_track-1)]+[gene_track_height]},
                             **fig_args)
        plt.subplots_adjust(hspace=region_hspace)

    if region_titles is None:
        region_titles=["" for i in range(n_plot) ]
    #elif mode=="m":
    #    n_plot = len(objects)
#
    #    fig_args["figsize"] = [10,subplot_height*n_plot]
    #    fig, axes = plt.subplots(n_plot, 1, sharex=True, 
    #                         gridspec_kw={'height_ratios': [1 for i in range(n_plot-1)]},
    #                         **fig_args)
    #    plt.subplots_adjust(hspace=region_hspace)
#
    #elif mode=="mqq":
    #    n_plot = len(objects)
#
    #    fig_args["figsize"] = [10,subplot_height*n_plot]
    #    fig, axes = plt.subplots(n_plot, 1, sharex=True, 
    #                         gridspec_kw={'height_ratios': [1 for i in range(n_plot-1)]},
    #                         **fig_args)
    #    plt.subplots_adjust(hspace=region_hspace)

    mqq_args_for_each_plot = _sort_args(mqq_args, n_plot)
    
    # load sumstats
    sumstats_list = objects
    
    lead_variants_is={}
    # plot manhattan plot
    for index,sumstats in enumerate(sumstats_list):
        if index==0:
            # plot last m and gene track 
            fig,log,lead_i,lead_i2 = mqqplot(sumstats,
                            chrom="CHR",
                            pos="POS",
                            p="P",
                            region=region,
                            mlog10p="MLOG10P",
                            snpid="SNPID",
                            vcf_path=vcfs[index],
                            region_title = region_titles[index],
                            mode=mode,
                            region_lead_grid=False,
                            gtf_path="default",
                            rr_ylabel=False,
                            figax=(fig,axes[index],axes[-1]),
                            _get_region_lead=True,
                            build=build,verbose=verbose,
                            **mqq_args_for_each_plot[index]
                            )  
            lead_variants_is[index] = (lead_i,lead_i2)
        else:
            fig,log,lead_i,lead_i2 = mqqplot(sumstats,
                            chrom="CHR",
                            pos="POS",
                            p="P",
                            region=region,
                            mlog10p="MLOG10P",
                            snpid="SNPID",
                            region_title = region_titles[index],
                            vcf_path=vcfs[index],
                            region_lead_grid=False,
                            mode=mode,
                            rr_ylabel=False,
                            region_ld_legend=False,
                            gtf_path=None,
                            figax=(fig,axes[index],axes[-1]),
                            _get_region_lead=True,
                            build=build,verbose=verbose,
                            **mqq_args_for_each_plot[index]
                            )
            lead_variants_is[index] = (lead_i,lead_i2)
    
    # adjust labels
    # drop labels for each plot 
    # set a common laebl for all plots
    for index in range(n_plot):
        axes[index].set_ylabel("")
    
    for index, sig_is in lead_variants_is.items():
        for sig_i in sig_is:
            if sig_i is not None:
                for each_axis_index in range(n_plot + 1):    
                    axes[each_axis_index].axvline(x=sig_i, zorder=2,**region_lead_grid_line)


    gene_track_height_ratio = gene_track_height/(gene_track_height + n_plot*subplot_height)
    ylabel_height = (1 - gene_track_height_ratio)*0.5 + gene_track_height_ratio
    fig.text(0.08, ylabel_height , "$-log_{10}(P)$", va='center', rotation='vertical')
    fig.text(0.93, ylabel_height, "Recombination rate(cM/Mb)", va='center', rotation=-90)

    save_figure(fig = fig, save = save, keyword= "stacked_" + mode, save_args=save_args, log = log, verbose=verbose)

    return fig


def _sort_args(mqq_args, n_plot):

    mqq_args_for_each_plot={i:{} for i in range(n_plot)}
    
    for key, value in mqq_args.items():
        if key[-1].isnumeric():
            mqq_args_for_each_plot[int(key[-1])-1][key[:-1]]=value
        else:
            # for both top and bottom
            for i in range(n_plot):
                mqq_args_for_each_plot[i][key] = value
    return mqq_args_for_each_plot