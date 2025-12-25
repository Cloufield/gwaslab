import gc as garbage_collect
import numpy as np
import pandas as pd
import scipy as sp
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from adjustText import adjust_text
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between

import gwaslab as gl
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.io.io_gtf import get_gtf
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_recombination_rate
from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats
from gwaslab.io.io_to_pickle import load_data_from_pickle
from gwaslab.io.io_to_pickle import load_pickle
from gwaslab.util.util_in_get_sig import _anno_gene
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.viz.viz_aux_annotate_plot import annotate_pair
from gwaslab.viz.viz_aux_chromatin import _plot_chromatin_state
from gwaslab.viz.viz_aux_quickfix import (
    _get_largenumber,
    _quick_add_tchrpos,
    _quick_assign_highlight_hue_pair,
    _quick_assign_i,
    _quick_assign_i_with_rank,
    _quick_assign_marker_relative_size,
    _quick_extract_snp_in_region,
    _quick_fix,
)
from gwaslab.viz.viz_aux_reposition_text import adjust_text_position
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.viz.viz_plot_mqqplot import _mqqplot


def plot_stacked_mqq(   objects,
                        vcfs=None,
                        mode="r",
                        mqqratio=3,
                        region=None,
                        region_chromatin_height=0.1,
                        region_chromatin_files = None,
                        region_chromatin_labels= None,
                        titles= None,
                        title_pos=None,
                        title_kwargs=None,
                        #title_box = None,
                        gtf=None,
                        mqq_height=1,
                        cs_height=0.5,
                        gene_track_height=0.5,
                        fig_kwargs=None,
                        region_hspace=0.07,
                        subplot_height=4,
                        region_lead_grids = None,
                        region_ld_legends = None,
                        fontsize=9,
                        font_family="Arial",
                        common_ylabel=True,
                        build="99", 
                        save=None,
                        save_kwargs=None,
                        verbose=True,
                        pm=None,
                        log=Log(),
                        **mqq_kwargs
                        ):
    
    log.write("Start to create stacked mqq plot by iteratively calling plot_mqq:",verbose=verbose)
    # load sumstats
    log.write(" -Stacked plot mode:{}...".format(mode),verbose=verbose)

    ##########################################################################################################################################
    # a list of modes for each panel
    if pm is None:
        pm=[]

    sumstats_list = [] 
    for each_object in objects:
        if type(each_object) is Sumstats:
            if "P" in each_object.data.columns or "MLOG10P" in each_object.data.columns:
                sumstats_list.append(each_object.data)
                pm.append("m")
        if type(each_object) is pd.DataFrame:
            if "PIP" in each_object.columns:
                sumstats_list.append(each_object)
                pm.append("pip")
                common_ylabel=False
            else:
                sumstats_list.append(each_object)
                pm.append("m")

    log.write(" -Panel mode:{}...".format(pm),verbose=verbose)
    ##########################################################################################################################################

    if common_ylabel==True:
        rr_ylabel=False
    else:
        rr_ylabel=True
    
    style = set_plot_style(
        plot="plot_stacked_mqq",
        fig_kwargs=fig_kwargs if fig_kwargs is not None else fig_kwargs,
        save_kwargs=save_kwargs if save_kwargs is not None else save_kwargs,
        save=save,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", style.get("fig_kwargs", {}))
    save_kwargs = style.get("save_kwargs", style.get("save_kwargs", {}))
    fontsize = style["fontsize"]
    font_family = style["font_family"]
    # title kwargs: merge priority -> explicit param -> mqq_kwargs -> style defaults
    base_title = style.get("title_kwargs", style.get("title_kwargs", {}))
    explicit_title_kwargs = title_kwargs
    title_kwargs = dict(base_title)
    if title_kwargs is not None:
        title_kwargs.update(title_kwargs)
    if explicit_title_kwargs is not None:
        title_kwargs.update(explicit_title_kwargs)
    extra = mqq_kwargs.get("title_kwargs") or mqq_kwargs.get("title_kwargs")
    if extra is not None:
        title_kwargs.update(extra)
    #title_pos = mqq_kwargs.get("title_pos") or style.get("title_pos", None)

    if region_chromatin_files is None:
        region_chromatin_files = []
        region_chromatin_height = len(region_chromatin_files) * region_chromatin_height
    if region_chromatin_labels is None:
        region_chromatin_labels = []
    if region_ld_legends is None:
        region_ld_legends = [0]
    if not title_kwargs:
        title_kwargs = {"family": font_family}
    else:
        if "family" not in title_kwargs.keys():
            title_kwargs["family"] = font_family

    from gwaslab.viz.viz_aux_style_options import figure_kwargs_for_vector_plot
    fig_kwargs, mqq_kwargs_scatter = figure_kwargs_for_vector_plot(save, fig_kwargs, mqq_kwargs.get("scatter_kwargs"))
    if mqq_kwargs_scatter is not None:
        mqq_kwargs["scatter_kwargs"] = mqq_kwargs_scatter
    # create figure and axes ##################################################################################################################
    #
    # subplot_height : subplot height
    # figsize : Width, height in inches
    if mode=="r":
    ##########################################################################################################################################   
        if not (len(vcfs)==1 or len(vcfs)==len(sumstats_list)):
            raise ValueError("Please make sure VCFs match Objects!")
    
        if len(vcfs)==1:
            vcfs = vcfs *len(sumstats_list)
    ##########################################################################################################################################    
        height_ratios=[]
        for index, i in enumerate(pm):
            if i =="m":
                height_ratios.append(mqq_height)
            elif i=="pip":
                height_ratios.append(cs_height)
                vcfs[index] = "NA"

        log.write(" -VCFs:{}...".format(vcfs),verbose=verbose)

        # n: sumstats
        # +1 : gene track
        n_plot = len(sumstats_list)
        n_plot_plus_gene_track = n_plot + 1
        
        if len(region_chromatin_files)>0 and mode=="r":
            #height_ratios = [1 for i in range(n_plot_plus_gene_track-1)]+[region_chromatin_height]+[gene_track_height]
            height_ratios += [region_chromatin_height]+[gene_track_height]
            # +1 : region_chromatin_files
            n_plot_plus_gene_track +=1
        else:
            #height_ratios = [1 for i in range(n_plot_plus_gene_track-1)]+[gene_track_height]
            height_ratios += [gene_track_height]
        
        if "figsize" not in fig_kwargs.keys():
            fig_kwargs["figsize"] = [16,subplot_height*n_plot_plus_gene_track]

        fig, axes = plt.subplots(n_plot_plus_gene_track, 1, sharex=True, 
                             gridspec_kw={'height_ratios': height_ratios},
                             **fig_kwargs)
        
        plt.subplots_adjust(hspace=region_hspace)
    elif mode=="m":
        n_plot = len(sumstats_list)
        if "figsize" not in fig_kwargs.keys():
            fig_kwargs["figsize"] = [10,subplot_height*n_plot]
        fig, axes = plt.subplots(n_plot, 1, sharex=True, 
                             gridspec_kw={'height_ratios': [1 for i in range(n_plot)]},
                             **fig_kwargs)
        plt.subplots_adjust(hspace=region_hspace)
        vcfs = [None for i in range(n_plot)]
    elif mode=="mqq":
        n_plot = len(objects)
        if "figsize" not in fig_kwargs.keys():
            fig_kwargs["figsize"] = [10,subplot_height*n_plot]
        fig, axes = plt.subplots(n_plot, 2, sharex=True, 
                             gridspec_kw={'height_ratios': [1 for i in range(n_plot)],
                                          'width_ratios':[mqqratio,1]},
                                          **fig_kwargs)
        plt.subplots_adjust(hspace=region_hspace)
        vcfs = [None for i in range(n_plot)]
    if region_lead_grids is None:
        region_lead_grids = [i for i in range(len(axes))]
    ##########################################################################################################################################
    mqq_kwargs_for_each_plot = _sort_kwargs(mqq_kwargs, n_plot)

    ##########################################################################################################################################
    # get x axis dict
    if mode=="m" or mode=="r":
        _posdiccul = _get_chrom_dic(sumstats_list,chrom="CHR",pos="POS",chrpad=0.02)
    else:
        _posdiccul=None
    
    ##########################################################################################################################################
    # a dict to store lead variants of each plot
    lead_variants_is={}
    lead_variants_is_color ={}

    ##########################################################################################################################################
    # plot manhattan plot
    for index,sumstats in enumerate(sumstats_list):

        #################################################################
        if mode=="m" or mode=="r":
            figax = (fig,axes[index],axes[-1])
        elif mode=="mqq":
            figax = (fig,axes[index,0],axes[index,1])
        if index in region_ld_legends:
            region_ld_legend = True
        else:
            region_ld_legend = False
        #################################################################
        if index==0:
            # plot last m and gene track 
            # Remove region_ld_legend and region_lead_grid from kwargs since they're set explicitly
            plot_kwargs = dict(mqq_kwargs_for_each_plot[index])
            plot_kwargs.pop("region_ld_legend", None)
            plot_kwargs.pop("region_lead_grid", None)
            fig,log,lead_snp_is,lead_snp_is_color = _mqqplot(sumstats,
                            chrom="CHR",
                            pos="POS",
                            p="P",
                            region=region,
                            mlog10p="MLOG10P",
                            snpid="SNPID",
                            vcf_path=vcfs[index],
                            mode=mode,
                            fontsize=fontsize,
                            font_family=font_family,
                            region_lead_grid=False,
                            region_ld_legend=region_ld_legend,
                            gtf_path="default",
                            rr_ylabel=rr_ylabel,
                            figax=figax,
                            _get_region_lead=True,
                            _if_quick_qc=False,
                            _posdiccul=_posdiccul,
                            build=build,
                            verbose=verbose,
                            log=log,
                            **plot_kwargs
                            )  
            lead_variants_is[index] = lead_snp_is
            lead_variants_is_color[index] = lead_snp_is_color
        else:
            if pm[index]=="m":
                # plot only the scatter plot
                # Remove region_ld_legend and region_lead_grid from kwargs since they're set explicitly
                plot_kwargs = dict(mqq_kwargs_for_each_plot[index])
                plot_kwargs.pop("region_ld_legend", None)
                plot_kwargs.pop("region_lead_grid", None)
                fig,log,lead_snp_is,lead_snp_is_color = _mqqplot(sumstats,
                                chrom="CHR",
                                pos="POS",
                                p="P",
                                region=region,
                                mlog10p="MLOG10P",
                                snpid="SNPID",
                                vcf_path=vcfs[index],
                                region_lead_grid=False,
                                fontsize=fontsize,
                                font_family=font_family,
                                mode=mode,
                                rr_ylabel=rr_ylabel,
                                region_ld_legend=region_ld_legend,
                                gtf_path=None,
                                figax=figax,
                                _get_region_lead=True,
                                _if_quick_qc=False,
                                _posdiccul=_posdiccul,
                                build=build,
                                verbose=verbose,
                                log=log,
                                **plot_kwargs
                                )
                lead_variants_is[index] = lead_snp_is
                lead_variants_is_color[index] = lead_snp_is_color
            elif pm[index]=="pip":
                fig,log =_plot_cs(sumstats,
                                  region=region,
                                  _posdiccul=_posdiccul,
                                  figax=figax,
                                  log=log,
                                  verbose=verbose,
                                  **mqq_kwargs_for_each_plot[index])
    if len(region_chromatin_files)>0 and mode=="r":
        xlim_i = axes[-1].get_xlim()
        fig = _plot_chromatin_state(     region_chromatin_files = region_chromatin_files, 
                                         region_chromatin_labels = region_chromatin_labels,
                                         region = region, 
                                         fig = fig, 
                                         ax = axes[-2],
                                         xlim_i=xlim_i,
                                         log=log,
                                         verbose=verbose,
                                         fontsize = fontsize,
                                         font_family = font_family)    
    # adjust labels
    # drop labels for each plot 
    # set a common laebl for all plots
    #if title_box is None:
    #    title_box = dict(boxstyle='square', facecolor='white', alpha=1.0, edgecolor="black")
    #    title_box = {}

    #if title_kwargs is None:
    #    title_kwargs = {}   
    #if titles is not None and mode=="r":    
    #    if title_pos is None:
    #        title_pos = [0.01,0.99]
    #    for index,title in enumerate(titles):
    #        
    #        current_text = axes[index].text(title_pos[0], title_pos[1] , title, transform=axes[index].transAxes,ha="left", va='top',zorder=999999, **title_kwargs)
    #        r = fig.canvas.get_renderer()
    #        bb = current_text.get_window_extent(renderer=r).transformed(axes[index].transAxes.inverted())
    #        width = bb.width
    #        height = bb.height
#
    #        rect = patches.Rectangle((0.0,1.0 - height),
    #                        height=height + 0.02*2,
    #                        width=width + 0.01*2, 
    #                        transform=axes[index].transAxes,
    #                        linewidth=1, 
    #                        edgecolor='black', 
    #                        facecolor='white',
    #                        alpha=1.0,
    #                        zorder=99998)
    #        axes[index].add_patch(rect)
    #        rect.set(zorder=99998) 
    #else:
    if title_pos is None:
        title_pos = [0.01,0.97]

    if mode=="mqq":
        axes = [axes[i,0] for i in range(len(objects))]

    if titles is not None:
        for index, title in enumerate(titles):
            axes[index].text(
                title_pos[0],
                title_pos[1],
                title,
                transform=axes[index].transAxes,
                ha="left",
                va='top',
                zorder=999999,
                **(title_kwargs or {"family": font_family, "fontsize": fontsize})
            )
    
    ##########################################################################################################################################
    # draw the line for lead variants
    region_lead_grid_line = mqq_kwargs.get("region_lead_grid_line", {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"})
    _draw_grid_line_for_lead_variants(mode, lead_variants_is,lead_variants_is_color, n_plot, axes, region_lead_grid_line,region_chromatin_files,region_lead_grids)
    
    ##########################################################################################################################################  
    if common_ylabel==True:
        _drop_old_y_labels(axes, n_plot)
        
        _add_new_y_label(mode, fig, gene_track_height,n_plot,subplot_height ,fontsize,font_family)
    
    ##########################################################################################################################################
    save_figure(fig = fig, save = save, keyword= "stacked_" + mode, save_kwargs=save_kwargs, log = log, verbose=verbose)
    
    log.write("Finished creating stacked mqq plot by iteratively calling plot_mqq.",verbose=verbose)
    
    return fig, log

def _drop_old_y_labels(axes, n_plot):
    for index in range(n_plot):
        axes[index].set_ylabel("")

def _draw_grid_line_for_lead_variants(mode, lead_variants_is,lead_variants_is_color, n_plot, axes, region_lead_grid_line,region_chromatin_files,region_lead_grids):
    if len(region_chromatin_files)>0:
        n_plot_and_track = n_plot+2
    else:
        n_plot_and_track = n_plot+1
    
    plotted=[None]
    if mode=="r":
        for index, sig_is in lead_variants_is.items():
            if index in region_lead_grids:
                for j, sig_i in enumerate(sig_is):
                    try:
                        region_lead_grid_line["color"] = lead_variants_is_color[index][j]
                    except:
                        pass
                    if sig_i not in plotted:
                        for each_axis_index in range(n_plot_and_track):  
                            axes[each_axis_index].axvline(x=sig_i, zorder=2,**region_lead_grid_line)

def _add_new_y_label(mode, fig, gene_track_height,n_plot,subplot_height ,fontsize,font_family):
    gene_track_height_ratio = gene_track_height/(gene_track_height + n_plot*subplot_height)
    ylabel_height = (1 - gene_track_height_ratio)*0.5 + gene_track_height_ratio
    if mode=="r":
        fig.text(0.08, ylabel_height , "$\mathregular{-log_{10}(P)}$", va='center', rotation='vertical',
                 fontsize=fontsize,
                 family=font_family)
        
        fig.text(0.93, ylabel_height, "Recombination rate(cM/Mb)", va='center', rotation=-90,fontsize=fontsize,family=font_family)
    elif mode=="m":
        fig.text(0.08, ylabel_height , "$\mathregular{-log_{10}(P)}$", va='center', rotation='vertical',
                 fontsize=fontsize,
                 family=font_family)    

def _sort_kwargs(mqq_kwargs, n_plot):
    mqq_kwargs_for_each_plot={i:{} for i in range(n_plot)}
    
    for key, value in mqq_kwargs.items():
        if key[-1].isnumeric():
            mqq_kwargs_for_each_plot[int(key[-1])-1][key[:-1]]=value
        else:
            # for both top and bottom
            for i in range(n_plot):
                mqq_kwargs_for_each_plot[i][key] = value
    return mqq_kwargs_for_each_plot

def _get_chrom_dic(sumstats_list,chrom="CHR",pos="POS",chrpad=0.02):
    posdiccul = {}
    max_chr = 0 
    max_pos = 0
    for sumstats in sumstats_list:
        posdic = sumstats.groupby(chrom)[pos].max()
        if sumstats[chrom].max() > max_chr:
            max_chr = sumstats[chrom].max()
        if sumstats[pos].max() > max_chr:
            max_pos = sumstats[pos].max()
        # convert to dictionary
        posdic = dict(posdic)
            
        # fill empty chr with 0
        for i in posdic.keys():
            if i in posdiccul.keys():
                if posdic[i] > posdiccul[i]:
                    posdiccul[i] = posdic[i]
            else:
                posdiccul[i] = posdic[i]

    for i in range(0,max_chr+1):
        if i in posdiccul: 
            continue
        else: 
            posdiccul[i]=0
            # cumulative sum dictionary
    posdiccul_raw = posdiccul.copy()
    for i in range(1,sumstats[chrom].max()+1):
        posdiccul[i]= posdiccul[i-1] + posdiccul[i] + max_pos*chrpad
    return posdiccul
