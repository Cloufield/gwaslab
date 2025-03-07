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
from gwaslab.viz_aux_chromatin import _plot_chromatin_state
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
from gwaslab.viz_plot_credible_sets import _plot_cs
import matplotlib.patches as patches

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
                        title_args=None,
                        #title_box = None,
                        gtf=None,
                        gene_track_height=0.5,
                        fig_args=None,
                        region_hspace=0.07,
                        subplot_height=4,
                        region_lead_grids = None,
                        region_lead_grid_line=None,
                        region_ld_legends = None,
                        fontsize=9,
                        font_family="Arial",
                        common_ylabel=True,
                        build="99", 
                        save=None,
                        save_args=None,
                        verbose=True,
                        pm=None,
                        log=Log(),
                        **mqq_args
                        ):
    
    log.write("Start to create stacked mqq plot by iteratively calling plot_mqq:",verbose=verbose)
    # load sumstats
    
    ##########################################################################################################################################
    if pm is None:
        pm=[]

    sumstats_list = [] 
    for each_object in objects:
        if type(each_object) is Sumstats:
            if "P" in each_object.data.columns or "MLOG10P" in each_object.data.columns:
                sumstats_list.append(each_object.data)
                pm.append("m")
        else:
            if "PIP" in each_object.columns:
                sumstats_list.append(each_object)
                pm.append("pip")
                common_ylabel=False
    
    if common_ylabel==True:
        rr_ylabel=False
    else:
        rr_ylabel=True

    log.write(" -Panel mode:{}...".format(pm),verbose=verbose)
    
    if fig_args is None:
        fig_args = {"dpi":200}
    if save_args is None:
        save_args = {"dpi":400,"facecolor":"white"}
    if region_lead_grid_line is None:
        region_lead_grid_line = {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}
    if region_chromatin_files is None:
        region_chromatin_files = []
        region_chromatin_height = len(region_chromatin_files) * region_chromatin_height
    if region_chromatin_labels is None:
        region_chromatin_labels = []
    if region_ld_legends is None:
        region_ld_legends = [0]
    if title_args is None:
        title_args = {"family":font_family}
    else:
        if "family" not in title_args.keys():
            title_args["family"] = font_family

    if save is not None:
        if type(save) is not bool:
            if len(save)>3:
                if save[-3:]=="pdf" or save[-3:]=="svg":
                    log.write(" -Adjusting options for saving as pdf/svg...",verbose=verbose)
                    fig_args["dpi"]=72
                    
                    if "scatter_args" not in  mqq_args.keys():
                        mqq_args["scatter_args"]={"rasterized":True}
                    else:
                        mqq_args["scatter_args"]["rasterized"] = True

                    if mode=="r":
                        if "scatter_args" not in  mqq_args.keys():
                            mqq_args["scatter_args"]={"rasterized":False}
                        else:
                            mqq_args["scatter_args"]["rasterized"] = False
                else:
                    fig_args["dpi"] = save_args["dpi"]
    # create figure and axes ##################################################################################################################
    #
    # subplot_height : subplot height
    # figsize : Width, height in inches

    if mode=="r":
        if len(vcfs)==1:
            vcfs = vcfs *len(sumstats_list)
        n_plot = len(sumstats_list)
        n_plot_plus_gene_track = n_plot + 1

        if len(region_chromatin_files)>0 and mode=="r":
            height_ratios = [1 for i in range(n_plot_plus_gene_track-1)]+[region_chromatin_height]+[gene_track_height]
            n_plot_plus_gene_track +=1
        else:
            height_ratios = [1 for i in range(n_plot_plus_gene_track-1)]+[gene_track_height]
        
        if "figsize" not in fig_args.keys():
            fig_args["figsize"] = [16,subplot_height*n_plot_plus_gene_track]

        fig, axes = plt.subplots(n_plot_plus_gene_track, 1, sharex=True, 
                             gridspec_kw={'height_ratios': height_ratios},
                             **fig_args)
        plt.subplots_adjust(hspace=region_hspace)
    elif mode=="m":
        n_plot = len(sumstats_list)
        if "figsize" not in fig_args.keys():
            fig_args["figsize"] = [10,subplot_height*n_plot]
        fig, axes = plt.subplots(n_plot, 1, sharex=True, 
                             gridspec_kw={'height_ratios': [1 for i in range(n_plot)]},
                             **fig_args)
        plt.subplots_adjust(hspace=region_hspace)
        vcfs = [None for i in range(n_plot)]
    elif mode=="mqq":
        n_plot = len(objects)
        if "figsize" not in fig_args.keys():
            fig_args["figsize"] = [10,subplot_height*n_plot]
        fig, axes = plt.subplots(n_plot, 2, sharex=True, 
                             gridspec_kw={'height_ratios': [1 for i in range(n_plot-1)],
                                          'width_ratios':[mqqratio,1]},
                                          **fig_args)
        plt.subplots_adjust(hspace=region_hspace)
    
    if region_lead_grids is None:
        region_lead_grids = [i for i in range(len(axes))]
    ##########################################################################################################################################
    mqq_args_for_each_plot = _sort_args(mqq_args, n_plot)

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
            fig,log,lead_snp_is,lead_snp_is_color = mqqplot(sumstats,
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
                            **mqq_args_for_each_plot[index]
                            )  
            lead_variants_is[index] = lead_snp_is
            lead_variants_is_color[index] = lead_snp_is_color
        else:
            if pm[index]=="m":
                # plot only the scatter plot
                fig,log,lead_snp_is,lead_snp_is_color = mqqplot(sumstats,
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
                                **mqq_args_for_each_plot[index]
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
                                  **mqq_args_for_each_plot[index])
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

    #if title_args is None:
    #    title_args = {}   
    #if titles is not None and mode=="r":    
    #    if title_pos is None:
    #        title_pos = [0.01,0.99]
    #    for index,title in enumerate(titles):
    #        
    #        current_text = axes[index].text(title_pos[0], title_pos[1] , title, transform=axes[index].transAxes,ha="left", va='top',zorder=999999, **title_args)
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
    if titles is not None:
        for index,title in enumerate(titles):
            axes[index].text(title_pos[0], title_pos[1] , title, transform=axes[index].transAxes,ha="left", va='top',zorder=999999, **title_args)
    
    ##########################################################################################################################################
    # draw the line for lead variants
    _draw_grid_line_for_lead_variants(mode, lead_variants_is,lead_variants_is_color, n_plot, axes, region_lead_grid_line,region_chromatin_files,region_lead_grids)
    
    ##########################################################################################################################################  
    if common_ylabel==True:
        _drop_old_y_labels(axes, n_plot)
        
        _add_new_y_label(mode, fig, gene_track_height,n_plot,subplot_height ,fontsize,font_family)
    
    ##########################################################################################################################################
    save_figure(fig = fig, save = save, keyword= "stacked_" + mode, save_args=save_args, log = log, verbose=verbose)
    
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
