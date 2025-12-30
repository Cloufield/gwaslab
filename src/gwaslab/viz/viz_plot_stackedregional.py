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
from gwaslab.util.util_in_filter_value import _filter_region
from gwaslab.qc.qc_normalize_args import _normalize_region
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
    """
    Create stacked Manhattan/QQ/regional plots by iteratively calling plot_mqq.
    
    Parameters
    ----------
    objects : list
        List of Sumstats objects or DataFrames to plot
    vcfs : list, optional
        List of VCF file paths for LD calculation. Default is None.
    mode : str, default="r"
        Plot mode: "r" for regional, "m" for Manhattan, "mqq" for Manhattan+QQ
    region_ld_legends : list, optional, default=None
        List of panel indices to show LD legends (cbar). 
        If None, automatically enables legends for all mqq panels.
        After plotting, if all panels with legends share the same legend (same lead variant),
        duplicate legends are automatically removed (keeping only the first panel's legend).
        If panels have different legends, all legends are kept.
    region_lead_grids : list, optional, default=None
        List of panel indices to draw vertical grid lines for lead variants.
        If None, draws grid lines for all panels.
    **mqq_kwargs
        Additional keyword arguments passed to plot_mqq for each panel.
        Panel-specific arguments can be specified with numeric suffix (e.g., region_ref1, region_ref2).
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib figure object.
    log : gwaslab.Log
        Updated logging object with operation records.
    
    Notes
    -----
    - By default, LD legends are automatically enabled for all mqq panels
    - Duplicate legends are removed only when ALL panels with legends share the same lead variant
    - If any panel has a different lead variant, all legends are kept to show the differences
    """
    
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
        # By default, enable LD legends for all mqq panels (will check for duplicates later)
        region_ld_legends = None  # Will be set after we know which panels are mqq
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
    # Auto-enable LD legends for all mqq panels by default
    if region_ld_legends is None:
        mqq_panel_indices = [i for i in range(n_plot) if pm[i] == "m"]
        region_ld_legends = mqq_panel_indices if len(mqq_panel_indices) > 0 else [0]
    
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
    
    ##########################################################################################################################################
    # Check for duplicate LD legends (cbar) and remove them ONLY when ALL mqq panels share the same legend
    # Only check panels that actually have legends enabled (are in region_ld_legends)
    if mode == "r":
        mqq_panel_indices = [i for i in range(n_plot) if pm[i] == "m"]
        # Get mqq panels that have legends enabled
        mqq_panels_with_legends = [idx for idx in mqq_panel_indices if idx in region_ld_legends]
        if len(mqq_panels_with_legends) > 1:
            # Check for duplicates only among panels that have legends enabled
            _remove_duplicate_ld_legends(fig, axes, mqq_panels_with_legends, lead_variants_is, log, verbose)
    
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
    print(lead_variants_is, lead_variants_is_color)
    #{0: [2793163258], 1: [2793301317]} {0: ['#FF0000'], 1: ['#FF0000']}
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

def _remove_duplicate_ld_legends(fig, axes, mqq_panel_indices, lead_variants_is, log, verbose):
    """
    Remove duplicate LD legends (cbar) ONLY when ALL mqq panels share the same legend.
    
    Logic: 
    - If ANY panel has a different legend, keep ALL legends (do NOT remove any)
    - Only when ALL panels have the EXACT SAME legend, remove duplicates (keep first one)
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure containing all axes
    axes : list
        List of axes for all panels
    mqq_panel_indices : list
        List of indices for mqq panels
    lead_variants_is : dict
        Dictionary mapping panel index to list of lead variant positions (i values)
    log : gwaslab.Log
        Logging object
    verbose : bool
        Whether to show verbose messages
    """
    # Step 1: Collect lead variant positions for all mqq panels
    lead_variant_positions = {}
    for idx in mqq_panel_indices:
        if idx in lead_variants_is and len(lead_variants_is[idx]) > 0:
            # Convert to tuple for comparison (sorted to handle multiple lead variants)
            lead_variant_positions[idx] = tuple(sorted(lead_variants_is[idx]))
        else:
            lead_variant_positions[idx] = None
    
    # Step 2: Check if ALL panels have the same lead variant positions
    # Only proceed if ALL panels have valid positions and they are ALL identical
    valid_positions = {k: v for k, v in lead_variant_positions.items() if v is not None}
    
    # Condition 1: All panels must have valid positions
    # Condition 2: All positions must be identical (only one unique value)
    # If ANY panel has a different position, we keep ALL legends
    if len(valid_positions) == len(mqq_panel_indices):
        unique_positions = set(valid_positions.values())
        # Only remove duplicates if ALL panels have the EXACT SAME legend
        if len(unique_positions) == 1:
            # ALL panels have the same lead variant positions - safe to remove duplicates
            # Keep the first panel's legend, remove from all others
            first_idx = list(valid_positions.keys())[0]
            if verbose:
                log.write(f" -All panels share the same lead variant positions: {valid_positions[first_idx]}. Removing duplicate legends.", verbose=verbose)
            for idx in list(valid_positions.keys())[1:]:
                _remove_cbar_from_axis(fig, axes[idx], idx, first_idx, log, verbose)
        else:
            # Panels have different legends - keep ALL legends
            if verbose:
                log.write(f" -Panels have different lead variant positions. Keeping all legends.", verbose=verbose)
    else:
        # Not all panels have valid positions
        if verbose:
            log.write(f" -Not all panels have valid lead variant positions. Keeping all legends.", verbose=verbose)

def _remove_cbar_from_axis(fig, ax, panel_idx, first_panel_idx, log, verbose):
    """
    Helper function to find and remove cbar (LD legend) from a specific axis.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure containing all axes
    ax : matplotlib.axes.Axes
        The axis to remove cbar from
    panel_idx : int
        Index of the panel (for logging)
    first_panel_idx : int
        Index of the first panel that keeps its legend (for logging)
    log : gwaslab.Log
        Logging object
    verbose : bool
        Whether to show verbose messages
    """
    cbar_found = False
    ax_pos = ax.get_position()
    
    # Method 1: Check all axes in the figure to find inset axes (cbar)
    # Check characteristics first, then verify position
    for fig_ax in list(fig.axes):
        if fig_ax == ax:
            continue
        try:
            # First check if it has LD legend characteristics
            yticklabels = [label.get_text() for label in fig_ax.get_yticklabels()]
            patches = fig_ax.patches
            
            # LD legend cbar has yticks (region_ref names) and patches (rectangles for LD thresholds)
            # Check if it has the characteristics of an LD legend
            if len(yticklabels) > 0 and len(patches) >= 3:  # LD legend has many patches
                # Check if it's positioned within or overlapping the parent axes (inset axes)
                pos = fig_ax.get_position()
                # Check for overlap with parent axes (very lenient - any overlap)
                overlaps = (pos.x0 < ax_pos.x1 and pos.x1 > ax_pos.x0 and 
                           pos.y0 < ax_pos.y1 and pos.y1 > ax_pos.y0)
                
                # Also check if it's smaller than parent (inset axes are typically much smaller)
                pos_width = pos.x1 - pos.x0
                pos_height = pos.y1 - pos.y0
                ax_width = ax_pos.x1 - ax_pos.x0
                ax_height = ax_pos.y1 - ax_pos.y0
                is_smaller = (pos_width < ax_width * 0.6 and pos_height < ax_height * 0.6)
                
                # If it has yticks, many patches, and (overlaps with parent OR is smaller), it's likely the cbar
                if overlaps or is_smaller:
                    # Additional check: white facecolor (LD legend has white background)
                    try:
                        facecolor = fig_ax.get_facecolor()
                        # Handle different facecolor formats
                        if isinstance(facecolor, (tuple, list, np.ndarray)):
                            if len(facecolor) >= 3:
                                is_white = (abs(float(facecolor[0]) - 1.0) < 0.1 and 
                                           abs(float(facecolor[1]) - 1.0) < 0.1 and 
                                           abs(float(facecolor[2]) - 1.0) < 0.1)
                            else:
                                is_white = False
                        else:
                            is_white = False
                    except Exception:
                        is_white = False
                    
                    # Strong indicators: white facecolor OR (many patches AND overlaps/smaller)
                    if is_white or (len(patches) >= 3):
                        # This is likely the LD legend cbar - remove it
                        fig_ax.remove()
                        log.write(f" -Removed duplicate LD legend (cbar) from panel {panel_idx} (all panels share the same legend).", verbose=verbose)
                        cbar_found = True
                        break
        except Exception as e:
            if verbose:
                log.write(f" -Warning: Error checking axis for cbar in panel {panel_idx}: {e}", verbose=verbose)
            continue
    
    # Method 2: Check children of ax directly (inset axes might be stored as children)
    if not cbar_found:
        for child in list(ax.get_children()):
            if isinstance(child, plt.Axes) and child != ax:
                try:
                    yticklabels = [label.get_text() for label in child.get_yticklabels()]
                    patches = child.patches
                    # Check if this looks like an LD legend cbar
                    if len(yticklabels) > 0 and len(patches) > 3:  # LD legend has many patches
                        child.remove()
                        log.write(f" -Removed duplicate LD legend (cbar) from panel {panel_idx} (all panels share the same legend).", verbose=verbose)
                        cbar_found = True
                        break
                except Exception as e:
                    if verbose:
                        log.write(f" -Warning: Error checking child axis for cbar: {e}", verbose=verbose)
                    continue
    
    if not cbar_found and verbose:
        log.write(f" -Warning: Could not find cbar to remove from panel {panel_idx}. Total axes in figure: {len(fig.axes)}.", verbose=verbose)

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
