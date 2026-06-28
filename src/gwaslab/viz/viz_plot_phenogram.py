"""Phenogram plot entry point.
"""

from typing import Optional, Union, Dict, Any, List
import pandas as pd
import matplotlib.pyplot as plt
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_aux_phenogram import (
    DEFAULT_MARKER_COLORS,
    DEFAULT_MARKER_SHAPES,
    _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    _DEFAULT_MARKER_MAX_PER_ROW,
    _TELOMERE_LENGTH,
    _build_leads_with_y,
    _bump_phenogram_ylim_from_artists,
    _chr_number_label_position,
    _chr_panel_bottom_from_row_top,
    _compute_phenogram_row_bands,
    _compute_phenogram_xlim,
    _estimate_global_annotation_xmax,
    _label_eligible_mask,
    _load_phenogram_cytobands,
    _phenogram_global_xmin,
    _phenogram_panel_ceiling_y,
    _phenogram_panel_floor_y,
    _phenogram_provisional_xmax,
    _plot_chr_annotations,
    _plot_chr_body,
    _plot_phenogram_marker_legends,
    _prepare_marker_label_layout,
    _prepare_phenogram_leads_dict,
    _prepare_phenogram_annotation_column,
    _prepare_text_label_layout,
    _resolve_marker_maps,
    _resolve_phenogram_cytoband_path,
    _resolve_phenogram_chr_width,
    _summarize_phenogram_cytobands,
)


def _plot_phenogram(
    insumstats,
    snpid: str = "SNPID",
    chrom: str = "CHR",
    pos: str = "POS",
    p: str = "P",
    mlog10p: str = "MLOG10P",
    windowsizekb: int = 500,
    sig_level: float = 5e-8,
    cytoband_path: Optional[str] = None,
    build: str = "19",
    include_sex_chr: bool = False,
    only_anno_chr: bool = False,
    ncols: int = 11,
    figsize: tuple = (20, 48),
    dpi: int = 400,
    chr_width: float = 0.35,
    chr_x: float = 0.0,
    anno_x_pad: float = 0.18,
    chr_label_pad: float = 0.06,
    use_lead_extraction: bool = True,
    anno: Optional[Union[bool, str]] = None,
    anno_set: Optional[List[str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
    anno_source: str = "ensembl",
    anno_gtf_path: Optional[str] = None,
    anno_sig_level: float = 5e-8,
    anno_max_rows: int = 200,
    anno_kwargs: Optional[Dict[str, Any]] = None,
    anno_kwargs_single: Optional[Dict[str, Any]] = None,
    arrow_kwargs: Optional[Dict[str, Any]] = None,
    anno_group: Optional[str] = None,
    anno_shape: Optional[str] = None,
    anno_color: Optional[str] = None,
    marker_shapes: Optional[List[str]] = None,
    marker_colors: Optional[List[str]] = None,
    marker_shape_map: Optional[Dict[str, str]] = None,
    marker_color_map: Optional[Dict[str, str]] = None,
    marker_size: float = 81,
    marker_gap_pt: float = 3,
    marker_label_gap_pt: float = 2.75,
    marker_max_per_row: int = _DEFAULT_MARKER_MAX_PER_ROW,
    marker_row_gap_pt: float = 2.5,
    marker_label_align: str = "center",
    marker_fontsize: float = 11,
    marker_linewidth: float = 0.6,
    marker_label_bbox: bool = True,
    group_min_vertical_gap_pt: float = 3.5,
    group_marker_to_marker_gap_pt: float = 3,
    group_label_box_pad_pt: float = 1.5,
    show_legend: bool = True,
    legend_ncol: int = 6,
    legend_kwargs: Optional[Dict[str, Any]] = None,
    anno_style: str = "expand",
    repel_force: float = 0.5,
    anno_max_iter: int = 300,
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
    save: Union[bool, str] = False,
    save_kwargs: Optional[Dict[str, Any]] = None,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    verbose: bool = True,
    log: Log = Log(),
) -> plt.Figure:
    """Create a karyotype-style phenogram with cytobands and lead variants.

    Registry-aligned parameters are on ``Sumstats.plot_phenogram()``.

Returns
-------
matplotlib.figure.Figure
    The phenogram figure.
"""
    
    log.write("Start to create phenogram plot...", verbose=verbose)
    
    chr_width = _resolve_phenogram_chr_width(chr_width, ncols)
    
    marker_mode = anno_group is not None
    if legend_kwargs is None:
        legend_kwargs = {}
    if anno_set is None:
        anno_set = []
    if anno_kwargs is None:
        anno_kwargs = {}
    if anno_kwargs_single is None:
        anno_kwargs_single = {}
    if arrow_kwargs is None:
        arrow_kwargs = {}
    if anno_alias is None:
        anno_alias = {}
    
    # Extract dataframe if Sumstats object is passed
    if hasattr(insumstats, 'data') and not isinstance(insumstats, pd.DataFrame):
        insumstats = insumstats.data
    
    # Create working copy to preserve original (plotting functions should not modify input)
    sumstats = insumstats.copy()
    
    # Get leads from sumstats or use input table directly
    if use_lead_extraction:
        log.write(" -Extracting lead variants...", verbose=verbose)
        leads = _get_sig(
            insumstats_or_dataframe=sumstats,
            variant_id=snpid,
            chrom=chrom,
            pos=pos,
            p=p,
            mlog10p=mlog10p,
            windowsizekb=windowsizekb,
            sig_level=anno_sig_level,
            log=log,
            verbose=verbose
        )
    else:
        log.write(" -Using input table directly (no lead extraction)...", verbose=verbose)
        leads = sumstats.copy()
    
    if leads is None or len(leads) == 0:
        log.write(" -No variants to annotate. Plotting chromosomes without annotations.", verbose=verbose)
        leads = pd.DataFrame()
    else:
        log.write(" -Found {} variant row(s) to annotate.".format(len(leads)), verbose=verbose)
    
    if marker_mode:
        resolved_shape_map, resolved_color_map = _resolve_marker_maps(
            data=leads if len(leads) > 0 else sumstats,
            anno_shape=anno_shape,
            anno_color=anno_color,
            marker_shapes=marker_shapes,
            marker_colors=marker_colors,
            marker_shape_map=marker_shape_map,
            marker_color_map=marker_color_map,
            log=log,
            verbose=verbose,
        )
        log.write(
            " -Marker mode: {} shape categories, {} color categories.".format(
                len(resolved_shape_map), len(resolved_color_map)
            ),
            verbose=verbose,
        )
    else:
        resolved_shape_map, resolved_color_map = {}, {}
        log.write(" -Text annotation mode.", verbose=verbose)
        if anno is not None:
            if anno is True:
                log.write(" -Label source: Chr:POS", verbose=verbose)
            elif anno == "GENENAME":
                log.write(" -Label source: GENENAME", verbose=verbose)
            else:
                log.write(" -Label source: column '{}'".format(anno), verbose=verbose)
        else:
            log.write(" -Label source: none (bars only).", verbose=verbose)

    label_eligible = pd.Series(dtype=bool)
    if len(leads) > 0 and anno is not None:
        leads = _prepare_phenogram_annotation_column(
            leads=leads,
            anno=anno,
            snpid=snpid,
            chrom=chrom,
            pos=pos,
            build=build,
            anno_source=anno_source,
            anno_gtf_path=anno_gtf_path,
            log=log,
            verbose=verbose,
        )
        if not marker_mode:
            label_eligible = _label_eligible_mask(
                leads=leads,
                snpid=snpid,
                anno_set=anno_set,
                anno_max_rows=anno_max_rows,
                p=p,
                mlog10p=mlog10p,
                log=log,
                verbose=verbose,
            )
    
    cytoband_path = _resolve_phenogram_cytoband_path(
        build, cytoband_path, log=log, verbose=verbose
    )
    log.write(" -Loading cytoband data from: {}".format(cytoband_path), verbose=verbose)

    try:
        cytobands = _load_phenogram_cytobands(str(cytoband_path))
    except Exception as e:
        log.write(" -Error loading cytoband data: {}".format(e), verbose=verbose)
        raise ValueError("Could not load cytoband data from {}".format(cytoband_path))

    chr_list, chr_sizes, centromere_by_chr, cytobands_by_chr = (
        _summarize_phenogram_cytobands(cytobands, include_sex_chr=include_sex_chr)
    )
    leads_dict = _prepare_phenogram_leads_dict(
        leads=leads,
        chrom=chrom,
        pos=pos,
        snpid=snpid,
        chr_list=chr_list,
        marker_mode=marker_mode,
        anno=anno,
        label_eligible=label_eligible,
        anno_alias=anno_alias,
        anno_group=anno_group,
        anno_shape=anno_shape,
        anno_color=anno_color,
    )

    if only_anno_chr:
        annotated_chr_list = [chr_name for chr_name in chr_list if leads_dict.get(chr_name)]
        if annotated_chr_list:
            chr_list = annotated_chr_list
            log.write(
                " -only_anno_chr: plotting {} chromosome(s) with annotations.".format(
                    len(chr_list)
                ),
                verbose=verbose,
            )
        else:
            log.write(
                " -only_anno_chr: no chromosomes with annotations; plotting all chromosomes.",
                verbose=verbose,
            )

    n_chr = len(chr_list)

    log.write(" -Plotting {} chromosomes...".format(n_chr), verbose=verbose)

    # Set up style
    style = set_plot_style(
        plot="plot_phenogram",
        fig_kwargs=fig_kwargs if fig_kwargs is not None else {},
        save_kwargs=save_kwargs if save_kwargs is not None else {},
        save=save,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})

    # Create figure
    if "figsize" not in fig_kwargs:
        fig_kwargs["figsize"] = figsize
    if "dpi" not in fig_kwargs:
        fig_kwargs["dpi"] = dpi

    fig, axes = plt.subplots(nrows=1, ncols=ncols, figsize=fig_kwargs["figsize"],
                             dpi=fig_kwargs["dpi"])

    # Handle single axis case
    if ncols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    # Calculate offsets for each chromosome (row pitch from tallest chr per row).
    max_chr_size = max(chr_sizes[chr_name] for chr_name in chr_list)
    chr_number_fontsize = 10.0
    chr_number_pad_pt = 4.0
    chr_number_ref_x = chr_x + chr_width / 2
    global_xmin = _phenogram_global_xmin(chr_x, chr_width)
    fixed_global_xmax = _phenogram_provisional_xmax(chr_x, chr_width, anno_x_pad)
    provisional_ymax = 3.0
    for col in range(ncols):
        axes[col].set_xlim(global_xmin, fixed_global_xmax)
        axes[col].set_ylim(-0.05, provisional_ymax)

    (
        n_grid_rows,
        row_ideogram_top,
        row_band_bottom,
        telemere_full_length,
        label_above_ref,
        global_ymin,
        global_ymax,
    ) = _compute_phenogram_row_bands(
        chr_list,
        chr_sizes,
        max_chr_size,
        ncols,
        axes[0],
        ref_x=chr_number_ref_x,
        chr_number_fontsize=chr_number_fontsize,
        chr_number_pad_pt=chr_number_pad_pt,
    )
    
    # Pass 1: chromosome geometry only (cytobands + red bars in text mode).
    chr_specs: List[Dict[str, Any]] = []
    for i, chr_name in enumerate(chr_list):
        chr_cytobands = cytobands_by_chr[chr_name]

        visual_row = i // ncols
        chr_size = chr_sizes[chr_name]
        chr_extent = chr_size / max_chr_size
        chr_offset = _chr_panel_bottom_from_row_top(
            row_ideogram_top[visual_row],
            chr_extent,
            telemere_full_length,
        )

        chr_centromere_u, chr_centromere_l = centromere_by_chr[chr_name]

        spec = _plot_chr_body(
            axes[i % ncols],
            chr_size=chr_size,
            chr_centromere_u=chr_centromere_u,
            chr_centromere_l=chr_centromere_l,
            max_chr_size=max_chr_size,
            offset=chr_offset,
            chr_cytobands=chr_cytobands,
            chr_name=chr_name,
            leads=leads_dict.get(chr_name, []),
            chr_width=chr_width,
            chr_x=chr_x,
            anno_x_pad=anno_x_pad,
            log=log,
            verbose=verbose,
        )
        spec["visual_row"] = visual_row
        spec["chr_offset"] = chr_offset
        spec["row_ideogram_top"] = row_ideogram_top[visual_row]
        spec["panel_ceiling_y"] = _phenogram_panel_ceiling_y(
            visual_row,
            row_band_bottom,
            has_chr_above=(visual_row > 0),
        )
        spec["panel_floor_y"] = _phenogram_panel_floor_y(
            visual_row,
            row_band_bottom,
            has_chr_below=(visual_row < n_grid_rows - 1),
        )
        chr_specs.append(spec)

    global_ymax = global_ymax + chr_label_pad
    provisional_xmax = fixed_global_xmax

    for col in range(ncols):
        axes[col].set_xlim(global_xmin, provisional_xmax)
        axes[col].set_ylim(global_ymin, global_ymax)

    has_annotation_layout = (
        any(spec.get("leads_with_y") for spec in chr_specs)
        if marker_mode
        else any(
            item["lead"].get("label_text")
            for spec in chr_specs
            for item in (spec.get("leads_with_y") or [])
        )
    )

    if has_annotation_layout:
        fig.canvas.draw()
        estimate_renderer = fig.canvas.get_renderer()
        estimated_xmax = _estimate_global_annotation_xmax(
            chr_specs,
            global_xmin=global_xmin,
            provisional_xmax=provisional_xmax,
            marker_mode=marker_mode,
            marker_size=marker_size,
            marker_gap_pt=marker_gap_pt,
            marker_linewidth=marker_linewidth,
            marker_max_per_row=marker_max_per_row,
            marker_row_gap_pt=marker_row_gap_pt,
            marker_fontsize=marker_fontsize,
            marker_label_gap_pt=marker_label_gap_pt,
            group_label_box_pad_pt=group_label_box_pad_pt,
            anno_kwargs=anno_kwargs,
            renderer=estimate_renderer,
            anno=anno,
            anno_alias=anno_alias,
            marker_label_align=marker_label_align,
            anno_wrap=anno_wrap,
            anno_wrap_width_pt=anno_wrap_width_pt,
            anno_wrap_chars_per_line=anno_wrap_chars_per_line,
            anno_max_len=anno_max_len,
        )
    else:
        estimated_xmax = provisional_xmax

    global_xmax = _compute_phenogram_xlim(
        global_xmin,
        chr_x,
        chr_width,
        estimated_xmax,
        anno_x_pad=anno_x_pad,
    )

    for col in range(ncols):
        axes[col].set_xlim(global_xmin, global_xmax)

    if has_annotation_layout:
        fig.canvas.draw()
        layout_renderer = fig.canvas.get_renderer()

        if marker_mode:
            for spec in chr_specs:
                layout_ymax = _prepare_marker_label_layout(
                    spec,
                    marker_size=marker_size,
                    marker_gap_pt=marker_gap_pt,
                    marker_label_gap_pt=marker_label_gap_pt,
                    marker_max_per_row=marker_max_per_row,
                    marker_row_gap_pt=marker_row_gap_pt,
                    marker_fontsize=marker_fontsize,
                    marker_linewidth=marker_linewidth,
                    group_min_vertical_gap_pt=group_min_vertical_gap_pt,
                    group_marker_to_marker_gap_pt=group_marker_to_marker_gap_pt,
                    group_label_box_pad_pt=group_label_box_pad_pt,
                    anno_max_iter=anno_max_iter,
                    repel_force=repel_force,
                    marker_label_align=marker_label_align,
                    anno=anno,
                    anno_alias=anno_alias,
                    anno_wrap=anno_wrap,
                    anno_wrap_width_pt=anno_wrap_width_pt,
                    anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                    anno_max_len=anno_max_len,
                    renderer=layout_renderer,
                )
                global_ymax = max(global_ymax, layout_ymax)
        else:
            for spec in chr_specs:
                layout_ymax = _prepare_text_label_layout(
                    spec,
                    anno_kwargs=anno_kwargs,
                    group_min_vertical_gap_pt=group_min_vertical_gap_pt,
                    group_marker_to_marker_gap_pt=group_marker_to_marker_gap_pt,
                    group_label_box_pad_pt=group_label_box_pad_pt,
                    marker_label_gap_pt=marker_label_gap_pt,
                    repel_force=repel_force,
                    anno_max_iter=anno_max_iter,
                    renderer=layout_renderer,
                    global_xmax=global_xmax,
                    anno_wrap=anno_wrap,
                    anno_wrap_width_pt=anno_wrap_width_pt,
                    anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                    anno_max_len=anno_max_len,
                )
                global_ymax = max(global_ymax, layout_ymax)

        for col in range(ncols):
            axes[col].set_ylim(global_ymin, global_ymax)

        fig.canvas.draw()

    if show_legend and marker_mode:
        same_legend_column = (
            anno_shape is not None
            and anno_color is not None
            and anno_shape == anno_color
        )
        legend_top_margin = 0.04 if same_legend_column else 0.075
        fig.tight_layout(rect=[0, 0, 1, 1 - legend_top_margin])
    else:
        fig.tight_layout()

    if has_annotation_layout:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    else:
        renderer = None

    # Pass 2: annotations on frozen axis limits.
    for spec in chr_specs:
        _plot_chr_annotations(
            spec,
            anno_kwargs=anno_kwargs,
            anno_kwargs_single=anno_kwargs_single,
            arrow_kwargs=arrow_kwargs,
            marker_shape_map=resolved_shape_map,
            marker_color_map=resolved_color_map,
            marker_size=marker_size,
            marker_gap_pt=marker_gap_pt,
            marker_label_gap_pt=marker_label_gap_pt,
            marker_max_per_row=marker_max_per_row,
            marker_row_gap_pt=marker_row_gap_pt,
            marker_label_align=marker_label_align,
            marker_fontsize=marker_fontsize,
            marker_linewidth=marker_linewidth,
            marker_label_bbox=marker_label_bbox,
            group_min_vertical_gap_pt=group_min_vertical_gap_pt,
            group_marker_to_marker_gap_pt=group_marker_to_marker_gap_pt,
            group_label_box_pad_pt=group_label_box_pad_pt,
            marker_mode=marker_mode,
            repel_force=repel_force,
            anno_max_iter=anno_max_iter,
            renderer=renderer,
            anno=anno,
            anno_alias=anno_alias,
            anno_wrap=anno_wrap,
            anno_wrap_width_pt=anno_wrap_width_pt,
            anno_wrap_chars_per_line=anno_wrap_chars_per_line,
            anno_max_len=anno_max_len,
            global_xmax=global_xmax,
        )

    measured_ymax, measured_ymin = _bump_phenogram_ylim_from_artists(
        fig, axes, ncols, global_ymax, global_ymin
    )
    if measured_ymax > global_ymax or measured_ymin < global_ymin:
        global_ymax = max(global_ymax, measured_ymax)
        global_ymin = min(global_ymin, measured_ymin)
        for col in range(ncols):
            axes[col].set_ylim(global_ymin, global_ymax)
        fig.canvas.draw()

    if not marker_mode:
        measured_xmax = global_xmax
        for col in range(ncols):
            ax_xmax = getattr(axes[col], "_phenogram_xmax", global_xmax)
            measured_xmax = max(measured_xmax, ax_xmax)
        if measured_xmax > global_xmax:
            global_xmax = measured_xmax
            for col in range(ncols):
                axes[col].set_xlim(global_xmin, global_xmax)
            fig.canvas.draw()

    # Axis cosmetics (labels drawn after layout so transforms match the final axes).
    for i in range(n_chr):
        col_ax = axes[i % ncols]
        col_ax.spines['top'].set_visible(False)
        col_ax.spines['right'].set_visible(False)
        col_ax.spines['bottom'].set_visible(False)
        col_ax.spines['left'].set_visible(False)

    # Hide unused axes
    for i in range(n_chr, len(axes)):
        axes[i].set_visible(False)

    if show_legend and marker_mode:
        _plot_phenogram_marker_legends(
            fig,
            axes,
            ncols,
            n_chr,
            resolved_shape_map,
            resolved_color_map,
            legend_ncol,
            legend_kwargs,
            anno_shape=anno_shape,
            anno_color=anno_color,
        )

    for i, spec in enumerate(chr_specs):
        col_ax = spec["ax"]
        chr_name = spec["chr_name"]
        chr_num = chr_name[3:] if chr_name.startswith("chr") else chr_name
        body_top_y = spec["body_top_y"]
        label_y, label_va = _chr_number_label_position(
            col_ax,
            body_top_y,
            chr_number_ref_x,
            fontsize=chr_number_fontsize,
            pad_pt=chr_number_pad_pt,
        )
        col_ax.text(
            chr_number_ref_x,
            label_y,
            chr_num,
            ha="center",
            va=label_va,
            fontsize=chr_number_fontsize,
            zorder=300,
        )
    
    # Save figure
    save_figure(fig=fig, save=save, keyword="phenogram", save_kwargs=save_kwargs, 
                log=log, verbose=verbose)
    
    log.write("Finished creating phenogram plot.", verbose=verbose)
    
    return fig


