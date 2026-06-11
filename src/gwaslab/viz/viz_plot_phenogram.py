"""Phenogram plot entry point."""

from pathlib import Path
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
    _lead_anno_field,
    _phenogram_global_xmin,
    _phenogram_panel_ceiling_y,
    _phenogram_panel_floor_y,
    _phenogram_provisional_xmax,
    _plot_chr_annotations,
    _plot_chr_body,
    _plot_phenogram_marker_legends,
    _prepare_marker_label_layout,
    _prepare_phenogram_annotation_column,
    _prepare_text_label_layout,
    _resolve_marker_maps,
    _resolve_phenogram_chr_width,
    _resolve_phenogram_label,
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
    """
    Create a phenogram plot showing cytobands for each chromosome with lead SNPs annotated.
    
    This function generates a karyotype-like visualization (phenogram) displaying:
    - Cytogenetic bands for each chromosome
    - Lead SNPs from the sumstats data
    - Optional annotations for SNP IDs
    
    Parameters
    ----------
    mysumstats : Sumstats
        Sumstats object containing the GWAS data
    snpid : str, default="SNPID"
        Column name for SNP identifier
    chrom : str, default="CHR"
        Column name for chromosome
    pos : str, default="POS"
        Column name for position
    p : str, default="P"
        Column name for p-value
    mlog10p : str, default="MLOG10P"
        Column name for -log10(p-value)
    windowsizekb : int, default=500
        Window size in kilobases for lead variant identification
    sig_level : float, default=5e-8
        Significance threshold for variant selection
    cytoband_path : str, optional
        Path to cytoband file. If None, uses default hg19 cytoband file
    build : str, default="19"
        Genome build version ("19" for hg19, "38" for hg38)
    ncols : int, default=11
        Number of columns for arranging chromosomes
    figsize : tuple, default=(20, 48)
        Figure size in inches
    dpi : int, default=400
        Resolution of the figure
    chr_width : float, default=0.35
        Chromosome width in data coordinates
    chr_x : float, default=0.0
        Left x boundary of the chromosome in data coordinates
    anno_x_pad : float, default=0.18
        Extra horizontal gap between chromosome right edge and annotation area (data coords)
    chr_label_pad : float, default=0.06
        Vertical space reserved below each chromosome for the chromosome number label
    use_lead_extraction : bool, default=True
        If True, extract lead variants via ``_get_sig()``. If False, use every input row.
    anno : None, bool, or str, default=None
        Text label source (same as MQQ ``plot_mqq``):

        - ``None``: draw lead bars only, no text labels
        - ``True``: label as ``Chr{CHR}:{POS}``
        - ``"GENENAME"``: nearest gene via reference annotation
        - column name (e.g. ``"SNPID"``): label from that column
    anno_set : list of str, optional
        Restrict text labels to these SNPID values. All extracted leads still get bars.
    anno_alias : dict, optional
        Map SNPID to custom label text.
    anno_source : str, default="ensembl"
        Gene annotation source when ``anno="GENENAME"``.
    anno_gtf_path : str, optional
        Custom GTF path for ``anno="GENENAME"``.
    anno_sig_level : float, default=5e-8
        Significance threshold for lead extraction.
    anno_max_rows : int, default=200
        Maximum number of text labels (top by significance).
    anno_kwargs : dict, optional
        Default matplotlib text/connector styling for all labels (``fontsize``,
        ``fontstyle``, ``arrow_shaft``, ``arrow_pad``, ``arrow_shrink_b``, etc.).
    anno_kwargs_single : dict, optional
        Per-SNPID overrides merged on top of ``anno_kwargs`` (MQQ-compatible).
    arrow_kwargs : dict, optional
        Matplotlib arrow/line styling for text-mode connectors (``color``,
        ``lw``, ``arrowstyle``, etc.), merged into connector drawing.
    anno_group : str, optional
        Column for marker-mode grouping. If set, enables marker annotation mode.
    anno_shape : str, optional
        Column for marker shape; unique values auto-mapped to ``marker_shapes`` pool.
    anno_color : str, optional
        Column for marker color; unique values auto-mapped to ``marker_colors`` pool.
    marker_shapes : list of str, optional
        Pool of matplotlib marker codes for auto shape mapping.
    marker_colors : list of str, optional
        Pool of colors for auto color mapping.
    marker_shape_map : dict, optional
        Explicit override for shape value -> marker code.
    marker_color_map : dict, optional
        Explicit override for color value -> color.
    marker_size : float, default=81
        Scatter marker size for annotation markers.
    marker_gap_pt : float, default=3
        Fixed edge-to-edge gap between markers within a row, in points. Center spacing
        is ``2 * marker_radius + marker_gap_pt``.
    marker_label_gap_pt : float, default=2.75
        Vertical gap between the bottom marker row and the group text label, in
        points. Also sets the default spacing between adjacent annotation units
        (marker cluster + label).
    marker_max_per_row : int, default=4
        Maximum markers on one row before wrapping to the next row within a group.
    marker_row_gap_pt : float, default=2.5
        Vertical gap between wrapped marker rows inside a group, in points.
    marker_label_align : str, default="center"
        Deprecated for group label placement: group labels are always
        left-aligned at the fixed annotation column (``anno_col_x``).
        Retained for API compatibility when ``anno_col_x`` is not passed.
    marker_fontsize : float, default=11
        Font size for marker-mode text labels.
    marker_linewidth : float, default=0.6
        Edge linewidth for marker scatter points.
    marker_label_bbox : bool, default=True
        If True, draw marker-mode labels with a white text outline (no filled box).
    group_min_vertical_gap_pt : float, default=3.5
        Optional floor for vertical spacing between annotation units, in points
        (default inter-unit gap is ``marker_label_gap_pt``).
    group_marker_to_marker_gap_pt : float, default=3
        Optional floor for vertical spacing between adjacent units, in points.
    group_label_box_pad_pt : float, default=1.5
        Extra padding around label text for overlap calculations, in points.
    show_legend : bool, default=True
        If True, draw two single-line figure legends above the plot in marker
        mode: one row ``Shape …`` and one row ``Color …``.
    legend_ncol : int, default=6
        Deprecated for marker legends (each row is kept on one line). Retained
        for API compatibility.
    legend_kwargs : dict, optional
        Extra keyword arguments passed to ``fig.legend()``.
    anno_style : str, default="expand"
        Annotation layout style for lead SNP labels
    repel_force : float, default=0.5
        Scales inter-unit vertical gap (based on ``marker_label_gap_pt``) for
        crowded loci. Display-space spreading uses ``anno_max_iter`` iterations.
    anno_max_iter : int, default=300
        Maximum iterations for annotation repulsion / marker group spreading
    anno_wrap : bool, default=True
        If True, word-wrap long text-mode labels instead of truncating.
    anno_wrap_width_pt : float, optional
        Maximum label line width in points when ``anno_wrap_chars_per_line`` is
        None. If both are None, width is derived from the annotation column.
    anno_wrap_chars_per_line : int, optional, default=9
        Maximum characters per line for text-mode labels. Set to ``None`` to
        use width-based wrapping only.
    anno_max_len : int, optional
        Optional hard character limit applied after wrapping (text mode) or to
        group labels (marker mode). When ``anno_wrap=False``, defaults to 20.
    save : bool or str, default=False
        If True or str, save the figure
    save_kwargs : dict, optional
        Additional keyword arguments for saving
    fig_kwargs : dict, optional
        Additional keyword arguments for figure creation
    verbose : bool, default=True
        If True, print progress messages
    log : Log, default=Log()
        Logging object
        
    Returns
    -------
    matplotlib.figure.Figure
        The created matplotlib figure object
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
    
    # Load cytoband data
    if cytoband_path is None:
        # Use default cytoband file path
        data_dir = Path(__file__).parent.parent / "data" / "cytoband"
        if build == "19":
            cytoband_path = data_dir / "cytoBand_hg19.txt.gz"
        elif build == "38":
            cytoband_path = data_dir / "cytoBand_hg38.txt.gz"
        else:
            # Default to hg19
            cytoband_path = data_dir / "cytoBand_hg19.txt.gz"
            log.write(" -Unknown build '{}', using hg19 cytoband data.".format(build), verbose=verbose)
    
    log.write(" -Loading cytoband data from: {}".format(cytoband_path), verbose=verbose)
    
    try:
        cytobands = pd.read_csv(cytoband_path, sep="\s+", header=None, compression='gzip')
        cytobands.columns = ["CHR", "START", "END", "ARM", "STAIN"]
    except Exception as e:
        log.write(" -Error loading cytoband data: {}".format(e), verbose=verbose)
        raise ValueError("Could not load cytoband data from {}".format(cytoband_path))
    
    # Color dictionary for cytobands
    color_dict = {
        "gpos100": (100/255, 100/255, 100/255),
        "gpos": (100/255, 100/255, 100/255),
        "gpos75": (110/255, 110/255, 110/255),
        "gpos66": (130/255, 130/255, 130/255),
        "gpos50": (160/255, 160/255, 160/255),
        "gpos33": (180/255, 180/255, 180/255),
        "gpos25": (200/255, 200/255, 200/255),
        "gvar": (220/255, 220/255, 220/255),
        "gneg": (255/255, 255/255, 255/255),
        "acen": (217/255, 217/255, 217/255),
        "stalk": (100/255, 127/255, 164/255)
    }
    cytobands["COLOR"] = cytobands["STAIN"].map(color_dict)
    
    # Get chromosome sizes from cytoband data
    chr_sizes = {}
    for chr_name in cytobands["CHR"].unique():
        chr_bands = cytobands[cytobands["CHR"] == chr_name]
        chr_sizes[chr_name] = chr_bands["END"].max()
    
    # Get numeric chromosome list (1-22, X, Y, MT)
    chr_list = []
    for chr_name in sorted(cytobands["CHR"].unique()):
        # Extract number from chr name (e.g., "chr1" -> 1)
        if chr_name.startswith("chr"):
            chr_num_str = chr_name[3:]
            try:
                if chr_num_str in ["X", "Y", "MT"]:
                    chr_list.append(chr_name)
                else:
                    chr_num = int(chr_num_str)
                    if 1 <= chr_num <= 22:
                        chr_list.append(chr_name)
            except ValueError:
                continue
    
    # Sort chromosomes numerically (1-22, then X, Y, MT)
    def chr_sort_key(chr_name):
        if not chr_name.startswith("chr"):
            return (999, chr_name)
        chr_num_str = chr_name[3:]
        if chr_num_str.isdigit():
            return (int(chr_num_str), chr_name)
        elif chr_num_str == "X":
            return (23, chr_name)
        elif chr_num_str == "Y":
            return (24, chr_name)
        elif chr_num_str == "MT":
            return (25, chr_name)
        else:
            return (999, chr_name)
    
    chr_list = sorted(chr_list, key=chr_sort_key)
    
    # Limit to autosomes (1-22) for now - can be extended later
    chr_list = [c for c in chr_list if c.startswith("chr") and c[3:].isdigit() and 1 <= int(c[3:]) <= 22]
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
    max_chr_size = max(chr_sizes.values())
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
    
    # Prepare leads data for annotation
    leads_dict = {}
    if len(leads) > 0:
        # Convert chromosome format in leads to match cytoband format
        leads_copy = leads.copy()
        if chrom in leads_copy.columns and pos in leads_copy.columns:
            # Convert numeric chromosomes to chr format
            for idx, row in leads_copy.iterrows():
                chr_val = row[chrom]
                pos_val = row[pos]
                if pd.notna(chr_val) and pd.notna(pos_val):
                    # Convert to chr format
                    if isinstance(chr_val, (int, float)):
                        chr_str = "chr{}".format(int(chr_val))
                    else:
                        chr_str = str(chr_val)
                        if not chr_str.startswith("chr"):
                            # Try to extract number
                            try:
                                chr_num = int(chr_str)
                                chr_str = "chr{}".format(chr_num)
                            except (ValueError, TypeError):
                                # If it's already a string like "X", "Y", etc., add "chr" prefix
                                chr_str = "chr{}".format(chr_str)
                    
                    # Only include if it's a valid chromosome (1-22, X, Y, MT)
                    if chr_str in chr_list or (chr_str.startswith("chr") and 
                                               (chr_str[3:].isdigit() or chr_str[3:] in ["X", "Y", "MT"])):
                        chrom_num = int(chr_val) if isinstance(chr_val, (int, float)) else chr_str.replace("chr", "")
                        annotation_val = None
                        if "Annotation" in leads_copy.columns and pd.notna(row.get("Annotation")):
                            annotation_val = row["Annotation"]
                        lead_entry = {
                            'pos': int(pos_val) if pd.notna(pos_val) else None,
                            'chrom': chrom_num,
                            'snpid': row[snpid] if snpid in leads_copy.columns and pd.notna(row[snpid]) else None,
                            'index': idx,
                            'annotation': annotation_val,
                            'anno_group': _lead_anno_field(row, anno_group, leads_copy.columns),
                            'anno_shape': _lead_anno_field(row, anno_shape, leads_copy.columns),
                            'anno_color': _lead_anno_field(row, anno_color, leads_copy.columns),
                        }
                        if not marker_mode and anno is not None and label_eligible.get(idx, False):
                            lead_entry['label_text'] = _resolve_phenogram_label(
                                lead_entry, anno, anno_alias
                            )
                        else:
                            lead_entry['label_text'] = None
                        leads_dict.setdefault(chr_str, []).append(lead_entry)
    
    # Pass 1: chromosome geometry only (cytobands + red bars in text mode).
    chr_specs: List[Dict[str, Any]] = []
    for i, chr_name in enumerate(chr_list):
        chr_cytobands = cytobands.loc[cytobands["CHR"] == chr_name, :].copy()

        visual_row = i // ncols
        chr_size = chr_sizes[chr_name]
        chr_extent = chr_size / max_chr_size
        chr_offset = _chr_panel_bottom_from_row_top(
            row_ideogram_top[visual_row],
            chr_extent,
            telemere_full_length,
        )

        acen_bands = chr_cytobands.loc[chr_cytobands["STAIN"] == "acen", :]
        if len(acen_bands) > 0:
            chr_centromere_u = acen_bands["START"].min()
            chr_centromere_l = acen_bands["END"].max()
        else:
            chr_centromere_u = chr_size * 0.44
            chr_centromere_l = chr_size * 0.46

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

    fig.canvas.draw()

    global_xmax = _estimate_global_annotation_xmax(
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
        renderer=fig.canvas.get_renderer(),
        anno=anno,
        anno_alias=anno_alias,
        marker_label_align=marker_label_align,
        anno_wrap=anno_wrap,
        anno_wrap_width_pt=anno_wrap_width_pt,
        anno_wrap_chars_per_line=anno_wrap_chars_per_line,
        anno_max_len=anno_max_len,
    )
    global_xmax = _compute_phenogram_xlim(
        global_xmin,
        chr_x,
        chr_width,
        global_xmax,
        anno_x_pad=anno_x_pad,
    )

    for col in range(ncols):
        axes[col].set_xlim(global_xmin, global_xmax)

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

        for col in range(ncols):
            axes[col].set_ylim(global_ymin, global_ymax)

        fig.canvas.draw()
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
        legend_top_margin = 0.075
        fig.tight_layout(rect=[0, 0, 1, 1 - legend_top_margin])
    else:
        fig.tight_layout()

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

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


