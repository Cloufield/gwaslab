import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.collections import PatchCollection
from pathlib import Path
from typing import Optional, Union, Dict, Any, List, Tuple
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.bd.bd_common_data import get_chr_to_number, get_number_to_chr
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_aux_reposition_text import adjust_text_position


DEFAULT_MARKER_SHAPES: List[str] = ["o", "D", "^", "s", "v", "P", "X", "*"]
DEFAULT_MARKER_COLORS: List[str] = [
    "black",
    "#78A6D6",
    "#B7354A",
    "#E6A23A",
    "#70AD47",
    "#7F7F7F",
    "#9467BD",
    "#8C564B",
]

_TEXT_SKIP_KEYS = frozenset({
    "arrow_pad", "arrow_shaft", "arrow_shrink_b", "arrowprops",
})


def _points_to_data_delta(ax, x, y, dx_points=0.0, dy_points=0.0):
    """Convert a display offset in matplotlib points to a data-coordinate delta at (x, y)."""
    scale = ax.figure.dpi / 72.0
    x_disp, y_disp = ax.transData.transform((x, y))
    x_data, y_data = ax.transData.inverted().transform(
        (x_disp + dx_points * scale, y_disp + dy_points * scale)
    )
    return x_data - x, y_data - y


def _resolve_anno_arrow_pt(
    explicit: Optional[float],
    kwargs_key: str,
    annotation_kwargs: Dict[str, Any],
    default: float,
) -> float:
    """Resolve arrow spacing in points: explicit param > annotation_kwargs > default."""
    if explicit is not None:
        return float(explicit)
    if kwargs_key in annotation_kwargs:
        return float(annotation_kwargs[kwargs_key])
    return float(default)


def _lead_anno_field(row: pd.Series, col_name: Optional[str], columns: pd.Index) -> Any:
    if col_name is None or col_name not in columns:
        return None
    val = row[col_name]
    if pd.isna(val):
        return None
    return val


def _get_group_key(lead: Dict[str, Any]) -> str:
    if lead.get("anno_group") is not None:
        return str(lead["anno_group"])
    if lead.get("anno_text") is not None:
        return str(lead["anno_text"])
    if lead.get("snpid") is not None:
        return str(lead["snpid"])
    return str(lead.get("index", "unknown"))


def _get_text_mode_label(lead: Dict[str, Any]) -> str:
    if lead.get("anno_text") is not None:
        return str(lead["anno_text"])
    if lead.get("snpid") is not None:
        return str(lead["snpid"])
    return ""


def _get_marker_group_label(group_key: str, group_items: List[Dict[str, Any]]) -> str:
    first_lead = group_items[0]["lead"]
    if first_lead.get("anno_text") is not None:
        return str(first_lead["anno_text"])
    if first_lead.get("anno_group") is not None:
        return str(first_lead["anno_group"])
    if first_lead.get("snpid") is not None:
        return str(first_lead["snpid"])
    return str(group_key)


def _resolve_marker_maps(
    data: pd.DataFrame,
    anno_shape: Optional[str],
    anno_color: Optional[str],
    marker_shapes: Optional[List[str]] = None,
    marker_colors: Optional[List[str]] = None,
    marker_shape_map: Optional[Dict[str, str]] = None,
    marker_color_map: Optional[Dict[str, str]] = None,
    log: Log = Log(),
    verbose: bool = True,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Build value -> marker shape/color maps from data unique values or user overrides."""
    shapes = list(marker_shapes) if marker_shapes is not None else list(DEFAULT_MARKER_SHAPES)
    colors = list(marker_colors) if marker_colors is not None else list(DEFAULT_MARKER_COLORS)

    if marker_shape_map is not None:
        resolved_shape_map = dict(marker_shape_map)
    elif anno_shape is not None and anno_shape in data.columns:
        shape_values = sorted([v for v in data[anno_shape].dropna().unique()], key=str)
        if len(shape_values) > len(shapes):
            log.write(
                " -Warning: more unique marker shape values than available marker shapes; shapes will be reused.",
                verbose=verbose,
            )
        resolved_shape_map = {
            value: shapes[i % len(shapes)]
            for i, value in enumerate(shape_values)
        }
    else:
        resolved_shape_map = {}

    if marker_color_map is not None:
        resolved_color_map = dict(marker_color_map)
    elif anno_color is not None and anno_color in data.columns:
        color_values = sorted([v for v in data[anno_color].dropna().unique()], key=str)
        if len(color_values) > len(colors):
            log.write(
                " -Warning: more unique marker color values than available marker colors; colors will be reused.",
                verbose=verbose,
            )
        resolved_color_map = {
            value: colors[i % len(colors)]
            for i, value in enumerate(color_values)
        }
    else:
        resolved_color_map = {}

    return resolved_shape_map, resolved_color_map


def _get_marker_style(
    lead: Dict[str, Any],
    marker_shape_map: Dict[str, str],
    marker_color_map: Dict[str, str],
) -> Tuple[str, str]:
    shape_value = lead.get("anno_shape")
    color_value = lead.get("anno_color")
    marker = marker_shape_map.get(shape_value, "o") if shape_value is not None else "o"
    color = marker_color_map.get(color_value, "black") if color_value is not None else "black"
    return marker, color


def _build_text_mode_kwargs(annotation_kwargs: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    text_kwargs = {
        "fontsize": 8,
        "ha": "left",
        "va": "center",
        "fontweight": "normal",
    }
    if annotation_kwargs:
        clean = {k: v for k, v in annotation_kwargs.items() if k not in _TEXT_SKIP_KEYS}
        text_kwargs.update(clean)
    text_kwargs["ha"] = "left"
    text_kwargs["va"] = "center"
    return text_kwargs


def _data_point_below(ax, x: float, y: float, dy_points: float) -> Tuple[float, float]:
    """Return data coordinates dy_points below (x, y) in display space."""
    scale = ax.figure.dpi / 72.0
    x_disp, y_disp = ax.transData.transform((x, y))
    return ax.transData.inverted().transform((x_disp, y_disp + dy_points * scale))


def _marker_radius_pt(marker_size: float) -> float:
    """Matplotlib scatter ``s`` is area in points²; return marker radius in points."""
    return np.sqrt(max(marker_size, 1.0) / np.pi)


def _marker_block_extents_pt(
    marker_size: float,
    marker_label_gap_pt: float,
    marker_fontsize: float,
    group_label_box_pad_pt: float,
    bbox_pad: float = 0.15,
) -> Tuple[float, float]:
    """
    Return (above_pt, below_pt) vertical extent of a marker-mode block from marker_y.

    Marker is centered at marker_y; label top sits below the marker row bottom edge.
    """
    marker_radius_pt = _marker_radius_pt(marker_size)
    text_height_pt = marker_fontsize * 1.25 + bbox_pad * 2.0
    above_pt = marker_radius_pt
    below_pt = (
        marker_radius_pt
        + marker_label_gap_pt
        + text_height_pt
        + group_label_box_pad_pt * 2.0
    )
    return above_pt, below_pt


def _marker_label_position(
    ax,
    label_x: float,
    marker_y: float,
    marker_size: float,
    marker_label_gap_pt: float,
) -> Tuple[float, float]:
    """Return (label_x, label_y) with va='top' anchor below the marker row bottom edge."""
    scale = ax.figure.dpi / 72.0
    marker_radius_pt = _marker_radius_pt(marker_size)
    x_disp, y_disp = ax.transData.transform((label_x, marker_y))
    label_top_disp = y_disp + (marker_radius_pt + marker_label_gap_pt) * scale
    return ax.transData.inverted().transform((x_disp, label_top_disp))


def _marker_row_x_positions(
    ax,
    start_x: float,
    start_y: float,
    n_markers: int,
    gap_pt: float,
) -> List[float]:
    """Return evenly spaced marker center x positions using fixed display-pixel gaps."""
    if n_markers <= 0:
        return []
    scale = ax.figure.dpi / 72.0
    gap_px = gap_pt * scale
    x0_disp, y_disp = ax.transData.transform((start_x, start_y))
    xs_data: List[float] = []
    for j in range(n_markers):
        x_disp = x0_disp + j * gap_px
        x_data, _ = ax.transData.inverted().transform((x_disp, y_disp))
        xs_data.append(x_data)
    return xs_data


def _marker_block_bottom_display(
    ax,
    ref_x: float,
    marker_y: float,
    below_pt: float,
) -> float:
    """Return display-y of the bottom edge of a group annotation block."""
    _, y_disp = ax.transData.transform((ref_x, marker_y))
    scale = ax.figure.dpi / 72.0
    return y_disp + below_pt * scale


def _marker_block_height_pt(
    marker_size: float,
    marker_label_gap_pt: float,
    marker_fontsize: float,
    group_label_box_pad_pt: float,
) -> float:
    """Approximate total vertical extent of a marker-mode annotation block in points."""
    above_pt, below_pt = _marker_block_extents_pt(
        marker_size, marker_label_gap_pt, marker_fontsize, group_label_box_pad_pt
    )
    return above_pt + below_pt


def _spread_group_blocks_display(
    ax,
    y_data_list: List[float],
    ref_x: float,
    above_pt: float,
    below_pt: float,
    gap_pt: float = 4.0,
    max_iter: int = 300,
) -> np.ndarray:
    """
    Spread group marker_y positions in display space so full blocks do not overlap.

    Each block spans from (marker_y - above) to (marker_y + below) in display
    coordinates. Groups are processed top-to-bottom; when a block would overlap
    the previous block (including its label below the markers), it is pushed down
    sequentially to leave enough room.
    """
    ys_data = np.array(y_data_list, dtype=float)
    if len(ys_data) <= 1:
        return ys_data.copy()

    scale = ax.figure.dpi / 72.0
    above_px = above_pt * scale
    below_px = below_pt * scale
    gap_px = gap_pt * scale

    ys_disp = np.array([ax.transData.transform((ref_x, y))[1] for y in ys_data])
    order = np.argsort(ys_disp)
    ys_sorted = ys_disp[order].copy()

    for _ in range(max_iter):
        changed = False
        for i in range(1, len(ys_sorted)):
            prev_bottom = ys_sorted[i - 1] + below_px
            min_marker_disp = prev_bottom + gap_px + above_px
            if ys_sorted[i] < min_marker_disp:
                ys_sorted[i] = min_marker_disp
                changed = True
        if not changed:
            break

    result = ys_data.copy()
    for rank, orig_idx in enumerate(order):
        _, y_data = ax.transData.inverted().transform((ref_x, ys_sorted[rank]))
        result[orig_idx] = y_data
    return result


def _spread_group_positions_display(
    ax,
    y_data_list: List[float],
    ref_x: float,
    min_gap_pt: float = 18.0,
    block_height_pt: Optional[float] = None,
    max_iter: int = 300,
) -> np.ndarray:
    """Legacy wrapper: spread by total block height with symmetric padding."""
    block_pt = block_height_pt if block_height_pt is not None else min_gap_pt
    half = block_pt / 2.0
    return _spread_group_blocks_display(
        ax,
        y_data_list,
        ref_x=ref_x,
        above_pt=half,
        below_pt=half,
        gap_pt=min_gap_pt,
        max_iter=max_iter,
    )


def _build_marker_text_kwargs(
    annotation_kwargs: Optional[Dict[str, Any]],
    marker_fontsize: float = 10,
    marker_label_bbox: bool = True,
) -> Dict[str, Any]:
    text_kwargs: Dict[str, Any] = {
        "fontsize": marker_fontsize,
        "ha": "center",
        "va": "top",
        "fontstyle": "italic",
        "fontweight": "bold",
    }
    if marker_label_bbox:
        text_kwargs["path_effects"] = [
            pe.withStroke(linewidth=2.5, foreground="white"),
        ]
    if annotation_kwargs:
        clean = {k: v for k, v in annotation_kwargs.items() if k not in _TEXT_SKIP_KEYS}
        text_kwargs.update(clean)
    for key in ("arrow_pad", "arrow_shaft", "arrow_shrink_b", "arrowprops"):
        text_kwargs.pop(key, None)
    return text_kwargs


def _make_marker_legend_handles(
    marker_shape_map: Dict[str, str],
    marker_color_map: Dict[str, str],
) -> List[Line2D]:
    handles: List[Line2D] = []
    for label, marker in marker_shape_map.items():
        handles.append(
            Line2D(
                [0], [0],
                marker=marker,
                linestyle="None",
                markerfacecolor="white",
                markeredgecolor="black",
                markersize=6,
                label=str(label),
            )
        )
    for label, color in marker_color_map.items():
        handles.append(
            Line2D(
                [0], [0],
                marker="s",
                linestyle="None",
                markerfacecolor=color,
                markeredgecolor="black",
                markersize=6,
                label=str(label),
            )
        )
    return handles


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
    figsize: tuple = (20, 40),
    dpi: int = 400,
    chr_width: float = 0.35,
    chr_x: float = 0.0,
    anno_x_pad: float = 0.14,
    use_lead_extraction: bool = True,
    annotate_snps: bool = True,
    annotation_kwargs: Optional[Dict[str, Any]] = None,
    anno_arrow_shaft: Optional[float] = 18,
    anno_arrow_pad: Optional[float] = 10,
    anno_arrow_shrink_b: Optional[float] = 4,
    anno_text: Optional[str] = None,
    anno_group: Optional[str] = None,
    anno_shape: Optional[str] = None,
    anno_color: Optional[str] = None,
    marker_shapes: Optional[List[str]] = None,
    marker_colors: Optional[List[str]] = None,
    marker_shape_map: Optional[Dict[str, str]] = None,
    marker_color_map: Optional[Dict[str, str]] = None,
    marker_size: float = 42,
    marker_gap_pt: float = 14,
    marker_label_gap_pt: float = 6,
    marker_label_align: str = "center",
    marker_fontsize: float = 10,
    marker_linewidth: float = 0.6,
    marker_label_bbox: bool = True,
    group_min_vertical_gap_pt: float = 18,
    group_marker_to_marker_gap_pt: float = 16,
    group_label_box_pad_pt: float = 1.5,
    show_legend: bool = True,
    legend_ncol: int = 6,
    legend_kwargs: Optional[Dict[str, Any]] = None,
    anno_style: str = "expand",
    repel_force: float = 0.5,
    anno_max_iter: int = 300,
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
    figsize : tuple, default=(20, 40)
        Figure size in inches
    dpi : int, default=400
        Resolution of the figure
    chr_width : float, default=0.35
        Chromosome width in data coordinates
    chr_x : float, default=0.0
        Left x boundary of the chromosome in data coordinates
    anno_x_pad : float, default=0.14
        Extra horizontal gap between chromosome right edge and annotation area (data coords)
    use_lead_extraction : bool, default=True
        If True, extract lead variants via ``_get_sig()``. If False, use every input row.
    annotate_snps : bool, default=True
        If True, annotate lead SNP positions
    annotation_kwargs : dict, optional
        Additional keyword arguments for text/annotate styling.
    anno_arrow_shaft : float, optional, default=18
        Arrow shaft length in points (text mode and marker connector offset).
    anno_arrow_pad : float, optional, default=10
        Gap between arrow end and text box (``shrinkA``), text mode only.
    anno_arrow_shrink_b : float, optional, default=4
        Gap between arrow head and chromosome (``shrinkB``), text mode only.
    anno_text : str, optional
        Column for annotation text labels (e.g. gene name).
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
    marker_size : float, default=42
        Scatter marker size for annotation markers.
    marker_gap_pt : float, default=14
        Horizontal gap between markers within a group, in points.
    marker_label_gap_pt : float, default=6
        Vertical gap between marker row and text label below, in points.
    marker_label_align : str, default="center"
        Horizontal alignment of text label relative to marker row (``"center"``,
        ``"left"``, or ``"right"``).
    marker_fontsize : float, default=10
        Font size for marker-mode text labels.
    marker_linewidth : float, default=0.6
        Edge linewidth for marker scatter points.
    marker_label_bbox : bool, default=True
        If True, draw marker-mode labels with a white text outline (no filled box).
    group_min_vertical_gap_pt : float, default=18
        Minimum vertical spacing between group annotation blocks, in points.
    group_marker_to_marker_gap_pt : float, default=16
        Minimum spacing between marker rows of adjacent groups, in points.
    group_label_box_pad_pt : float, default=1.5
        Extra padding around label text for overlap calculations, in points.
    show_legend : bool, default=True
        If True, draw a figure-level marker legend below the plot (marker mode only).
    legend_ncol : int, default=6
        Number of columns in the figure legend.
    legend_kwargs : dict, optional
        Extra keyword arguments passed to ``fig.legend()``.
    anno_style : str, default="expand"
        Annotation layout style for lead SNP labels
    repel_force : float, default=0.5
        Repulsion force for text-mode label separation. Marker mode uses display
        spreading with ``anno_max_iter`` and always applies the full adjustment.
    anno_max_iter : int, default=300
        Maximum iterations for annotation repulsion / marker group spreading
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
    
    marker_mode = anno_group is not None
    if legend_kwargs is None:
        legend_kwargs = {}
    
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
            sig_level=sig_level,
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
    
    # Convert chromosome names to numbers for sorting
    chr_to_num = get_chr_to_number(out_chr=True, xymt=["X", "Y", "MT"])
    num_to_chr = get_number_to_chr()
    
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
    
    # Calculate offsets for each chromosome
    max_chr_size = max(chr_sizes.values())
    offset = [0 for _ in range(ncols)]
    max_row_offset = 0
    row_n = 0
    
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
                        leads_dict.setdefault(chr_str, []).append({
                            'pos': int(pos_val) if pd.notna(pos_val) else None,
                            'snpid': row[snpid] if snpid in leads_copy.columns and pd.notna(row[snpid]) else None,
                            'index': idx,
                            'anno_text': _lead_anno_field(row, anno_text, leads_copy.columns),
                            'anno_group': _lead_anno_field(row, anno_group, leads_copy.columns),
                            'anno_shape': _lead_anno_field(row, anno_shape, leads_copy.columns),
                            'anno_color': _lead_anno_field(row, anno_color, leads_copy.columns),
                        })
    
    # Plot each chromosome
    chr_bottom_positions = {}  # Store bottom position for each chromosome for label placement
    for i, chr_name in enumerate(chr_list):
        chr_cytobands = cytobands.loc[cytobands["CHR"] == chr_name, :].copy()
        
        if i // ncols > row_n:
            row_n += 1
            max_row_offset += 1.2
        
        chr_size = chr_sizes[chr_name]
        
        # Get the offset for this specific chromosome
        chr_offset = max_row_offset
        
        # Calculate bottom position of chromosome for label placement
        # The chromosome extends from offset to offset + chr_size/max_chr_size
        # Plus telomere at bottom (0.02)
        # Since y-axis is inverted, bottom is at higher y value
        telemere_full_length = 0.02
        # Calculate the actual bottom of the chromosome including telomere
        # The telomere is centered at offset + chr_size/max_chr_size
        # and extends telemere_full_length/2 below that center
        # So the bottom edge is at: offset + chr_size/max_chr_size + telemere_full_length/2
        chr_bottom_with_telomere = chr_size / max_chr_size + chr_offset + telemere_full_length / 2
        chr_bottom_positions[i] = chr_bottom_with_telomere
        
        # Get centromere boundaries
        acen_bands = chr_cytobands.loc[chr_cytobands["STAIN"] == "acen", :]
        if len(acen_bands) > 0:
            chr_centromere_u = acen_bands["START"].min()
            chr_centromere_l = acen_bands["END"].max()
        else:
            # Default centromere position if not found
            chr_centromere_u = chr_size * 0.44
            chr_centromere_l = chr_size * 0.46
        
        # Plot chromosome
        _plot_chr(
            axes[i % ncols],
            chr_size=chr_size,
            chr_centromere_u=chr_centromere_u,
            chr_centromere_l=chr_centromere_l,
            max_chr_size=max_chr_size,
            offset=max_row_offset,
            chr_cytobands=chr_cytobands,
            chr_name=chr_name,
            leads=leads_dict.get(chr_name, []),
            chr_width=chr_width,
            chr_x=chr_x,
            anno_x_pad=anno_x_pad,
            annotate_snps=annotate_snps,
            annotation_kwargs=annotation_kwargs if annotation_kwargs else {},
            anno_arrow_shaft=anno_arrow_shaft,
            anno_arrow_pad=anno_arrow_pad,
            anno_arrow_shrink_b=anno_arrow_shrink_b,
            anno_text=anno_text,
            anno_group=anno_group,
            anno_shape=anno_shape,
            anno_color=anno_color,
            marker_shape_map=resolved_shape_map,
            marker_color_map=resolved_color_map,
            marker_size=marker_size,
            marker_gap_pt=marker_gap_pt,
            marker_label_gap_pt=marker_label_gap_pt,
            marker_label_align=marker_label_align,
            marker_fontsize=marker_fontsize,
            marker_linewidth=marker_linewidth,
            marker_label_bbox=marker_label_bbox,
            group_min_vertical_gap_pt=group_min_vertical_gap_pt,
            group_marker_to_marker_gap_pt=group_marker_to_marker_gap_pt,
            group_label_box_pad_pt=group_label_box_pad_pt,
            marker_mode=marker_mode,
            anno_style=anno_style,
            repel_force=repel_force,
            anno_max_iter=anno_max_iter,
            log=log,
            verbose=verbose
        )
    
    # All axes must share the same xlim. Otherwise, the same chr_width in data
    # coordinates will be rendered with different physical widths across subplots.
    global_xmin = chr_x - chr_width * 0.2
    default_xmax = chr_x + chr_width + anno_x_pad + 0.5
    auto_global_xmax = max(
        getattr(ax, "_phenogram_xmax", default_xmax)
        for ax in axes[:ncols]
    )
    fixed_global_xmax = chr_x + chr_width + anno_x_pad + 1.2
    global_xmax = max(auto_global_xmax, fixed_global_xmax)
    for ax in axes[:ncols]:
        ax.set_xlim(global_xmin, global_xmax)
    
    # Set axis properties
    for i in range(n_chr):
        axis_ymax = max(
            max_row_offset + 1,
            getattr(axes[i % ncols], "_phenogram_ymax", max_row_offset + 1),
        )
        axes[i % ncols].set_ylim(-0.1, axis_ymax)
        axes[i % ncols].invert_yaxis()
        axes[i % ncols].spines['top'].set_visible(False)
        axes[i % ncols].spines['right'].set_visible(False)
        axes[i % ncols].spines['bottom'].set_visible(False)
        axes[i % ncols].spines['left'].set_visible(False)
        # Add chromosome label - position right below the chromosome (avoiding telomere overlap)
        chr_num = chr_list[i][3:] if chr_list[i].startswith("chr") else chr_list[i]
        # Use stored bottom position for this chromosome
        if i in chr_bottom_positions:
            chr_bottom_with_telomere = chr_bottom_positions[i]
            # Position label below the telomere to avoid overlap
            # Since y-axis is inverted, "below" means higher y value
            # chr_bottom_with_telomere is already the bottom edge of telomere
            # Add padding (0.03) to ensure no overlap
            label_y = chr_bottom_with_telomere + 0.03
            axes[i % ncols].text(
                chr_x + chr_width / 2,
                label_y,
                chr_num,
                ha="center",
                va="bottom",
                fontsize=10,
                zorder=300,
            )
    
    # Hide unused axes
    for i in range(n_chr, len(axes)):
        axes[i].set_visible(False)
    
    if show_legend and marker_mode:
        legend_handles = _make_marker_legend_handles(resolved_shape_map, resolved_color_map)
        if legend_handles:
            legend_bottom_margin = 0.035
            fig.tight_layout(rect=[0, legend_bottom_margin, 1, 1])
            fig.canvas.draw()
            ax_bottom = min(
                axes[i % ncols].get_position().y0
                for i in range(n_chr)
            )
            fig.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, ax_bottom - 0.004),
                bbox_transform=fig.transFigure,
                ncol=legend_ncol,
                frameon=False,
                fontsize=9,
                **legend_kwargs,
            )
    else:
        fig.tight_layout()
    
    # Save figure
    save_figure(fig=fig, save=save, keyword="phenogram", save_kwargs=save_kwargs, 
                log=log, verbose=verbose)
    
    log.write("Finished creating phenogram plot.", verbose=verbose)
    
    return fig


def _plot_chr(
    ax,
    chr_size: float,
    chr_centromere_u: float,
    chr_centromere_l: float,
    max_chr_size: float,
    offset: float,
    chr_cytobands: pd.DataFrame,
    chr_name: str,
    leads: List[Dict[str, Any]],
    chr_width: float = 0.35,
    chr_x: float = 0.0,
    anno_x_pad: float = 0.14,
    annotate_snps: bool = True,
    annotation_kwargs: Dict[str, Any] = None,
    anno_arrow_shaft: Optional[float] = 18,
    anno_arrow_pad: Optional[float] = 10,
    anno_arrow_shrink_b: Optional[float] = 4,
    anno_text: Optional[str] = None,
    anno_group: Optional[str] = None,
    anno_shape: Optional[str] = None,
    anno_color: Optional[str] = None,
    marker_shape_map: Optional[Dict[str, str]] = None,
    marker_color_map: Optional[Dict[str, str]] = None,
    marker_size: float = 42,
    marker_gap_pt: float = 14,
    marker_label_gap_pt: float = 6,
    marker_label_align: str = "center",
    marker_fontsize: float = 10,
    marker_linewidth: float = 0.6,
    marker_label_bbox: bool = True,
    group_min_vertical_gap_pt: float = 18,
    group_marker_to_marker_gap_pt: float = 16,
    group_label_box_pad_pt: float = 1.5,
    marker_mode: bool = False,
    anno_style: str = "expand",
    repel_force: float = 0.5,
    anno_max_iter: int = 300,
    log: Log = Log(),
    verbose: bool = True
):
    """
    Plot a single chromosome with cytobands and optional SNP annotations.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to plot on
    chr_size : float
        Size of the chromosome
    chr_centromere_u : float
        Upper boundary of centromere
    chr_centromere_l : float
        Lower boundary of centromere
    max_chr_size : float
        Maximum chromosome size (for normalization)
    offset : float
        Y-axis offset for this chromosome
    chr_cytobands : pd.DataFrame
        Cytoband data for this chromosome
    chr_name : str
        Chromosome name
    leads : list
        List of lead SNP dictionaries with 'pos', 'snpid', 'index' keys
    chr_width : float, default=0.35
        Chromosome width in data coordinates
    chr_x : float, default=0.0
        Left x boundary of the chromosome
    anno_x_pad : float, default=0.14
        Extra horizontal gap between chromosome and annotation text (data coords)
    annotate_snps : bool
        If True, annotate SNP positions
    annotation_kwargs : dict
        Additional keyword arguments for matplotlib annotate
    anno_arrow_shaft : float, optional
        Arrow shaft length in points (after anno_x_pad, before text anchor)
    anno_arrow_pad : float, optional
        Arrow-to-text gap in points (shrinkA)
    anno_arrow_shrink_b : float, optional
        Arrow-to-marker gap in points (shrinkB)
    log : Log
        Logging object
    verbose : bool
        If True, print progress messages
    """
    
    if annotation_kwargs is None:
        annotation_kwargs = {}
    if marker_shape_map is None:
        marker_shape_map = {}
    if marker_color_map is None:
        marker_color_map = {}
    
    chr_center_x = chr_x + chr_width / 2
    chr_right_x = chr_x + chr_width
    chr_inner_x = chr_x + chr_width * 0.05
    chr_inner_width = chr_width * 0.90
    
    max_xlim = chr_right_x + anno_x_pad + 0.5
    
    positions = [
        0 + offset,
        chr_centromere_u / max_chr_size + offset,
        chr_centromere_l / max_chr_size + offset,
        chr_size / max_chr_size + offset
    ]
    
    centromere_full_length = (chr_centromere_l - chr_centromere_u) / max_chr_size
    telemere_full_length = 0.02
    
    height_for_arm1 = positions[1] - positions[0]
    height_for_arm2 = positions[3] - positions[2]
    
    full_length = height_for_arm1 + height_for_arm2 + centromere_full_length
    
    # Draw chromosome arms
    arms = [
        Rectangle((chr_x, positions[0]), width=chr_width, height=height_for_arm1),
        Rectangle((chr_x, positions[1] + centromere_full_length), width=chr_width, height=height_for_arm2)
    ]
    
    # Draw centromeres
    centromeres = [
        Ellipse(
            xy=(chr_center_x, positions[1] + centromere_full_length),
            width=chr_width,
            height=centromere_full_length,
        ),
        Ellipse(
            xy=(chr_center_x, positions[1]),
            width=chr_width,
            height=centromere_full_length,
        ),
    ]
    
    # Draw telomeres
    telemeres = [
        Ellipse(
            xy=(chr_center_x, positions[0]),
            width=chr_width,
            height=telemere_full_length,
        ),
        Ellipse(
            xy=(chr_center_x, positions[0] + full_length),
            width=chr_width,
            height=telemere_full_length,
        ),
    ]
    
    # Create patch collections
    arms_pc = PatchCollection(arms, facecolor="white", edgecolor="grey", zorder=100)
    centromeres_pc = PatchCollection(centromeres, facecolor="grey", edgecolor="black", zorder=99)
    telemeres_pc = PatchCollection(telemeres, facecolor="white", edgecolor="grey", zorder=1)
    
    # Add collections to axes
    ax.add_collection(arms_pc)
    ax.add_collection(centromeres_pc)
    ax.add_collection(telemeres_pc)
    
    # Draw cytobands
    for index, row in chr_cytobands.iterrows():
        if row["END"] <= chr_centromere_u:
            # Upper arm
            band_start = row["START"] / max_chr_size + offset
            band_end = row["END"] / max_chr_size + offset
        elif row["START"] >= chr_centromere_l:
            # Lower arm
            band_start = (row["START"] - chr_centromere_l) / max_chr_size + height_for_arm1 + offset + centromere_full_length
            band_end = (row["END"] - chr_centromere_l) / max_chr_size + height_for_arm1 + offset + centromere_full_length
        else:
            # Skip centromere bands (already drawn)
            continue
        
        band_height = band_end - band_start
        
        if band_height > 0:
            band = Rectangle(
                (chr_inner_x, band_start),
                width=chr_inner_width,
                height=band_height,
            )
            facecolor = row["COLOR"]
            
            if row["STAIN"] == "stalk":
                bands_pc = PatchCollection([band], facecolor=facecolor, edgecolor=None, 
                                           linewidths=0, zorder=102)
            else:
                bands_pc = PatchCollection([band], facecolor=facecolor, edgecolor=None, 
                                         linewidths=0, zorder=101)
            ax.add_collection(bands_pc)
    
    # Lead SNP annotations (text mode or marker mode)
    if annotate_snps and len(leads) > 0:
        leads_sorted = sorted(
            [l for l in leads if l["pos"] is not None and pd.notna(l["pos"])],
            key=lambda x: x["pos"],
        )
        
        if len(leads_sorted) > 0:
            leads_with_y: List[Dict[str, Any]] = []
            for lead in leads_sorted:
                pos = lead["pos"]
                
                if pos <= chr_centromere_u:
                    y_pos_chr = pos / max_chr_size + offset
                elif pos >= chr_centromere_l:
                    y_pos_chr = (
                        (pos - chr_centromere_l) / max_chr_size
                        + height_for_arm1 + offset + centromere_full_length
                    )
                else:
                    continue
                
                leads_with_y.append({
                    "lead": lead,
                    "y_pos_chr": y_pos_chr,
                    "original_y": y_pos_chr,
                })
            
            if len(leads_with_y) > 0:
                if marker_mode:
                    text_kwargs = _build_marker_text_kwargs(
                        annotation_kwargs,
                        marker_fontsize=marker_fontsize,
                        marker_label_bbox=marker_label_bbox,
                    )
                    anno_fontsize = float(text_kwargs.get("fontsize", marker_fontsize))
                    
                    grouped: Dict[str, List[Dict[str, Any]]] = {}
                    for item in leads_with_y:
                        group_key = _get_group_key(item["lead"])
                        grouped.setdefault(group_key, []).append(item)
                    
                    group_representatives: List[Dict[str, Any]] = []
                    for group_key, group_items in grouped.items():
                        group_representatives.append({
                            "group_key": group_key,
                            "group_items": group_items,
                            "original_y": group_items[0]["original_y"],
                        })
                    
                    block_above_pt, block_below_pt = _marker_block_extents_pt(
                        marker_size,
                        marker_label_gap_pt,
                        marker_fontsize,
                        group_label_box_pad_pt,
                    )
                    
                    orig_y = np.array([g["original_y"] for g in group_representatives])
                    spread_y = _spread_group_blocks_display(
                        ax,
                        orig_y.tolist(),
                        ref_x=chr_right_x,
                        above_pt=block_above_pt,
                        below_pt=block_below_pt,
                        gap_pt=group_min_vertical_gap_pt,
                        max_iter=anno_max_iter,
                    )
                    marker_ys = spread_y
                    
                    for i, g in enumerate(group_representatives):
                        g["marker_y"] = float(marker_ys[i])
                    
                    arrow_shaft_pt = _resolve_anno_arrow_pt(
                        anno_arrow_shaft,
                        "arrow_shaft",
                        annotation_kwargs,
                        18.0,
                    )
                    
                    autoscale_x_on = ax.get_autoscalex_on()
                    ax.set_autoscalex_on(False)
                    
                    for g in group_representatives:
                        group_key = g["group_key"]
                        group_items = g["group_items"]
                        marker_y = float(g["marker_y"])
                        label_text = _get_marker_group_label(group_key, group_items)
                        if len(str(label_text)) > 20:
                            label_text = str(label_text)[:17] + "..."
                        
                        main_item = group_items[0]
                        main_y = main_item["original_y"]
                        
                        dx_shaft, _ = _points_to_data_delta(
                            ax, chr_right_x, marker_y, dx_points=arrow_shaft_pt
                        )
                        marker_start_x = chr_right_x + anno_x_pad + dx_shaft
                        
                        marker_xs = _marker_row_x_positions(
                            ax,
                            marker_start_x,
                            marker_y,
                            len(group_items),
                            marker_gap_pt,
                        )
                        for j, item in enumerate(group_items):
                            marker, color = _get_marker_style(
                                item["lead"], marker_shape_map, marker_color_map
                            )
                            marker_x = marker_xs[j]
                            ax.scatter(
                                marker_x,
                                marker_y,
                                marker=marker,
                                s=marker_size,
                                facecolor=color,
                                edgecolor="black",
                                linewidth=marker_linewidth,
                                clip_on=False,
                                zorder=235,
                            )
                        
                        if not marker_xs:
                            continue
                        
                        marker_center_x = (marker_xs[0] + marker_xs[-1]) / 2
                        if marker_label_align == "left":
                            label_x = marker_xs[0]
                        elif marker_label_align == "right":
                            label_x = marker_xs[-1]
                        else:
                            label_x = marker_center_x
                        label_x, label_y = _marker_label_position(
                            ax,
                            label_x,
                            marker_y,
                            marker_size,
                            marker_label_gap_pt,
                        )
                        
                        ax.plot(
                            [chr_right_x, marker_xs[0]],
                            [main_y, marker_y],
                            color="black",
                            linewidth=0.7,
                            clip_on=False,
                            zorder=190,
                        )
                        
                        ax.text(
                            label_x,
                            label_y,
                            str(label_text),
                            clip_on=False,
                            zorder=230,
                            **text_kwargs,
                        )
                        
                        max_xlim = max(max_xlim, marker_xs[-1])
                        text_width_dx, _ = _points_to_data_delta(
                            ax,
                            label_x,
                            label_y,
                            dx_points=anno_fontsize * len(str(label_text)) * 0.55,
                        )
                        ha = text_kwargs.get("ha", "center")
                        if ha == "center":
                            text_right_x = label_x + text_width_dx / 2
                        elif ha == "left":
                            text_right_x = label_x + text_width_dx
                        else:
                            text_right_x = label_x
                        max_xlim = max(max_xlim, text_right_x)
                        
                        _, label_bottom_y = _data_point_below(
                            ax, chr_right_x, marker_y, block_below_pt + 4.0
                        )
                        ax._phenogram_ymax = max(
                            getattr(ax, "_phenogram_ymax", 0),
                            label_bottom_y + 0.01,
                        )
                    
                    ax.set_autoscalex_on(autoscale_x_on)
                else:
                    text_kwargs = _build_text_mode_kwargs(annotation_kwargs)
                    anno_fontsize = float(text_kwargs.get("fontsize", 8))
                    
                    y_positions = np.array([l["y_pos_chr"] for l in leads_with_y])
                    if len(y_positions) > 1:
                        y_span = max(y_positions) - min(y_positions)
                        if y_span < 0.01:
                            y_span = 0.1
                    else:
                        y_span = 0.1
                    
                    adjusted_y = adjust_text_position(
                        y_positions.copy(),
                        y_span,
                        repel_force=repel_force,
                        max_iter=anno_max_iter,
                        amode="float",
                        log=log,
                        verbose=verbose,
                    )
                    
                    for i, l in enumerate(leads_with_y):
                        l["y_pos_chr"] = adjusted_y[i]
                    
                    arrow_pad_pt = _resolve_anno_arrow_pt(
                        anno_arrow_pad,
                        "arrow_pad",
                        annotation_kwargs,
                        10.0,
                    )
                    arrow_shaft_pt = _resolve_anno_arrow_pt(
                        anno_arrow_shaft,
                        "arrow_shaft",
                        annotation_kwargs,
                        18.0,
                    )
                    shrink_b_pt = _resolve_anno_arrow_pt(
                        anno_arrow_shrink_b,
                        "arrow_shrink_b",
                        annotation_kwargs,
                        4.0,
                    )
                    
                    for l in leads_with_y:
                        lead = l["lead"]
                        original_y = l["original_y"]
                        text_y = l["y_pos_chr"]
                        
                        ax.plot(
                            [chr_x, chr_right_x],
                            [original_y, original_y],
                            color="red",
                            linewidth=2,
                            clip_on=False,
                            zorder=200,
                        )
                        
                        label_text = _get_text_mode_label(lead)
                        if not label_text:
                            continue
                        if len(str(label_text)) > 20:
                            label_text = str(label_text)[:17] + "..."
                        
                        chr_point_x = chr_right_x
                        chr_point_y = original_y
                        
                        dx_shaft, _ = _points_to_data_delta(
                            ax, chr_point_x, text_y, dx_points=arrow_shaft_pt
                        )
                        text_x = chr_point_x + anno_x_pad + dx_shaft
                        
                        text_width_dx, _ = _points_to_data_delta(
                            ax,
                            text_x,
                            text_y,
                            dx_points=anno_fontsize * len(str(label_text)) * 0.55,
                        )
                        max_xlim = max(max_xlim, text_x + text_width_dx)
                        
                        ax.annotate(
                            str(label_text),
                            xy=(chr_point_x, chr_point_y),
                            xytext=(text_x, text_y),
                            arrowprops=dict(
                                arrowstyle="-|>",
                                color="gray",
                                lw=0.8,
                                relpos=(0, 0.5),
                                shrinkA=arrow_pad_pt,
                                shrinkB=shrink_b_pt,
                            ),
                            clip_on=False,
                            zorder=201,
                            **text_kwargs,
                        )
    
    ax.set_xticks(ticks=[])
    ax.set_yticks(ticks=[])
    # Record required right boundary for this axis; xlim is set globally in _plot_phenogram().
    xlim_default = chr_right_x + anno_x_pad + 0.5
    xlim_right = max(
        getattr(ax, "_phenogram_xmax", xlim_default),
        max_xlim * 1.05,
    )
    ax._phenogram_xmax = xlim_right
    # Do NOT call ax.set_xlim() here — all axes share one xlim after all chromosomes are drawn.
