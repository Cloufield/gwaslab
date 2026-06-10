import math
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
from gwaslab.util.util_in_get_sig import _get_sig, _anno_gene
from gwaslab.bd.bd_common_data import get_chr_to_number, get_number_to_chr
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style


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
    kwargs_key: str,
    anno_kwargs: Dict[str, Any],
    default: float,
) -> float:
    """Resolve arrow spacing in points from ``anno_kwargs`` or default."""
    if kwargs_key in anno_kwargs:
        return float(anno_kwargs[kwargs_key])
    return float(default)


def _merge_lead_anno_kwargs(
    anno_kwargs: Dict[str, Any],
    anno_kwargs_single: Optional[Dict[str, Any]],
    snpid: Optional[str],
) -> Dict[str, Any]:
    """Merge global and per-SNP annotation kwargs (MQQ-compatible)."""
    merged = dict(anno_kwargs)
    if anno_kwargs_single and snpid is not None and snpid in anno_kwargs_single:
        merged.update(anno_kwargs_single[snpid])
    return merged


def _build_lead_text_kwargs(
    base_kwargs: Dict[str, Any],
    anno_kwargs_single: Optional[Dict[str, Any]],
    lead: Dict[str, Any],
    snpid_key: str = "snpid",
    force_ha: Optional[str] = "left",
    force_va: Optional[str] = "center",
) -> Dict[str, Any]:
    """Apply ``anno_kwargs_single`` overrides for one lead."""
    kwargs = dict(base_kwargs)
    snpid = lead.get(snpid_key)
    if anno_kwargs_single and snpid is not None and snpid in anno_kwargs_single:
        single = {
            k: v
            for k, v in anno_kwargs_single[snpid].items()
            if k not in _TEXT_SKIP_KEYS
        }
        kwargs.update(single)
    if force_ha is not None:
        kwargs["ha"] = force_ha
    if force_va is not None:
        kwargs["va"] = force_va
    return kwargs


def _resolve_lead_connector_style(
    anno_kwargs: Dict[str, Any],
    arrow_kwargs: Optional[Dict[str, Any]],
    anno_kwargs_single: Optional[Dict[str, Any]],
    snpid: Optional[str],
) -> Tuple[float, float, float, Dict[str, Any], Dict[str, Any]]:
    """Resolve connector point lengths and matplotlib arrow/line styling."""
    merged = _merge_lead_anno_kwargs(anno_kwargs, anno_kwargs_single, snpid)
    arrow_pad_pt = _resolve_anno_arrow_pt("arrow_pad", merged, 10.0)
    arrow_shaft_pt = _resolve_anno_arrow_pt("arrow_shaft", merged, 18.0)
    shrink_b_pt = _resolve_anno_arrow_pt("arrow_shrink_b", merged, 4.0)

    arrowprops: Dict[str, Any] = {
        "arrowstyle": "-|>",
        "color": "gray",
        "lw": 0.8,
        "shrinkA": arrow_pad_pt,
        "shrinkB": 0,
    }
    line_kwargs: Dict[str, Any] = {
        "color": "gray",
        "linewidth": 0.8,
    }
    if arrow_kwargs:
        arrowprops.update(arrow_kwargs)
        if "color" in arrow_kwargs:
            line_kwargs["color"] = arrow_kwargs["color"]
        lw = arrow_kwargs.get("lw", arrow_kwargs.get("linewidth"))
        if lw is not None:
            line_kwargs["linewidth"] = lw
    arrowprops["shrinkA"] = _resolve_anno_arrow_pt("arrow_pad", merged, arrow_pad_pt)
    return arrow_shaft_pt, arrow_pad_pt, shrink_b_pt, arrowprops, line_kwargs


def _significance_sort_key(
    data: pd.DataFrame,
    p: str,
    mlog10p: str,
) -> Tuple[Optional[str], bool]:
    """Return (sort_column, ascending) for ranking variants by significance."""
    if mlog10p in data.columns:
        return mlog10p, False
    if p in data.columns:
        return p, True
    return None, True


def _label_eligible_mask(
    leads: pd.DataFrame,
    snpid: str,
    anno_set: Optional[List[str]],
    anno_max_rows: int,
    p: str,
    mlog10p: str,
    log: Log = Log(),
    verbose: bool = True,
) -> pd.Series:
    """Boolean mask of lead rows that may receive text labels."""
    if len(leads) == 0:
        return pd.Series(dtype=bool)

    mask = pd.Series(True, index=leads.index)
    if anno_set is not None and len(anno_set) > 0:
        if snpid in leads.columns:
            mask = leads[snpid].isin(anno_set)
        else:
            mask = pd.Series(False, index=leads.index)
        if mask.any():
            log.write(
                " -Found {} variant(s) in anno_set for labeling.".format(int(mask.sum())),
                verbose=verbose,
            )

    eligible = leads.index[mask]
    if len(eligible) > anno_max_rows:
        log.write(
            " -Limiting labels to top {} variant(s) by significance.".format(anno_max_rows),
            verbose=verbose,
        )
        subset = leads.loc[eligible]
        sort_column, ascending = _significance_sort_key(subset, p, mlog10p)
        if sort_column is not None:
            subset = subset.sort_values(by=sort_column, ascending=ascending).head(anno_max_rows)
        else:
            subset = subset.head(anno_max_rows)
        mask = pd.Series(leads.index.isin(subset.index), index=leads.index)
    return mask


def _prepare_phenogram_annotation_column(
    leads: pd.DataFrame,
    anno: Optional[Union[bool, str]],
    snpid: str,
    chrom: str,
    pos: str,
    build: str,
    anno_source: str,
    anno_gtf_path: Optional[str],
    log: Log = Log(),
    verbose: bool = True,
) -> pd.DataFrame:
    """Add an ``Annotation`` column for column-based or GENENAME labels."""
    if len(leads) == 0 or anno is None or anno is True:
        return leads

    if anno == "GENENAME":
        log.write(" -Annotating lead variants with nearest gene names...", verbose=verbose)
        annotated = _anno_gene(
            leads,
            id=snpid,
            chrom=chrom,
            pos=pos,
            log=log,
            build=build,
            source=anno_source,
            gtf_path=anno_gtf_path,
            verbose=verbose,
        )
        annotated = annotated.copy()
        annotated["Annotation"] = annotated["GENE"]
        return annotated

    if isinstance(anno, str):
        if anno not in leads.columns:
            log.warning(" -Column '{}' not found for annotation.".format(anno), verbose=verbose)
            out = leads.copy()
            out["Annotation"] = pd.NA
            return out
        out = leads.copy()
        out["Annotation"] = leads[anno].astype("string")
        return out

    return leads


def _resolve_phenogram_label(
    lead: Dict[str, Any],
    anno: Optional[Union[bool, str]],
    anno_alias: Optional[Dict[str, str]] = None,
) -> Optional[str]:
    """Resolve display label for a phenogram lead (MQQ-compatible ``anno`` rules)."""
    if anno is None:
        return None

    snpid_val = lead.get("snpid")
    alias = anno_alias if anno_alias is not None else {}
    if snpid_val is not None and snpid_val in alias:
        return str(alias[snpid_val])

    if anno is True:
        chrom_val = lead.get("chrom")
        pos_val = lead.get("pos")
        if chrom_val is None or pos_val is None or pd.isna(pos_val):
            return None
        return "Chr{}:{}".format(chrom_val, int(pos_val))

    if isinstance(anno, str):
        ann = lead.get("annotation")
        if ann is None or (isinstance(ann, float) and pd.isna(ann)):
            return None
        return str(ann)

    return None


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
    if lead.get("annotation") is not None:
        return str(lead["annotation"])
    if lead.get("snpid") is not None:
        return str(lead["snpid"])
    return str(lead.get("index", "unknown"))


def _get_marker_group_label(group_key: str, group_items: List[Dict[str, Any]]) -> str:
    first_lead = group_items[0]["lead"]
    if first_lead.get("annotation") is not None:
        return str(first_lead["annotation"])
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


def _build_text_mode_kwargs(anno_kwargs: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    text_kwargs = {
        "fontsize": 9,
        "ha": "left",
        "va": "center",
        "fontweight": "normal",
    }
    if anno_kwargs:
        clean = {k: v for k, v in anno_kwargs.items() if k not in _TEXT_SKIP_KEYS}
        text_kwargs.update(clean)
    text_kwargs["ha"] = "left"
    text_kwargs["va"] = "center"
    return text_kwargs


def _text_label_block_extents_pt(
    anno_fontsize: float,
    group_label_box_pad_pt: float,
) -> Tuple[float, float]:
    """Vertical half-extents of a horizontal text label (``va='center'``) in points."""
    half_text_pt = anno_fontsize * 0.65 + group_label_box_pad_pt
    return half_text_pt, half_text_pt


def _spread_phenogram_text_labels(
    ax,
    y_data_list: List[float],
    ref_x: float,
    anno_fontsize: float,
    gap_pt: float,
    group_label_box_pad_pt: float,
    max_iter: int,
) -> np.ndarray:
    """
    Spread text-label anchor rows in display space so label bounding boxes do not overlap.

    ``adjust_text_position`` uses data-coordinate steps derived from ``y_span`` that do
    not match rendered text height on the phenogram axes, so dense / neighboring variants
    still overlap. Marker mode already uses display-space spreading; text mode does the same.
    """
    above_pt, below_pt = _text_label_block_extents_pt(anno_fontsize, group_label_box_pad_pt)
    return _spread_group_blocks_display(
        ax,
        y_data_list,
        ref_x=ref_x,
        above_pt=above_pt,
        below_pt=below_pt,
        gap_pt=gap_pt,
        max_iter=max_iter,
    )


def _data_point_below(ax, x: float, y: float, dy_points: float) -> Tuple[float, float]:
    """Return data coordinates dy_points below (x, y) in display space."""
    scale = ax.figure.dpi / 72.0
    x_disp, y_disp = ax.transData.transform((x, y))
    return ax.transData.inverted().transform((x_disp, y_disp + dy_points * scale))


def _anno_column_x(chr_right_x: float, anno_x_pad: float) -> float:
    """Fixed left edge of the annotation column to the right of the chromosome."""
    return chr_right_x + anno_x_pad


def _plot_phenogram_connector(
    ax,
    chr_right_x: float,
    anno_col_x: float,
    locus_y: float,
    target_y: float,
    color: str = "black",
    linewidth: float = 0.7,
    zorder: int = 190,
) -> None:
    """Draw an elbow connector that stays outside the chromosome body."""
    ax.plot(
        [chr_right_x, anno_col_x],
        [locus_y, locus_y],
        color=color,
        linewidth=linewidth,
        clip_on=False,
        zorder=zorder,
    )
    if abs(target_y - locus_y) > 1e-12:
        ax.plot(
            [anno_col_x, anno_col_x],
            [locus_y, target_y],
            color=color,
            linewidth=linewidth,
            clip_on=False,
            zorder=zorder,
        )


def _cap_phenogram_arrow_shaft_pt(
    ax,
    chr_right_x: float,
    locus_y: float,
    text_x: float,
    arrow_shaft_pt: float,
    min_slant_gap_pt: float = 6.0,
) -> float:
    """Keep the fixed horizontal arm shorter than the chr-edge to text gap."""
    if text_x <= chr_right_x:
        return 0.0
    chr_disp = ax.transData.transform((chr_right_x, locus_y))
    text_disp = ax.transData.transform((text_x, locus_y))
    gap_pt = abs(text_disp[0] - chr_disp[0]) * 72.0 / ax.figure.dpi
    max_shaft_pt = max(2.0, gap_pt - min_slant_gap_pt)
    return min(float(arrow_shaft_pt), max_shaft_pt)


def _plot_phenogram_text_connector(
    ax,
    chr_right_x: float,
    locus_y: float,
    text_x: float,
    text_y: float,
    arrow_shaft_pt: float,
    arrow_pad_pt: float = 10.0,
    shrink_b_pt: float = 4.0,
    arrowprops: Optional[Dict[str, Any]] = None,
    line_kwargs: Optional[Dict[str, Any]] = None,
    zorder: int = 201,
) -> None:
    """
    Phenogram text annotation connector (vertical MQQ analogue).

    A fixed-length horizontal segment points at the locus on the chromosome edge.
    A single slanted segment links the right end of that arm to the left-middle
    anchor of the label (``ha='left'``, ``va='center'``).
    """
    shaft_pt = _cap_phenogram_arrow_shaft_pt(
        ax, chr_right_x, locus_y, text_x, arrow_shaft_pt
    )
    dx_arm, _ = _points_to_data_delta(
        ax, chr_right_x, locus_y, dx_points=shaft_pt
    )
    elbow_x = chr_right_x + dx_arm

    head_props = {
        "arrowstyle": "-|>",
        "color": "gray",
        "lw": 0.8,
        "shrinkA": arrow_pad_pt,
        "shrinkB": 0,
    }
    if arrowprops:
        head_props.update(arrowprops)
    head_props.setdefault("shrinkA", arrow_pad_pt)
    head_props.setdefault("shrinkB", 0)

    line_style = {"color": "gray", "linewidth": 0.8}
    if line_kwargs:
        line_style.update(line_kwargs)
    line_style.setdefault("color", head_props.get("color", "gray"))
    if "lw" in head_props and "linewidth" not in line_style:
        line_style["linewidth"] = head_props["lw"]

    ax.annotate(
        "",
        xy=(chr_right_x, locus_y),
        xytext=(elbow_x, locus_y),
        arrowprops=head_props,
        clip_on=False,
        zorder=zorder,
    )

    text_attach_x, text_attach_y = text_x, text_y
    if shrink_b_pt > 0:
        elbow_disp = np.array(ax.transData.transform((elbow_x, locus_y)))
        text_disp = np.array(ax.transData.transform((text_x, text_y)))
        vec = text_disp - elbow_disp
        norm = float(np.linalg.norm(vec))
        if norm > 1e-6:
            shrink_px = shrink_b_pt * ax.figure.dpi / 72.0
            attach_disp = text_disp - vec / norm * min(shrink_px, norm * 0.95)
            text_attach_x, text_attach_y = ax.transData.inverted().transform(
                attach_disp
            )

    if (
        abs(text_attach_x - elbow_x) > 1e-12
        or abs(text_attach_y - locus_y) > 1e-12
    ):
        ax.plot(
            [elbow_x, text_attach_x],
            [locus_y, text_attach_y],
            clip_on=False,
            zorder=zorder,
            **line_style,
        )


def _marker_radius_pt(marker_size: float) -> float:
    """Matplotlib scatter ``s`` is area in points²; return marker radius in points."""
    return np.sqrt(max(marker_size, 1.0) / np.pi)


def _effective_marker_radius_pt(
    marker_size: float,
    marker_linewidth: float = 0.0,
) -> float:
    """Scatter radius including edge linewidth for layout spacing."""
    return _marker_radius_pt(marker_size) + float(marker_linewidth) / 2.0


def _marker_group_block_extents_pt(
    marker_size: float,
    marker_label_gap_pt: float,
    marker_fontsize: float,
    group_label_box_pad_pt: float,
    n_marker_rows: int = 1,
    marker_row_gap_pt: float = 10.0,
    marker_linewidth: float = 0.6,
) -> Tuple[float, float]:
    """Vertical half-extents of a marker group block including wrapped marker rows."""
    marker_radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    text_height_pt = marker_fontsize * 1.25 + 0.15 * 2.0
    n_marker_rows = max(1, int(n_marker_rows))
    if n_marker_rows <= 1:
        row_stack_pt = 0.0
    else:
        row_pitch_pt = 2.0 * marker_radius_pt + float(marker_row_gap_pt)
        row_stack_pt = (n_marker_rows - 1) * row_pitch_pt
    above_pt = marker_radius_pt
    below_pt = (
        marker_radius_pt
        + row_stack_pt
        + marker_label_gap_pt
        + text_height_pt
        + group_label_box_pad_pt * 2.0
    )
    return above_pt, below_pt


def _marker_block_extents_pt(
    marker_size: float,
    marker_label_gap_pt: float,
    marker_fontsize: float,
    group_label_box_pad_pt: float,
    bbox_pad: float = 0.15,
) -> Tuple[float, float]:
    """Return (above_pt, below_pt) for a single-row marker group (legacy helper)."""
    return _marker_group_block_extents_pt(
        marker_size,
        marker_label_gap_pt,
        marker_fontsize,
        group_label_box_pad_pt,
        n_marker_rows=1,
    )


def _marker_center_step_pt(
    marker_size: float,
    marker_gap_pt: float,
    marker_linewidth: float = 0.6,
) -> float:
    """Fixed center-to-center spacing between markers in display points."""
    return (
        2.0 * _effective_marker_radius_pt(marker_size, marker_linewidth)
        + float(marker_gap_pt)
    )


def _compute_group_marker_positions(
    ax,
    anno_col_x: float,
    anchor_y: float,
    group_items: List[Dict[str, Any]],
    marker_size: float,
    marker_gap_pt: float,
    marker_max_per_row: int,
    marker_row_gap_pt: float,
    marker_linewidth: float = 0.6,
) -> Tuple[List[Tuple[float, float, Dict[str, Any]]], float]:
    """
    Lay out group markers on a fixed grid: ``marker_max_per_row`` per row, then wrap.

    Returns list of (x, y, item) and the y coordinate of the bottom marker row.
    """
    if marker_max_per_row < 1:
        marker_max_per_row = 1

    radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    step_pt = _marker_center_step_pt(marker_size, marker_gap_pt, marker_linewidth)
    row_pitch_pt = 2.0 * radius_pt + float(marker_row_gap_pt)

    dx_radius, _ = _points_to_data_delta(ax, anno_col_x, anchor_y, dx_points=radius_pt)
    start_center_x = anno_col_x + dx_radius

    positions: List[Tuple[float, float, Dict[str, Any]]] = []
    for idx, item in enumerate(group_items):
        row = idx // marker_max_per_row
        col = idx % marker_max_per_row
        _, dy_row = _points_to_data_delta(
            ax, anno_col_x, anchor_y, dy_points=row * row_pitch_pt
        )
        y = anchor_y + dy_row
        dx_col, _ = _points_to_data_delta(
            ax, start_center_x, y, dx_points=col * step_pt
        )
        x = start_center_x + dx_col
        positions.append((x, y, item))

    bottom_row_y = positions[-1][1] if positions else anchor_y
    return positions, bottom_row_y


def _marker_label_position(
    ax,
    label_x: float,
    marker_y: float,
    marker_size: float,
    marker_label_gap_pt: float,
    marker_linewidth: float = 0.6,
) -> Tuple[float, float]:
    """Return (label_x, label_y) with va='top' anchor below the marker row bottom edge."""
    scale = ax.figure.dpi / 72.0
    marker_radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    x_disp, y_disp = ax.transData.transform((label_x, marker_y))
    label_top_disp = y_disp + (marker_radius_pt + marker_label_gap_pt) * scale
    return ax.transData.inverted().transform((x_disp, label_top_disp))


def _marker_row_x_positions(
    ax,
    start_x: float,
    start_y: float,
    n_markers: int,
    gap_pt: float,
    marker_size: float = 42.0,
) -> List[float]:
    """Return evenly spaced marker center x positions with fixed center-to-center step."""
    if n_markers <= 0:
        return []
    step_pt = _marker_center_step_pt(marker_size, gap_pt)
    scale = ax.figure.dpi / 72.0
    step_px = step_pt * scale
    x0_disp, y_disp = ax.transData.transform((start_x, start_y))
    xs_data: List[float] = []
    for j in range(n_markers):
        x_disp = x0_disp + j * step_px
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
    anno_kwargs: Optional[Dict[str, Any]],
    marker_fontsize: float = 11,
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
    if anno_kwargs:
        clean = {k: v for k, v in anno_kwargs.items() if k not in _TEXT_SKIP_KEYS}
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


def _build_leads_with_y(
    leads: List[Dict[str, Any]],
    chr_centromere_u: float,
    chr_centromere_l: float,
    max_chr_size: float,
    offset: float,
    height_for_arm1: float,
    centromere_full_length: float,
) -> List[Dict[str, Any]]:
    """Map lead genomic positions to y coordinates on the chromosome drawing."""
    leads_sorted = sorted(
        [l for l in leads if l.get("pos") is not None and pd.notna(l["pos"])],
        key=lambda x: x["pos"],
    )
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
    return leads_with_y


def _estimate_global_annotation_xmax(
    chr_specs: List[Dict[str, Any]],
    global_xmin: float,
    provisional_xmax: float,
    marker_mode: bool,
    marker_size: float,
    marker_gap_pt: float,
    marker_linewidth: float,
    marker_max_per_row: int,
    marker_row_gap_pt: float,
    marker_fontsize: float,
    marker_label_gap_pt: float,
    group_label_box_pad_pt: float,
    anno_kwargs: Optional[Dict[str, Any]],
) -> float:
    """Estimate right xlim from marker/text width using frozen axis transforms."""
    global_xmax = provisional_xmax
    text_kwargs = (
        _build_marker_text_kwargs(
            anno_kwargs,
            marker_fontsize=marker_fontsize,
            marker_label_bbox=True,
        )
        if marker_mode
        else _build_text_mode_kwargs(anno_kwargs)
    )
    anno_fontsize = float(
        text_kwargs.get("fontsize", marker_fontsize if marker_mode else 9)
    )

    for spec in chr_specs:
        ax = spec["ax"]
        ax.set_xlim(global_xmin, provisional_xmax)
        leads_with_y = spec.get("leads_with_y") or []
        if not leads_with_y:
            continue
        anno_col_x = spec["anno_col_x"]
        chr_xmax = spec["chr_x"] + spec["chr_width"] + spec["anno_x_pad"] + 1.2
        local_xmax = chr_xmax

        if marker_mode:
            grouped: Dict[str, List[Dict[str, Any]]] = {}
            for item in leads_with_y:
                group_key = _get_group_key(item["lead"])
                grouped.setdefault(group_key, []).append(item)
            for group_key, group_items in grouped.items():
                anchor_y = float(group_items[0]["original_y"])
                label_text = _get_marker_group_label(group_key, group_items)
                if len(str(label_text)) > 20:
                    label_text = str(label_text)[:17] + "..."
                marker_positions, bottom_row_y = _compute_group_marker_positions(
                    ax,
                    anno_col_x,
                    anchor_y,
                    group_items,
                    marker_size,
                    marker_gap_pt,
                    marker_max_per_row,
                    marker_row_gap_pt,
                    marker_linewidth=marker_linewidth,
                )
                if marker_positions:
                    local_xmax = max(local_xmax, max(x for x, _, _ in marker_positions))
                label_x = anno_col_x
                _, label_y = _marker_label_position(
                    ax,
                    label_x,
                    bottom_row_y if marker_positions else anchor_y,
                    marker_size,
                    marker_label_gap_pt,
                    marker_linewidth=marker_linewidth,
                )
                text_width_dx, _ = _points_to_data_delta(
                    ax,
                    label_x,
                    label_y,
                    dx_points=anno_fontsize * len(str(label_text)) * 0.55,
                )
                local_xmax = max(local_xmax, label_x + text_width_dx)
        else:
            for item in leads_with_y:
                label_text = item["lead"].get("label_text")
                if not label_text:
                    continue
                if len(str(label_text)) > 20:
                    label_text = str(label_text)[:17] + "..."
                text_y = item["y_pos_chr"]
                text_width_dx, _ = _points_to_data_delta(
                    ax,
                    anno_col_x,
                    text_y,
                    dx_points=anno_fontsize * len(str(label_text)) * 0.55,
                )
                local_xmax = max(local_xmax, anno_col_x + text_width_dx)

        global_xmax = max(global_xmax, local_xmax * 1.05)

    return global_xmax


def _chr_ideogram_top_y(chr_offset: float, telemere_full_length: float = 0.02) -> float:
    """Data y of the top edge of the chromosome (including upper telomere cap)."""
    return chr_offset - telemere_full_length / 2.0


def _display_y_offset_data(
    ax,
    x: float,
    y: float,
    pt_offset: float,
    direction: str = "above",
) -> float:
    """
    Convert a display-point vertical offset at ``(x, y)`` to a data y coordinate.

    ``direction='above'`` moves toward the visual top of the axes (works with
    inverted y as well as normal y).
    """
    scale = ax.figure.dpi / 72.0
    x_disp, y_disp = ax.transData.transform((x, y))
    sign = -1.0 if direction == "above" else 1.0
    _, y_data = ax.transData.inverted().transform(
        (x_disp, y_disp + sign * pt_offset * scale)
    )
    return float(y_data)


def _chr_number_label_above_data(
    ax,
    ref_x: float,
    body_top_y: float,
    fontsize: float = 10.0,
    pad_pt: float = 4.0,
) -> float:
    """Return data-space height reserved above the ideogram top for its number."""
    text_pt = fontsize * 1.25
    total_pt = text_pt + pad_pt
    label_top_y = _display_y_offset_data(ax, ref_x, body_top_y, total_pt, direction="above")
    return abs(body_top_y - label_top_y)


def _chr_number_label_position(
    ax,
    body_top_y: float,
    ref_x: float,
    fontsize: float = 10.0,
    pad_pt: float = 4.0,
) -> Tuple[float, str]:
    """
    Return (y, va) with the label entirely above the ideogram top edge.

    Uses ``va='top'`` so upright text extends away from the body. With an
    inverted y-axis, ``va='bottom'`` would hang text downward into the chr body.
    """
    label_y = _display_y_offset_data(ax, ref_x, body_top_y, pad_pt, direction="above")
    return label_y, "top"


def _compute_phenogram_chr_offsets(
    chr_list: List[str],
    chr_sizes: Dict[str, float],
    max_chr_size: float,
    ncols: int,
    ax,
    ref_x: float,
    chr_number_fontsize: float = 10.0,
    chr_number_pad_pt: float = 4.0,
    inter_row_gap: float = 0.12,
) -> Tuple[int, List[float], float, float]:
    """
    Compute y offsets for stacked chromosomes on each column axis.

    Each visual row starts below the tallest chromosome in the row above (including
    its number label above the body), plus ``inter_row_gap``.
    """
    n_chr = len(chr_list)
    n_grid_rows = max(1, (n_chr + ncols - 1) // ncols)
    telemere_full_length = 0.02
    label_above_ref = _chr_number_label_above_data(
        ax,
        ref_x,
        _chr_ideogram_top_y(0.0, telemere_full_length),
        fontsize=chr_number_fontsize,
        pad_pt=chr_number_pad_pt,
    )

    row_max_extent = [0.0] * n_grid_rows
    for i, chr_name in enumerate(chr_list):
        grid_row = i // ncols
        visual_row = (n_grid_rows - 1) - grid_row
        chr_top_y = 0.0  # reference; label height weakly y-dependent
        body_top_y = _chr_ideogram_top_y(chr_top_y, telemere_full_length)
        label_above = _chr_number_label_above_data(
            ax,
            ref_x,
            body_top_y,
            fontsize=chr_number_fontsize,
            pad_pt=chr_number_pad_pt,
        )
        extent = (
            label_above
            + chr_sizes[chr_name] / max_chr_size
            + telemere_full_length / 2
        )
        row_max_extent[visual_row] = max(row_max_extent[visual_row], extent)

    offset_by_visual_row = [0.0]
    for visual_row in range(1, n_grid_rows):
        offset_by_visual_row.append(
            offset_by_visual_row[-1] + row_max_extent[visual_row - 1] + inter_row_gap
        )

    return n_grid_rows, offset_by_visual_row, telemere_full_length, label_above_ref


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
    marker_max_per_row: int = 5,
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
        Vertical gap between the bottom marker row and the group text label, in points.
    marker_max_per_row : int, default=5
        Maximum markers on one row before wrapping to the next row within a group.
    marker_row_gap_pt : float, default=2.5
        Vertical gap between wrapped marker rows inside a group, in points.
    marker_label_align : str, default="center"
        Horizontal alignment of text label relative to marker row (``"center"``,
        ``"left"``, or ``"right"``).
    marker_fontsize : float, default=11
        Font size for marker-mode text labels.
    marker_linewidth : float, default=0.6
        Edge linewidth for marker scatter points.
    marker_label_bbox : bool, default=True
        If True, draw marker-mode labels with a white text outline (no filled box).
    group_min_vertical_gap_pt : float, default=3.5
        Minimum vertical spacing between group annotation blocks, in points.
    group_marker_to_marker_gap_pt : float, default=3
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
        Scales minimum vertical gap between text labels in points
        (``group_min_vertical_gap_pt * repel_force``). Marker mode uses display
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
    if not marker_mode and len(leads) > 0 and anno is not None:
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
    
    # Calculate offsets for each chromosome (row pitch from tallest chr per row).
    max_chr_size = max(chr_sizes.values())
    chr_number_fontsize = 10.0
    chr_number_pad_pt = 4.0
    chr_number_ref_x = chr_x + chr_width / 2
    global_xmin = chr_x - chr_width * 0.2
    fixed_global_xmax = chr_x + chr_width + anno_x_pad + 1.2
    provisional_ymax = 3.0
    for col in range(ncols):
        axes[col].set_xlim(global_xmin, fixed_global_xmax)
        axes[col].set_ylim(-0.05, provisional_ymax)

    n_grid_rows, offset_by_visual_row, telemere_full_length, label_above_ref = (
        _compute_phenogram_chr_offsets(
            chr_list,
            chr_sizes,
            max_chr_size,
            ncols,
            axes[0],
            ref_x=chr_number_ref_x,
            chr_number_fontsize=chr_number_fontsize,
            chr_number_pad_pt=chr_number_pad_pt,
        )
    )
    max_row_offset = offset_by_visual_row[-1] if offset_by_visual_row else 0.0
    global_ymin = -label_above_ref - 0.03
    
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
    col_ymax = [global_ymin for _ in range(ncols)]
    chr_specs: List[Dict[str, Any]] = []
    for i, chr_name in enumerate(chr_list):
        chr_cytobands = cytobands.loc[cytobands["CHR"] == chr_name, :].copy()

        grid_row = i // ncols
        visual_row = (n_grid_rows - 1) - grid_row
        chr_offset = offset_by_visual_row[visual_row]

        chr_size = chr_sizes[chr_name]

        chr_bottom_with_telomere = (
            chr_size / max_chr_size + chr_offset + telemere_full_length / 2
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
        chr_specs.append(spec)

        col = i % ncols
        col_ymax[col] = max(
            col_ymax[col],
            chr_offset + 1,
            chr_bottom_with_telomere + chr_label_pad,
        )

    # All column axes must share the same ylim. Otherwise ``chr_size / max_chr_size``
    # in data coordinates maps to different pixel heights per column (e.g. chr11
    # looks larger than chr12 when column 0 stacks chr1+chr12 but column 10 stacks
    # chr11+chr22 with a shorter ymax).
    global_ymax = max(col_ymax) if col_ymax else max_row_offset + 1
    provisional_xmax = fixed_global_xmax

    for col in range(ncols):
        axes[col].set_xlim(global_xmin, provisional_xmax)
        axes[col].set_ylim(global_ymin, global_ymax)
        axes[col].invert_yaxis()

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
    )
    global_xmax = max(global_xmax, fixed_global_xmax)

    for col in range(ncols):
        axes[col].set_xlim(global_xmin, global_xmax)

    fig.canvas.draw()

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
        )

    global_ymax = max(
        global_ymax,
        max(
            getattr(axes[col], "_phenogram_ymax", global_ymax)
            for col in range(ncols)
        ),
    )
    for col in range(ncols):
        axes[col].set_ylim(global_ymin, global_ymax)

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

    fig.canvas.draw()

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


def _plot_chr_body(
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
    log: Log = Log(),
    verbose: bool = True,
) -> Dict[str, Any]:
    """Pass 1: draw chromosome geometry and build lead y-positions (no annotation layout)."""
    chr_center_x = chr_x + chr_width / 2
    chr_right_x = chr_x + chr_width
    chr_inner_x = chr_x + chr_width * 0.05
    chr_inner_width = chr_width * 0.90
    anno_col_x = _anno_column_x(chr_right_x, anno_x_pad)

    positions = [
        0 + offset,
        chr_centromere_u / max_chr_size + offset,
        chr_centromere_l / max_chr_size + offset,
        chr_size / max_chr_size + offset,
    ]

    centromere_full_length = (chr_centromere_l - chr_centromere_u) / max_chr_size
    telemere_full_length = 0.02

    height_for_arm1 = positions[1] - positions[0]
    height_for_arm2 = positions[3] - positions[2]
    full_length = height_for_arm1 + height_for_arm2 + centromere_full_length

    arms = [
        Rectangle((chr_x, positions[0]), width=chr_width, height=height_for_arm1),
        Rectangle(
            (chr_x, positions[1] + centromere_full_length),
            width=chr_width,
            height=height_for_arm2,
        ),
    ]
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

    ax.add_collection(PatchCollection(arms, facecolor="white", edgecolor="grey", zorder=100))
    ax.add_collection(PatchCollection(centromeres, facecolor="grey", edgecolor="black", zorder=99))
    ax.add_collection(PatchCollection(telemeres, facecolor="white", edgecolor="grey", zorder=1))

    for _, row in chr_cytobands.iterrows():
        if row["END"] <= chr_centromere_u:
            band_start = row["START"] / max_chr_size + offset
            band_end = row["END"] / max_chr_size + offset
        elif row["START"] >= chr_centromere_l:
            band_start = (
                (row["START"] - chr_centromere_l) / max_chr_size
                + height_for_arm1 + offset + centromere_full_length
            )
            band_end = (
                (row["END"] - chr_centromere_l) / max_chr_size
                + height_for_arm1 + offset + centromere_full_length
            )
        else:
            continue

        band_height = band_end - band_start
        if band_height <= 0:
            continue

        band = Rectangle(
            (chr_inner_x, band_start),
            width=chr_inner_width,
            height=band_height,
        )
        facecolor = row["COLOR"]
        if row["STAIN"] == "stalk":
            bands_pc = PatchCollection(
                [band], facecolor=facecolor, edgecolor=None, linewidths=0, zorder=102
            )
        else:
            bands_pc = PatchCollection(
                [band], facecolor=facecolor, edgecolor=None, linewidths=0, zorder=101
            )
        ax.add_collection(bands_pc)

    leads_with_y = _build_leads_with_y(
        leads,
        chr_centromere_u,
        chr_centromere_l,
        max_chr_size,
        offset,
        height_for_arm1,
        centromere_full_length,
    )

    ax.set_xticks(ticks=[])
    ax.set_yticks(ticks=[])

    return {
        "ax": ax,
        "chr_name": chr_name,
        "chr_x": chr_x,
        "chr_width": chr_width,
        "anno_x_pad": anno_x_pad,
        "chr_center_x": chr_center_x,
        "chr_right_x": chr_right_x,
        "anno_col_x": anno_col_x,
        "leads_with_y": leads_with_y,
        "body_top_y": positions[0] - telemere_full_length / 2.0,
    }


def _plot_chr_annotations(
    spec: Dict[str, Any],
    anno_kwargs: Optional[Dict[str, Any]] = None,
    anno_kwargs_single: Optional[Dict[str, Any]] = None,
    arrow_kwargs: Optional[Dict[str, Any]] = None,
    marker_shape_map: Optional[Dict[str, str]] = None,
    marker_color_map: Optional[Dict[str, str]] = None,
    marker_size: float = 81,
    marker_gap_pt: float = 3,
    marker_label_gap_pt: float = 2.75,
    marker_max_per_row: int = 5,
    marker_row_gap_pt: float = 2.5,
    marker_label_align: str = "center",
    marker_fontsize: float = 11,
    marker_linewidth: float = 0.6,
    marker_label_bbox: bool = True,
    group_min_vertical_gap_pt: float = 3.5,
    group_marker_to_marker_gap_pt: float = 3,
    group_label_box_pad_pt: float = 1.5,
    marker_mode: bool = False,
    repel_force: float = 0.5,
    anno_max_iter: int = 300,
) -> None:
    """Pass 2: layout and draw annotations after axis limits are frozen."""
    if anno_kwargs is None:
        anno_kwargs = {}
    if anno_kwargs_single is None:
        anno_kwargs_single = {}
    if arrow_kwargs is None:
        arrow_kwargs = {}
    if marker_shape_map is None:
        marker_shape_map = {}
    if marker_color_map is None:
        marker_color_map = {}

    ax = spec["ax"]
    chr_x = spec["chr_x"]
    chr_right_x = spec["chr_right_x"]
    anno_col_x = spec["anno_col_x"]
    leads_with_y = spec.get("leads_with_y") or []

    if not leads_with_y:
        return

    max_xlim = anno_col_x + 0.5

    if marker_mode:
        text_kwargs = _build_marker_text_kwargs(
            anno_kwargs,
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

        spread_above_pt = 0.0
        spread_below_pt = 0.0
        for g in group_representatives:
            n_rows = math.ceil(len(g["group_items"]) / max(1, marker_max_per_row))
            above_pt, below_pt = _marker_group_block_extents_pt(
                marker_size,
                marker_label_gap_pt,
                marker_fontsize,
                group_label_box_pad_pt,
                n_marker_rows=n_rows,
                marker_row_gap_pt=marker_row_gap_pt,
                marker_linewidth=marker_linewidth,
            )
            spread_above_pt = max(spread_above_pt, above_pt)
            spread_below_pt = max(spread_below_pt, below_pt)

        orig_y = np.array([g["original_y"] for g in group_representatives])
        spread_y = _spread_group_blocks_display(
            ax,
            orig_y.tolist(),
            ref_x=chr_right_x,
            above_pt=spread_above_pt,
            below_pt=spread_below_pt,
            gap_pt=max(group_min_vertical_gap_pt, group_marker_to_marker_gap_pt),
            max_iter=anno_max_iter,
        )

        for i, g in enumerate(group_representatives):
            g["marker_y"] = float(spread_y[i])

        for g in group_representatives:
            group_key = g["group_key"]
            group_items = g["group_items"]
            anchor_y = float(g["marker_y"])
            label_text = _get_marker_group_label(group_key, group_items)
            if len(str(label_text)) > 20:
                label_text = str(label_text)[:17] + "..."

            main_item = group_items[0]
            main_y = main_item["original_y"]

            marker_positions, bottom_row_y = _compute_group_marker_positions(
                ax,
                anno_col_x,
                anchor_y,
                group_items,
                marker_size,
                marker_gap_pt,
                marker_max_per_row,
                marker_row_gap_pt,
                marker_linewidth=marker_linewidth,
            )
            for marker_x, marker_y, item in marker_positions:
                marker, color = _get_marker_style(
                    item["lead"], marker_shape_map, marker_color_map
                )
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

            if not marker_positions:
                continue

            label_x = anno_col_x
            group_text_kwargs = _build_lead_text_kwargs(
                text_kwargs,
                anno_kwargs_single,
                main_item["lead"],
                force_ha="left",
                force_va="top",
            )
            label_x, label_y = _marker_label_position(
                ax,
                label_x,
                bottom_row_y,
                marker_size,
                marker_label_gap_pt,
                marker_linewidth=marker_linewidth,
            )

            _plot_phenogram_connector(
                ax,
                chr_right_x,
                anno_col_x,
                main_y,
                anchor_y,
                color="black",
                linewidth=0.7,
                zorder=190,
            )

            ax.text(
                label_x,
                label_y,
                str(label_text),
                clip_on=False,
                zorder=230,
                **group_text_kwargs,
            )

            max_marker_x = max(x for x, _, _ in marker_positions)
            max_xlim = max(max_xlim, max_marker_x)
            group_anno_fontsize = float(
                group_text_kwargs.get("fontsize", marker_fontsize)
            )
            text_width_dx, _ = _points_to_data_delta(
                ax,
                label_x,
                label_y,
                dx_points=group_anno_fontsize * len(str(label_text)) * 0.55,
            )
            ha = group_text_kwargs.get("ha", "left")
            if ha == "center":
                text_right_x = label_x + text_width_dx / 2
            elif ha == "left":
                text_right_x = label_x + text_width_dx
            else:
                text_right_x = label_x
            max_xlim = max(max_xlim, text_right_x, max_marker_x)

            _, block_below_pt = _marker_group_block_extents_pt(
                marker_size,
                marker_label_gap_pt,
                marker_fontsize,
                group_label_box_pad_pt,
                n_marker_rows=math.ceil(
                    len(group_items) / max(1, marker_max_per_row)
                ),
                marker_row_gap_pt=marker_row_gap_pt,
                marker_linewidth=marker_linewidth,
            )
            _, label_bottom_y = _data_point_below(
                ax, chr_right_x, bottom_row_y, block_below_pt + 4.0
            )
            ax._phenogram_ymax = max(
                getattr(ax, "_phenogram_ymax", 0),
                label_bottom_y + 0.01,
            )
    else:
        text_kwargs = _build_text_mode_kwargs(anno_kwargs)
        anno_fontsize = float(text_kwargs.get("fontsize", 9))

        for l in leads_with_y:
            original_y = l["original_y"]
            ax.plot(
                [chr_x, chr_right_x],
                [original_y, original_y],
                color="red",
                linewidth=2,
                clip_on=False,
                zorder=200,
            )

        labelled_items = sorted(
            [l for l in leads_with_y if l["lead"].get("label_text")],
            key=lambda item: item["y_pos_chr"],
        )
        if len(labelled_items) > 0:
            orig_y = [l["y_pos_chr"] for l in labelled_items]
            label_gap_pt = max(
                4.0,
                group_min_vertical_gap_pt * max(float(repel_force), 0.05),
            )
            spread_y = _spread_phenogram_text_labels(
                ax,
                orig_y,
                ref_x=chr_right_x,
                anno_fontsize=anno_fontsize,
                gap_pt=label_gap_pt,
                group_label_box_pad_pt=group_label_box_pad_pt,
                max_iter=anno_max_iter,
            )
            for i, l in enumerate(labelled_items):
                l["text_y"] = float(spread_y[i])

            _, below_pt = _text_label_block_extents_pt(
                anno_fontsize, group_label_box_pad_pt
            )
            for l in labelled_items:
                _, label_bottom_y = _data_point_below(
                    ax,
                    chr_right_x,
                    l["text_y"],
                    below_pt + 2.0,
                )
                ax._phenogram_ymax = max(
                    getattr(ax, "_phenogram_ymax", 0),
                    label_bottom_y + 0.01,
                )

        for l in leads_with_y:
            lead = l["lead"]
            original_y = l["original_y"]
            label_text = lead.get("label_text")
            if not label_text:
                continue
            if len(str(label_text)) > 20:
                label_text = str(label_text)[:17] + "..."

            text_y = l.get("text_y", l["y_pos_chr"])
            text_x = anno_col_x
            lead_text_kwargs = _build_lead_text_kwargs(
                text_kwargs,
                anno_kwargs_single,
                lead,
            )
            lead_fontsize = float(lead_text_kwargs.get("fontsize", anno_fontsize))
            (
                lead_arrow_shaft_pt,
                lead_arrow_pad_pt,
                lead_shrink_b_pt,
                lead_arrowprops,
                lead_line_kwargs,
            ) = _resolve_lead_connector_style(
                anno_kwargs,
                arrow_kwargs,
                anno_kwargs_single,
                lead.get("snpid"),
            )

            text_width_dx, _ = _points_to_data_delta(
                ax,
                text_x,
                text_y,
                dx_points=lead_fontsize * len(str(label_text)) * 0.55,
            )
            max_xlim = max(max_xlim, text_x + text_width_dx)

            _plot_phenogram_text_connector(
                ax,
                chr_right_x,
                original_y,
                text_x,
                text_y,
                arrow_shaft_pt=lead_arrow_shaft_pt,
                arrow_pad_pt=lead_arrow_pad_pt,
                shrink_b_pt=lead_shrink_b_pt,
                arrowprops=lead_arrowprops,
                line_kwargs=lead_line_kwargs,
            )
            ax.text(
                text_x,
                text_y,
                str(label_text),
                clip_on=False,
                zorder=202,
                **lead_text_kwargs,
            )

    xlim_default = chr_right_x + spec["anno_x_pad"] + 0.5
    ax._phenogram_xmax = max(
        getattr(ax, "_phenogram_xmax", xlim_default),
        max_xlim * 1.05,
    )
