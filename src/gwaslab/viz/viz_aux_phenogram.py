"""Phenogram plot helpers (ideogram, layout, labels, markers, drawing).
"""

import math
import textwrap
from functools import lru_cache
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.text import Text as MplText
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.collections import PatchCollection
from typing import Optional, Union, Dict, Any, List, Tuple, Sequence
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _anno_gene


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

_DEFAULT_CHR_WIDTH = 0.35
_MIN_LABEL_RESERVE = 0.35
_MIN_CHR_X_FRACTION = 0.25
_CHR_OUTLINE_LW = 0.8
_TELOMERE_LENGTH = 0.02
_CHR_OUTLINE_COLOR = "grey"
_CHR_OUTLINE_ARC_STEPS = 48
_DEFAULT_ANNO_WRAP_CHARS_PER_LINE = 9
_DEFAULT_ANNO_FONTSIZE = 11
_DEFAULT_MARKER_MAX_PER_ROW = 4
_LEGEND_MARKER_SIZE = 11.0
_LEGEND_FONTSIZE = 12.0
_PHENOGRAM_INTER_ROW_GAP = 0.12
_PANEL_FLOOR_PAD_FRACTION = 0.35
_ANNOTATION_Y_MARGIN_DATA = 0.03
_ANNOTATION_PATH_PAD_PT = 3.0
_CYTOBAND_COLORS = {
    "gpos100": (100 / 255, 100 / 255, 100 / 255),
    "gpos": (100 / 255, 100 / 255, 100 / 255),
    "gpos75": (110 / 255, 110 / 255, 110 / 255),
    "gpos66": (130 / 255, 130 / 255, 130 / 255),
    "gpos50": (160 / 255, 160 / 255, 160 / 255),
    "gpos33": (180 / 255, 180 / 255, 180 / 255),
    "gpos25": (200 / 255, 200 / 255, 200 / 255),
    "gvar": (220 / 255, 220 / 255, 220 / 255),
    "gneg": (255 / 255, 255 / 255, 255 / 255),
    "acen": (217 / 255, 217 / 255, 217 / 255),
    "stalk": (100 / 255, 127 / 255, 164 / 255),
}


def _resolve_phenogram_cytoband_path(
    build: str,
    cytoband_path: Optional[str],
    log: Log = Log(),
    verbose: bool = True,
) -> Path:
    """Resolve user/default cytoband path for phenogram plotting.
"""
    if cytoband_path is not None:
        return Path(cytoband_path)

    data_dir = Path(__file__).parent.parent / "data" / "cytoband"
    if build == "38":
        return data_dir / "cytoBand_hg38.txt.gz"
    if build != "19":
        log.write(" -Unknown build '{}', using hg19 cytoband data.".format(build), verbose=verbose)
    return data_dir / "cytoBand_hg19.txt.gz"


@lru_cache(maxsize=8)
def _load_phenogram_cytobands(cytoband_path: str) -> pd.DataFrame:
    """Load and color-map cytobands; cached because packaged files are static.
"""
    cytobands = pd.read_csv(cytoband_path, sep=r"\s+", header=None, compression="infer")
    cytobands.columns = ["CHR", "START", "END", "ARM", "STAIN"]
    cytobands["COLOR"] = cytobands["STAIN"].map(_CYTOBAND_COLORS)
    return cytobands


def _phenogram_chr_sort_key(chr_name: str) -> Tuple[int, str]:
    if not str(chr_name).startswith("chr"):
        return (999, str(chr_name))
    chr_num_str = str(chr_name)[3:]
    if chr_num_str.isdigit():
        return (int(chr_num_str), str(chr_name))
    if chr_num_str == "X":
        return (23, str(chr_name))
    if chr_num_str == "Y":
        return (24, str(chr_name))
    if chr_num_str == "MT":
        return (25, str(chr_name))
    return (999, str(chr_name))


def _summarize_phenogram_cytobands(
    cytobands: pd.DataFrame,
    include_sex_chr: bool = False,
) -> Tuple[List[str], Dict[str, int], Dict[str, Tuple[float, float]], Dict[str, pd.DataFrame]]:
    """Precompute chromosome order, sizes, centromeres, and per-chromosome bands.
"""
    chr_sizes = cytobands.groupby("CHR", sort=False)["END"].max().to_dict()
    chr_list = [
        chr_name
        for chr_name in sorted(chr_sizes, key=_phenogram_chr_sort_key)
        if str(chr_name).startswith("chr")
        and str(chr_name)[3:].isdigit()
        and 1 <= int(str(chr_name)[3:]) <= 22
    ]
    if include_sex_chr:
        chr_list.extend(chr_name for chr_name in ("chrX", "chrY") if chr_name in chr_sizes)
    cytobands_by_chr = {
        chr_name: chr_bands
        for chr_name, chr_bands in cytobands.groupby("CHR", sort=False)
        if chr_name in chr_list
    }
    centromere_by_chr: Dict[str, Tuple[float, float]] = {}
    for chr_name in chr_list:
        chr_bands = cytobands_by_chr[chr_name]
        acen_bands = chr_bands.loc[chr_bands["STAIN"] == "acen", :]
        if len(acen_bands) > 0:
            centromere_by_chr[chr_name] = (
                float(acen_bands["START"].min()),
                float(acen_bands["END"].max()),
            )
        else:
            chr_size = float(chr_sizes[chr_name])
            centromere_by_chr[chr_name] = (chr_size * 0.44, chr_size * 0.46)
    return chr_list, chr_sizes, centromere_by_chr, cytobands_by_chr


# --- Constants and ideogram geometry ---


@dataclass(frozen=True)
class _ChrIdeogramGeometry:
    """Single source of truth for chromosome ideogram layout (data coordinates).
"""

    left_x: float
    right_x: float
    width: float
    center_x: float
    y_arm1_top: float
    y_cent_top: float
    y_cent_bottom: float
    y_arm2_bottom: float
    telomere_h: float
    centromere_h: float
    centromere_marker_h: float

    @property
    def arm1_h(self) -> float:
        return self.y_cent_top - self.y_arm1_top

    @property
    def arm2_h(self) -> float:
        return self.y_arm2_bottom - self.y_cent_bottom

    @property
    def ideogram_top(self) -> float:
        return self.y_arm1_top + self.telomere_h / 2.0

    @property
    def ideogram_bottom(self) -> float:
        return self.y_arm2_bottom - self.telomere_h / 2.0

    @property
    def y_cent_mid(self) -> float:
        return (self.y_cent_top + self.y_cent_bottom) / 2.0


def _compute_chr_ideogram_geometry(
    chr_x: float,
    chr_width: float,
    chr_size: float,
    chr_centromere_u: float,
    chr_centromere_l: float,
    max_chr_size: float,
    offset: float,
    telomere_h: float = _TELOMERE_LENGTH,
    centromere_marker_h: Optional[float] = None,
) -> _ChrIdeogramGeometry:
    """Derive all ideogram bounds from one geometry definition (y increases upward).
"""
    chr_extent = chr_size / max_chr_size
    y_arm2_bottom = offset + telomere_h / 2.0
    y_arm1_top = offset + chr_extent + telomere_h / 2.0
    centromere_mid_bp = (float(chr_centromere_u) + float(chr_centromere_l)) / 2.0
    y_cent_mid = y_arm1_top - centromere_mid_bp / max_chr_size
    marker_h = float(centromere_marker_h if centromere_marker_h is not None else telomere_h)
    y_cent_top = y_cent_mid + marker_h / 2.0
    y_cent_bottom = y_cent_mid - marker_h / 2.0
    return _ChrIdeogramGeometry(
        left_x=chr_x,
        right_x=chr_x + chr_width,
        width=chr_width,
        center_x=chr_x + chr_width / 2.0,
        y_arm1_top=y_arm1_top,
        y_cent_top=y_cent_top,
        y_cent_bottom=y_cent_bottom,
        y_arm2_bottom=y_arm2_bottom,
        telomere_h=telomere_h,
        centromere_h=marker_h,
        centromere_marker_h=marker_h,
    )


def _semicircle_arc(
    center_x: float,
    base_y: float,
    half_width: float,
    half_height: float,
    upper: bool,
    n_steps: int = _CHR_OUTLINE_ARC_STEPS,
) -> Tuple[np.ndarray, np.ndarray]:
    """Sample a semicircular arc with flat edge at ``base_y`` (y increases upward).
"""
    if upper:
        theta = np.linspace(0.0, np.pi, n_steps)
        ys = base_y + half_height * np.sin(theta)
    else:
        theta = np.linspace(np.pi, 2.0 * np.pi, n_steps)
        ys = base_y + half_height * np.sin(theta)
    xs = center_x + half_width * np.cos(theta)
    return xs, ys


def _chr_ideogram_body_clip_path(geom: _ChrIdeogramGeometry) -> MplPath:
    """Closed outer boundary used for clipping chromosome fills.
"""
    hw = geom.width / 2.0
    hh = geom.telomere_h / 2.0

    top_x, top_y = _semicircle_arc(
        geom.center_x, geom.y_arm1_top, hw, hh, upper=True
    )
    bot_x, bot_y = _semicircle_arc(
        geom.center_x, geom.y_arm2_bottom, hw, hh, upper=False
    )

    verts: List[Tuple[float, float]] = []
    codes: List[int] = []

    for i, (x, y) in enumerate(zip(top_x, top_y)):
        verts.append((float(x), float(y)))
        codes.append(MplPath.MOVETO if i == 0 else MplPath.LINETO)

    for x, y in (
        (geom.right_x, geom.y_arm1_top),
        (geom.right_x, geom.y_arm2_bottom),
    ):
        verts.append((float(x), float(y)))
        codes.append(MplPath.LINETO)

    for x, y in zip(bot_x, bot_y):
        verts.append((float(x), float(y)))
        codes.append(MplPath.LINETO)

    for x, y in (
        (geom.left_x, geom.y_arm2_bottom),
        (geom.left_x, geom.y_arm1_top),
    ):
        verts.append((float(x), float(y)))
        codes.append(MplPath.LINETO)

    codes.append(MplPath.CLOSEPOLY)
    verts.append((0.0, 0.0))

    return MplPath(verts, codes)


def _chr_ideogram_outline_path(geom: _ChrIdeogramGeometry) -> MplPath:
    """Visible ideogram outline with side gaps at the centromere marker.
"""
    hw = geom.width / 2.0
    hh = geom.telomere_h / 2.0

    top_x, top_y = _semicircle_arc(
        geom.center_x, geom.y_arm1_top, hw, hh, upper=True
    )
    bot_x, bot_y = _semicircle_arc(
        geom.center_x, geom.y_arm2_bottom, hw, hh, upper=False
    )

    verts: List[Tuple[float, float]] = []
    codes: List[int] = []

    for i, (x, y) in enumerate(zip(top_x, top_y)):
        verts.append((float(x), float(y)))
        codes.append(MplPath.MOVETO if i == 0 else MplPath.LINETO)
    verts.append((geom.left_x, geom.y_cent_top))
    codes.append(MplPath.LINETO)

    verts.append((geom.left_x, geom.y_cent_bottom))
    codes.append(MplPath.MOVETO)
    verts.append((geom.left_x, geom.y_arm2_bottom))
    codes.append(MplPath.LINETO)

    for x, y in zip(bot_x, bot_y):
        verts.append((float(x), float(y)))
        codes.append(MplPath.LINETO)

    verts.append((geom.right_x, geom.y_cent_bottom))
    codes.append(MplPath.LINETO)

    verts.append((geom.right_x, geom.y_cent_top))
    codes.append(MplPath.MOVETO)
    verts.append((geom.right_x, geom.y_arm1_top))
    codes.append(MplPath.LINETO)

    return MplPath(verts, codes)


def _chr_ideogram_fill_patches(
    geom: _ChrIdeogramGeometry,
) -> Tuple[List[Any], List[Any], List[Any]]:
    """Return arm, centromere, and telomere fill patches (no stroke).
"""
    arms = [
        Rectangle(
            (geom.left_x, geom.y_arm2_bottom),
            geom.width,
            geom.y_arm1_top - geom.y_arm2_bottom,
        ),
    ]
    telemeres: List[Any] = []
    return arms, [], telemeres


def _add_chr_fill_collection(
    ax,
    patches: List[Any],
    facecolor: Any,
    zorder: int,
) -> PatchCollection:
    """Add ideogram fill patches with no edge stroke.
"""
    pc = PatchCollection(
        patches,
        facecolor=facecolor,
        edgecolor="none",
        linewidths=0.0,
        antialiased=False,
        zorder=zorder,
    )
    pc.set_clip_on(False)
    ax.add_collection(pc)
    return pc


def _plot_chr_ideogram_outline(
    ax,
    geom: _ChrIdeogramGeometry,
    color: str = _CHR_OUTLINE_COLOR,
    linewidth: float = _CHR_OUTLINE_LW,
    zorder: int = 104,
) -> PathPatch:
    """Draw one unified outer stroke from the shared geometry path.
"""
    outline = PathPatch(
        _chr_ideogram_outline_path(geom),
        facecolor="none",
        edgecolor=color,
        linewidth=linewidth,
        antialiased=True,
        clip_on=False,
        zorder=zorder,
        joinstyle="miter",
        capstyle="butt",
    )
    ax.add_patch(outline)
    return outline


def _telomere_cap_path(
    geom: _ChrIdeogramGeometry,
    upper: bool,
) -> MplPath:
    """Filled semicircle cap with flat edge on the arm base.
"""
    hw = geom.width / 2.0
    hh = geom.telomere_h / 2.0
    base_y = geom.y_arm1_top if upper else geom.y_arm2_bottom
    xs, ys = _semicircle_arc(geom.center_x, base_y, hw, hh, upper=upper)
    verts = [(float(x), float(y)) for x, y in zip(xs, ys)]
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(verts) - 1)
    codes.append(MplPath.CLOSEPOLY)
    verts.append((0.0, 0.0))
    return MplPath(verts, codes)


def _add_telomere_cap(
    ax,
    geom: _ChrIdeogramGeometry,
    upper: bool,
    facecolor: str,
    zorder: int,
) -> PathPatch:
    cap = PathPatch(
        _telomere_cap_path(geom, upper=upper),
        facecolor=facecolor,
        edgecolor="none",
        linewidth=0.0,
        antialiased=False,
        clip_on=False,
        zorder=zorder,
    )
    ax.add_patch(cap)
    return cap


def _centromere_hourglass_path(geom: _ChrIdeogramGeometry) -> MplPath:
    """Fixed-size hourglass centromere decoration.

    The underlying chromosome body remains continuous and position-bearing; this
    path only marks the centromere location visually.
"""
    mid_y = geom.y_cent_mid
    verts = [
        (geom.left_x, geom.y_cent_top),
        (geom.right_x, geom.y_cent_top),
        (geom.center_x, mid_y),
        (geom.left_x, geom.y_cent_top),
        (geom.left_x, geom.y_cent_bottom),
        (geom.right_x, geom.y_cent_bottom),
        (geom.center_x, mid_y),
        (geom.left_x, geom.y_cent_bottom),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.LINETO,
        MplPath.CLOSEPOLY,
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.LINETO,
        MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def _centromere_fill_patches(geom: _ChrIdeogramGeometry) -> List[Rectangle]:
    """Centromere is rendered by fixed triangles, not a proportional band.
"""
    del geom
    return []


def _add_centromere_hourglass(
    ax,
    geom: _ChrIdeogramGeometry,
    facecolor: str,
    zorder: int,
    edgecolor: str = "grey",
    linewidth: float = 0.4,
) -> PathPatch:
    hourglass = PathPatch(
        _centromere_hourglass_path(geom),
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        antialiased=False,
        clip_on=False,
        zorder=zorder,
    )
    ax.add_patch(hourglass)
    return hourglass


def _add_centromere_decoration(
    ax,
    geom: _ChrIdeogramGeometry,
    facecolor: str = "grey",
    zorder: int = 103,
) -> PathPatch:
    """Draw one hourglass centromere decoration patch.
"""
    return _add_centromere_hourglass(
        ax, geom, facecolor=facecolor, zorder=zorder
    )
    return top, bottom


def _cytoband_y_range(
    row: pd.Series,
    geom: _ChrIdeogramGeometry,
    chr_centromere_u: float,
    chr_centromere_l: float,
    max_chr_size: float,
    offset: float,
) -> Optional[Tuple[float, float]]:
    """Map a cytoband row to data-y start/end; acen is shown by fixed triangles.
"""
    del offset, chr_centromere_u, chr_centromere_l
    if row["STAIN"] == "acen":
        return None
    band_start = geom.y_arm1_top - row["END"] / max_chr_size
    band_end = geom.y_arm1_top - row["START"] / max_chr_size
    if band_end <= band_start:
        return None
    return band_start, band_end


def _points_to_data_delta(ax, x, y, dx_points=0.0, dy_points=0.0):
    """Convert a display offset in matplotlib points to a data-coordinate delta at (x, y).
"""
    scale = ax.figure.dpi / 72.0
    x_disp, y_disp = ax.transData.transform((x, y))
    x_data, y_data = ax.transData.inverted().transform(
        (x_disp + dx_points * scale, y_disp + dy_points * scale)
    )
    return x_data - x, y_data - y


# --- Annotation kwargs and label resolution ---


def _resolve_anno_arrow_pt(
    kwargs_key: str,
    anno_kwargs: Dict[str, Any],
    default: float,
) -> float:
    """Resolve arrow spacing in points from ``anno_kwargs`` or default.
"""
    if kwargs_key in anno_kwargs:
        return float(anno_kwargs[kwargs_key])
    return float(default)


def _merge_lead_anno_kwargs(
    anno_kwargs: Dict[str, Any],
    anno_kwargs_single: Optional[Dict[str, Any]],
    snpid: Optional[str],
) -> Dict[str, Any]:
    """Merge global and per-SNP annotation kwargs (MQQ-compatible).
"""
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
    """Apply ``anno_kwargs_single`` overrides for one lead.
"""
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
    """Resolve connector point lengths and matplotlib arrow/line styling.
"""
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
    """Return (sort_column, ascending) for ranking variants by significance.
"""
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
    """Boolean mask of lead rows that may receive text labels.
"""
    if len(leads) == 0:
        return pd.Series(dtype=bool)

    mask = pd.Series(True, index=leads.index)
    if anno_set is not None and len(anno_set) > 0:
        if snpid in leads.columns:
            mask = leads[snpid].isin(anno_set)
        else:
            mask = pd.Series(False, index=leads.index)
        n_matched = int(mask.sum())
        n_requested = len(anno_set)
        if n_matched == 0:
            log.warning(
                " -anno_set: 0 of {} SNPID(s) matched extracted leads. "
                "Verify SNPID strings (not gene names); all leads still get bars.".format(
                    n_requested
                ),
                verbose=verbose,
            )
        elif n_matched < n_requested:
            log.write(
                " -anno_set: found {}/{} SNPID(s) for labeling.".format(
                    n_matched, n_requested
                ),
                verbose=verbose,
            )
        else:
            log.write(
                " -anno_set: found {}/{} SNPID(s) for labeling.".format(
                    n_matched, n_requested
                ),
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
    """Add an ``Annotation`` column for column-based or GENENAME labels.
"""
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
    """Resolve display label for a phenogram lead (MQQ-compatible ``anno`` rules).
"""
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


def _phenogram_chr_string(chr_val: Any) -> Optional[str]:
    if pd.isna(chr_val):
        return None
    if isinstance(chr_val, (int, float)):
        if int(chr_val) == 23:
            return "chrX"
        if int(chr_val) == 24:
            return "chrY"
        return "chr{}".format(int(chr_val))
    chr_str = str(chr_val)
    if chr_str.startswith("chr"):
        return chr_str
    try:
        return "chr{}".format(int(chr_str))
    except (ValueError, TypeError):
        return "chr{}".format(chr_str)


def _phenogram_chr_display_value(chr_val: Any, chr_str: str) -> Any:
    if isinstance(chr_val, (int, float)):
        return int(chr_val)
    return chr_str.replace("chr", "")


def _notna_or_none(value: Any) -> Any:
    if pd.isna(value):
        return None
    return value


def _prepare_phenogram_leads_dict(
    leads: pd.DataFrame,
    chrom: str,
    pos: str,
    snpid: str,
    chr_list: List[str],
    marker_mode: bool,
    anno: Optional[Union[bool, str]],
    label_eligible: pd.Series,
    anno_alias: Optional[Dict[str, str]],
    anno_group: Optional[str],
    anno_shape: Optional[str],
    anno_color: Optional[str],
) -> Dict[str, List[Dict[str, Any]]]:
    """Normalize lead rows into per-chromosome dictionaries used by plot layout.
"""
    if len(leads) == 0 or chrom not in leads.columns or pos not in leads.columns:
        return {}

    valid_chr_values = set(chr_list)
    valid_chr_values.update({"chrX", "chrY", "chrMT"})

    chr_strings = leads[chrom].map(_phenogram_chr_string)
    valid_mask = chr_strings.notna() & leads[pos].notna()
    valid_mask &= chr_strings.isin(valid_chr_values)
    if not valid_mask.any():
        return {}

    prepared = pd.DataFrame(
        {
            "lead_index": leads.index,
            "chr_str": chr_strings,
            "chr_val": leads[chrom],
            "pos_value": leads[pos],
            "snpid_value": leads[snpid] if snpid in leads.columns else None,
            "annotation_value": leads["Annotation"] if "Annotation" in leads.columns else None,
            "anno_group_value": leads[anno_group] if anno_group in leads.columns else None,
            "anno_shape_value": leads[anno_shape] if anno_shape in leads.columns else None,
            "anno_color_value": leads[anno_color] if anno_color in leads.columns else None,
        },
        index=leads.index,
    ).loc[valid_mask]

    label_lookup = set(label_eligible.index[label_eligible]) if len(label_eligible) else set()
    leads_dict: Dict[str, List[Dict[str, Any]]] = {}
    for row in prepared.itertuples(index=False):
        lead_entry = {
            "pos": int(row.pos_value),
            "chrom": _phenogram_chr_display_value(row.chr_val, row.chr_str),
            "snpid": _notna_or_none(row.snpid_value),
            "index": row.lead_index,
            "annotation": _notna_or_none(row.annotation_value),
            "anno_group": _notna_or_none(row.anno_group_value),
            "anno_shape": _notna_or_none(row.anno_shape_value),
            "anno_color": _notna_or_none(row.anno_color_value),
        }
        if not marker_mode and anno is not None and row.lead_index in label_lookup:
            lead_entry["label_text"] = _resolve_phenogram_label(
                lead_entry, anno, anno_alias
            )
        else:
            lead_entry["label_text"] = None
        leads_dict.setdefault(row.chr_str, []).append(lead_entry)
    return leads_dict


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


def _resolve_marker_group_label(
    group_key: str,
    group_items: List[Dict[str, Any]],
    anno: Optional[Union[bool, str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
) -> str:
    """Resolve marker-mode group label text, honoring ``anno`` when set.
"""
    if anno is not None:
        main_lead = group_items[0]["lead"]
        resolved = _resolve_phenogram_label(main_lead, anno, anno_alias)
        if resolved is not None:
            return resolved
    return _get_marker_group_label(group_key, group_items)


def _prepare_marker_group_display_label(
    group_key: str,
    group_items: List[Dict[str, Any]],
    *,
    anno: Optional[Union[bool, str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
    fontproperties: Optional[FontProperties] = None,
    renderer=None,
) -> Tuple[str, int]:
    """Wrap or truncate a marker-group label; return display text and line count.
"""
    raw = _resolve_marker_group_label(
        group_key, group_items, anno=anno, anno_alias=anno_alias
    )
    if fontproperties is None:
        fontproperties = FontProperties(size=11)
    if renderer is not None:
        display, n_lines, _ = _prepare_phenogram_display_label(
            raw,
            anno_wrap=anno_wrap,
            anno_wrap_width_pt=anno_wrap_width_pt,
            anno_wrap_chars_per_line=anno_wrap_chars_per_line,
            anno_max_len=anno_max_len,
            fontproperties=fontproperties,
            renderer=renderer,
        )
        return display, n_lines
    if anno_wrap and anno_wrap_chars_per_line and anno_wrap_chars_per_line > 0:
        display = _wrap_phenogram_label_chars(raw, anno_wrap_chars_per_line)
    else:
        display = _format_phenogram_label(
            raw,
            max_len=anno_max_len if anno_max_len is not None else 20,
        )
    return display, max(1, display.count("\n") + 1)


# --- Marker mode ---


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
    """Build value -> marker shape/color maps from data unique values or user overrides.
"""
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
        "fontsize": _DEFAULT_ANNO_FONTSIZE,
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
    n_lines: int = 1,
) -> Tuple[float, float]:
    """Vertical half-extents of a text label block (``va='center'``) in points.
"""
    half_text_pt = (
        _phenogram_line_height_pt(anno_fontsize) * max(1, int(n_lines)) / 2.0
        + group_label_box_pad_pt
    )
    return half_text_pt, half_text_pt


@dataclass
# --- Annotation layout solver ---


class _AnnotationLayoutBlock:
    """One rigid annotation unit for bounded 1D vertical layout (data coordinates).
"""

    block_id: int
    target_y: float
    height: float
    y: float
    min_gap: float
    panel_min: float
    panel_max: float
    original_y: float = 0.0
    n_marker_rows: int = 1
    n_text_lines: int = 1
    group_key: str = ""
    group_items: Optional[List[Dict[str, Any]]] = None
    display_label: str = ""
    above_pt: float = 0.0
    below_pt: float = 0.0
    is_marker: bool = False
    lead_item: Optional[Dict[str, Any]] = None


def _pt_to_data_height(ax, ref_x: float, ref_y: float, pt: float) -> float:
    """Convert a vertical extent in points to data-y delta at ``(ref_x, ref_y)``.
"""
    _, dy = _points_to_data_delta(ax, ref_x, ref_y, dy_points=float(pt))
    return abs(float(dy))


def _marker_block_center_from_anchor(
    ax,
    ref_x: float,
    anchor_y: float,
    above_pt: float,
    below_pt: float,
) -> Tuple[float, float]:
    """Return (block_center_y, block_height_data) for a marker unit anchored at the top marker row.
"""
    bias_pt = (float(below_pt) - float(above_pt)) / 2.0
    _, bias_dy = _points_to_data_delta(ax, ref_x, anchor_y, dy_points=-bias_pt)
    height = _pt_to_data_height(ax, ref_x, anchor_y, above_pt + below_pt)
    return float(anchor_y) + float(bias_dy), height


def _marker_anchor_from_block_center(
    ax,
    ref_x: float,
    center_y: float,
    above_pt: float,
    below_pt: float,
) -> float:
    """Top marker-row y from a solved block center.
"""
    bias_pt = (float(below_pt) - float(above_pt)) / 2.0
    _, bias_dy = _points_to_data_delta(ax, ref_x, center_y, dy_points=-bias_pt)
    return float(center_y) - float(bias_dy)


def _annotation_panel_limits(
    spec: Dict[str, Any],
    ideogram_geom: Optional[_ChrIdeogramGeometry] = None,
) -> Tuple[float, float]:
    """Return (panel_min, panel_max) data-y bounds for annotation blocks on one chromosome.
"""
    floor_y = spec.get("panel_floor_y")
    ceiling_y = spec.get("panel_ceiling_y")
    if ideogram_geom is not None:
        panel_min = float(floor_y) if floor_y is not None else ideogram_geom.ideogram_bottom
        if ceiling_y is not None:
            panel_max = float(ceiling_y)
        else:
            row_top = spec.get("row_ideogram_top")
            panel_max = (
                float(row_top)
                if row_top is not None
                else ideogram_geom.ideogram_top
            )
    else:
        panel_min = float(floor_y) if floor_y is not None else -1e9
        panel_max = float(ceiling_y) if ceiling_y is not None else 1e9
    if panel_min > panel_max:
        panel_min, panel_max = panel_max, panel_min
    return panel_min, panel_max


def _block_top(block: _AnnotationLayoutBlock) -> float:
    return float(block.y) + float(block.height) / 2.0


def _block_bottom(block: _AnnotationLayoutBlock) -> float:
    return float(block.y) - float(block.height) / 2.0


def _block_bottom_from_center(
    ax,
    ref_x: float,
    center_y: float,
    above_pt: float,
    below_pt: float,
) -> float:
    """Data-y of the bottom edge for a symmetric block centered at ``center_y``.
"""
    half_h = _pt_to_data_height(
        ax, ref_x, center_y, (float(above_pt) + float(below_pt)) / 2.0
    )
    return float(center_y) - half_h


def _block_top_from_center(
    ax,
    ref_x: float,
    center_y: float,
    above_pt: float,
    below_pt: float,
) -> float:
    """Data-y of the top edge for a symmetric block centered at ``center_y``.
"""
    half_h = _pt_to_data_height(
        ax, ref_x, center_y, (float(above_pt) + float(below_pt)) / 2.0
    )
    return float(center_y) + half_h


def _blocks_overlap(
    blocks: List[_AnnotationLayoutBlock],
    gap_tol: float = 1e-9,
) -> bool:
    for j in range(1, len(blocks)):
        if _block_bottom(blocks[j]) < _block_top(blocks[j - 1]) + blocks[j].min_gap - gap_tol:
            return True
    return False


def _blocks_inside_panel(blocks: List[_AnnotationLayoutBlock], tol: float = 1e-9) -> bool:
    if not blocks:
        return True
    panel_min = blocks[0].panel_min
    panel_max = blocks[0].panel_max
    for block in blocks:
        if _block_bottom(block) < panel_min - tol or _block_top(block) > panel_max + tol:
            return False
    return True


def _shift_blocks_y(blocks: List[_AnnotationLayoutBlock], delta: float) -> None:
    for block in blocks:
        block.y += float(delta)


def _solve_bounded_stack_layout(
    blocks: List[_AnnotationLayoutBlock],
    max_iterations: int = 300,
    min_gap_floor: float = 0.0,
) -> Tuple[List[_AnnotationLayoutBlock], bool]:
    """Lay out rigid annotation blocks on a vertical axis without overlap.

    Returns solved blocks and whether layout fits inside panel bounds without overlap.
"""
    if not blocks:
        return blocks, True

    panel_min = blocks[0].panel_min
    panel_max = blocks[0].panel_max
    for block in blocks:
        block.y = float(block.target_y)

    blocks.sort(key=lambda b: b.target_y)

    for _ in range(max_iterations):
        changed = False

        for j in range(1, len(blocks)):
            prev_blk = blocks[j - 1]
            curr_blk = blocks[j]
            overlap = (
                _block_top(prev_blk)
                + curr_blk.min_gap
                - _block_bottom(curr_blk)
            )
            if overlap > 0:
                shift = float(overlap)
                new_bottom = _block_bottom(curr_blk) + shift
                if new_bottom < panel_min:
                    curr_blk.y += shift
                    remainder = panel_min - _block_bottom(curr_blk)
                    if remainder > 0:
                        prev_blk.y += remainder
                else:
                    curr_blk.y += shift
                changed = True

        if blocks:
            top_excess = _block_top(blocks[-1]) - panel_max
            if top_excess > 0:
                _shift_blocks_y(blocks, -top_excess)
                changed = True
            bottom_excess = panel_min - _block_bottom(blocks[0])
            if bottom_excess > 0:
                _shift_blocks_y(blocks, bottom_excess)
                changed = True

        for block in blocks:
            old_y = block.y
            if block.y > block.target_y:
                block.y = block.target_y
                if _blocks_overlap(blocks) or not _blocks_inside_panel(blocks):
                    block.y = old_y
                elif block.y != old_y:
                    changed = True

        for j in range(1, len(blocks)):
            prev_blk = blocks[j - 1]
            curr_blk = blocks[j]
            overlap = (
                _block_top(prev_blk)
                + curr_blk.min_gap
                - _block_bottom(curr_blk)
            )
            if overlap > 0:
                curr_blk.y += overlap
                changed = True

        if blocks:
            top_excess = _block_top(blocks[-1]) - panel_max
            if top_excess > 0:
                _shift_blocks_y(blocks, -top_excess)
                changed = True
            bottom_excess = panel_min - _block_bottom(blocks[0])
            if bottom_excess > 0:
                _shift_blocks_y(blocks, bottom_excess)
                changed = True

        if not changed and not _blocks_overlap(blocks) and _blocks_inside_panel(blocks):
            break

    fits = not _blocks_overlap(blocks) and _blocks_inside_panel(blocks)
    if not fits and min_gap_floor > 0:
        reduced = max(float(min_gap_floor), 0.0)
        for block in blocks:
            block.min_gap = min(block.min_gap, reduced)
        return _solve_bounded_stack_layout(
            blocks,
            max_iterations=max_iterations,
            min_gap_floor=0.0,
        )
    return blocks, fits


def _gap_pt_to_data(
    ax,
    ref_x: float,
    ref_y: float,
    gap_pt: float,
) -> float:
    return _pt_to_data_height(ax, ref_x, ref_y, gap_pt)


def _place_marker_unit_at_center(
    ax,
    anno_col_x: float,
    block_center_y: float,
    group_items: List[Dict[str, Any]],
    marker_size: float,
    marker_gap_pt: float,
    marker_label_gap_pt: float,
    marker_max_per_row: int,
    marker_row_gap_pt: float,
    marker_linewidth: float,
    marker_fontsize: float,
    n_text_lines: int,
    group_label_box_pad_pt: float,
    marker_label_align: str,
) -> Dict[str, Any]:
    """Place marker cluster + label as one rigid unit from solved block center.
"""
    n_rows = math.ceil(len(group_items) / max(1, marker_max_per_row))
    above_pt, below_pt = _marker_group_block_extents_pt(
        marker_size,
        marker_label_gap_pt,
        marker_fontsize,
        group_label_box_pad_pt,
        n_marker_rows=n_rows,
        marker_row_gap_pt=marker_row_gap_pt,
        marker_linewidth=marker_linewidth,
        n_text_lines=n_text_lines,
    )
    anchor_y = _marker_anchor_from_block_center(
        ax, anno_col_x, block_center_y, above_pt, below_pt
    )
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
        row_left, row_center, row_right = _marker_row_x_bounds(
            ax, marker_positions, marker_size, marker_linewidth
        )
    else:
        row_left = row_center = row_right = anno_col_x
    label_x, label_ha = _marker_label_x_and_ha(
        row_left, row_center, row_right, marker_label_align, anno_col_x
    )
    label_x, label_y = _marker_label_position(
        ax,
        label_x,
        bottom_row_y,
        marker_size,
        marker_label_gap_pt,
        marker_linewidth=marker_linewidth,
        n_text_lines=n_text_lines,
        marker_fontsize=marker_fontsize,
        group_label_box_pad_pt=group_label_box_pad_pt,
    )
    return {
        "marker_positions": marker_positions,
        "label_x": label_x,
        "label_y": label_y,
        "label_ha": label_ha,
        "anchor_y": anchor_y,
        "above_pt": above_pt,
        "below_pt": below_pt,
        "row_right": row_right,
    }


def _layout_blocks_from_extents(
    ax,
    ref_x: float,
    target_ys: List[float],
    above_pt_list: List[float],
    below_pt_list: List[float],
    gap_pt: float,
    panel_min: float,
    panel_max: float,
    max_iter: int,
    min_gap_floor: float = 0.0,
) -> np.ndarray:
    """Build symmetric blocks from pt extents and run the bounded 1D stack solver.
"""
    if not target_ys:
        return np.array([], dtype=float)
    gap_data = _gap_pt_to_data(ax, ref_x, float(np.mean(target_ys)), gap_pt)
    blocks: List[_AnnotationLayoutBlock] = []
    for i, target_y in enumerate(target_ys):
        above_pt = above_pt_list[i]
        below_pt = below_pt_list[i]
        height = _pt_to_data_height(ax, ref_x, target_y, above_pt + below_pt)
        blocks.append(
            _AnnotationLayoutBlock(
                block_id=i,
                target_y=float(target_y),
                height=height,
                y=float(target_y),
                min_gap=gap_data,
                panel_min=panel_min,
                panel_max=panel_max,
                original_y=float(target_y),
                above_pt=above_pt,
                below_pt=below_pt,
            )
        )
    solved, _ = _solve_bounded_stack_layout(
        blocks,
        max_iterations=max_iter,
        min_gap_floor=min_gap_floor,
    )
    result = np.empty(len(target_ys), dtype=float)
    for block in solved:
        result[block.block_id] = block.y
    return result


def _spread_phenogram_text_labels(
    ax,
    y_data_list: List[float],
    ref_x: float,
    anno_fontsize: float,
    gap_pt: float,
    group_label_box_pad_pt: float,
    max_iter: int,
    n_lines: int = 1,
    n_lines_list: Optional[List[int]] = None,
    panel_floor_y: Optional[float] = None,
    panel_ceiling_y: Optional[float] = None,
    prefer_down: bool = True,
    min_gap_floor: float = 0.0,
) -> np.ndarray:
    """Lay out text-label units with the bounded 1D stack solver (data coordinates).
"""
    del prefer_down
    if n_lines_list is None:
        n_lines_list = [n_lines] * len(y_data_list)
    above_list: List[float] = []
    below_list: List[float] = []
    for nl in n_lines_list:
        above_pt, below_pt = _text_label_block_extents_pt(
            anno_fontsize, group_label_box_pad_pt, n_lines=nl
        )
        above_list.append(above_pt)
        below_list.append(below_pt)
    panel_min = float(panel_floor_y) if panel_floor_y is not None else -1e9
    panel_max = float(panel_ceiling_y) if panel_ceiling_y is not None else 1e9
    return _layout_blocks_from_extents(
        ax,
        ref_x,
        y_data_list,
        above_list,
        below_list,
        gap_pt,
        panel_min,
        panel_max,
        max_iter,
        min_gap_floor=min_gap_floor,
    )


def _phenogram_inter_unit_gap_pt(
    marker_mode: bool,
    marker_label_gap_pt: float,
    anno_fontsize: float,
    group_min_vertical_gap_pt: float,
    group_marker_to_marker_gap_pt: float,
    repel_force: float = 0.5,
) -> float:
    """Vertical gap between adjacent annotation units (marker cluster + label, or text).

    Defaults to ``marker_label_gap_pt`` so inter-unit spacing matches marker-to-label
    spacing; explicit group gap kwargs remain as optional floors.
"""
    if marker_mode:
        base_gap = float(marker_label_gap_pt)
    else:
        base_gap = max(
            float(marker_label_gap_pt),
            float(anno_fontsize) * 0.25,
        )
    scaled = base_gap * max(float(repel_force), 0.2)
    return max(
        scaled,
        float(group_min_vertical_gap_pt),
        float(group_marker_to_marker_gap_pt),
    )


def _phenogram_panel_ceiling_y(
    visual_row: int,
    row_band_bottom: List[float],
    has_chr_above: bool,
    inter_row_gap: float = _PHENOGRAM_INTER_ROW_GAP,
) -> Optional[float]:
    """Upper data-y bound for annotation anchors (larger y = visual top).

    Prevents labels from entering the stacked chromosome panel above in the same column.
    ``visual_row`` 0 is the top row (chr1–chr11 with default ordering).
"""
    if not has_chr_above:
        return None
    above_row = int(visual_row) - 1
    return (
        float(row_band_bottom[above_row])
        - float(inter_row_gap) * _PANEL_FLOOR_PAD_FRACTION
    )


def _phenogram_panel_floor_y(
    visual_row: int,
    row_band_bottom: List[float],
    has_chr_below: bool,
    inter_row_gap: float = _PHENOGRAM_INTER_ROW_GAP,
) -> Optional[float]:
    """Lower data-y bound for annotation anchors (smaller y = visual bottom).

    Prevents labels from entering the stacked chromosome panel below in the same
    column (larger ``visual_row`` index = lower on the figure).
"""
    if not has_chr_below:
        return None
    return (
        float(row_band_bottom[int(visual_row)])
        + float(inter_row_gap) * _PANEL_FLOOR_PAD_FRACTION
    )


def _bump_phenogram_ylim_from_artists(
    fig,
    axes: List[Any],
    ncols: int,
    global_ymax: float,
    global_ymin: float,
    zorders: Tuple[int, ...] = (230,),
) -> Tuple[float, float]:
    """Expand ymax/ymin if rendered annotation text extends outside the current ylim.

    Catches estimate vs render mismatch (path effects, wrapping, tight_layout).
"""
    renderer = fig.canvas.get_renderer()
    max_data_y = float(global_ymax)
    min_data_y = float(global_ymin)
    for col in range(ncols):
        ax = axes[col]
        inv = ax.transData.inverted()
        for artist in ax.get_children():
            if not isinstance(artist, MplText) or not artist.get_visible():
                continue
            if int(artist.get_zorder()) not in zorders:
                continue
            try:
                bbox = artist.get_window_extent(renderer=renderer)
            except Exception:
                continue
            for x_disp in (bbox.x0, bbox.x1):
                for y_disp in (bbox.y0, bbox.y1):
                    _, y_data = inv.transform((x_disp, y_disp))
                    max_data_y = max(max_data_y, float(y_data))
                    min_data_y = min(min_data_y, float(y_data))
    return (
        max_data_y + _ANNOTATION_Y_MARGIN_DATA,
        min_data_y - _ANNOTATION_Y_MARGIN_DATA,
    )


def _unit_bottom_data_y(ax, x: float, anchor_y: float, below_pt: float) -> float:
    """Return data-y of the visual bottom edge of a spread annotation unit.
"""
    _, dy = _points_to_data_delta(ax, x, anchor_y, dy_points=-float(below_pt))
    return float(anchor_y) + dy




def _anno_column_x(chr_right_x: float, anno_x_pad: float) -> float:
    """Fixed left edge of the annotation column to the right of the chromosome.
"""
    return chr_right_x + anno_x_pad


def _format_phenogram_label(
    label_text: Any,
    max_len: Optional[int] = 20,
) -> str:
    """Truncate long labels when ``max_len`` is set.
"""
    text = str(label_text)
    if max_len is not None and len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


def _phenogram_text_fontproperties(text_kwargs: Dict[str, Any]) -> FontProperties:
    """Build font properties from matplotlib text kwargs.
"""
    fp = FontProperties()
    if "fontsize" in text_kwargs:
        fp.set_size(text_kwargs["fontsize"])
    if "fontweight" in text_kwargs:
        fp.set_weight(text_kwargs["fontweight"])
    if "fontstyle" in text_kwargs:
        fp.set_style(text_kwargs["fontstyle"])
    if "fontfamily" in text_kwargs:
        fp.set_family(text_kwargs["fontfamily"])
    return fp


def _phenogram_line_height_pt(fontsize: float) -> float:
    return float(fontsize) * 1.25


def _measure_text_line_width_pt(
    renderer,
    line: str,
    fontproperties: FontProperties,
) -> float:
    if not line:
        return 0.0
    ismath = "$" in line
    try:
        width, _, _ = renderer.get_text_width_height_descent(
            line, fontproperties, ismath=ismath
        )
        return float(width) * 72.0 / renderer.figure.dpi
    except Exception:
        return len(line) * fontproperties.get_size_in_points() * 0.6


def _wrap_phenogram_label_chars(text: str, max_chars: int) -> str:
    """Wrap label text to at most ``max_chars`` characters per line.
"""
    if max_chars < 1 or not text:
        return text
    return textwrap.fill(
        str(text),
        width=int(max_chars),
        break_long_words=True,
        break_on_hyphens=False,
    )


def _wrap_phenogram_label(
    text: str,
    max_width_pt: float,
    fontproperties: FontProperties,
    renderer,
) -> str:
    """Word-wrap label text to fit ``max_width_pt`` using renderer metrics.
"""
    if max_width_pt <= 0 or not text:
        return text
    words = text.split()
    if not words:
        return text
    lines: List[str] = []
    current = words[0]
    for word in words[1:]:
        candidate = current + " " + word
        if _measure_text_line_width_pt(renderer, candidate, fontproperties) <= max_width_pt:
            current = candidate
        else:
            lines.append(current)
            current = word
    lines.append(current)
    wrapped_lines: List[str] = []
    for line in lines:
        if _measure_text_line_width_pt(renderer, line, fontproperties) <= max_width_pt:
            wrapped_lines.append(line)
            continue
        chunk = ""
        for char in line:
            candidate = chunk + char
            if (
                chunk
                and _measure_text_line_width_pt(renderer, candidate, fontproperties)
                > max_width_pt
            ):
                wrapped_lines.append(chunk)
                chunk = char
            else:
                chunk = candidate
        if chunk:
            wrapped_lines.append(chunk)
    return "\n".join(wrapped_lines)


def _measure_phenogram_label(
    text: str,
    fontproperties: FontProperties,
    renderer,
) -> Tuple[float, float, int]:
    """Return (max_line_width_pt, block_height_pt, n_lines) for wrapped text.
"""
    lines = str(text).split("\n") if text else [""]
    n_lines = max(1, len(lines))
    fontsize = fontproperties.get_size_in_points()
    line_height_pt = _phenogram_line_height_pt(fontsize)
    max_width_pt = max(
        (_measure_text_line_width_pt(renderer, line, fontproperties) for line in lines),
        default=0.0,
    )
    block_height_pt = n_lines * line_height_pt
    return max_width_pt, block_height_pt, n_lines


def _default_anno_wrap_width_pt(
    ax,
    anno_col_x: float,
    global_xmax: float,
    margin_pt: float = 8.0,
) -> float:
    """Available annotation column width in points.
"""
    right_disp = ax.transData.transform((global_xmax, anno_col_x))[0]
    left_disp = ax.transData.transform((anno_col_x, anno_col_x))[0]
    width_px = max(0.0, right_disp - left_disp)
    return width_px * 72.0 / ax.figure.dpi - margin_pt


def _prepare_phenogram_display_label(
    label_text: Any,
    *,
    anno_wrap: bool,
    anno_wrap_width_pt: Optional[float],
    anno_wrap_chars_per_line: Optional[int],
    anno_max_len: Optional[int],
    fontproperties: FontProperties,
    renderer,
) -> Tuple[str, int, float]:
    """Return wrapped/truncated display string, line count, and max line width (pt).
"""
    text = str(label_text)
    if anno_wrap:
        if anno_wrap_chars_per_line is not None and anno_wrap_chars_per_line > 0:
            display = _wrap_phenogram_label_chars(text, anno_wrap_chars_per_line)
        elif (
            anno_wrap_width_pt is not None
            and anno_wrap_width_pt > 0
            and renderer is not None
        ):
            display = _wrap_phenogram_label(
                text, anno_wrap_width_pt, fontproperties, renderer
            )
        else:
            display = text
    else:
        trunc_len = anno_max_len if anno_max_len is not None else 20
        display = _format_phenogram_label(text, max_len=trunc_len)

    if anno_max_len is not None and len(display.replace("\n", "")) > anno_max_len:
        display = _format_phenogram_label(display.replace("\n", " "), max_len=anno_max_len)

    if renderer is not None:
        width_pt, _, n_lines = _measure_phenogram_label(
            display, fontproperties, renderer
        )
    else:
        lines = display.split("\n") if display else [""]
        n_lines = max(1, len(lines))
        fontsize = fontproperties.get_size_in_points()
        width_pt = max(len(line) for line in lines) * fontsize * 0.55
    return display, n_lines, width_pt


def _label_width_data_dx(
    ax,
    x: float,
    y: float,
    width_pt: float,
    ha: str = "left",
) -> Tuple[float, float]:
    """Return (left_x, right_x) in data coords for a label of ``width_pt``.
"""
    dx, _ = _points_to_data_delta(ax, x, y, dx_points=width_pt)
    if ha == "center":
        return x - dx / 2.0, x + dx / 2.0
    if ha == "right":
        return x - dx, x
    return x, x + dx


def _resolve_phenogram_chr_width(chr_width: float, ncols: int) -> float:
    """Widen default ideogram when fewer columns leave more horizontal space.
"""
    if chr_width != _DEFAULT_CHR_WIDTH:
        return chr_width
    return max(_DEFAULT_CHR_WIDTH, _DEFAULT_CHR_WIDTH + (11 - ncols) * 0.02)


def _phenogram_global_xmin(chr_x: float, chr_width: float) -> float:
    return chr_x - chr_width * 0.25


def _phenogram_provisional_xmax(
    chr_x: float,
    chr_width: float,
    anno_x_pad: float,
    min_label_reserve: float = _MIN_LABEL_RESERVE,
) -> float:
    return chr_x + chr_width + anno_x_pad + min_label_reserve


def _compute_phenogram_xlim(
    global_xmin: float,
    chr_x: float,
    chr_width: float,
    estimated_label_xmax: float,
    anno_x_pad: float = 0.18,
    min_chr_fraction: float = _MIN_CHR_X_FRACTION,
) -> float:
    """Ensure the ideogram keeps a minimum share of the horizontal axis.
"""
    global_xmax = max(
        estimated_label_xmax,
        _phenogram_provisional_xmax(chr_x, chr_width, anno_x_pad),
    )
    x_span = global_xmax - global_xmin
    if x_span <= 0:
        return global_xmax
    if chr_width / x_span < min_chr_fraction:
        target_span = chr_width / min_chr_fraction
        global_xmax = global_xmin + target_span
    return global_xmax


def _text_bbox_attach_point(
    ax,
    text_artist,
    renderer,
    side: str = "left",
) -> Tuple[float, float]:
    """Return data coords of the text bbox attach point (left-center by default).
"""
    bbox = text_artist.get_window_extent(renderer=renderer)
    if side == "left":
        x_disp = bbox.x0
    elif side == "right":
        x_disp = bbox.x1
    else:
        x_disp = (bbox.x0 + bbox.x1) / 2.0
    y_disp = (bbox.y0 + bbox.y1) / 2.0
    return ax.transData.inverted().transform((x_disp, y_disp))




def _marker_row_x_bounds(
    ax,
    marker_positions: List[Tuple[float, float, Dict[str, Any]]],
    marker_size: float,
    marker_linewidth: float,
) -> Tuple[float, float, float]:
    """Return (row_left, row_center, row_right) in data coordinates.
"""
    if not marker_positions:
        return 0.0, 0.0, 0.0
    radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    xs = [pos[0] for pos in marker_positions]
    ref_x, ref_y = xs[0], marker_positions[0][1]
    dx_radius, _ = _points_to_data_delta(ax, ref_x, ref_y, dx_points=radius_pt)
    row_left = min(xs) - dx_radius
    row_right = max(xs) + dx_radius
    row_center = (min(xs) + max(xs)) / 2.0
    return row_left, row_center, row_right


def _plot_locus_bar(
    ax,
    left_x: float,
    right_x: float,
    y: float,
    color: str = "red",
    linewidth: float = 2.0,
    zorder: int = 200,
) -> None:
    """Draw a flat locus indicator on the chromosome body (no round line caps).
"""
    ax.plot(
        [left_x, right_x],
        [y, y],
        color=color,
        linewidth=linewidth,
        clip_on=False,
        zorder=zorder,
        solid_capstyle="butt",
        solid_joinstyle="miter",
    )


def _marker_label_x_and_ha(
    row_left: float,
    row_center: float,
    row_right: float,
    marker_label_align: str,
    anno_col_x: Optional[float] = None,
) -> Tuple[float, str]:
    """Group label anchor: fixed left column at ``anno_col_x`` when provided.

    Marker-mode labels sit below the marker row, left-aligned at a fixed
    horizontal offset from the chromosome (``anno_col_x``), not centered on
    the marker row.
"""
    if anno_col_x is not None:
        return float(anno_col_x), "left"
    align = (marker_label_align or "center").lower()
    if align == "left":
        return row_left, "left"
    if align == "right":
        return row_right, "right"
    return row_center, "center"


def _cap_phenogram_arrow_shaft_pt(
    ax,
    chr_right_x: float,
    locus_y: float,
    text_x: float,
    arrow_shaft_pt: float,
    min_slant_gap_pt: float = 6.0,
) -> float:
    """Keep the fixed horizontal arm shorter than the chr-edge to text gap.
"""
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
    text_artist=None,
    renderer=None,
) -> None:
    """Phenogram text annotation connector (vertical MQQ analogue).

    1. A fixed-length horizontal arm at ``locus_y`` from the elbow (right) toward
       the chromosome edge (left), with an arrowhead pointing into the body.
    2. A slanted segment from the elbow to the label left-middle anchor
       (``ha='left'``, ``va='center'`` → ``(text_x, text_y)``).
"""
    del shrink_b_pt  # anchor point is the label left-middle; no shrink offset
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
        "shrinkA": 0,
        "shrinkB": arrow_pad_pt,
    }
    if arrowprops:
        head_props.update(arrowprops)
    head_props.setdefault("shrinkA", 0)
    head_props.setdefault("shrinkB", arrow_pad_pt)

    line_style = {"color": "gray", "linewidth": 0.8}
    if line_kwargs:
        line_style.update(line_kwargs)
    line_style.setdefault("color", head_props.get("color", "gray"))
    if "lw" in head_props and "linewidth" not in line_style:
        line_style["linewidth"] = head_props["lw"]

    ax.plot(
        [elbow_x, chr_right_x],
        [locus_y, locus_y],
        clip_on=False,
        zorder=zorder,
        **line_style,
    )
    head_dx, _ = _points_to_data_delta(ax, chr_right_x, locus_y, dx_points=4.0)
    ax.annotate(
        "",
        xy=(chr_right_x, locus_y),
        xytext=(chr_right_x + head_dx, locus_y),
        arrowprops=head_props,
        clip_on=False,
        zorder=zorder + 1,
    )

    if text_artist is not None and renderer is not None:
        text_attach_x, text_attach_y = _text_bbox_attach_point(
            ax, text_artist, renderer, side="left"
        )
    else:
        text_attach_x, text_attach_y = text_x, text_y

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



def _prepare_lead_display_label_simple(
    label_text: Any,
    *,
    anno_wrap: bool,
    anno_wrap_width_pt: Optional[float],
    anno_wrap_chars_per_line: Optional[int],
    anno_max_len: Optional[int],
    fontproperties: FontProperties,
    renderer=None,
) -> Tuple[str, int]:
    """Return display string and line count (renderer path or char-wrap fallback).
"""
    if renderer is not None:
        display, n_lines, _ = _prepare_phenogram_display_label(
            label_text,
            anno_wrap=anno_wrap,
            anno_wrap_width_pt=anno_wrap_width_pt,
            anno_wrap_chars_per_line=anno_wrap_chars_per_line,
            anno_max_len=anno_max_len,
            fontproperties=fontproperties,
            renderer=renderer,
        )
        return display, n_lines
    if anno_wrap and anno_wrap_chars_per_line and anno_wrap_chars_per_line > 0:
        display = _wrap_phenogram_label_chars(label_text, anno_wrap_chars_per_line)
    else:
        display = _format_phenogram_label(
            label_text,
            max_len=anno_max_len if anno_max_len is not None else 20,
        )
    return display, max(1, display.count("\n") + 1)


def _measure_label_right_x(
    ax,
    x: float,
    y: float,
    label_text: Any,
    *,
    anno_wrap: bool,
    anno_wrap_width_pt: Optional[float],
    anno_wrap_chars_per_line: Optional[int],
    anno_max_len: Optional[int],
    fontproperties: FontProperties,
    renderer=None,
    ha: str = "left",
    fallback_fontsize: Optional[float] = None,
) -> float:
    """Return right x data coordinate for a label.
"""
    if renderer is not None:
        _, _, width_pt = _prepare_phenogram_display_label(
            label_text,
            anno_wrap=anno_wrap,
            anno_wrap_width_pt=anno_wrap_width_pt,
            anno_wrap_chars_per_line=anno_wrap_chars_per_line,
            anno_max_len=anno_max_len,
            fontproperties=fontproperties,
            renderer=renderer,
        )
    else:
        display, _ = _prepare_lead_display_label_simple(
            label_text,
            anno_wrap=anno_wrap,
            anno_wrap_width_pt=anno_wrap_width_pt,
            anno_wrap_chars_per_line=anno_wrap_chars_per_line,
            anno_max_len=anno_max_len,
            fontproperties=fontproperties,
            renderer=None,
        )
        fontsize = fallback_fontsize or fontproperties.get_size_in_points()
        lines = display.split("\n")
        width_pt = fontsize * max(len(line) for line in lines) * 0.55
    _, label_right = _label_width_data_dx(ax, x, y, width_pt, ha=ha)
    return label_right


def _draw_pending_connectors(
    ax,
    chr_right_x: float,
    pending_items: List[Dict[str, Any]],
    max_xlim: float,
) -> float:
    """Flush text artists and draw connectors; return updated max_xlim.
"""
    if not pending_items:
        return max_xlim
    ax.figure.canvas.draw()
    label_renderer = ax.figure.canvas.get_renderer()
    for item in pending_items:
        text_artist = item["text_artist"]
        bbox = text_artist.get_window_extent(renderer=label_renderer)
        right_disp = bbox.x1
        right_x, _ = ax.transData.inverted().transform((right_disp, bbox.y0))
        max_xlim = max(max_xlim, right_x)
        _plot_phenogram_text_connector(
            ax,
            chr_right_x,
            item["original_y"],
            item["text_x"],
            item["text_y"],
            arrow_shaft_pt=item["lead_arrow_shaft_pt"],
            arrow_pad_pt=item["lead_arrow_pad_pt"],
            shrink_b_pt=item["lead_shrink_b_pt"],
            arrowprops=item["lead_arrowprops"],
            line_kwargs=item["lead_line_kwargs"],
            text_artist=text_artist,
            renderer=label_renderer,
        )
    return max_xlim


def _marker_radius_pt(marker_size: float) -> float:
    """Matplotlib scatter ``s`` is area in points²; return marker radius in points.
"""
    return np.sqrt(max(marker_size, 1.0) / np.pi)


def _effective_marker_radius_pt(
    marker_size: float,
    marker_linewidth: float = 0.0,
) -> float:
    """Scatter radius including edge linewidth for layout spacing.
"""
    return _marker_radius_pt(marker_size) + float(marker_linewidth) / 2.0


def _marker_group_block_extents_pt(
    marker_size: float,
    marker_label_gap_pt: float,
    marker_fontsize: float,
    group_label_box_pad_pt: float,
    n_marker_rows: int = 1,
    marker_row_gap_pt: float = 10.0,
    marker_linewidth: float = 0.6,
    n_text_lines: int = 1,
) -> Tuple[float, float]:
    """Vertical half-extents of a marker group block including wrapped marker rows.
"""
    marker_radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    text_height_pt = (
        _phenogram_line_height_pt(marker_fontsize) * max(1, int(n_text_lines))
        + 0.15 * 2.0
    )
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
    """Return (above_pt, below_pt) for a single-row marker group (legacy helper).
"""
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
    """Fixed center-to-center spacing between markers in display points.
"""
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
    """Lay out group markers on a fixed grid: ``marker_max_per_row`` per row, then wrap.

    Returns list of (x, y, item) and the y coordinate of the bottom marker row.
"""
    if marker_max_per_row < 1:
        marker_max_per_row = 1

    radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    step_pt = _marker_center_step_pt(marker_size, marker_gap_pt, marker_linewidth)
    row_pitch_pt = 2.0 * radius_pt + float(marker_row_gap_pt)

    dx_radius, _ = _points_to_data_delta(ax, anno_col_x, anchor_y, dx_points=radius_pt)
    start_center_x = anno_col_x + dx_radius
    dx_step, _ = _points_to_data_delta(ax, start_center_x, anchor_y, dx_points=step_pt)
    _, dy_row_step = _points_to_data_delta(
        ax, anno_col_x, anchor_y, dy_points=-row_pitch_pt
    )

    positions: List[Tuple[float, float, Dict[str, Any]]] = []
    for idx, item in enumerate(group_items):
        row = idx // marker_max_per_row
        col = idx % marker_max_per_row
        y = anchor_y + row * dy_row_step
        x = start_center_x + col * dx_step
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
    n_text_lines: int = 1,
    marker_fontsize: float = 11.0,
    group_label_box_pad_pt: float = 1.5,
) -> Tuple[float, float]:
    """Return (label_x, label_y) for ``va='bottom'`` with glyphs fully below ``marker_y``.

    ``marker_y`` is the bottom marker-row center. Offset matches
    ``_marker_group_block_extents_pt`` below the bottom row (radius, gap, text,
    label box pad).
"""
    marker_radius_pt = _effective_marker_radius_pt(marker_size, marker_linewidth)
    text_height_pt = (
        _phenogram_line_height_pt(marker_fontsize) * max(1, int(n_text_lines))
        + 0.15 * 2.0
    )
    offset_pt = (
        marker_radius_pt
        + float(marker_label_gap_pt)
        + text_height_pt
        + float(group_label_box_pad_pt) * 2.0
    )
    label_y = _display_y_offset_data(
        ax, label_x, marker_y, offset_pt, direction="below"
    )
    return label_x, label_y








def _spread_group_blocks_display(
    ax,
    y_data_list: List[float],
    ref_x: float,
    above_pt: Union[float, Sequence[float]],
    below_pt: Union[float, Sequence[float]],
    gap_pt: float = 4.0,
    max_iter: int = 300,
    panel_floor_y: Optional[float] = None,
    panel_ceiling_y: Optional[float] = None,
    prefer_down: bool = True,
) -> np.ndarray:
    """Lay out rigid annotation-unit centers via the bounded 1D stack solver.

    ``above_pt`` / ``below_pt`` are half-extents in points; block height is their sum.
"""
    del prefer_down
    n = len(y_data_list)
    if n == 0:
        return np.array([], dtype=float)
    if isinstance(above_pt, (int, float)):
        above_list = [float(above_pt)] * n
    else:
        above_list = [float(v) for v in above_pt]
    if isinstance(below_pt, (int, float)):
        below_list = [float(below_pt)] * n
    else:
        below_list = [float(v) for v in below_pt]
    panel_min = float(panel_floor_y) if panel_floor_y is not None else -1e9
    panel_max = float(panel_ceiling_y) if panel_ceiling_y is not None else 1e9
    return _layout_blocks_from_extents(
        ax,
        ref_x,
        list(y_data_list),
        above_list,
        below_list,
        gap_pt,
        panel_min,
        panel_max,
        max_iter,
    )




def _build_marker_text_kwargs(
    anno_kwargs: Optional[Dict[str, Any]],
    marker_fontsize: float = 11,
    marker_label_bbox: bool = True,
) -> Dict[str, Any]:
    text_kwargs: Dict[str, Any] = {
        "fontsize": marker_fontsize,
        "ha": "center",
        "va": "bottom",
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


def _make_legend_prefix_handle(text: str) -> Line2D:
    """Invisible legend entry used as an inline row label (e.g. ``Shape``).
"""
    return Line2D(
        [],
        [],
        linestyle="None",
        marker="None",
        markersize=0.01,
        label=str(text),
    )


def _marker_legend_row_handles(
    prefix: str,
    item_handles: List[Line2D],
) -> List[Line2D]:
    """Return ``[prefix, ...items]`` for a single horizontal legend row.
"""
    if not item_handles:
        return []
    return [_make_legend_prefix_handle(prefix)] + item_handles


def _make_marker_shape_legend_handles(
    marker_shape_map: Dict[str, str],
    markersize: float = _LEGEND_MARKER_SIZE,
) -> List[Line2D]:
    """Legend handles for ``anno_shape`` categories.
"""
    handles: List[Line2D] = []
    for label, marker in marker_shape_map.items():
        handles.append(
            Line2D(
                [0], [0],
                marker=marker,
                linestyle="None",
                markerfacecolor="white",
                markeredgecolor="black",
                markeredgewidth=0.8,
                markersize=markersize,
                label=str(label),
            )
        )
    return handles


def _make_marker_color_legend_handles(
    marker_color_map: Dict[str, str],
    markersize: float = _LEGEND_MARKER_SIZE,
) -> List[Line2D]:
    """Legend handles for ``anno_color`` categories.
"""
    handles: List[Line2D] = []
    for label, color in marker_color_map.items():
        handles.append(
            Line2D(
                [0], [0],
                marker="s",
                linestyle="None",
                markerfacecolor=color,
                markeredgecolor="black",
                markeredgewidth=0.8,
                markersize=markersize,
                label=str(label),
            )
        )
    return handles


def _make_marker_combined_legend_handles(
    marker_shape_map: Dict[str, str],
    marker_color_map: Dict[str, str],
    markersize: float = _LEGEND_MARKER_SIZE,
) -> List[Line2D]:
    """Legend handles showing both shape and fill color per category.
"""
    labels = sorted(
        set(marker_shape_map) | set(marker_color_map),
        key=str,
    )
    handles: List[Line2D] = []
    for label in labels:
        marker = marker_shape_map.get(label, "o")
        color = marker_color_map.get(label, "black")
        handles.append(
            Line2D(
                [0], [0],
                marker=marker,
                linestyle="None",
                markerfacecolor=color,
                markeredgecolor="black",
                markeredgewidth=0.8,
                markersize=markersize,
                label=str(label),
            )
        )
    return handles


def _plot_phenogram_marker_legends(
    fig,
    axes: List[Any],
    ncols: int,
    n_chr: int,
    marker_shape_map: Dict[str, str],
    marker_color_map: Dict[str, str],
    legend_ncol: int,
    legend_kwargs: Optional[Dict[str, Any]] = None,
    legend_fontsize: float = _LEGEND_FONTSIZE,
    legend_markersize: float = _LEGEND_MARKER_SIZE,
    anno_shape: Optional[str] = None,
    anno_color: Optional[str] = None,
) -> None:
    """Draw marker legend rows above the ideogram grid.

    When ``anno_shape`` and ``anno_color`` refer to the same column, draw one
    combined row (``Marker …``). Otherwise draw separate ``Shape …`` and/or
    ``Color …`` rows.
"""
    del legend_ncol  # each row is one line: prefix + all entries
    extra = dict(legend_kwargs) if legend_kwargs else {}
    legend_fontsize = float(extra.pop("fontsize", legend_fontsize))
    legend_markersize = float(extra.pop("markersize", legend_markersize))

    same_column = (
        anno_shape is not None
        and anno_color is not None
        and anno_shape == anno_color
    )

    ax_top = max(
        axes[i % ncols].get_position().y1
        for i in range(n_chr)
    )
    common: Dict[str, Any] = {
        "loc": "lower center",
        "bbox_transform": fig.transFigure,
        "frameon": False,
        "fontsize": legend_fontsize,
        "handlelength": 1.8,
        "handletextpad": 0.55,
        "columnspacing": 1.0,
        "borderaxespad": 0.0,
    }
    common.update(extra)

    if same_column:
        combined_items = _make_marker_combined_legend_handles(
            marker_shape_map,
            marker_color_map,
            markersize=legend_markersize,
        )
        combined_row = _marker_legend_row_handles("Marker", combined_items)
        if not combined_row:
            return
        fig.legend(
            handles=combined_row,
            bbox_to_anchor=(0.5, ax_top + 0.020),
            ncol=len(combined_row),
            **common,
        )
        return

    shape_items = _make_marker_shape_legend_handles(
        marker_shape_map, markersize=legend_markersize
    )
    color_items = _make_marker_color_legend_handles(
        marker_color_map, markersize=legend_markersize
    )
    shape_row = _marker_legend_row_handles("Shape", shape_items)
    color_row = _marker_legend_row_handles("Color", color_items)
    if not shape_row and not color_row:
        return

    # Above the ideogram grid (annotations extend right / downward, not here).
    shape_anchor = (0.5, ax_top + 0.006)
    color_anchor = (0.5, ax_top + 0.034)

    if shape_row and color_row:
        color_legend = fig.legend(
            handles=color_row,
            bbox_to_anchor=color_anchor,
            ncol=len(color_row),
            **common,
        )
        fig.add_artist(color_legend)
        fig.legend(
            handles=shape_row,
            bbox_to_anchor=shape_anchor,
            ncol=len(shape_row),
            **common,
        )
    elif shape_row:
        fig.legend(
            handles=shape_row,
            bbox_to_anchor=shape_anchor,
            ncol=len(shape_row),
            **common,
        )
    else:
        fig.legend(
            handles=color_row,
            bbox_to_anchor=color_anchor,
            ncol=len(color_row),
            **common,
        )


def _build_leads_with_y(
    leads: List[Dict[str, Any]],
    chr_centromere_u: float,
    chr_centromere_l: float,
    max_chr_size: float,
    y_arm1_top: float,
    y_cent_bottom: float,
) -> List[Dict[str, Any]]:
    """Map lead genomic positions to y coordinates on the chromosome drawing.
"""
    del chr_centromere_u, chr_centromere_l, y_cent_bottom
    leads_sorted = sorted(
        [l for l in leads if l.get("pos") is not None and pd.notna(l["pos"])],
        key=lambda x: x["pos"],
    )
    leads_with_y: List[Dict[str, Any]] = []
    for lead in leads_sorted:
        pos = lead["pos"]
        y_pos_chr = y_arm1_top - pos / max_chr_size
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
    renderer=None,
    anno: Optional[Union[bool, str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
    marker_label_align: str = "center",
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
) -> float:
    """Estimate right xlim from marker/text width using frozen axis transforms.
"""
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
        text_kwargs.get(
            "fontsize",
            marker_fontsize if marker_mode else _DEFAULT_ANNO_FONTSIZE,
        )
    )

    for spec in chr_specs:
        ax = spec["ax"]
        ax.set_xlim(global_xmin, provisional_xmax)
        leads_with_y = spec.get("leads_with_y") or []
        if not leads_with_y:
            continue
        anno_col_x = spec["anno_col_x"]
        chr_xmax = _phenogram_provisional_xmax(
            spec["chr_x"], spec["chr_width"], spec["anno_x_pad"]
        )
        local_xmax = chr_xmax

        if marker_mode:
            grouped: Dict[str, List[Dict[str, Any]]] = {}
            for item in leads_with_y:
                group_key = _get_group_key(item["lead"])
                grouped.setdefault(group_key, []).append(item)
            group_display_cache = spec.setdefault("_marker_group_display_cache", {})
            for group_key, group_items in grouped.items():
                anchor_y = float(group_items[0]["original_y"])
                display_label, n_lines = _prepare_marker_group_display_label(
                    group_key,
                    group_items,
                    anno=anno,
                    anno_alias=anno_alias,
                    anno_wrap=anno_wrap,
                    anno_wrap_width_pt=anno_wrap_width_pt,
                    anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                    anno_max_len=anno_max_len,
                    fontproperties=_phenogram_text_fontproperties(text_kwargs),
                    renderer=renderer,
                )
                group_display_cache[group_key] = (display_label, n_lines)
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
                    row_left, row_center, row_right = _marker_row_x_bounds(
                        ax, marker_positions, marker_size, marker_linewidth
                    )
                else:
                    row_left, row_center, row_right = anno_col_x, anno_col_x, anno_col_x
                label_x, label_ha = _marker_label_x_and_ha(
                    row_left, row_center, row_right, marker_label_align, anno_col_x
                )
                _, label_y = _marker_label_position(
                    ax,
                    label_x,
                    bottom_row_y if marker_positions else anchor_y,
                    marker_size,
                    marker_label_gap_pt,
                    marker_linewidth=marker_linewidth,
                    n_text_lines=n_lines,
                    marker_fontsize=anno_fontsize,
                    group_label_box_pad_pt=group_label_box_pad_pt,
                )
                label_right = _measure_label_right_x(
                    ax,
                    label_x,
                    label_y,
                    display_label,
                    anno_wrap=False,
                    anno_wrap_width_pt=None,
                    anno_wrap_chars_per_line=None,
                    anno_max_len=anno_max_len,
                    fontproperties=_phenogram_text_fontproperties(text_kwargs),
                    renderer=renderer,
                    ha=label_ha,
                    fallback_fontsize=anno_fontsize,
                )
                local_xmax = max(local_xmax, label_right)
        else:
            fontproperties = _phenogram_text_fontproperties(text_kwargs)
            wrap_width_pt = anno_wrap_width_pt
            if (
                wrap_width_pt is None
                and renderer is not None
                and not (anno_wrap_chars_per_line and anno_wrap_chars_per_line > 0)
            ):
                wrap_width_pt = _default_anno_wrap_width_pt(
                    ax, anno_col_x, provisional_xmax
                )
            for item in leads_with_y:
                label_text = item["lead"].get("label_text")
                if not label_text:
                    continue
                text_y = item["y_pos_chr"]
                display, n_lines, width_pt = _prepare_phenogram_display_label(
                    label_text,
                    anno_wrap=anno_wrap,
                    anno_wrap_width_pt=wrap_width_pt,
                    anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                    anno_max_len=anno_max_len,
                    fontproperties=fontproperties,
                    renderer=renderer,
                )
                item["display_label"] = display
                item["n_lines"] = n_lines
                _, label_right = _label_width_data_dx(
                    ax, anno_col_x, text_y, width_pt, ha="left"
                )
                local_xmax = max(local_xmax, label_right)

        global_xmax = max(global_xmax, local_xmax * 1.05)

    return global_xmax


def _prepare_text_label_layout(
    spec: Dict[str, Any],
    anno_kwargs: Optional[Dict[str, Any]] = None,
    group_min_vertical_gap_pt: float = 3.5,
    group_marker_to_marker_gap_pt: float = 3,
    group_label_box_pad_pt: float = 1.5,
    marker_label_gap_pt: float = 2.75,
    repel_force: float = 0.5,
    anno_max_iter: int = 300,
    renderer=None,
    global_xmax: Optional[float] = None,
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
) -> float:
    """Compute spread text_y positions and bottom extent without drawing annotations.

    Returns the maximum data-y needed below the ideogram for label overflow.
"""
    if anno_kwargs is None:
        anno_kwargs = {}

    ax = spec["ax"]
    chr_right_x = spec["chr_right_x"]
    leads_with_y = spec.get("leads_with_y") or []
    if not leads_with_y:
        return float(getattr(ax, "_phenogram_ymax", 0) or 0)

    text_kwargs = _build_text_mode_kwargs(anno_kwargs)
    anno_fontsize = float(text_kwargs.get("fontsize", _DEFAULT_ANNO_FONTSIZE))
    fontproperties = _phenogram_text_fontproperties(text_kwargs)
    wrap_width_pt = anno_wrap_width_pt
    if (
        wrap_width_pt is None
        and global_xmax is not None
        and renderer is not None
        and not (anno_wrap_chars_per_line and anno_wrap_chars_per_line > 0)
    ):
        wrap_width_pt = _default_anno_wrap_width_pt(
            ax, spec["anno_col_x"], global_xmax
        )

    labelled_items = sorted(
        [l for l in leads_with_y if l["lead"].get("label_text")],
        key=lambda item: item["y_pos_chr"],
    )
    layout_ymax = float(getattr(ax, "_phenogram_ymax", 0) or 0)
    if not labelled_items:
        spec["label_layout_done"] = True
        return layout_ymax

    for l in labelled_items:
        if l.get("display_label") and l.get("n_lines"):
            display = l["display_label"]
            n_lines = l["n_lines"]
        else:
            display, n_lines = _prepare_lead_display_label_simple(
                l["lead"].get("label_text"),
                anno_wrap=anno_wrap,
                anno_wrap_width_pt=wrap_width_pt,
                anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                anno_max_len=anno_max_len,
                fontproperties=fontproperties,
                renderer=renderer,
            )
        l["display_label"] = display
        l["n_lines"] = n_lines

    orig_y = [l["y_pos_chr"] for l in labelled_items]
    n_lines_list = [l.get("n_lines", 1) for l in labelled_items]
    label_gap_pt = _phenogram_inter_unit_gap_pt(
        marker_mode=False,
        marker_label_gap_pt=marker_label_gap_pt,
        anno_fontsize=anno_fontsize,
        group_min_vertical_gap_pt=group_min_vertical_gap_pt,
        group_marker_to_marker_gap_pt=group_marker_to_marker_gap_pt,
        repel_force=repel_force,
    )
    ideogram_geom = spec.get("ideogram_geom")
    panel_min, panel_max = _annotation_panel_limits(spec, ideogram_geom)
    min_gap_floor = _gap_pt_to_data(
        ax, spec["anno_col_x"], 0.5, group_min_vertical_gap_pt
    )
    above_list = [
        _text_label_block_extents_pt(
            anno_fontsize, group_label_box_pad_pt, n_lines=nl
        )
        for nl in n_lines_list
    ]
    spread_y = _spread_phenogram_text_labels(
        ax,
        orig_y,
        ref_x=spec["anno_col_x"],
        anno_fontsize=anno_fontsize,
        gap_pt=label_gap_pt,
        group_label_box_pad_pt=group_label_box_pad_pt,
        max_iter=anno_max_iter,
        n_lines_list=n_lines_list,
        panel_floor_y=panel_min,
        panel_ceiling_y=panel_max,
        min_gap_floor=min_gap_floor,
    )
    for i, l in enumerate(labelled_items):
        l["text_y"] = float(spread_y[i])

    for i, l in enumerate(labelled_items):
        above_pt, below_pt = above_list[i]
        half_h = _pt_to_data_height(
            ax,
            spec["anno_col_x"],
            l["text_y"],
            (above_pt + below_pt) / 2.0,
        )
        label_bottom_y = float(l["text_y"]) - half_h
        layout_ymax = max(
            layout_ymax, label_bottom_y + _ANNOTATION_Y_MARGIN_DATA
        )

    spec["label_layout_done"] = True
    return layout_ymax


def _compute_marker_group_layout(
    ax,
    leads_with_y: List[Dict[str, Any]],
    chr_right_x: float,
    marker_size: float,
    marker_gap_pt: float,
    marker_label_gap_pt: float,
    marker_max_per_row: int,
    marker_row_gap_pt: float,
    marker_fontsize: float,
    marker_linewidth: float,
    group_label_box_pad_pt: float,
    group_min_vertical_gap_pt: float,
    group_marker_to_marker_gap_pt: float,
    anno_max_iter: int,
    repel_force: float = 0.5,
    panel_floor_y: Optional[float] = None,
    panel_ceiling_y: Optional[float] = None,
    anno: Optional[Union[bool, str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
    fontproperties: Optional[FontProperties] = None,
    renderer=None,
    anno_col_x: Optional[float] = None,
    ideogram_geom: Optional[_ChrIdeogramGeometry] = None,
    spec: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """Lay out marker groups as rigid units; return representatives with ``block_y``.
"""
    ref_x = float(anno_col_x if anno_col_x is not None else chr_right_x)
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for item in leads_with_y:
        group_key = _get_group_key(item["lead"])
        grouped.setdefault(group_key, []).append(item)

    group_representatives: List[Dict[str, Any]] = []
    blocks: List[_AnnotationLayoutBlock] = []
    gap_pt = _phenogram_inter_unit_gap_pt(
        marker_mode=True,
        marker_label_gap_pt=marker_label_gap_pt,
        anno_fontsize=marker_fontsize,
        group_min_vertical_gap_pt=group_min_vertical_gap_pt,
        group_marker_to_marker_gap_pt=group_marker_to_marker_gap_pt,
        repel_force=repel_force,
    )
    if spec is not None:
        panel_min, panel_max = _annotation_panel_limits(spec, ideogram_geom)
    else:
        panel_min = float(panel_floor_y) if panel_floor_y is not None else -1e9
        panel_max = float(panel_ceiling_y) if panel_ceiling_y is not None else 1e9
    gap_data = _gap_pt_to_data(ax, ref_x, 0.5, gap_pt)
    min_gap_floor = _gap_pt_to_data(ax, ref_x, 0.5, group_min_vertical_gap_pt)

    display_cache = spec.get("_marker_group_display_cache", {}) if spec is not None else {}
    for block_id, (group_key, group_items) in enumerate(grouped.items()):
        cached_display = display_cache.get(group_key)
        if cached_display is not None:
            display_label, n_lines = cached_display
        else:
            display_label, n_lines = _prepare_marker_group_display_label(
                group_key,
                group_items,
                anno=anno,
                anno_alias=anno_alias,
                anno_wrap=anno_wrap,
                anno_wrap_width_pt=anno_wrap_width_pt,
                anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                anno_max_len=anno_max_len,
                fontproperties=fontproperties,
                renderer=renderer,
            )
        locus_y = float(group_items[0]["original_y"])
        n_rows = math.ceil(len(group_items) / max(1, marker_max_per_row))
        above_pt, below_pt = _marker_group_block_extents_pt(
            marker_size,
            marker_label_gap_pt,
            marker_fontsize,
            group_label_box_pad_pt,
            n_marker_rows=n_rows,
            marker_row_gap_pt=marker_row_gap_pt,
            marker_linewidth=marker_linewidth,
            n_text_lines=n_lines,
        )
        center_y, height = _marker_block_center_from_anchor(
            ax, ref_x, locus_y, above_pt, below_pt
        )
        block = _AnnotationLayoutBlock(
            block_id=block_id,
            target_y=center_y,
            height=height,
            y=center_y,
            min_gap=gap_data,
            panel_min=panel_min,
            panel_max=panel_max,
            original_y=locus_y,
            n_marker_rows=n_rows,
            n_text_lines=n_lines,
            group_key=group_key,
            group_items=group_items,
            display_label=display_label,
            above_pt=above_pt,
            below_pt=below_pt,
            is_marker=True,
        )
        blocks.append(block)
        group_representatives.append({
            "group_key": group_key,
            "group_items": group_items,
            "original_y": locus_y,
            "display_label": display_label,
            "n_lines": n_lines,
            "block_id": block_id,
        })

    solved, fits = _solve_bounded_stack_layout(
        blocks,
        max_iterations=anno_max_iter,
        min_gap_floor=min_gap_floor,
    )
    if spec is not None and not fits:
        spec["_annotation_layout_warning"] = True

    block_y_by_id = {block.block_id: block.y for block in solved}
    for g in group_representatives:
        block_y = float(block_y_by_id[g["block_id"]])
        g["block_y"] = block_y
        g["marker_y"] = block_y
    return group_representatives


def _prepare_marker_label_layout(
    spec: Dict[str, Any],
    marker_size: float = 81,
    marker_gap_pt: float = 3,
    marker_label_gap_pt: float = 2.75,
    marker_max_per_row: int = _DEFAULT_MARKER_MAX_PER_ROW,
    marker_row_gap_pt: float = 2.5,
    marker_fontsize: float = 11,
    marker_linewidth: float = 0.6,
    group_min_vertical_gap_pt: float = 3.5,
    group_marker_to_marker_gap_pt: float = 3,
    group_label_box_pad_pt: float = 1.5,
    anno_max_iter: int = 300,
    repel_force: float = 0.5,
    marker_label_align: str = "center",
    anno: Optional[Union[bool, str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
    renderer=None,
) -> float:
    """Compute spread marker-group rows and bottom extent without drawing annotations.

    Returns the maximum data-y needed below the ideogram for marker blocks.
"""
    ax = spec["ax"]
    chr_right_x = spec["chr_right_x"]
    leads_with_y = spec.get("leads_with_y") or []
    layout_ymax = float(getattr(ax, "_phenogram_ymax", 0) or 0)
    if not leads_with_y:
        spec["marker_group_layout"] = []
        spec["label_layout_done"] = True
        return layout_ymax

    fontproperties = _phenogram_text_fontproperties(
        _build_marker_text_kwargs(
            None,
            marker_fontsize=marker_fontsize,
            marker_label_bbox=True,
        )
    )

    group_representatives = _compute_marker_group_layout(
        ax,
        leads_with_y,
        chr_right_x,
        marker_size,
        marker_gap_pt,
        marker_label_gap_pt,
        marker_max_per_row,
        marker_row_gap_pt,
        marker_fontsize,
        marker_linewidth,
        group_label_box_pad_pt,
        group_min_vertical_gap_pt,
        group_marker_to_marker_gap_pt,
        anno_max_iter,
        repel_force=repel_force,
        anno=anno,
        anno_alias=anno_alias,
        anno_wrap=anno_wrap,
        anno_wrap_width_pt=anno_wrap_width_pt,
        anno_wrap_chars_per_line=anno_wrap_chars_per_line,
        anno_max_len=anno_max_len,
        fontproperties=fontproperties,
        renderer=renderer,
        anno_col_x=spec["anno_col_x"],
        ideogram_geom=spec.get("ideogram_geom"),
        spec=spec,
    )

    for g in group_representatives:
        group_items = g["group_items"]
        n_lines = int(g.get("n_lines", 1))
        block_y = float(g["block_y"])
        placed = _place_marker_unit_at_center(
            ax,
            spec["anno_col_x"],
            block_y,
            group_items,
            marker_size,
            marker_gap_pt,
            marker_label_gap_pt,
            marker_max_per_row,
            marker_row_gap_pt,
            marker_linewidth,
            marker_fontsize,
            n_lines,
            group_label_box_pad_pt,
            marker_label_align,
        )
        half_h = _pt_to_data_height(
            ax,
            spec["anno_col_x"],
            block_y,
            (placed["above_pt"] + placed["below_pt"]) / 2.0,
        )
        layout_ymax = max(
            layout_ymax,
            float(block_y) + half_h + _ANNOTATION_Y_MARGIN_DATA,
        )

    spec["marker_group_layout"] = group_representatives
    spec["label_layout_done"] = True
    return layout_ymax


def _chr_ideogram_top_y(
    panel_bottom: float,
    telemere_full_length: float = _TELOMERE_LENGTH,
    chr_extent: float = 0.0,
) -> float:
    """Data y of the top edge of the chromosome (including upper telomere cap).
"""
    return panel_bottom + chr_extent + telemere_full_length


def _chr_panel_bottom_from_row_top(
    row_ideogram_top: float,
    chr_extent: float,
    telemere_full_length: float = _TELOMERE_LENGTH,
) -> float:
    """Panel bottom for a top-aligned chromosome in a row band.
"""
    return float(row_ideogram_top) - float(chr_extent) - float(telemere_full_length)


def _display_y_offset_data(
    ax,
    x: float,
    y: float,
    pt_offset: float,
    direction: str = "above",
) -> float:
    """Convert a display-point vertical offset at ``(x, y)`` to a data y coordinate.

    ``direction='above'`` moves toward the visual top of the axes.
"""
    scale = ax.figure.dpi / 72.0
    x_disp, y_disp = ax.transData.transform((x, y))
    sign = 1.0 if direction == "above" else -1.0
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
    """Return data-space height reserved above the ideogram top for its number.
"""
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
    """Return (y, va) with the label entirely above the ideogram top edge.

    Uses ``va='bottom'`` so ``y`` is the text baseline side closest to the body;
    only ``pad_pt`` of clear gap sits between the ideogram cap and the glyphs.
"""
    label_y = _display_y_offset_data(ax, ref_x, body_top_y, pad_pt, direction="above")
    return label_y, "bottom"


def _compute_phenogram_row_bands(
    chr_list: List[str],
    chr_sizes: Dict[str, float],
    max_chr_size: float,
    ncols: int,
    ax,
    ref_x: float,
    chr_number_fontsize: float = 10.0,
    chr_number_pad_pt: float = 4.0,
    inter_row_gap: float = _PHENOGRAM_INTER_ROW_GAP,
    bottom_margin: float = 0.03,
) -> Tuple[int, List[float], List[float], float, float, float, float]:
    """Compute top-down row bands for stacked chromosomes on each column axis.

    ``visual_row`` 0 is the top figure row (chr1–chr11 with default ``ncols``).
    Ideogram tops within a row share ``row_ideogram_top[vr]``; shorter chromosomes
    hang downward from that line.
"""
    n_chr = len(chr_list)
    n_grid_rows = max(1, (n_chr + ncols - 1) // ncols)
    telemere_full_length = _TELOMERE_LENGTH

    row_body_extent = [0.0] * n_grid_rows
    for i, chr_name in enumerate(chr_list):
        visual_row = i // ncols
        chr_extent = chr_sizes[chr_name] / max_chr_size
        row_body_extent[visual_row] = max(
            row_body_extent[visual_row],
            chr_extent + telemere_full_length,
        )

    ref_ideogram_top = (
        row_body_extent[0]
        if row_body_extent
        else telemere_full_length
    )
    label_above_ref = _chr_number_label_above_data(
        ax,
        ref_x,
        ref_ideogram_top,
        fontsize=chr_number_fontsize,
        pad_pt=chr_number_pad_pt,
    )

    total_body = sum(row_body_extent) + inter_row_gap * max(0, n_grid_rows - 1)
    row_ideogram_top = [0.0] * n_grid_rows
    row_band_bottom = [0.0] * n_grid_rows

    row_ideogram_top[0] = float(total_body)
    row_band_bottom[0] = row_ideogram_top[0] - row_body_extent[0]
    for visual_row in range(1, n_grid_rows):
        row_ideogram_top[visual_row] = (
            row_band_bottom[visual_row - 1] - inter_row_gap
        )
        row_band_bottom[visual_row] = (
            row_ideogram_top[visual_row] - row_body_extent[visual_row]
        )

    global_ymin = row_band_bottom[n_grid_rows - 1] - bottom_margin
    global_ymax = row_ideogram_top[0] + label_above_ref

    return (
        n_grid_rows,
        row_ideogram_top,
        row_band_bottom,
        telemere_full_length,
        label_above_ref,
        global_ymin,
        global_ymax,
    )


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
    """Pass 1: draw chromosome geometry from a unified layout model.
"""
    del log, verbose
    geom = _compute_chr_ideogram_geometry(
        chr_x,
        chr_width,
        chr_size,
        chr_centromere_u,
        chr_centromere_l,
        max_chr_size,
        offset,
        telomere_h=_TELOMERE_LENGTH,
    )
    anno_col_x = _anno_column_x(geom.right_x, anno_x_pad)

    arms, _, _ = _chr_ideogram_fill_patches(geom)
    _add_chr_fill_collection(ax, arms, facecolor="white", zorder=100)
    _add_telomere_cap(ax, geom, upper=True, facecolor="white", zorder=101)
    _add_telomere_cap(ax, geom, upper=False, facecolor="white", zorder=101)

    body_clip = PathPatch(
        _chr_ideogram_body_clip_path(geom),
        facecolor="none",
        edgecolor="none",
        linewidth=0.0,
        clip_on=False,
        zorder=50,
    )
    ax.add_patch(body_clip)

    band_groups: Dict[Tuple[Any, int], List[Rectangle]] = {}
    for _, row in chr_cytobands.iterrows():
        y_range = _cytoband_y_range(
            row,
            geom,
            chr_centromere_u,
            chr_centromere_l,
            max_chr_size,
            offset,
        )
        if y_range is None:
            continue
        band_start, band_end = y_range
        band = Rectangle(
            (geom.left_x, band_start),
            geom.width,
            band_end - band_start,
        )
        facecolor = row["COLOR"]
        zorder = 102 if row["STAIN"] == "stalk" else 101
        band_groups.setdefault((facecolor, zorder), []).append(band)

    for (facecolor, zorder), patches in band_groups.items():
        band_pc = PatchCollection(
            patches,
            facecolor=facecolor,
            edgecolor="none",
            linewidths=0.0,
            antialiased=False,
            zorder=zorder,
        )
        band_pc.set_clip_path(body_clip)
        ax.add_collection(band_pc)

    _add_centromere_decoration(ax, geom, facecolor="grey", zorder=103)
    _plot_chr_ideogram_outline(ax, geom)

    leads_with_y = _build_leads_with_y(
        leads,
        chr_centromere_u,
        chr_centromere_l,
        max_chr_size,
        geom.y_arm1_top,
        geom.y_cent_bottom,
    )

    ax.set_xticks(ticks=[])
    ax.set_yticks(ticks=[])

    return {
        "ax": ax,
        "chr_name": chr_name,
        "chr_x": geom.left_x,
        "chr_width": geom.width,
        "anno_x_pad": anno_x_pad,
        "chr_center_x": geom.center_x,
        "chr_right_x": geom.right_x,
        "anno_col_x": anno_col_x,
        "leads_with_y": leads_with_y,
        "body_top_y": geom.ideogram_top,
        "ideogram_geom": geom,
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
    marker_max_per_row: int = _DEFAULT_MARKER_MAX_PER_ROW,
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
    renderer=None,
    anno: Optional[Union[bool, str]] = None,
    anno_alias: Optional[Dict[str, str]] = None,
    anno_wrap: bool = True,
    anno_wrap_width_pt: Optional[float] = None,
    anno_wrap_chars_per_line: Optional[int] = _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    anno_max_len: Optional[int] = None,
    global_xmax: Optional[float] = None,
) -> None:
    """Pass 2: layout and draw annotations after axis limits are frozen.
"""
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

        group_representatives = spec["marker_group_layout"]
        marker_fontproperties = _phenogram_text_fontproperties(text_kwargs)

        pending_marker_connectors: List[Dict[str, Any]] = []

        for g in group_representatives:
            group_key = g["group_key"]
            group_items = g["group_items"]
            block_y = float(g.get("block_y", g["marker_y"]))
            n_lines = int(g.get("n_lines", 1))
            if g.get("display_label"):
                label_text = g["display_label"]
            else:
                label_text, n_lines = _prepare_marker_group_display_label(
                    group_key,
                    group_items,
                    anno=anno,
                    anno_alias=anno_alias,
                    anno_wrap=anno_wrap,
                    anno_wrap_width_pt=anno_wrap_width_pt,
                    anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                    anno_max_len=anno_max_len,
                    fontproperties=marker_fontproperties,
                    renderer=renderer,
                )

            main_item = group_items[0]
            main_y = main_item["original_y"]

            _plot_locus_bar(ax, chr_x, chr_right_x, main_y)

            placed = _place_marker_unit_at_center(
                ax,
                anno_col_x,
                block_y,
                group_items,
                marker_size,
                marker_gap_pt,
                marker_label_gap_pt,
                marker_max_per_row,
                marker_row_gap_pt,
                marker_linewidth,
                marker_fontsize,
                n_lines,
                group_label_box_pad_pt,
                marker_label_align,
            )
            marker_positions = placed["marker_positions"]
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

            label_x = placed["label_x"]
            label_y = placed["label_y"]
            group_text_kwargs = _build_lead_text_kwargs(
                text_kwargs,
                anno_kwargs_single,
                main_item["lead"],
                force_ha="left",
                force_va="bottom",
            )
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
                main_item["lead"].get("snpid"),
            )

            text_artist = ax.text(
                label_x,
                label_y,
                str(label_text),
                clip_on=False,
                zorder=230,
                **group_text_kwargs,
            )
            pending_marker_connectors.append({
                "original_y": main_y,
                "text_x": label_x,
                "text_y": label_y,
                "text_artist": text_artist,
                "lead_arrow_shaft_pt": lead_arrow_shaft_pt,
                "lead_arrow_pad_pt": lead_arrow_pad_pt,
                "lead_shrink_b_pt": lead_shrink_b_pt,
                "lead_arrowprops": lead_arrowprops,
                "lead_line_kwargs": lead_line_kwargs,
            })

            max_marker_x = max(x for x, _, _ in marker_positions)
            max_xlim = max(max_xlim, placed["row_right"])
            label_right = _measure_label_right_x(
                ax,
                anno_col_x,
                label_y,
                label_text,
                anno_wrap=False,
                anno_wrap_width_pt=None,
                anno_wrap_chars_per_line=None,
                anno_max_len=anno_max_len,
                fontproperties=marker_fontproperties,
                renderer=renderer,
                ha="left",
                fallback_fontsize=float(
                    group_text_kwargs.get("fontsize", marker_fontsize)
                ),
            )
            max_xlim = max(max_xlim, label_right, max_marker_x)

        max_xlim = _draw_pending_connectors(
            ax, chr_right_x, pending_marker_connectors, max_xlim
        )
    else:
        text_kwargs = _build_text_mode_kwargs(anno_kwargs)
        anno_fontsize = float(text_kwargs.get("fontsize", _DEFAULT_ANNO_FONTSIZE))
        fontproperties = _phenogram_text_fontproperties(text_kwargs)
        wrap_width_pt = anno_wrap_width_pt
        if (
            wrap_width_pt is None
            and global_xmax is not None
            and renderer is not None
            and not (anno_wrap_chars_per_line and anno_wrap_chars_per_line > 0)
        ):
            wrap_width_pt = _default_anno_wrap_width_pt(
                ax, anno_col_x, global_xmax
            )

        for l in leads_with_y:
            original_y = l["original_y"]
            _plot_locus_bar(ax, chr_x, chr_right_x, original_y)

        pending_connectors: List[Dict[str, Any]] = []
        for l in leads_with_y:
            lead = l["lead"]
            original_y = l["original_y"]
            label_text = lead.get("label_text")
            if not label_text:
                continue
            if not l.get("display_label"):
                label_text, _ = _prepare_lead_display_label_simple(
                    label_text,
                    anno_wrap=anno_wrap,
                    anno_wrap_width_pt=wrap_width_pt,
                    anno_wrap_chars_per_line=anno_wrap_chars_per_line,
                    anno_max_len=anno_max_len,
                    fontproperties=fontproperties,
                    renderer=renderer,
                )
            else:
                label_text = l["display_label"]

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

            text_artist = ax.text(
                text_x,
                text_y,
                str(label_text),
                clip_on=False,
                zorder=202,
                **lead_text_kwargs,
            )
            pending_connectors.append({
                "original_y": original_y,
                "text_x": text_x,
                "text_y": text_y,
                "text_artist": text_artist,
                "lead_arrow_shaft_pt": lead_arrow_shaft_pt,
                "lead_arrow_pad_pt": lead_arrow_pad_pt,
                "lead_shrink_b_pt": lead_shrink_b_pt,
                "lead_arrowprops": lead_arrowprops,
                "lead_line_kwargs": lead_line_kwargs,
            })

        max_xlim = _draw_pending_connectors(
            ax, chr_right_x, pending_connectors, max_xlim
        )

    xlim_default = chr_right_x + spec["anno_x_pad"] + 0.5
    ax._phenogram_xmax = max(
        getattr(ax, "_phenogram_xmax", xlim_default),
        max_xlim * 1.05,
    )
