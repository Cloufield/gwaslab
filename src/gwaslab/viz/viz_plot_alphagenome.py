"""
Plot GWASLab track bundles (from AlphaGenome or other sources) on matplotlib axes.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection

from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_track_bundle import (
    ContactBundle,
    JunctionBundle,
    OverlayBundle,
    TrackBundle,
    x_positions,
)

_DEFAULT_TRACK_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]
_OVERLAY_COLORS = {"REF": "#555555", "ALT": "#E51819"}
_VARIANT_LINE_COLOR = "#333333"
_DEFAULT_AG_YLABEL_FONTSIZE = 9
DEFAULT_AG_YLABEL_TEMPLATE = "{biosample_name}\n{ontology_curie}\n{name}"


def _clean_label_part(val) -> Optional[str]:
    if val is None:
        return None
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    s = str(val).strip()
    if not s or s in (".", "nan", "None"):
        return None
    return s


def _modality_from_row(row: dict) -> Optional[str]:
    """Assay/modality label without redundant ontology or biosample tokens."""
    for key in ("Assay title", "assay"):
        part = _clean_label_part(row.get(key))
        if part:
            return part

    name = _clean_label_part(row.get("name"))
    if not name:
        return None

    ontology_curie = _clean_label_part(row.get("ontology_curie"))
    if ontology_curie:
        prefix = f"{ontology_curie} "
        if name.startswith(prefix):
            name = name[len(prefix):].strip()
        elif name == ontology_curie:
            name = ""

    gtex_tissue = _clean_label_part(row.get("gtex_tissue"))
    if gtex_tissue:
        for prefix in (f"gtex {gtex_tissue}", f"GTEx {gtex_tissue}", gtex_tissue):
            if name.lower().startswith(prefix.lower()):
                name = name[len(prefix):].strip()
                break

    if name.lower().startswith("gtex "):
        name = name[5:].strip()

    biosample_name = _clean_label_part(row.get("biosample_name"))
    if biosample_name and name.lower().startswith(biosample_name.lower()):
        name = name[len(biosample_name):].strip()

    return _clean_label_part(name)


def _default_ag_ylabel(row: dict, modality_suffix: str = "") -> str:
    """Three-line y-label: ontology (biosample), code (CURIE), modality (assay)."""
    lines: List[str] = []
    for key in ("biosample_name", "ontology_curie"):
        part = _clean_label_part(row.get(key))
        if part:
            lines.append(part)
    modality = _modality_from_row(row)
    if modality:
        if modality_suffix:
            modality = f"{modality} {modality_suffix}"
        lines.append(modality)
    elif modality_suffix:
        lines.append(modality_suffix.strip(" ()"))
    if lines:
        return "\n".join(lines)
    return _clean_label_part(row.get("name")) or ""


def _ylabel_from_metadata(
    metadata,
    axis_index: int,
    template: Optional[str],
    modality_suffix: str = "",
) -> str:
    if template is None:
        return ""
    row = metadata.iloc[axis_index].to_dict()
    if template == DEFAULT_AG_YLABEL_TEMPLATE:
        return _default_ag_ylabel(row, modality_suffix=modality_suffix)
    try:
        return template.format(**row)
    except (KeyError, IndexError):
        return str(row.get("name", ""))


def _set_ag_ylabel(ax: plt.Axes, ylab: str, fontsize: float) -> None:
    if ylab:
        ax.set_ylabel(ylab, fontsize=fontsize, rotation=0, ha="right", va="center")


def _style_ag_panel_spines(ax: plt.Axes) -> None:
    """Match regional panel spines (full box on all sides)."""
    ax.spines["top"].set_visible(True)
    ax.spines["top"].set_zorder(1)
    ax.spines["right"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)


def finalize_ag_panel_spines(ag_axes: Sequence[plt.Axes]) -> None:
    """Re-apply regional-style spines after stacked layout adjustments."""
    for ax in ag_axes:
        _style_ag_panel_spines(ax)


def _ag_ylabel_with_allele(
    metadata,
    axis_index: int,
    template: Optional[str],
    allele: str,
    modality_suffix: str = "",
) -> str:
    ylab = _ylabel_from_metadata(metadata, axis_index, template, modality_suffix=modality_suffix)
    return f"{allele}\n{ylab}" if ylab else allele


def _plot_overlay_allele_axis(
    ax: plt.Axes,
    tb: TrackBundle,
    track_index: int,
    *,
    region: Optional[Tuple[int, int, int]],
    color: str,
    allele: str,
    linewidth: float,
    alpha: float,
    ylabel_template: Optional[str],
    fontsize: float,
    variant_pos: Optional[int],
) -> None:
    x = x_positions(tb)
    vals = tb.values[:, track_index]
    if region is not None:
        mask = (x >= region[1]) & (x <= region[2])
        x, vals = x[mask], vals[mask]
    ax.plot(x, vals, color=color, linewidth=linewidth, alpha=alpha)
    ylab = _ag_ylabel_with_allele(tb.metadata, track_index, ylabel_template, allele)
    _set_ag_ylabel(ax, ylab, fontsize)
    if variant_pos is not None:
        ax.axvline(
            variant_pos,
            color=_VARIANT_LINE_COLOR,
            linestyle=":",
            alpha=0.75,
            linewidth=1,
            zorder=1,
        )
    _style_ag_panel_spines(ax)


def plot_ag_tracks(
    bundle: TrackBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    filled: bool = False,
    color: Optional[Union[str, Sequence[str]]] = None,
    ylabel_template: Optional[str] = DEFAULT_AG_YLABEL_TEMPLATE,
    linewidth: float = 1.0,
    verbose: bool = True,
    log: Log = Log(),
    **kwargs: Any,
) -> List[plt.Axes]:
    if len(axes) < bundle.num_axes:
        raise ValueError(
            f"plot_ag_tracks needs {bundle.num_axes} axes, got {len(axes)}"
        )
    x = x_positions(bundle)
    if region is not None:
        mask = (x >= region[1]) & (x <= region[2])
        x = x[mask]
        values = bundle.values[mask, :]
    else:
        values = bundle.values

    for i in range(bundle.num_axes):
        ax = axes[i]
        c = color
        if isinstance(color, (list, tuple)):
            c = color[i % len(color)]
        elif color is None:
            c = _DEFAULT_TRACK_COLORS[i % len(_DEFAULT_TRACK_COLORS)]
        arr = values[:, i]
        if filled:
            ax.fill_between(x, arr, color=c, alpha=kwargs.get("alpha", 0.8), linewidth=0)
        ax.plot(x, arr, color=c, linewidth=linewidth)
        ylab = _ylabel_from_metadata(bundle.metadata, i, ylabel_template)
        _set_ag_ylabel(ax, ylab, kwargs.get("fontsize", _DEFAULT_AG_YLABEL_FONTSIZE))
        _style_ag_panel_spines(ax)
    log.write(f" -Plotted ag_tracks ({bundle.num_axes} track axes)", verbose=verbose)
    return list(axes[: bundle.num_axes])


def plot_ag_overlay(
    bundle: OverlayBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    colors: Optional[Dict[str, str]] = None,
    alpha: float = 0.85,
    ylabel_template: Optional[str] = DEFAULT_AG_YLABEL_TEMPLATE,
    linewidth: float = 1.0,
    verbose: bool = True,
    log: Log = Log(),
    **kwargs: Any,
) -> List[plt.Axes]:
    """Plot REF/ALT as interleaved axes per modality (half height each)."""
    if colors is None:
        colors = _OVERLAY_COLORS.copy()
    n_tracks = bundle.num_axes
    n_axes = 2 * n_tracks
    if len(axes) < n_axes:
        raise ValueError(f"plot_ag_overlay needs {n_axes} axes, got {len(axes)}")

    ref_key = "REF" if "REF" in bundle.tracks else next(iter(bundle.tracks))
    alt_key = "ALT" if "ALT" in bundle.tracks else None
    if alt_key is None:
        keys = list(bundle.tracks.keys())
        alt_key = keys[1] if len(keys) > 1 else ref_key

    ref_tb = bundle.tracks[ref_key]
    alt_tb = bundle.tracks[alt_key]
    fontsize = kwargs.get("fontsize", _DEFAULT_AG_YLABEL_FONTSIZE)
    plot_kw = dict(
        region=region,
        linewidth=linewidth,
        alpha=alpha,
        ylabel_template=ylabel_template,
        fontsize=fontsize,
        variant_pos=bundle.variant_pos,
    )

    for i in range(n_tracks):
        _plot_overlay_allele_axis(
            axes[2 * i],
            ref_tb,
            i,
            color=colors.get("REF", colors.get(ref_key, "#555555")),
            allele="REF",
            **plot_kw,
        )
        _plot_overlay_allele_axis(
            axes[2 * i + 1],
            alt_tb,
            i,
            color=colors.get("ALT", colors.get(alt_key, "#E51819")),
            allele="ALT",
            **plot_kw,
        )

    log.write(
        f" -Plotted ag_overlay ({n_tracks} REF + {n_tracks} ALT track axes)",
        verbose=verbose,
    )
    return list(axes[:n_axes])


def plot_ag_contact(
    bundle: ContactBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    cmap: str = "autumn_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    ylabel_template: Optional[str] = DEFAULT_AG_YLABEL_TEMPLATE,
    verbose: bool = True,
    log: Log = Log(),
    **kwargs: Any,
) -> List[plt.Axes]:
    if len(axes) < bundle.num_axes:
        raise ValueError(
            f"plot_ag_contact needs {bundle.num_axes} axes, got {len(axes)}"
        )
    start = bundle.region[1]
    n_bins = bundle.values.shape[0]
    pos = np.arange(n_bins) * bundle.resolution + start
    if region is not None:
        i0 = max(0, int((region[1] - start) // bundle.resolution))
        i1 = min(n_bins, int((region[2] - start) // bundle.resolution) + 1)
        pos = pos[i0:i1]
        data = bundle.values[i0:i1, i0:i1, :]
    else:
        data = bundle.values

    for i in range(bundle.num_axes):
        ax = axes[i]
        arr = data[:, :, i]
        if vmin is None:
            vmin = float(np.nanmin(arr))
        if vmax is None:
            vmax = float(np.nanmax(arr))
        im = ax.pcolormesh(pos, pos, arr, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_aspect("equal")
        ylab = _ylabel_from_metadata(bundle.metadata, i, ylabel_template)
        if not ylab:
            ylab = str(bundle.metadata.iloc[i].get("name", f"contact_{i}"))
        _set_ag_ylabel(ax, ylab, kwargs.get("fontsize", _DEFAULT_AG_YLABEL_FONTSIZE))
        _style_ag_panel_spines(ax)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    log.write(f" -Plotted ag_contact ({bundle.num_axes} maps)", verbose=verbose)
    return list(axes[: bundle.num_axes])


def plot_ag_sashimi(
    bundle: JunctionBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    color: str = "#1f77b4",
    filter_threshold: Optional[float] = None,
    ylabel_template: Optional[str] = DEFAULT_AG_YLABEL_TEMPLATE,
    verbose: bool = True,
    log: Log = Log(),
    **kwargs: Any,
) -> List[plt.Axes]:
    """Arc-style splice junction plot (simplified sashimi)."""
    if len(axes) < bundle.num_axes:
        raise ValueError(
            f"plot_ag_sashimi needs {bundle.num_axes} axes, got {len(axes)}"
        )
    jdf = bundle.junctions.copy()
    if region is not None:
        jdf = jdf[(jdf["donor"] >= region[1]) & (jdf["acceptor"] <= region[2])]
    if filter_threshold is None and len(bundle.values):
        filter_threshold = 0.05 * float(np.max(bundle.values))
    axis_idx = 0
    for track_i in range(bundle.num_tracks):
        for strand in bundle.strands:
            ax = axes[axis_idx]
            sub = jdf
            if "strand" in jdf.columns:
                sub = jdf[jdf["strand"] == strand]
            segments = []
            widths = []
            for _, row in sub.iterrows():
                donor, acceptor = int(row["donor"]), int(row["acceptor"])
                w = float(row.get("value", row.get("count", 1.0)))
                if filter_threshold is not None and w < filter_threshold:
                    continue
                segments.append([(donor, 0), (donor, 1), (acceptor, 1), (acceptor, 0)])
                widths.append(w)
            if segments:
                lc = LineCollection(
                    segments,
                    linewidths=np.array(widths) * 2 + 0.5,
                    colors=color,
                    alpha=0.7,
                )
                ax.add_collection(lc)
                ax.autoscale()
                ax.set_ylim(-0.5, 1.5)
            ylab = _ylabel_from_metadata(
                bundle.metadata,
                track_i,
                ylabel_template,
                modality_suffix=f"({strand})",
            )
            if not ylab:
                name = bundle.metadata.iloc[track_i].get("name", f"junction_{track_i}")
                ylab = f"{name} ({strand})"
            _set_ag_ylabel(ax, ylab, kwargs.get("fontsize", _DEFAULT_AG_YLABEL_FONTSIZE))
            _style_ag_panel_spines(ax)
            axis_idx += 1
    log.write(f" -Plotted ag_sashimi ({bundle.num_axes} axes)", verbose=verbose)
    return list(axes[: bundle.num_axes])
