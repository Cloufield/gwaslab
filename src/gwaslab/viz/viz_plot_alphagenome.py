"""
Plot GWASLab track bundles (from AlphaGenome or other sources) on matplotlib axes.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
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
_OVERLAY_COLORS = {"REF": "dimgrey", "ALT": "#E51819"}


def _ylabel_from_metadata(metadata, axis_index: int, template: Optional[str]) -> str:
    if template is None:
        return ""
    row = metadata.iloc[axis_index]
    try:
        return template.format(**row.to_dict())
    except (KeyError, IndexError):
        return str(row.get("name", ""))


def plot_ag_tracks(
    bundle: TrackBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    filled: bool = False,
    color: Optional[Union[str, Sequence[str]]] = None,
    ylabel_template: Optional[str] = "{biosample_name}\n({strand}) {name}",
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
        if ylab:
            ax.set_ylabel(ylab, fontsize=kwargs.get("fontsize", 8), rotation=0, ha="right", va="center")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    log.write(f" -Plotted ag_tracks ({bundle.num_axes} track axes)", verbose=verbose)
    return list(axes[: bundle.num_axes])


def plot_ag_overlay(
    bundle: OverlayBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    colors: Optional[Dict[str, str]] = None,
    alpha: float = 0.85,
    ylabel_template: Optional[str] = "{biosample_name}\n({strand}) {name}",
    linewidth: float = 1.0,
    verbose: bool = True,
    log: Log = Log(),
    **kwargs: Any,
) -> List[plt.Axes]:
    if colors is None:
        colors = _OVERLAY_COLORS.copy()
    ref_key = "REF" if "REF" in bundle.tracks else next(iter(bundle.tracks))
    meta = bundle.tracks[ref_key].metadata
    n = bundle.num_axes
    if len(axes) < n:
        raise ValueError(f"plot_ag_overlay needs {n} axes, got {len(axes)}")

    for i in range(n):
        ax = axes[i]
        for label, tb in bundle.tracks.items():
            x = x_positions(tb)
            vals = tb.values[:, i]
            if region is not None:
                mask = (x >= region[1]) & (x <= region[2])
                x, vals = x[mask], vals[mask]
            ax.plot(x, vals, color=colors.get(label, None), alpha=alpha, linewidth=linewidth, label=label)
        if i == 0 and len(bundle.tracks) > 1:
            ax.legend(loc="upper right", fontsize=7, frameon=False)
        ylab = _ylabel_from_metadata(meta, i, ylabel_template)
        if ylab:
            ax.set_ylabel(ylab, fontsize=8, rotation=0, ha="right", va="center")
        if bundle.variant_pos is not None:
            ax.axvline(bundle.variant_pos, color="#FF0000", linestyle="--", alpha=0.6, linewidth=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    log.write(f" -Plotted ag_overlay ({n} track axes)", verbose=verbose)
    return list(axes[:n])


def plot_ag_contact(
    bundle: ContactBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    cmap: str = "autumn_r",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
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
        name = bundle.metadata.iloc[i].get("name", f"contact_{i}")
        ax.set_ylabel(str(name), fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    log.write(f" -Plotted ag_contact ({bundle.num_axes} maps)", verbose=verbose)
    return list(axes[: bundle.num_axes])


def plot_ag_sashimi(
    bundle: JunctionBundle,
    axes: Sequence[plt.Axes],
    region: Optional[Tuple[int, int, int]] = None,
    color: str = "#1f77b4",
    filter_threshold: Optional[float] = None,
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
            name = bundle.metadata.iloc[track_i].get("name", f"junction_{track_i}")
            ax.set_ylabel(f"{name} ({strand})", fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            axis_idx += 1
    log.write(f" -Plotted ag_sashimi ({bundle.num_axes} axes)", verbose=verbose)
    return list(axes[: bundle.num_axes])
