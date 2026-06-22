"""
Sankey / alluvial plot for categorical sumstats columns.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path

from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_sankey import (
    DEFAULT_BETA_BINS,
    GWASLAB_CATEGORICAL,
    NEUTRAL_NODE_COLOR,
    assign_sankey_colors,
    prepare_sankey_data,
)
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style


def _ribbon_bezier_path(
    x0: float,
    x1: float,
    sy0: float,
    sy1: float,
    ty0: float,
    ty1: float,
    node_width: float,
    curvature: float,
) -> Path:
    """Build a smooth Sankey ribbon with cubic Bezier top/bottom edges."""
    left = x0 + node_width / 2.0
    right = x1 - node_width / 2.0
    if right <= left:
        right = left + 1e-6

    span = right - left
    cx = left + span * float(np.clip(curvature, 0.05, 0.95))

    verts = [
        (left, sy0),
        (left, sy1),
        (cx, sy1),
        (cx, ty1),
        (right, ty1),
        (right, ty0),
        (cx, ty0),
        (cx, sy0),
        (left, sy0),
    ]
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
    ]
    return Path(verts, codes)


def _draw_ribbon(
    ax: plt.Axes,
    x0: float,
    x1: float,
    sy0: float,
    sy1: float,
    ty0: float,
    ty1: float,
    color: str,
    alpha: float,
    node_width: float,
    curvature: float = 0.5,
) -> None:
    path = _ribbon_bezier_path(x0, x1, sy0, sy1, ty0, ty1, node_width, curvature)
    ax.add_patch(
        PathPatch(
            path,
            facecolor=color,
            edgecolor="none",
            alpha=alpha,
            linewidth=0,
            zorder=1,
        )
    )


def _resolve_link_color(
    link: pd.Series,
    color_by: str,
    flow_colors: Dict[str, str],
    node_color_overrides: Dict[str, str],
) -> str:
    color_key = str(link["color_key"])
    if color_key in flow_colors:
        return flow_colors[color_key]
    if link["source"] in node_color_overrides:
        return node_color_overrides[link["source"]]
    return GWASLAB_CATEGORICAL[0]


def _node_label_xy(
    x: float,
    y0: float,
    y1: float,
    stage_idx: int,
    n_stages: int,
    node_width: float,
) -> Tuple[float, float, str]:
    y_mid = (y0 + y1) / 2.0
    pad = node_width * 0.75
    if stage_idx == 0:
        return x - node_width / 2.0 - pad, y_mid, "right"
    if stage_idx == n_stages - 1:
        return x + node_width / 2.0 + pad, y_mid, "left"
    return x + node_width / 2.0 + pad, y_mid, "left"


def plot_sankey(
    data: Union[pd.DataFrame, Any],
    columns: Sequence[str],
    column_map: Optional[Dict[str, str]] = None,
    weight: Union[str, float] = "count",
    dropna: bool = True,
    stage_labels: Optional[Sequence[str]] = None,
    colors: Optional[Dict[str, str]] = None,
    color_by: str = "first",
    palette: str = "auto",
    node_color_mode: str = "stacked",
    beta_bins: Sequence[float] = DEFAULT_BETA_BINS,
    link_alpha: float = 0.55,
    node_width: float = 0.025,
    ribbon_curvature: float = 0.5,
    gap_frac: float = 0.02,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    save: Union[bool, str, None] = None,
    save_kwargs: Optional[Dict[str, Any]] = None,
    fontsize: int = 12,
    font_family: str = "Arial",
    title: Optional[str] = None,
    verbose: bool = True,
    log: Log = Log(),
) -> Tuple[Optional[plt.Figure], Optional[plt.Axes], Dict[str, Any]]:
    """
    Draw a Sankey / alluvial diagram across ordered categorical stages.

    Parameters
    ----------
    data : pandas.DataFrame or Sumstats
        Input table. Sumstats objects use ``.data``.
    columns : sequence of str
        Ordered stage names. Each name may be an existing column or a preset
        (``MAF``, ``P``, ``BETA``, ``STATUS``).
    column_map : dict, optional
        Map preset names to source columns, e.g. ``{"MAF": "EAF"}``.
    weight : str, default ``"count"``
        ``"count"`` to count rows, or a numeric column name to sum.
    dropna : bool, default True
        Drop rows missing any stage category.
    stage_labels : sequence of str, optional
        Labels shown under each stage column.
    colors : dict, optional
        Override colors keyed by first-stage category or full ``node_id``.
    color_by : str, default ``"first"``
        Ribbon coloring: ``"first"`` (track stage-0 category), ``"source"``,
        or ``"target"``.
    palette : str, default ``"auto"``
        ``"auto"`` uses semantic palettes for MAF/P/BETA presets on stage 0.
    node_color_mode : str, default ``"stacked"``
        ``"stacked"`` colored bands inside nodes, or ``"neutral"`` gray nodes.
    beta_bins : sequence of float, optional
        |BETA| thresholds for the BETA preset.
    link_alpha : float, default 0.55
        Ribbon transparency.
    node_width : float, default 0.025
        Node bar width in x-axis data units.
    ribbon_curvature : float, default 0.5
        Horizontal position of Bezier control points as a fraction of link
        span (0.5 = classic symmetric Sankey curves).
    gap_frac : float, default 0.02
        Vertical gap between nodes within a stage.
    fig_kwargs, save, save_kwargs, fontsize, font_family, title, verbose, log
        Standard GWASLab plotting options.

    Returns
    -------
    tuple
        ``(fig, ax, tables)`` where ``tables`` contains ``nodes``, ``links``,
        ``node_bands``, ``flow_colors``, and ``stages``.
    """
    log.write("Start to create Sankey plot...", verbose=verbose)
    log.write(f" -Stages: {list(columns)}", verbose=verbose)
    log.write(f" -Color by: {color_by}", verbose=verbose)

    if node_color_mode not in {"stacked", "neutral"}:
        raise ValueError("node_color_mode must be 'stacked' or 'neutral'")

    nodes_df, links_df, node_bands_df, stage_names, work_df, _ = prepare_sankey_data(
        data=data,
        columns=columns,
        column_map=column_map,
        weight=weight,
        dropna=dropna,
        beta_bins=beta_bins,
        color_by=color_by,
        gap_frac=gap_frac,
    )

    flow_colors, node_color_overrides = assign_sankey_colors(
        stage_names,
        work_df,
        palette=palette,
        colors=colors,
        beta_bins=beta_bins,
    )

    if not node_bands_df.empty:
        node_bands_df = node_bands_df.copy()
        node_bands_df["color"] = node_bands_df["color_key"].map(
            lambda k: flow_colors.get(str(k), GWASLAB_CATEGORICAL[0])
        )
        for node_id, override in node_color_overrides.items():
            node_bands_df.loc[node_bands_df["node_id"] == node_id, "color"] = override

    tables: Dict[str, Any] = {
        "nodes": nodes_df,
        "links": links_df,
        "node_bands": node_bands_df,
        "flow_colors": flow_colors,
        "stages": stage_names,
        "work": work_df,
    }

    if nodes_df is None or nodes_df.empty:
        log.warning("No rows available to plot after filtering.")
        return None, None, tables

    labels = list(stage_labels) if stage_labels is not None else list(stage_names)
    if len(labels) != len(stage_names):
        raise ValueError("stage_labels must have the same length as columns")

    n_stages = len(stage_names)
    max_nodes = int(nodes_df.groupby("stage_idx").size().max())
    auto_width = max(8.0, 3.0 * n_stages)
    auto_height = max(6.0, 1.2 * max_nodes)

    if fig_kwargs is None:
        fig_kwargs = {}
    if "figsize" not in fig_kwargs:
        fig_kwargs = {**fig_kwargs, "figsize": (auto_width, auto_height)}

    style = set_plot_style(
        plot="plot_sankey",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", fig_kwargs)
    save_kwargs = style.get("save_kwargs", save_kwargs)
    fontsize = style.get("fontsize", fontsize)
    font_family = style.get("font_family", font_family)

    fig, ax = plt.subplots(**fig_kwargs)

    for _, link in links_df.iterrows():
        ribbon_color = _resolve_link_color(link, color_by, flow_colors, node_color_overrides)
        _draw_ribbon(
            ax=ax,
            x0=float(nodes_df.loc[nodes_df["node_id"] == link["source"], "x"].iloc[0]),
            x1=float(nodes_df.loc[nodes_df["node_id"] == link["target"], "x"].iloc[0]),
            sy0=float(link["sy0"]),
            sy1=float(link["sy1"]),
            ty0=float(link["ty0"]),
            ty1=float(link["ty1"]),
            color=ribbon_color,
            alpha=link_alpha,
            node_width=node_width,
            curvature=ribbon_curvature,
        )

    for _, node in nodes_df.iterrows():
        x = float(node["x"])
        y0 = float(node["y0"])
        y1 = float(node["y1"])
        height = y1 - y0
        if height <= 0:
            continue
        node_id = node["node_id"]

        if node_color_mode == "neutral":
            bar = FancyBboxPatch(
                (x - node_width / 2.0, y0),
                node_width,
                height,
                boxstyle="round,pad=0,rounding_size=0.004",
                facecolor=NEUTRAL_NODE_COLOR,
                edgecolor="white",
                linewidth=0.6,
                zorder=2,
            )
            ax.add_patch(bar)
        else:
            bands = node_bands_df[node_bands_df["node_id"] == node_id]
            if bands.empty:
                bar_color = node_color_overrides.get(
                    node_id, flow_colors.get(str(node["category"]), GWASLAB_CATEGORICAL[0])
                )
                ax.add_patch(
                    FancyBboxPatch(
                        (x - node_width / 2.0, y0),
                        node_width,
                        height,
                        boxstyle="round,pad=0,rounding_size=0.004",
                        facecolor=bar_color,
                        edgecolor="white",
                        linewidth=0.6,
                        zorder=2,
                    )
                )
            else:
                for _, band in bands.iterrows():
                    band_h = float(band["y1"]) - float(band["y0"])
                    if band_h <= 0:
                        continue
                    ax.add_patch(
                        Rectangle(
                            (x - node_width / 2.0, float(band["y0"])),
                            node_width,
                            band_h,
                            facecolor=band["color"],
                            edgecolor="white",
                            linewidth=0.3,
                            zorder=2,
                        )
                    )

        if height > 0.02:
            lx, ly, ha = _node_label_xy(
                x, y0, y1, int(node["stage_idx"]), n_stages, node_width
            )
            ax.text(
                lx,
                ly,
                str(node["category"]),
                ha=ha,
                va="center",
                fontsize=max(fontsize - 2, 7),
                fontfamily=font_family,
                color="#333333",
                zorder=3,
            )

    x_positions = [i / max(n_stages - 1, 1) for i in range(n_stages)]
    for x_pos, label in zip(x_positions, labels):
        ax.text(
            x_pos,
            -0.05,
            label,
            ha="center",
            va="top",
            fontsize=fontsize,
            fontfamily=font_family,
            transform=ax.get_xaxis_transform(),
        )

    ax.set_xlim(-0.22, 1.22)
    ax.set_ylim(-0.12, 1.02)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontfamily=font_family, pad=12)

    save_figure(fig, save, keyword="sankey", save_kwargs=save_kwargs, log=log, verbose=verbose)
    log.write("Finished creating Sankey plot.", verbose=verbose)
    return fig, ax, tables
