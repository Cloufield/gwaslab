"""Helpers for Sankey / alluvial diagrams from categorical sumstats columns.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

SANKEY_PRESETS = frozenset({"MAF", "P", "BETA", "STATUS"})

MAF_CATEGORY_ORDER = [
    "Common (MAF>=0.05)",
    "Low_frequency (0.01<MAF<=0.05)",
    "Rare (0.001<MAF<=0.01)",
    "Ultra Rare (MAF<=0.001)",
]

P_CATEGORY_ORDER = [
    "P<5e-8",
    "P<5e-6",
    "Non-significant",
]

DEFAULT_BETA_BINS = (0.3, 0.1, 0.05)

MAF_PALETTE = ["#08519C", "#3182BD", "#6BAED6", "#DEEBF7"]
P_PALETTE = ["#CB181D", "#FDAE6B", "#E0E0E0"]
BETA_PALETTE = ["#542788", "#807DBA", "#BCBDDC", "#F0F0F0"]

GWASLAB_CATEGORICAL = [
    "#597FBD",
    "#74BAD3",
    "#E8A838",
    "#D95F5F",
    "#6BA368",
    "#9B6DB6",
    "#4C9F9F",
    "#C77EB5",
    "#8C8C8C",
    "#E07B39",
]

NEUTRAL_NODE_COLOR = "#E6E6E6"


def _extract_dataframe(data: Any) -> pd.DataFrame:
    if hasattr(data, "data") and not isinstance(data, pd.DataFrame):
        return data.data.copy()
    if isinstance(data, pd.DataFrame):
        return data.copy()
    raise TypeError(f"data must be a pandas DataFrame or Sumstats object, got {type(data)}")


def _mapped_col(name: str, column_map: Optional[Dict[str, str]], default: str) -> str:
    if column_map and name in column_map:
        return column_map[name]
    if column_map:
        for key, value in column_map.items():
            if key.upper() == name.upper():
                return value
    return default


def _categorize_maf(eaf: pd.Series) -> pd.Series:
    maf = pd.to_numeric(eaf, errors="coerce")
    maf = maf.where(maf <= 0.5, 1.0 - maf)
    out = pd.Series(np.nan, index=eaf.index, dtype=object)
    out[maf >= 0.05] = MAF_CATEGORY_ORDER[0]
    out[(maf > 0.01) & (maf <= 0.05)] = MAF_CATEGORY_ORDER[1]
    out[(maf > 0.001) & (maf <= 0.01)] = MAF_CATEGORY_ORDER[2]
    out[maf <= 0.001] = MAF_CATEGORY_ORDER[3]
    return out


def _categorize_p(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=p_values.index, dtype=object)
    out[p < 5e-8] = P_CATEGORY_ORDER[0]
    out[(p >= 5e-8) & (p < 5e-6)] = P_CATEGORY_ORDER[1]
    out[p >= 5e-6] = P_CATEGORY_ORDER[2]
    return out


def _categorize_beta(
    beta: pd.Series,
    thresholds: Sequence[float] = DEFAULT_BETA_BINS,
) -> pd.Series:
    abs_beta = pd.to_numeric(beta, errors="coerce").abs()
    out = pd.Series("Small", index=beta.index, dtype=object)
    for threshold in sorted(thresholds, reverse=True):
        out[abs_beta > threshold] = f"|BETA|>{threshold}"
    out[abs_beta.isna()] = np.nan
    return out


def _beta_category_order(thresholds: Sequence[float] = DEFAULT_BETA_BINS) -> List[str]:
    return [f"|BETA|>{t}" for t in sorted(thresholds, reverse=True)] + ["Small"]


def _resolve_stage(
    df: pd.DataFrame,
    name: str,
    column_map: Optional[Dict[str, str]] = None,
    beta_bins: Sequence[float] = DEFAULT_BETA_BINS,
) -> Tuple[pd.Series, str]:
    preset = name.upper()
    if preset in SANKEY_PRESETS:
        pass
    elif name in df.columns:
        series = df[name].astype("string").replace({"nan": np.nan, "None": np.nan})
        return series, name
    else:
        available = ", ".join(sorted(df.columns))
        presets = ", ".join(sorted(SANKEY_PRESETS))
        raise ValueError(
            f"Unknown stage '{name}'. Use an existing column ({available}) "
            f"or a preset ({presets})."
        )

    if preset == "MAF":
        col = _mapped_col("MAF", column_map, "EAF")
        if col not in df.columns:
            raise ValueError(
                f"MAF preset requires column '{col}'. Pass column_map={{'MAF': '<eaf_col>'}} "
                f"or add the column first."
            )
        return _categorize_maf(df[col]), name

    if preset == "P":
        col = _mapped_col("P", column_map, "P")
        if col not in df.columns:
            raise ValueError(
                f"P preset requires column '{col}'. Pass column_map={{'P': '<p_col>'}} "
                f"or add the column first."
            )
        return _categorize_p(df[col]), name

    if preset == "BETA":
        col = _mapped_col("BETA", column_map, "BETA")
        if col not in df.columns:
            raise ValueError(
                f"BETA preset requires column '{col}'. Pass column_map={{'BETA': '<beta_col>'}} "
                f"or add the column first."
            )
        return _categorize_beta(df[col], thresholds=beta_bins), name

    col = _mapped_col("STATUS", column_map, "STATUS")
    if col not in df.columns:
        raise ValueError(
            f"STATUS preset requires column '{col}'. Pass column_map={{'STATUS': '<status_col>'}} "
            f"or add the column first."
        )
    series = df[col].astype("string").replace({"nan": np.nan, "None": np.nan})
    return series, name


def _category_order(stage_name: str, categories: Sequence[str], beta_bins: Sequence[float]) -> List[str]:
    preset = stage_name.upper()
    if preset == "MAF":
        order = [c for c in MAF_CATEGORY_ORDER if c in categories]
    elif preset == "P":
        order = [c for c in P_CATEGORY_ORDER if c in categories]
    elif preset == "BETA":
        order = [c for c in _beta_category_order(beta_bins) if c in categories]
    else:
        seen = []
        for cat in categories:
            if cat not in seen:
                seen.append(cat)
        order = seen
    remaining = [c for c in categories if c not in order]
    return order + remaining


def _palette_for_ordered_categories(
    ordered_categories: Sequence[str],
    preset_palette: Sequence[str],
    preset_order: Sequence[str],
) -> List[str]:
    order_index = {cat: i for i, cat in enumerate(preset_order)}
    palette = []
    for cat in ordered_categories:
        idx = order_index.get(cat, len(preset_palette) - 1)
        palette.append(preset_palette[min(idx, len(preset_palette) - 1)])
    return palette


def _categorical_palette(categories: Sequence[str]) -> List[str]:
    return [GWASLAB_CATEGORICAL[i % len(GWASLAB_CATEGORICAL)] for i in range(len(categories))]


def assign_sankey_colors(
    stage_names: Sequence[str],
    work_df: pd.DataFrame,
    palette: str = "auto",
    colors: Optional[Dict[str, str]] = None,
    beta_bins: Sequence[float] = DEFAULT_BETA_BINS,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Assign flow colors keyed by first-stage category and optional node_id overrides.
"""
    if work_df.empty or "_stage_0" not in work_df.columns:
        return {}, {}

    flow_keys = _category_order(
        stage_names[0],
        work_df["_stage_0"].dropna().unique().tolist(),
        beta_bins,
    )

    if palette == "auto":
        preset = stage_names[0].upper()
        if preset == "MAF":
            hex_list = _palette_for_ordered_categories(flow_keys, MAF_PALETTE, MAF_CATEGORY_ORDER)
        elif preset == "P":
            hex_list = _palette_for_ordered_categories(flow_keys, P_PALETTE, P_CATEGORY_ORDER)
        elif preset == "BETA":
            hex_list = _palette_for_ordered_categories(
                flow_keys, BETA_PALETTE, _beta_category_order(beta_bins)
            )
        else:
            hex_list = _categorical_palette(flow_keys)
    else:
        hex_list = _categorical_palette(flow_keys)

    flow_colors = {key: hex_list[i] for i, key in enumerate(flow_keys)}
    node_color_overrides: Dict[str, str] = {}

    if colors:
        for key, value in colors.items():
            if "|" in key:
                node_color_overrides[key] = value
            else:
                flow_colors[key] = value

    return flow_colors, node_color_overrides


def _link_color_key(
    color_by: str,
    stage_idx: int,
    row: pd.Series,
    stage_names: Sequence[str],
) -> str:
    if color_by == "first":
        return str(row["_stage_0"])
    if color_by == "source":
        return str(row[f"_stage_{stage_idx}"])
    if color_by == "target":
        return str(row[f"_stage_{stage_idx + 1}"])
    raise ValueError("color_by must be one of: 'first', 'source', 'target'")


def _link_groupby_cols(color_by: str, stage_idx: int, n_stages: int) -> List[str]:
    if color_by == "first":
        cols = ["_stage_0", f"_stage_{stage_idx}", f"_stage_{stage_idx + 1}"]
        deduped: List[str] = []
        for col in cols:
            if col not in deduped:
                deduped.append(col)
        return deduped
    return [f"_stage_{stage_idx}", f"_stage_{stage_idx + 1}"]


def prepare_sankey_data(
    data: Any,
    columns: Sequence[str],
    column_map: Optional[Dict[str, str]] = None,
    weight: Union[str, float] = "count",
    dropna: bool = True,
    beta_bins: Sequence[float] = DEFAULT_BETA_BINS,
    color_by: str = "first",
    gap_frac: float = 0.02,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], pd.DataFrame, Dict[str, str]]:
    if len(columns) < 2:
        raise ValueError("columns must contain at least 2 stages")
    if color_by not in {"first", "source", "target"}:
        raise ValueError("color_by must be one of: 'first', 'source', 'target'")

    df = _extract_dataframe(data)
    stage_series: List[pd.Series] = []
    stage_names: List[str] = []

    for name in columns:
        series, stage_label = _resolve_stage(df, name, column_map=column_map, beta_bins=beta_bins)
        stage_series.append(series)
        stage_names.append(stage_label)

    work = pd.DataFrame({f"_stage_{i}": s for i, s in enumerate(stage_series)})
    if weight == "count":
        work["_weight"] = 1.0
    else:
        if weight not in df.columns:
            raise ValueError(f"weight column '{weight}' not found in data")
        work["_weight"] = pd.to_numeric(df[weight], errors="coerce").fillna(0.0)

    stage_cols = [f"_stage_{i}" for i in range(len(columns))]
    if dropna:
        work = work.dropna(subset=stage_cols)
    work = work[work["_weight"] > 0]

    empty_nodes = pd.DataFrame(
        columns=["node_id", "stage_idx", "stage_name", "category", "value", "x", "y0", "y1"]
    )
    empty_links = pd.DataFrame(
        columns=[
            "source",
            "target",
            "value",
            "color_key",
            "flow_id",
            "stage_from",
            "stage_to",
            "sy0",
            "sy1",
            "ty0",
            "ty1",
        ]
    )
    empty_bands = pd.DataFrame(
        columns=["node_id", "color_key", "y0", "y1", "color"]
    )

    if work.empty:
        return empty_nodes, empty_links, empty_bands, stage_names, work, {}

    flow_colors, _ = assign_sankey_colors(
        stage_names, work, palette="auto", colors=None, beta_bins=beta_bins
    )

    nodes_records: List[Dict[str, Any]] = []
    for stage_idx, stage_name in enumerate(stage_names):
        grouped = work.groupby(f"_stage_{stage_idx}", sort=False)["_weight"].sum()
        categories = _category_order(stage_name, list(grouped.index), beta_bins)
        for category in categories:
            if category not in grouped.index:
                continue
            nodes_records.append(
                {
                    "node_id": f"{stage_name}|{category}",
                    "stage_idx": stage_idx,
                    "stage_name": stage_name,
                    "category": category,
                    "value": float(grouped[category]),
                }
            )

    nodes_df = pd.DataFrame(nodes_records)
    link_records: List[Dict[str, Any]] = []
    n_stages = len(stage_names)

    for stage_idx in range(n_stages - 1):
        group_cols = _link_groupby_cols(color_by, stage_idx, n_stages)
        grouped = work.groupby(group_cols, sort=False)["_weight"].sum().reset_index()
        for _, row in grouped.iterrows():
            source_cat = row[f"_stage_{stage_idx}"]
            target_cat = row[f"_stage_{stage_idx + 1}"]
            color_key = _link_color_key(color_by, stage_idx, row, stage_names)
            source_id = f"{stage_names[stage_idx]}|{source_cat}"
            target_id = f"{stage_names[stage_idx + 1]}|{target_cat}"
            link_records.append(
                {
                    "source": source_id,
                    "target": target_id,
                    "value": float(row["_weight"]),
                    "color_key": color_key,
                    "flow_id": f"{color_key}|{source_id}|{target_id}",
                    "stage_from": stage_idx,
                    "stage_to": stage_idx + 1,
                }
            )

    links_df = pd.DataFrame(link_records)
    nodes_df, links_df, node_bands_df = layout_sankey(
        nodes_df,
        links_df,
        n_stages,
        gap_frac=gap_frac,
        flow_colors=flow_colors,
        color_by=color_by,
    )
    return nodes_df, links_df, node_bands_df, stage_names, work, flow_colors


def _build_node_bands(
    nodes_df: pd.DataFrame,
    links_df: pd.DataFrame,
    flow_colors: Dict[str, str],
    color_by: str,
) -> pd.DataFrame:
    records: List[Dict[str, Any]] = []
    for _, node in nodes_df.iterrows():
        node_id = node["node_id"]
        stage_idx = int(node["stage_idx"])
        if stage_idx == 0 and color_by == "first":
            touching = links_df[links_df["source"] == node_id]
            y0_col, y1_col = "sy0", "sy1"
        else:
            touching = links_df[links_df["target"] == node_id]
            if touching.empty and stage_idx == 0:
                touching = links_df[links_df["source"] == node_id]
                y0_col, y1_col = "sy0", "sy1"
            else:
                y0_col, y1_col = "ty0", "ty1"

        if touching.empty:
            continue

        for color_key, grp in touching.groupby("color_key", sort=False):
            records.append(
                {
                    "node_id": node_id,
                    "color_key": color_key,
                    "y0": float(grp[y0_col].min()),
                    "y1": float(grp[y1_col].max()),
                    "color": flow_colors.get(str(color_key), GWASLAB_CATEGORICAL[0]),
                }
            )

    return pd.DataFrame(records)


def layout_sankey(
    nodes_df: pd.DataFrame,
    links_df: pd.DataFrame,
    n_stages: int,
    gap_frac: float = 0.02,
    flow_colors: Optional[Dict[str, str]] = None,
    color_by: str = "first",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if nodes_df.empty:
        return nodes_df, links_df, pd.DataFrame(columns=["node_id", "color_key", "y0", "y1", "color"])

    if flow_colors is None:
        flow_colors = {}

    nodes_df = nodes_df.copy()
    links_df = links_df.copy()
    nodes_df["x"] = nodes_df["stage_idx"] / max(n_stages - 1, 1)
    nodes_df["y0"] = np.nan
    nodes_df["y1"] = np.nan

    node_pos: Dict[str, Tuple[float, float]] = {}
    for stage_idx in range(n_stages):
        stage_nodes = nodes_df[nodes_df["stage_idx"] == stage_idx].copy()
        stage_nodes = stage_nodes.sort_values("category")
        total = stage_nodes["value"].sum()
        n_nodes = len(stage_nodes)
        usable = 1.0 - gap_frac * max(n_nodes - 1, 0)
        y = 0.0
        for node_id, row in stage_nodes.set_index("node_id").iterrows():
            height = (row["value"] / total) * usable if total > 0 else 0.0
            nodes_df.loc[nodes_df["node_id"] == node_id, ["y0", "y1"]] = [y, y + height]
            node_pos[node_id] = (y, y + height)
            y += height + gap_frac

    for col in ("sy0", "sy1", "ty0", "ty1"):
        links_df[col] = np.nan

    for stage_idx in range(n_stages - 1):
        stage_links = links_df[links_df["stage_from"] == stage_idx].copy()
        if stage_links.empty:
            continue

        for source_id, source_links in stage_links.groupby("source", sort=False):
            source_links = source_links.sort_values(["color_key", "target"])
            sy0, sy1 = node_pos[source_id]
            source_total = source_links["value"].sum()
            source_height = sy1 - sy0
            y = sy0
            for idx, link in source_links.iterrows():
                height = (link["value"] / source_total) * source_height if source_total > 0 else 0.0
                links_df.loc[idx, "sy0"] = y
                links_df.loc[idx, "sy1"] = y + height
                y += height

        target_links = links_df[links_df["stage_from"] == stage_idx].copy()
        for target_id, incoming in target_links.groupby("target", sort=False):
            incoming = incoming.sort_values(["color_key", "source"])
            ty0, ty1 = node_pos[target_id]
            target_total = incoming["value"].sum()
            target_height = ty1 - ty0
            y = ty0
            for idx, link in incoming.iterrows():
                height = (link["value"] / target_total) * target_height if target_total > 0 else 0.0
                links_df.loc[idx, "ty0"] = y
                links_df.loc[idx, "ty1"] = y + height
                y += height

    node_bands_df = _build_node_bands(nodes_df, links_df, flow_colors, color_by)
    return nodes_df, links_df, node_bands_df


def maf_bin_counts(eaf: pd.Series) -> Dict[str, int]:
    categorized = _categorize_maf(eaf).dropna()
    return {label: int((categorized == label).sum()) for label in MAF_CATEGORY_ORDER}
