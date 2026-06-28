# Sankey plot

!!! info "Available since v4.2.0"

## Description

GWASLab can draw **Sankey / alluvial diagrams** to show how variants flow across ordered categorical stages—for example MAF bin → P-value significance → effect-size magnitude, or custom labels such as overall vs subtype GWAS signal.

**Key features:**

- **Preset stages** (`MAF`, `P`, `BETA`, `STATUS`) auto-bin from standard sumstats columns
- **Custom categorical columns** from your DataFrame or `Sumstats.data`
- Stacked node coloring or neutral nodes with colored ribbons
- Returns node/link tables for downstream analysis

## API

```python
# Top-level
fig, ax, tables = gl.plot_sankey(data, columns=["MAF", "P", "BETA"])

# Sumstats method (uses mysumstats.data)
fig, ax, tables = mysumstats.plot_sankey(columns=["MAF", "P"])
```

`data` may be a `pandas.DataFrame` or a `gl.Sumstats` object.

## Stage columns

Each entry in `columns` is one left-to-right stage.

### Preset stages

| Preset | Source column (default) | Categories |
|--------|-------------------------|------------|
| `MAF` | `EAF` (via `column_map`) | Common, low-frequency, rare, ultra-rare bins |
| `P` | `P` | `P<5e-8`, `P<5e-6`, non-significant |
| `BETA` | `BETA` | `\|BETA\|>0.3`, `>0.1`, `>0.05`, Small |
| `STATUS` | `STATUS` | Values from the STATUS column |

Remap source columns with `column_map`, e.g. `column_map={"MAF": "EAF", "P": "P_OVERALL"}`.

### Custom columns

Any column already present in the data is used as-is (string categories). Example: `["overall_signal", "subtype_signal", "MAF"]`.

## Coloring

| Option | Description |
|--------|-------------|
| `color_by="first"` | Ribbon colors follow the first stage category (default) |
| `color_by="source"` | Colors follow the source node category at each hop |
| `color_by="target"` | Colors follow the target node category at each hop |
| `palette="auto"` | GWASLab palettes for MAF/P/BETA presets; categorical palette otherwise |
| `colors` | Dict mapping category label → hex color (overrides auto palette) |
| `node_color_mode="stacked"` | Nodes show stacked color bands (default) |
| `node_color_mode="neutral"` | Gray nodes; color on ribbons only |

## Returns

`fig, ax, tables = gl.plot_sankey(...)`

| Key in `tables` | Description |
|-----------------|-------------|
| `nodes` | Node positions and categories per stage |
| `links` | Ribbon source/target, flow values, ribbon geometry |
| `node_bands` | Stacked color bands within nodes |
| `flow_colors` | Category → color mapping |
| `stages` | Resolved stage names |
| `work` | Working DataFrame with resolved stage categories |

If no rows remain after filtering, `(None, None, tables)` is returned.

## Options

### Data and stages

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `data` | `DataFrame` or `Sumstats` | Input table | Required |
| `columns` | `list` of `str` | Ordered stage names (presets or column names) | Required |
| `column_map` | `dict` | Map preset name → source column | `None` |
| `weight` | `str` or `float` | Column name for flow weight, or `"count"` | `"count"` |
| `dropna` | `bool` | Drop rows with missing stage values | `True` |
| `stage_labels` | `list` of `str` | Axis labels for each stage | Same as `columns` |
| `beta_bins` | `tuple` | `\|BETA\|` thresholds for `BETA` preset | `(0.3, 0.1, 0.05)` |

### Appearance

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `colors` | `dict` | Category → color overrides | `None` |
| `color_by` | `str` | `"first"`, `"source"`, or `"target"` | `"first"` |
| `palette` | `str` | Color palette mode | `"auto"` |
| `node_color_mode` | `str` | `"stacked"` or `"neutral"` | `"stacked"` |
| `link_alpha` | `float` | Ribbon transparency | `0.55` |
| `node_width` | `float` | Node bar width (axes fraction) | `0.025` |
| `ribbon_curvature` | `float` | Bezier curvature for ribbons | `0.5` |
| `gap_frac` | `float` | Vertical gap between node segments | `0.02` |
| `fig_kwargs` | `dict` | Figure kwargs (`figsize`, `dpi`; auto-sized if omitted) | `{"figsize": (10, 6), "dpi": 300}` |
| `fontsize` | `int` | Label font size | `12` |
| `font_family` | `str` | Label font family | `"Arial"` |
| `title` | `str` | Figure title | `None` |

### Output

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `save` | `bool` or `str` | Save figure | `False` |
| `save_kwargs` | `dict` | `savefig` kwargs | `{"dpi": 300, "facecolor": "white"}` |
| `verbose` | `bool` | Print progress messages | `True` |

## Examples

!!! example "MAF → P → BETA preset pipeline"
    ```python
    import gwaslab as gl

    mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")

    fig, ax, tables = mysumstats.plot_sankey(
        columns=["MAF", "P", "BETA"],
        title="Variant flow by MAF, significance, and effect size",
    )
    ```

!!! example "Custom categorical stages"
    ```python
    df = mysumstats.data.copy()
    df["cohort"] = df["STATUS"].map({"significant": "GW sig", "other": "Not sig"})

    fig, ax, tables = gl.plot_sankey(
        df,
        columns=["cohort", "MAF", "P"],
        colors={"GW sig": "#E51819", "Not sig": "#BBBBBB"},
        color_by="first",
    )
    ```

!!! example "Save figure"
    ```python
    gl.plot_sankey(
        mysumstats,
        columns=["MAF", "P"],
        save="sankey_maf_p.png",
        verbose=False,
    )
    ```

Runnable step-by-step examples with figures: [Sankey workflow](visualization_sankey.md).

API reference: [plot_sankey](api/plotting.md#plot_sankey).
