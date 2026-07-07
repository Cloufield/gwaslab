# Variant set plot

!!! info "In development"

GWASLab can compare **effect sizes for a predefined set of variants** across multiple GWAS summary-statistics studies. The workflow uses `gl.SumstatsSet` to extract matching rows from each study, then `plot_effect()` to render a forest-style figure with optional EAF and SNPR2 side panels.

**Key features:**

- Extract variants by **rsID**, **`[CHR, POS]`**, or **`chr:pos:ea:nea`** strings
- Load studies from a **dictionary** of `Sumstats` objects or a **glob pattern**
- Forest-style layout grouped by variant with study-level error bars
- Optional **EAF** and **SNPR2** side panels

!!! note "Not the same as compare_effect or plot_forest"
    - [`compare_effect()`](EffectSize.md) — scatter plot for **two** studies with optional auto lead extraction
    - [`plot_forest()`](ForestPlot.md) — meta-analysis forest plot with a **combined** estimate
    - **Variant set plot** — fixed SNP list pulled from **two or more** sumstats files, no meta-analysis line

## `gl.SumstatsSet()`

```python
sset = gl.SumstatsSet(
    sumstats_dic={"StudyA": gl1, "StudyB": gl2},
    variant_set=["rs12345", [1, 100000], "7:500000:G:A"],
    build="19",
)
```

Or load from files:

```python
sset = gl.SumstatsSet(
    "./data/study_*.tsv.gz",
    variant_set=["rs12345"],
    fmt="auto",
    build="19",
)
```

## `SumstatsSet.plot_effect()`

```python
sset.plot_effect(
    y="STUDY",
    group=["CHR", "POS"],
    y_sort=["CHR", "POS", "STUDY"],
    hue="STUDY",
    save="variant_set.png",
)
```

The same method exists on `gl.Sumstats` when `self.data` already contains multiple studies (e.g. after manual merge), but **`SumstatsSet` is the intended API** for cross-study extraction.

## Prerequisites

- At least **two** loaded `gl.Sumstats` objects (or two files matched by glob)
- **CHR**, **POS**, **SNPID**, **BETA**, and **SE** on each study
- **EAF** / **SNPR2** optional (enable side panels when present)
- **Maximum 100 rows** in the extracted table (plot is skipped above this limit)

## Variant ID formats

| Format | Example | Matching logic |
|--------|---------|----------------|
| rsID string | `"rs12345"` | Exact match on `SNPID` |
| Coordinate pair | `[1, 100000]` | Match `CHR` and `POS` |
| Variant ID string | `"1:100000:A:G"`, `"chr1-100000-A-G"` | Parse separators (`:`, `_`, `-`) and match `CHR` + `POS` |

## Returns

`fig = sset.plot_effect(...)`

| Output | Description |
|--------|-------------|
| `fig` | matplotlib figure (main effect panel + optional EAF/SNPR2 panels) |
| `sset.data` | Long-format DataFrame with a **STUDY** column |

## Options

### SumstatsSet construction

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `sumstats_dic` | `dict` or `str` | Dict of `{study_name: Sumstats}` or glob pattern | Required |
| `variant_set` | `list` | Variants to extract (see formats above) | Required |
| `build` | `str` | Genome build metadata | `"99"` |
| `set` | `str` | Name for this variant set | `"set1"` |
| `**readargs` | — | Passed to `Sumstats()` when loading from files | — |

### plot_effect layout

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `x` | `str` | Effect-size column | `"BETA"` |
| `se` | `str` | Standard-error column | `"SE"` |
| `y` | `str` or `list` | Y-axis tick label column(s) | auto |
| `group` | `list` | Columns defining variant groups | `["CHR", "POS"] + y_sort` |
| `y_sort` | `list` | Sort order within groups | `["CHR", "POS", "STUDY"]` |
| `gap` | `float` | Vertical spacing between variant groups | `0.3` |
| `hue` | `str` | Color points by column (e.g. `"STUDY"`) | `None` |
| `ylabel` | `str` | Y-axis title | `"Variant"` |
| `title` | `str` | Figure title | `None` |
| `effect_label` | `str` | X-axis label override | `None` |

### Side panels

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `eaf` | `str` | EAF column name | `"EAF"` |
| `eaf_panel` | `bool` | Show EAF bar panel | `True` |
| `snpr2` | `str` | SNPR2 column name | `"SNPR2"` |
| `snpvar_panel` | `bool` | Show SNPR2 bar panel | `True` |
| `xlim_eaf` | `tuple` | X limits for EAF panel | `None` |
| `xlim_snpr2` | `tuple` | X limits for SNPR2 panel | `None` |

### Figure styling and output

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `fig_kwargs` | `dict` | matplotlib figure kwargs | `{"figsize": [8, 8], "dpi": 200}` |
| `scatter_kwargs` | `dict` | seaborn scatter kwargs | `{"s": 20}` |
| `err_kwargs` | `dict` | Error-bar kwargs | grey bars |
| `eaf_kwargs` | `dict` | EAF bar kwargs | blue bars |
| `snpr2_kwargs` | `dict` | SNPR2 bar kwargs | blue bars |
| `save` | `bool` or `str` | Save figure path | `None` |
| `save_kwargs` | `dict` | `savefig` kwargs | `{}` |
| `fontsize` | `int` | Font size | `12` |
| `font_family` | `str` | Font family | `"Arial"` |
| `verbose` | `bool` | Print progress messages | `True` |

## Examples

!!! example "Three studies, mixed variant IDs"
    ```python
    import gwaslab as gl

    sset = gl.SumstatsSet(
        {"EUR": gl_eur, "EAS": gl_eas, "AFR": gl_afr},
        variant_set=["rs12345", [7, 156938803]],
        build="19",
    )
    sset.plot_effect(
        y="STUDY",
        group=["CHR", "POS"],
        hue="STUDY",
        effect_label="Per-allele effect size",
        save="variant_set.png",
        verbose=False,
    )
    ```

!!! example "Glob pattern loading"
    ```python
    sset = gl.SumstatsSet(
        "./meta/cohort_*.tsv.gz",
        variant_set=["rs12345"],
        fmt="auto",
        build="19",
    )
    sset.plot_effect(y="STUDY", verbose=False)
    ```

Runnable step-by-step examples with figures: [Variant set plot workflow](visualization_variant_set.md).

Development notebook: `examples/_development/variant_set/visualization_variant_set.ipynb`.
