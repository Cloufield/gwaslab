# Lead overlap plot

!!! info "Available since v4.2.0"

## Description

GWASLab can visualize **overlap of annotated lead loci** across two or more studies. For each `Sumstats` object, lead variants are extracted with the same logic as [`get_lead()`](ExtractLead.md) (`_get_sig()`), then merged into shared locus groups when leads fall within `windowsizekb_for_overlap` on the same chromosome.

**Key features:**

- **Venn diagram** for 2–3 studies (`mode="venn"` or `mode="auto"`)
- **UpSet plot** for four or more studies (`mode="upset"` or `mode="auto"`)
- Optional gene labels on UpSet intersections (`show_genes=True`)
- Returns a tabular overlap summary plus the matplotlib figure

!!! note "Top-level API only"
    `plot_lead_overlap()` is exported as `gl.plot_lead_overlap()`. It is **not** a method on `Sumstats`.

## `gl.plot_lead_overlap()`

```python
gl.plot_lead_overlap(
    objects=[gl1, gl2, gl3],
    titles=["Study A", "Study B", "Study C"],
    mode="auto",
)
```

## Prerequisites

- At least **two** loaded `gl.Sumstats` objects
- **CHR**, **POS**, and **P** or **MLOG10P** on each object (same requirements as lead extraction)
- Study labels come from `titles=` or from `obj.meta["gwaslab"]["study_name"]`

## Plot modes

| `mode` | Behavior |
|--------|----------|
| `"auto"` | Venn for 2–3 studies; UpSet for 4+ (default) |
| `"venn"` | Venn diagram (exactly 2 or 3 studies) |
| `"upset"` | UpSet matrix (any count ≥ 2) |

Overlap clustering uses **`windowsizekb_for_overlap`** (default 1000 kb): leads on the same chromosome whose positions are within this window are treated as the same locus group. Lead extraction inside each study uses **`windowsizekb`** (default 500 kb), matching `get_lead()`.

See [Extract Lead Variants](ExtractLead.md) for lead-extraction details.

## Returns

`overlap_df, fig, log = gl.plot_lead_overlap(...)`

| Output | Description |
|--------|-------------|
| `overlap_df` | One row per overlapping locus group |
| `fig` | matplotlib figure (Venn or UpSet) |
| `log` | GWASLab log object |

**Key columns in `overlap_df`:**

| Column | Description |
|--------|-------------|
| `LOCUS_ID` | Locus identifier (`Locus_1`, …) |
| `CHR`, `START`, `END`, `POS` | Genomic span and representative position |
| `STUDIES` | Comma-separated study names with a lead in this group |
| `N_STUDIES` | Number of studies represented |
| `GENE` | Merged gene annotation (when `anno=True`) |
| `LEAD_SNPS` | Dict mapping study label → lead SNP IDs |
| `MEMBERSHIP_KEY` | Bit string per study order, e.g. `1\|0\|1` |
| `IN_<Study>` | Boolean membership flag per study |

For UpSet plots, `overlap_df.attrs["set_list"]` maps compact set IDs (`Set1`, …) to membership patterns.

## Options

### Input and lead extraction

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `objects` | `list` of `Sumstats` | GWAS summary-statistics objects to compare | Required |
| `titles` | `list` of `str` | Display names (one per object) | From `meta["gwaslab"]["study_name"]` or `Sumstats_N` |
| `sig_level` | `float` | Genome-wide significance threshold for lead extraction | `5e-8` |
| `windowsizekb` | `int` | Sliding window (kb) for lead extraction per study | `500` |
| `use_p` | `bool` | Use **P** instead of **MLOG10P** for lead ranking | `False` |
| `get_lead_kwargs` | `dict` | Extra kwargs forwarded to lead extraction | `None` |
| `anno` | `bool` | Annotate leads with nearest gene names | `True` |
| `build` | `str` or `list` | Genome build per object for annotation | object `build` or `"19"` |
| `source` | `str` | Gene annotation backend (`"ensembl"`, `"refseq"`) | `"ensembl"` |
| `gtf_path` | `str` | Custom GTF path for annotation | `None` |
| `wc_correction` | `bool` | Winner's Curse correction on effect sizes | `False` |

### Overlap clustering and display

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `mode` | `str` | `"auto"`, `"venn"`, or `"upset"` | `"auto"` |
| `windowsizekb_for_overlap` | `int` | Window (kb) for merging leads into shared loci | `1000` |
| `show_counts` | `bool` | Show intersection counts on the plot | `True` |
| `show_genes` | `bool` | Show gene names on UpSet intersections | `True` |
| `max_gene_labels` | `int` | Maximum gene labels on UpSet plot | `30` |
| `sort_by` | `str` | UpSet row sort key (`"count"`, …) | `"count"` |
| `title` | `str` | Figure title | `None` |

### Figure styling and output

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `fig_kwargs` | `dict` | matplotlib figure kwargs (`figsize`, `dpi`, …) | `{"figsize": (8, 6), "dpi": 200}` |
| `font_kwargs` | `dict` | Font kwargs for labels | `{"fontsize": 9}` |
| `legend_kwargs` | `dict` | Legend kwargs (UpSet) | `{}` |
| `venn_kwargs` | `dict` | Extra kwargs for Venn drawing | `None` |
| `upset_kwargs` | `dict` | Extra kwargs for UpSet drawing | `None` |
| `save` | `bool` or `str` | Save figure (`True` → default path; `str` → custom path) | `False` |
| `save_kwargs` | `dict` | `savefig` kwargs | `{}` |
| `verbose` | `bool` | Print progress messages | `True` |

## Examples

!!! example "Two-study Venn (auto mode)"
    ```python
    import gwaslab as gl

    gl1 = gl.Sumstats("study_a.txt.gz", fmt="plink2", study_name="EUR")
    gl2 = gl.Sumstats("study_b.txt.gz", fmt="plink2", study_name="EAS")

    overlap_df, fig, log = gl.plot_lead_overlap(
        objects=[gl1, gl2],
        mode="auto",
        anno=True,
    )
    ```

!!! example "Four studies — UpSet (auto mode)"
    ```python
    overlap_df, fig, log = gl.plot_lead_overlap(
        objects=[gl1, gl2, gl3, gl4],
        titles=["GWAS1", "GWAS2", "GWAS3", "GWAS4"],
        mode="auto",
        windowsizekb_for_overlap=1000,
    )
    # UpSet membership mapping:
    set_list = overlap_df.attrs.get("set_list")
    ```

!!! example "Save figure"
    ```python
    gl.plot_lead_overlap(
        objects=[gl1, gl2],
        save="lead_overlap.png",
        verbose=False,
    )
    ```

Runnable step-by-step examples with figures: [Lead overlap workflow](visualization_lead_overlap.md).

API reference: [plot_lead_overlap](api/plotting.md#plot_lead_overlap).
