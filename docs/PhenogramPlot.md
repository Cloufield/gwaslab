# Phenogram plot

!!! info "Available since v4.2.0"

## Description

GWASLab can create a **phenogram** (chromosome ideogram): a karyotype-style layout with UCSC cytoband coloring, one chromosome per panel, and genome-wide significant lead variants shown as red bars at their genomic positions. Optional labels appear to the right of each chromosome.

**Key features:**

- **Layout**: Autosomes 1–22 in a multi-column grid (default 11 columns); optional chrX/chrY
- **Cytobands**: Bundled hg19/hg38 cytoband files (auto-selected from `build`)
- **Leads**: By default, lead variants are extracted with the same logic as `get_lead()` before plotting
- **Annotations**: Text labels or **group mode** with categorized markers and a figure legend

## `.plot_phenogram()`

```python
mysumstats.plot_phenogram()
```

## Prerequisites

- **CHR**, **POS**, and **P** or **MLOG10P** on the Sumstats object
- `build` (`"19"` or `"38"`) should match your variant coordinates; cytobands are loaded from the packaged files unless `cytoband_path` is set
- With default `use_lead_extraction=True`, only lead variants at `anno_sig_level` (default `5e-8`) within `windowsizekb` (default 500 kb) are annotated

## Annotation modes

| Mode | Trigger | Behavior |
|------|---------|----------|
| **Text mode** | `anno` is set and `anno_group` is `None` | Text labels with connector lines; layout uses repel/wrap (`anno_style`, `repel_force`, `anno_wrap*`) |
| **Group mode** | `anno_group` is set to a column name | Markers grouped by category; shape and/or color mapped per group; figure legend when `show_legend=True` |

In group mode, set `anno_group` to a categorical column (e.g. `"GENE"` after `anno_gene()`, or `"LOCUS"`). Use `anno_shape` and/or `anno_color` for additional columns:

- **Same column** for both (e.g. `"TRAIT"`) — one combined legend row; shape and color encode the same category.
- **Different columns** (e.g. `anno_shape="SUBTYPE"`, `anno_color="ANCESTRY"`) — separate Shape and Color legend rows for dual encoding.

See [Multi-ancestry and multi-subtype (simulated)](visualization_phenogram.md#multi-ancestry-and-multi-subtype-simulated) for a runnable example.

## Options

### Data and lead extraction

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `snpid` | `str` | Variant ID column | `"SNPID"` |
| `chrom` | `str` | Chromosome column | `"CHR"` |
| `pos` | `str` | Position column | `"POS"` |
| `p` | `str` | P-value column | `"P"` |
| `mlog10p` | `str` | -log10(P) column | `"MLOG10P"` |
| `windowsizekb` | `int` | Sliding window (kb) for lead extraction | `500` |
| `sig_level` | `float` | Significance level (reserved; lead extraction uses `anno_sig_level`) | `5e-8` |
| `use_lead_extraction` | `bool` | If `True`, extract leads via `_get_sig()`; if `False`, annotate every row in the input table | `True` |
| `anno_sig_level` | `float` | P threshold for lead extraction / annotation eligibility | `5e-8` |

### Layout

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `build` | `str` | Genome build (`"19"`, `"38"`) for cytoband selection | object `build` or `"19"` |
| `cytoband_path` | `str` | Path to UCSC cytoband file (override bundled data) | `None` |
| `include_sex_chr` | `bool` | Include chrX and chrY after autosomes | `False` |
| `only_anno_chr` | `bool` | Plot only chromosomes that contain at least one lead/annotation | `False` |
| `ncols` | `int` | Number of chromosome columns per row | `11` |
| `figsize` | `tuple` | Figure size `(width, height)` in inches | `(20, 48)` |
| `dpi` | `int` | Figure DPI | `100` |
| `chr_width` | `float` | Chromosome body width in data coordinates | `0.35` |
| `chr_x` | `float` | Chromosome left x boundary in data coordinates | `0.0` |
| `anno_x_pad` | `float` | Extra horizontal gap between chromosome and annotation area | `0.18` |
| `chr_label_pad` | `float` | Vertical padding above chromosome number labels | `0.06` |

### Text annotations

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `anno` | `bool`, `str`, or `None` | Label source: `None` (bars only), `True` (Chr:POS), `"GENENAME"`, or column name | `None` |
| `anno_set` | `list` | Restrict annotation to these variant IDs | `[]` |
| `anno_alias` | `dict` | Map SNPID to custom label text | `{}` |
| `anno_source` | `str` | Gene annotation backend for `"GENENAME"` (`"ensembl"`, `"refseq"`) | `"ensembl"` |
| `anno_gtf_path` | `str` | Custom GTF path for gene lookup | `None` |
| `anno_max_rows` | `int` | Maximum number of variants to label (sorted by significance) | `200` |
| `anno_style` | `str` | Label placement style (e.g. `"expand"`) | `"expand"` |
| `repel_force` | `float` | Text repulsion strength | `0.5` |
| `anno_max_iter` | `int` | Maximum layout iterations for label repel | `300` |
| `anno_wrap` | `bool` | Wrap long labels onto multiple lines | `True` |
| `anno_wrap_width_pt` | `float` | Max label width in points for wrapping | `None` |
| `anno_wrap_chars_per_line` | `int` | Character-based wrap width when point width is unset | auto |
| `anno_max_len` | `int` | Truncate labels longer than this many characters | `None` |
| `anno_kwargs` | `dict` | Default matplotlib kwargs for text annotations | `None` |
| `anno_kwargs_single` | `dict` | Per-variant override kwargs (keyed by SNPID) | `{}` |
| `arrow_kwargs` | `dict` | Kwargs for connector arrows/lines | `{}` |

### Group mode

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `anno_group` | `str` | Column that defines marker groups (enables group mode) | `None` |
| `anno_shape` | `str` | Column for marker shape mapping | `None` |
| `anno_color` | `str` | Column for marker color mapping | `None` |
| `marker_shapes` | `list` | Pool of matplotlib marker codes for auto shape mapping | built-in pool |
| `marker_colors` | `list` | Pool of colors for auto color mapping | built-in pool |
| `marker_shape_map` | `dict` | Explicit map from group value to marker code | `None` |
| `marker_color_map` | `dict` | Explicit map from group value to color | `None` |
| `marker_size` | `float` | Marker scatter size | `81` |
| `marker_gap_pt` | `float` | Horizontal gap between markers within a group (points) | `3` |
| `marker_label_gap_pt` | `float` | Vertical gap between marker row and text label (points) | `2.75` |
| `marker_max_per_row` | `int` | Maximum markers per row within a group | `4` |
| `marker_row_gap_pt` | `float` | Vertical gap between marker rows (points) | `2.5` |
| `marker_label_align` | `str` | Label alignment relative to marker row (`center`, `left`, `right`) | `"center"` |
| `marker_fontsize` | `float` | Font size for group-mode text labels | `11` |
| `marker_linewidth` | `float` | Marker edge linewidth | `0.6` |
| `marker_label_bbox` | `bool` | Draw white text outline on marker labels | `True` |
| `group_min_vertical_gap_pt` | `float` | Minimum vertical spacing between group blocks (points) | `3.5` |
| `group_marker_to_marker_gap_pt` | `float` | Minimum spacing between marker rows of adjacent groups (points) | `3` |
| `group_label_box_pad_pt` | `float` | Extra padding around labels for overlap checks (points) | `1.5` |
| `show_legend` | `bool` | Show figure legend below the plot (group mode) | `True` |
| `legend_ncol` | `int` | Number of legend columns | `6` |
| `legend_kwargs` | `dict` | Extra kwargs passed to the legend | `{}` |

### Output

| Option | DataType | Description | Default |
|--------|----------|-------------|---------|
| `save` | `bool` or `str` | If `True`, save to default path; if string, save to that path | `False` |
| `save_kwargs` | `dict` | Arguments for `matplotlib.savefig()` | `{}` |
| `fig_kwargs` | `dict` | Extra figure creation kwargs (can set `figsize`, `dpi`) | `{}` |
| `verbose` | `bool` | Print progress messages | `True` |

## Examples

!!! example "Basic phenogram with Chr:POS labels"
    ```python
    import gwaslab as gl

    # Load BBJ T2D sumstats (hg19)
    mysumstats = gl.Sumstats(
        "examples/0_sample_data/t2d_bbj.txt.gz",
        snpid="SNP", chrom="CHR", pos="POS", p="P", build="19", verbose=False,
    )
    mysumstats.fix_chr(verbose=False)

    # Default: extract leads, label with Chr:POS
    mysumstats.plot_phenogram(anno=True, build="19")
    ```

!!! example "Compact layout (chromosomes with hits only)"
    ```python
    # Plot only chromosomes that contain at least one lead
    mysumstats.plot_phenogram(anno=True, only_anno_chr=True)
    ```

!!! example "Pre-extracted lead list"
    ```python
    # Pass a lead table directly (skip re-extraction inside plot_phenogram)
    lead_ss = mysumstats.get_lead(gls=True, verbose=False)
    lead_ss.plot_phenogram(use_lead_extraction=False, anno=True)
    ```

!!! example "Group mode with leads from three traits"
    ```python
    import pandas as pd

    # Shared column mapping for all three BBJ files (hg19)
    LOAD_KW = dict(
        snpid="SNP", chrom="CHR", pos="POS",
        ea="ALT", nea="REF", neaf="Frq", p="P",
        build="19", verbose=False,
    )

    def load_and_leads(path, trait):
        ss = gl.Sumstats(path, **LOAD_KW)
        ss.fix_chr(verbose=False)
        leads = ss.get_lead(verbose=False)
        leads["TRAIT"] = trait  # trait label for marker shape/color
        return leads

    # Merge leads from T2D and BMI (male/female)
    combined_leads = pd.concat([
        load_and_leads("examples/0_sample_data/t2d_bbj.txt.gz", "T2D"),
        load_and_leads("examples/0_sample_data/bmi_male_bbj.txt.gz", "BMI (male)"),
        load_and_leads("examples/0_sample_data/bmi_female_bbj.txt.gz", "BMI (female)"),
    ], ignore_index=True)

    lead_ss = gl.Sumstats(combined_leads, build="19", verbose=False)

    # Add nearest gene name column (required for anno_group="GENE")
    lead_ss.data = lead_ss.anno_gene(source="ensembl", verbose=False)

    lead_ss.plot_phenogram(
        use_lead_extraction=False,   # input is already a lead table
        anno_group="GENE",           # group markers by nearest gene
        anno_source="ensembl",
        anno_shape="TRAIT",          # marker shape per trait
        anno_color="TRAIT",          # marker color per trait
        marker_shapes=["o", "^", "s"],
        marker_colors=["#CB132D", "#597FBD", "#2ECC71"],
        show_legend=True,
        legend_ncol=3,
        figsize=(12, 20),
        dpi=100,
        verbose=False,
    )
    ```

!!! example "Group mode: multi-subtype and multi-ancestry (simulated)"
    ```python
    import pandas as pd

    def add_row(rows, rs_i, chrom, pos, locus, subtype, ancestry, p=1e-12):
        rows.append({
            "SNPID": f"rs{rs_i}", "CHR": chrom, "POS": int(pos), "P": p,
            "LOCUS": locus, "SUBTYPE": subtype, "ANCESTRY": ancestry,
        })
        return rs_i + 1

    rows, rs_i = [], 1

    # Shared loci (partial overlap — not every subtype × ancestry pair)
    for subtype, ancestry in [("TypeA", "EUR"), ("TypeB", "EAS"), ("TypeC", "AFR"), ("TypeA", "SAS")]:
        rs_i = add_row(rows, rs_i, 16, 53_800_000 + rs_i * 1_000, "FTO", subtype, ancestry)
    for subtype, ancestry in [("TypeA", "EUR"), ("TypeC", "SAS")]:
        rs_i = add_row(rows, rs_i, 10, 114_700_000 + rs_i * 1_000, "TCF7L2", subtype, ancestry)
    for ancestry in ["EUR", "AFR"]:
        rs_i = add_row(rows, rs_i, 19, 45_410_000 + rs_i * 1_000, "APOE", "TypeB", ancestry)

    # Solo loci: one lead per remaining subtype + ancestry pair
    for chrom, pos, locus, subtype, ancestry in [
        (7, 44_000_000, "GCK", "TypeA", "EAS"),
        (3, 123_000_000, "PPARG", "TypeB", "EUR"),
        (18, 60_000_000, "MC4R", "TypeA", "EUR"),
        # ... additional pair-specific loci
    ]:
        rs_i = add_row(rows, rs_i, chrom, pos, locus, subtype, ancestry)

    sim_ss = gl.Sumstats(pd.DataFrame(rows), build="19", verbose=False)
    sim_ss.plot_phenogram(
        use_lead_extraction=False,
        anno_group="LOCUS",
        anno_shape="SUBTYPE",
        anno_color="ANCESTRY",
        marker_shapes=["o", "^", "s"],
        marker_colors=["#CB132D", "#597FBD", "#2ECC71", "#F39C12"],
        show_legend=True,
        only_anno_chr=True,
    )
    ```

    Full list of solo loci: see [tutorial](visualization_phenogram.md#multi-ancestry-and-multi-subtype-simulated).

See runnable examples [here](visualization_phenogram.md).

For the full parameter list generated from the visualization registry, see the [API reference](api/sumstats/plot.md).
