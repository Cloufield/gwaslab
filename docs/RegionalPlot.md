# Regional plots

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

!!! note "Color issue"
    - gwaslab<=3.4.39 : the color assigned to each variant is actually the color for the lower LD r2 category. For example, variants with LD>0.8 will be colored with the color for 0.8>LD>0.6.
    - gwaslab v3.4.40 : the color for region_ref_second was assigned based on region_ref LD.
    - Solution: Update to new version (>=3.4.41) of gwaslab.

GWASLab provides functions for creating regional plots.

## .plot_mqq(mode="r")
```
.plot_mqq(mode="r",
          region = None,
          ...
          ):
```

GWASLab regional plot function is based on plot_mqq().
Most options are largely the same as [Manhattan plot](https://cloufield.github.io/gwaslab/Visualization/).


## Options

| Option                    | DataType     | Description                                                                                                                          | Default                                                                   |
|---------------------------|--------------|--------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `mode`                    | `r`          | specify regional plot mode                                                                                                           | -                                                                         |
| `region`                  | `tuple`      | a three elements tuple (chr, start, end); for example, (7,156538803,157538803)                                                       | -                                                                         |
| `vcf_path`                | `string`     | path to LD reference in VCF format: if None, LD information will not be plotted.                                                     | `None`                                                                    |
| `region_ref`              | `list`       | the SNPID or rsID for reference variants; if None, lead variants will be selected; support up to 7 reference markers (since v3.4.47) | [`None`]                                                                  |
| `region_grid`             | `boolean`    | If True, plot the grid line                                                                                                          | `False`                                                                   |
| `region_grid_line`        | `dict`       | parameters for the grid line                                                                                                         | `{"linewidth": 2,"linestyle":"--"}`                                       |
| `region_lead_grid`        | `boolean`    | If True, plot a line to show the reference variants                                                                                  | `True`                                                                    |
| `region_lead_grid_line`   | `dict`       | parameters for the line to show the reference variants                                                                               | `{"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}`       |
| `region_ld_threshold`     | `list`       | LD r2 categories                                                                                                                     | `[0.2,0.4,0.6,0.8]`                                                       |
| `region_ld_colors`        | `list`       | LD r2 categories colors for single reference marker                                                                                  | `["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"]` |
| `region_ld_colors_m`      | `list`       | list of colors used for multiple reference markers (since v3.4.47)                                                                   | `["#E51819","#367EB7","green","#F07818","#AD5691","yellow","purple"]`     |
| `region_marker_shapes`    | `list`       | list of shapes used for multiple reference markers (since v3.4.47)                                                                   | `['o', 's','^','D','*','P','X','h','8']`                                  |
| `region_chromatin_files`  | `list`       | list of paths of Roadmap `15_coreMarks_mnemonics.bed.gz` files                                                                       | []                                                                        |
| `region_chromatin_labels` | `list`       | list of labels for region_chromatin_files                                                                                            | []                                                                        |
| `region_hspace`           | `float`      | the space between the scatter plot and the gene track                                                                                | `0.02`                                                                    |
| `region_step`             | `int`        | number of X axis ticks                                                                                                               | `21`                                                                      |
| `region_recombination`    | `boolean`    |                                                                                                                                      | `True`                                                                    |
| `tabix`                   | `string`     | path to tabix; if None, GWASLab will search in environmental path; Note: if tabix is available, the speed is much faster!!!          | `None`                                                                    |
| `taf`                     | `list`       | a five-element list; number of gene track lanes, offset for gene track, font_ratio, exon_ratio, text_offset                          | `[4,0,0.95,1,1]`                                                          |
| `build`                   | `19` or `38` | reference genome build; `99` for unknown                                                                                             | `99`                                                                      |


!!! quote "Calculation of LD r2"
    The calculation is based on [Rogers and Huff r implemented in scikit-alle](https://scikit-allel.readthedocs.io/en/stable/stats/ld.html). Variants in reference vcf file should be biallelic format. Unphased data is acceptable. AF information is not needed. Variant ID is not required. Missing genotype is allowed.


## Examples

!!! example
    See [Regional plot](https://cloufield.github.io/gwaslab/visualization_regional/)


## gl.plot_stacked_mqq()

Creates stacked Manhattan-QQ plots or regional plots for multiple GWAS datasets, allowing side-by-side comparison of multiple studies or traits.

```
gl.plot_stacked_mqq(objects, **kwargs)
```

### Parameters

#### Required Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `objects` | `list` | List of `gl.Sumstats` objects or pandas DataFrames containing GWAS summary statistics | Required |

#### Plot Mode Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `mode` | `str` | Plot mode: `"r"` for regional plots, `"m"` for Manhattan plots, `"mqq"` for Manhattan-QQ plots | `"r"` |
| `pm` | `list` | List of panel modes for each object: `"m"` for Manhattan, `"pip"` for PIP/credible sets | `None` (auto-detected) |
| `region` | `tuple` | For regional plots: three-element tuple (chr, start, end), e.g., `(7, 156538803, 157538803)` | `None` |
| `vcfs` | `list` | List of VCF file paths for LD reference. If single VCF provided, it will be used for all panels. For regional plots, must match number of objects or be length 1 | `[]` |

#### Layout Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `titles` | `list` | List of titles for each panel | `None` |
| `title_pos` | `str` or `tuple` | Position of titles | `None` |
| `title_kwargs` | `dict` | Keyword arguments for title styling | `None` |
| `subplot_height` | `float` | Height of each subplot in inches | `4` |
| `region_hspace` | `float` | Space between subplots | `0.07` |
| `mqqratio` | `float` | Width ratio of Manhattan to QQ plot when `mode="mqq"` | `3` |
| `mqq_height` | `float` | Height ratio for Manhattan plot panels | `1` |
| `cs_height` | `float` | Height ratio for credible set (PIP) panels | `0.5` |
| `gene_track_height` | `float` | Height ratio for gene track (regional plots only) | `0.5` |
| `region_chromatin_height` | `float` | Height ratio for chromatin track (regional plots only) | `0.1` |
| `fig_kwargs` | `dict` | Keyword arguments for matplotlib figure creation | `None` |

#### Regional Plot Specific Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `region_chromatin_files` | `list` | List of paths to Roadmap `15_coreMarks_mnemonics.bed.gz` files | `[]` |
| `region_chromatin_labels` | `list` | List of labels for chromatin tracks | `[]` |
| `region_lead_grids` | `list` | List of panel indices to show lead variant grid lines | `None` (all panels) |
| `region_ld_legends` | `list` | List of panel indices to show LD legend | `[0]` |
| `gtf` | `str` | Path to GTF file for gene annotation. Use `"default"` for built-in GTF | `None` |
| `build` | `str` | Reference genome build: `"19"`, `"38"`, or `"99"` for unknown | `"99"` |

#### Styling Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `fontsize` | `float` | Font size for labels and text | `9` |
| `font_family` | `str` | Font family name | `"Arial"` |
| `common_ylabel` | `bool` | If `True`, use common y-axis label for all panels | `True` |

#### Output Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `save` | `str` | Path to save the figure | `None` |
| `save_kwargs` | `dict` | Keyword arguments for saving (dpi, bbox_inches, etc.) | `None` |
| `verbose` | `bool` | Print progress messages | `True` |

#### Additional Parameters

All parameters from `plot_mqq()` can be passed via `**mqq_kwargs` to customize individual panels. These include:

- `highlight`, `anno_set`, `pinpoint` for variant annotation
- `colors`, `scatter_kwargs` for styling
- `sig_line`, `suggestive_sig_line` for significance lines
- `region_ld_threshold`, `region_ld_colors` for LD coloring
- And many more (see [Manhattan plot documentation](https://cloufield.github.io/gwaslab/Visualization/))

### Plot Modes

#### Regional Mode (`mode="r"`)

Creates stacked regional plots with LD information and gene tracks. Each panel shows a regional plot for a specific genomic region. Requires:

- `region` parameter specifying the genomic region
- `vcfs` parameter for LD reference (optional but recommended)

#### Manhattan Mode (`mode="m"`)

Creates stacked Manhattan plots across the genome. All panels share the same x-axis (genomic position). Useful for comparing multiple traits or studies genome-wide.

#### Manhattan-QQ Mode (`mode="mqq"`)

Creates stacked panels with Manhattan plot on the left and QQ plot on the right for each dataset. Useful for quality control and comparing multiple studies.

### Panel Types

The function automatically detects panel types based on input data:

- **Manhattan panels** (`pm="m"`): For DataFrames with `P` or `MLOG10P` columns
- **Credible set panels** (`pm="pip"`): For DataFrames with `PIP` column (e.g., from fine-mapping)

### Examples

**Stacked regional plots for multiple studies:**

```
import gwaslab as gl

# Load multiple sumstats
sumstats1 = gl.Sumstats("study1.txt.gz")
sumstats2 = gl.Sumstats("study2.txt.gz")

# Create stacked regional plot
gl.plot_stacked_mqq(
    objects=[sumstats1, sumstats2],
    mode="r",
    region=(7, 156538803, 157538803),
    vcfs=["ld_ref.vcf.gz"],  # Single VCF used for all panels
    titles=["Study 1", "Study 2"],
    build="38"
)
```

**Stacked Manhattan plots:**

```
# Compare multiple traits genome-wide
gl.plot_stacked_mqq(
    objects=[trait1_sumstats, trait2_sumstats, trait3_sumstats],
    mode="m",
    titles=["Trait 1", "Trait 2", "Trait 3"],
    colors=["#1f77b4", "#ff7f0e", "#2ca02c"],
    sig_line=True,
    sig_level_plot=5e-8
)
```

**Stacked Manhattan-QQ plots:**

```
# Quality control comparison
gl.plot_stacked_mqq(
    objects=[study1, study2],
    mode="mqq",
    titles=["Study 1", "Study 2"],
    mqqratio=3
)
```

**Regional plot with credible sets:**

```
# Include fine-mapping results
sumstats = gl.Sumstats("gwas.txt.gz")
finemap_results = pd.read_csv("finemap_results.txt")  # Contains PIP column

gl.plot_stacked_mqq(
    objects=[sumstats, finemap_results],
    mode="r",
    region=(7, 156538803, 157538803),
    vcfs=["ld_ref.vcf.gz"],
    titles=["GWAS", "Fine-mapping"],
    build="38"
)
```

**Custom styling per panel:**

```
# Different colors and styles for each panel
gl.plot_stacked_mqq(
    objects=[sumstats1, sumstats2],
    mode="m",
    titles=["Panel 1", "Panel 2"],
    colors=[["#1f77b4"], ["#ff7f0e"]],  # Different colors per panel
    scatter_kwargs=[{"s": 10}, {"s": 20}]  # Different sizes per panel
)
```

### Notes

- When `mode="r"` (regional), the number of VCF files must match the number of objects, or a single VCF can be provided which will be used for all panels.
- For credible set panels (`pm="pip"`), VCF files are automatically set to `"NA"` as LD information is not applicable.
- The function automatically detects panel types based on column names (`P`/`MLOG10P` for Manhattan, `PIP` for credible sets).
- All parameters from `plot_mqq()` can be passed to customize individual panels. Use lists to specify different values for each panel.

### Notes and Troubleshooting

**Gene track limitations:**

The gene track only displays protein-coding genes from the reference GTF files. Non-coding genes, pseudogenes, and other gene types are excluded from the visualization.

**Missing exons in gene track:**

Very short exons may not appear in the gene track if they are too small to render at the current resolution. To display these exons, you can either increase the figure DPI (dots per inch) when saving the plot or reduce the length of the genomic region being plotted.

**LD calculation errors:**

LD (linkage disequilibrium) cannot be calculated when the reference variant in the VCF file is mono-allelic (i.e., has only one allele present in the reference panel). This will result in an error even if both variants are present in the reference VCF file. Ensure your reference VCF contains biallelic variants for accurate LD calculation.