# LD Block Plot

The `plot_ld_block()` function visualizes linkage disequilibrium (LD) matrices as 45°-rotated inverted triangles, providing an intuitive representation of LD structure across genomic regions.

## Overview

The LD block plot displays the upper triangle of an LD matrix (typically r² values) as a rotated inverted triangle, where:
- The x-axis represents genomic position
- The y-axis represents the rotated triangle structure
- Color intensity represents LD strength (r² values)
- Each cell shows the LD between a pair of variants

## Two Modes

### 1. Standalone Mode

Creates an independent figure with the LD block plot. Use this when you want to visualize LD structure separately.

```python
gl.plot_ld_block(
    ld=ld_matrix,  # or use vcf_path + region
    pos=positions,
    ...
)
```

### 2. Regional Mode

Integrates with `plot_mqq(mode="r")` to add an LD block visualization below the regional plot. The LD block automatically aligns with the regional plot's x-axis using the "i" coordinate system.

```python
mysumstats.plot_mqq(
    mode="r",
    region=(7, 156538803, 157538803),
    vcf_path="path/to/reference.vcf.gz",
    ld_block=True,  # Enable LD block visualization
    ...
)
```

## Function Signature

```python
gl.plot_ld_block(
    ld=None,
    pos=None,
    vcf_path=None,
    region=None,
    sumstats=None,
    pos_col="POS",
    nea_col="NEA",
    ea_col="EA",
    tabix=None,
    mapper=None,
    ax=None,
    ax_pos=None,
    mode=None,
    cmap=None,
    vmin=0.0,
    vmax=1.0,
    xlabel="Genomic position",
    title=None,
    cbar=True,
    cbar_label="LD (r²)",
    cbar_kwargs=None,
    fig_kwargs=None,
    save=None,
    save_kwargs=None,
    fontsize=10,
    font_family="Arial",
    region_step=21,
    lead_snp_is=None,
    lead_snp_is_color=None,
    anno_cell=False,
    anno_cell_fmt="{:.2f}",
    anno_cell_kwargs=None,
    ld_block_grid=False,
    ld_block_grid_kwargs=None,
    ld_block_anno=False,
    ld_block_anno_kwargs=None,
    ld_block_anno_set=None,
    ld_block_anno_max_rows=100,
    log=Log(),
    verbose=True,
    **kwargs
)
```

## Parameters

### Data Input Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `ld` | `np.ndarray` | (n, n) symmetric LD matrix (e.g., r² values). If None, will be extracted from VCF when `vcf_path` and `region` are provided | `None` |
| `pos` | `np.ndarray` | Length n array of genomic positions (bp or index). If None, uses `np.arange(n)` or extracts from sumstats when using VCF | `None` |
| `vcf_path` | `str` | Path to VCF file for LD calculation. If provided, `region` must also be provided | `None` |
| `region` | `tuple` | Region specification (chromosome, start, end) for extracting LD from VCF. Example: `(1, 1000000, 2000000)` for chr1:1000000-2000000 | `None` |
| `sumstats` | `pd.DataFrame` or `Sumstats` | Sumstats object or DataFrame for matching variants with VCF. Required when using `vcf_path`. Should contain POS, NEA, and EA columns | `None` |
| `pos_col` | `str` | Column name for position in sumstats | `"POS"` |
| `nea_col` | `str` | Column name for non-effect allele in sumstats | `"NEA"` |
| `ea_col` | `str` | Column name for effect allele in sumstats | `"EA"` |
| `tabix` | `bool` | Whether to use tabix indexing for VCF. If None, auto-detects | `None` |
| `mapper` | `ChromosomeMapper` | ChromosomeMapper instance for chromosome name conversion. If None, creates a default mapper with automatic format detection | `None` |

### Plotting Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `ax` | `matplotlib.axes.Axes` | Axes object to plot on. If None, creates a new figure (standalone mode). If provided with region/sumstats/i-coordinates, uses regional mode | `None` |
| `ax_pos` | `matplotlib.axes.Axes` | Axes object for position bar (optional). If provided, displays genomic positions above the LD block | `None` |
| `mode` | `str` | Plot mode: 'regional' or 'standalone'. If None, automatically determines based on whether ax is provided and i-coordinates are available | `None` (auto-detect) |
| `cmap` | `str` or `Colormap` | Colormap for LD values. Default: continuous gradient matching region plot LD categories (dark blue → light blue → green → orange → red) | `None` |
| `vmin` | `float` | Minimum value for colormap scaling | `0.0` |
| `vmax` | `float` | Maximum value for colormap scaling | `1.0` |
| `xlabel` | `str` | Label for x-axis | `"Genomic position"` |
| `title` | `str` | Title for the plot | `None` |
| `fontsize` | `float` | Font size for labels | `10` |
| `font_family` | `str` | Font family | `"Arial"` |
| `region_step` | `int` | Number of x-axis tick positions for regional mode | `21` |

### Colorbar Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `cbar` | `bool` | Whether to draw a colorbar | `True` |
| `cbar_label` | `str` | Label for the colorbar | `"LD (r²)"` |
| `cbar_kwargs` | `dict` | Additional arguments for colorbar | `None` |

### Annotation Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `lead_snp_is` | `list` | List of lead SNP i-coordinates for colored annotation lines | `None` |
| `lead_snp_is_color` | `list` | List of colors corresponding to lead SNPs | `None` |
| `anno_cell` | `bool` | Whether to annotate cells with LD r² values | `False` |
| `anno_cell_fmt` | `str` | Format string for cell annotations (e.g., `"{:.2f}"` for 2 decimal places) | `"{:.2f}"` |
| `anno_cell_kwargs` | `dict` | Additional keyword arguments for text annotations (e.g., `{'fontsize': 10, 'weight': 'bold'}`) | `None` |
| `ld_block_anno` | `bool` or `str` | Whether to add annotations on the left side of the LD block triangle. If `True`, uses 'chr:pos' format. If a string (e.g., 'SNPID', 'rsID'), uses that column for annotation text | `False` |
| `ld_block_anno_kwargs` | `dict` | Additional keyword arguments for left-side annotations | `None` |
| `ld_block_anno_set` | `list`, `set`, or `np.ndarray` | List of SNPIDs to annotate. If None, annotates all variants | `None` |
| `ld_block_anno_max_rows` | `int` | Maximum number of variants to annotate. If exceeded, annotations will be skipped | `100` |

### Grid Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `ld_block_grid` | `bool` | Whether to draw grid lines on the LD block | `False` |
| `ld_block_grid_kwargs` | `dict` | Additional keyword arguments for grid lines | `None` |

### Figure and Saving Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `fig_kwargs` | `dict` | Additional arguments for figure creation | `None` |
| `save` | `bool` or `str` | Whether to save the figure. If str, path to save file | `None` |
| `save_kwargs` | `dict` | Additional arguments for saving the figure | `None` |

### Other Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `log` | `Log` | Logger instance | `Log()` |
| `verbose` | `bool` | Whether to print log messages | `True` |
| `**kwargs` | `dict` | Additional arguments passed to `pcolormesh` | - |

## Examples

### Example 1: Standalone LD Block from VCF

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("my_sumstats.tsv")

# Plot LD block from VCF
gl.plot_ld_block(
    vcf_path="path/to/reference.vcf.gz",
    region=(7, 156538803, 157538803),
    sumstats=mysumstats,
    save="ld_block_plot.png"
)
```

### Example 2: Standalone LD Block from Pre-computed Matrix

```python
import numpy as np
import gwaslab as gl

# Pre-computed LD matrix (r² values)
ld_matrix = np.array([[1.0, 0.8, 0.5],
                      [0.8, 1.0, 0.7],
                      [0.5, 0.7, 1.0]])

# Genomic positions
positions = np.array([156538803, 156540000, 156541000])

# Plot LD block
gl.plot_ld_block(
    ld=ld_matrix,
    pos=positions,
    anno_cell=True,  # Show r² values in cells
    save="ld_block_standalone.png"
)
```

### Example 3: LD Block in Regional Plot

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("my_sumstats.tsv")

# Create regional plot with LD block
mysumstats.plot_mqq(
    mode="r",
    region=(7, 156538803, 157538803),
    vcf_path="path/to/reference.vcf.gz",
    ld_block=True,  # Enable LD block visualization
    ld_block_grid=True,  # Add grid lines
    anno_cell=True,  # Show r² values in cells
    save="regional_with_ld_block.png"
)
```

### Example 4: Customized LD Block with Annotations

```python
import gwaslab as gl

mysumstats = gl.Sumstats("my_sumstats.tsv")

gl.plot_ld_block(
    vcf_path="path/to/reference.vcf.gz",
    region=(7, 156538803, 157538803),
    sumstats=mysumstats,
    ld_block_anno="SNPID",  # Annotate with SNPID on left side
    ld_block_anno_set=["rs123", "rs456"],  # Only annotate specific variants
    lead_snp_is=[0, 5, 10],  # Highlight lead SNPs at indices 0, 5, 10
    lead_snp_is_color=["red", "blue", "green"],
    anno_cell=True,
    anno_cell_fmt="{:.3f}",  # 3 decimal places
    ld_block_grid=True,
    save="ld_block_annotated.png"
)
```

### Example 5: LD Block with Custom Colormap and Value Range

```python
import gwaslab as gl
import matplotlib.pyplot as plt

mysumstats = gl.Sumstats("my_sumstats.tsv")

gl.plot_ld_block(
    vcf_path="path/to/reference.vcf.gz",
    region=(7, 156538803, 157538803),
    sumstats=mysumstats,
    cmap="viridis",  # Use viridis colormap
    vmin=0.0,  # Minimum LD value to display
    vmax=0.8,  # Maximum LD value to display (focus on high LD)
    cbar=True,
    cbar_label="LD r² (0-0.8)",
    save="ld_block_custom_cmap.png"
)
```

### Example 6: LD Block without Colorbar

```python
import gwaslab as gl

mysumstats = gl.Sumstats("my_sumstats.tsv")

gl.plot_ld_block(
    vcf_path="path/to/reference.vcf.gz",
    region=(7, 156538803, 157538803),
    sumstats=mysumstats,
    cbar=False,  # Disable colorbar
    title="LD Block (no colorbar)",
    save="ld_block_no_cbar.png"
)
```

## Integration with plot_mqq()

The LD block can be seamlessly integrated into regional plots using the `ld_block=True` option in `plot_mqq(mode="r")`. When enabled:

1. **Automatic Layout**: The LD block is placed below the regional plot with proper spacing
2. **X-axis Alignment**: The LD block's x-axis automatically aligns with the regional plot using the "i" coordinate system
3. **Position Bar**: A position bar is added above the LD block showing genomic positions
4. **Annotation Lines**: Vertical lines connect variants in the position bar to their corresponding positions in the LD block

### Parameters for LD Block in plot_mqq()

When using `ld_block=True` in `plot_mqq()`, you can control the LD block appearance with these parameters:

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `ld_block` | `bool` | Enable/disable LD block visualization | `False` |
| `anno_cell` | `bool` | Show r² values in cells | `False` |
| `anno_cell_fmt` | `str` | Format string for cell annotations (e.g., `"{:.2f}"` for 2 decimal places) | `"{:.2f}"` |
| `anno_cell_kwargs` | `dict` | Additional arguments for cell text annotations | `None` |
| `ld_block_grid` | `bool` | Draw grid lines on the LD block | `False` |
| `ld_block_grid_kwargs` | `dict` | Additional arguments for grid lines (e.g., `{'linewidth': 1, 'alpha': 0.3}`) | `None` |
| `ld_block_anno` | `bool` or `str` | Add annotations on the left side. If `True`, uses 'chr:pos' format. If a string (e.g., 'SNPID', 'rsID'), uses that column | `False` |
| `ld_block_anno_kwargs` | `dict` | Additional arguments for left-side annotations | `None` |
| `ld_block_anno_set` | `list`, `set`, or `np.ndarray` | List of SNPIDs to annotate. If None, annotates all variants | `None` |
| `ld_block_anno_max_rows` | `int` | Maximum number of variants to annotate. If exceeded, annotations are skipped | `100` |

## Visual Features

### Position Bar
In standalone mode, a position bar is automatically created above the LD block showing genomic positions. Vertical lines connect variants in the position bar to their corresponding positions in the LD block triangle.

### Colorbar
By default, a horizontal colorbar is displayed in the lower right corner of the LD block, showing the LD (r²) value scale. The colorbar can be disabled by setting `cbar=False`.

### Grid Lines
Grid lines can be added to help visualize LD block boundaries using `ld_block_grid=True`. Customize grid appearance with `ld_block_grid_kwargs`.

### Cell Annotations
LD r² values can be displayed directly in each cell using `anno_cell=True`. Format the values with `anno_cell_fmt` (e.g., `"{:.2f}"` for 2 decimal places).

### Left-Side Annotations
Variant annotations can be added on the left side of the LD block triangle using `ld_block_anno`. Options include:
- `ld_block_anno=True`: Uses 'chr:pos' format
- `ld_block_anno="SNPID"`: Uses SNPID column from sumstats
- `ld_block_anno="rsID"`: Uses rsID column from sumstats

Use `ld_block_anno_set` to specify which variants to annotate, and `ld_block_anno_max_rows` to limit the number of annotations.

### Lead SNP Markers
Lead SNPs can be highlighted with colored vertical lines using `lead_snp_is` (list of i-coordinates) and `lead_snp_is_color` (corresponding colors).

## Notes

!!! note "LD Calculation"
    The LD calculation is based on [Rogers and Huff r implemented in scikit-allel](https://scikit-allel.readthedocs.io/en/stable/stats/ld.html). Variants in the reference VCF file should be in biallelic format. Unphased data is acceptable. AF information is not needed. Variant ID is not required. Missing genotypes are allowed.

!!! tip "Performance"
    Using tabix-indexed VCF files significantly speeds up LD extraction. GWASLab will auto-detect tabix availability, or you can specify it explicitly with the `tabix` parameter.

!!! warning "Coordinate System"
    In regional mode, the LD block uses the "i" coordinate system to align with the regional plot. This ensures that variants in both plots are correctly aligned even if there are gaps in the genomic positions. The "i" coordinate represents the ordered index of variants in the region, providing consistent spacing regardless of actual genomic distances.

!!! info "Standalone Mode Figure Structure"
    In standalone mode, the function creates a figure with 3 axes:
    1. Position bar axes (top) - displays genomic positions
    2. LD block axes (main) - displays the rotated triangle
    3. Colorbar axes (inset) - displays the color scale (if `cbar=True`)

