# Panel System

The Panel system in GWASLab provides a flexible framework for creating stacked multi-panel visualizations. It allows you to combine different types of genomic visualizations (tracks, arcs, regional plots, LD blocks, chromatin states, and more) into a single aligned figure with shared x-axis coordinates.

## Overview

The Panel system consists of two main components:

1. **`Panel` class**: Stores configuration for individual panels
2. **`plot_panels()` function**: Creates a stacked figure from multiple Panel objects

This system enables you to create publication-quality multi-panel figures that combine GWAS results with genomic annotations, chromatin states, contact maps, and other genomic data.

## Panel Class

The `Panel` class stores panel configuration information for stacked plots. Each panel has a type and associated parameters that are processed and validated before plotting.

### Creating a Panel

```python
import gwaslab as gl

# Create a panel with panel type and parameters
panel = gl.Panel(
    panel_type="track",  # Type of panel
    track_path="genes.gtf",  # Panel-specific parameters
    region=(1, 1000000, 2000000),
    color="#020080"
)
```

### Panel Types

GWASLab supports the following panel types:

| Panel Type | Description | Required Parameters |
|------------|-------------|---------------------|
| `track` | Genomic tracks from GTF, BED, bigWig, or bigBed files | `track_path`, `region` |
| `arc` | BEDPE arc plots for contact maps | `bedpe_path`, `region` |
| `ld_block` | LD block plots (45° rotated inverted triangle) | `vcf_path`, `region`, `sumstats` or `insumstats` |
| `region` | Regional plots with LD information and gene track | `sumstats` or `insumstats`, `region`, `vcf_path` |
| `chromatin` | Chromatin state tracks (Roadmap 15-state model) | `region_chromatin_files`, `region_chromatin_labels`, `region` |
| `pipcs` | PIP (Posterior Inclusion Probability) and Credible Sets plots | `pipcs_raw`, `region` |

### Panel Methods

The `Panel` class provides several methods for managing panel parameters:

```python
# Get panel type
panel_type = panel.get_type()

# Get all parameters
kwargs = panel.get_kwargs()

# Get a specific parameter
value = panel.get_kwarg("color", default=None)

# Set a parameter
panel.set_kwarg("color", "#FF0000")

# Update multiple parameters
panel.update_kwargs(color="#FF0000", alpha=0.5)
```

### Creating Panels from Sumstats Objects

If you're working with a `Sumstats` object, you can use the `.Panel()` method to automatically pass the sumstats data:

```python
mysumstats = gl.Sumstats("data.txt.gz")

# Create a region panel - sumstats data is automatically passed
panel1 = mysumstats.Panel(
    "region",
    region=(1, 1000000, 2000000),
    vcf_path="ld.vcf.gz",
    build="38"
)

# Create an LD block panel
panel2 = mysumstats.Panel(
    "ld_block",
    region=(1, 1000000, 2000000),
    vcf_path="ld.vcf.gz"
)
```

## plot_panels() Function

The `plot_panels()` function creates a stacked figure from a list of Panel objects with aligned x-axes.

### Basic Usage

```python
import gwaslab as gl

# Create multiple panels
panel1 = gl.Panel("track", track_path="genes.gtf", region=(1, 1000000, 2000000))
panel2 = gl.Panel("arc", bedpe_path="contacts.bedpe.gz", region=(1, 1000000, 2000000))
panel3 = gl.Panel("region", insumstats=sumstats, region=(1, 1000000, 2000000), vcf_path="ld.vcf.gz")

# Plot panels together
fig, axes = gl.plot_panels(
    [panel1, panel2, panel3],
    region=(1, 1000000, 2000000),
    titles=["Genes", "Contacts", "Regional Plot"],
    save="stacked_panels.png"
)
```

### Parameters

#### Required Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `panels` | `list` | List of `Panel` objects to plot | Required |

#### Optional Parameters

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `region` | `tuple` | Genomic region as (chromosome, start, end) in 1-based coordinates. If None, extracted from panels. | `None` |
| `height_ratios` | `list` | Height ratios for each panel. If None, defaults are used based on panel type. | `None` |
| `hspace` | `float` | Space between subplots (vertical spacing) | `0.07` |
| `subplot_height` | `float` | Height of each subplot in inches | `1.0` |
| `titles` | `list` | List of titles for each panel | `None` |
| `title_pos` | `str` or `tuple` | Position of titles: "left", "right", "center", or (x, y) tuple | `None` |
| `title_kwargs` | `dict` | Keyword arguments for title styling | `None` |
| `fig_kwargs` | `dict` | Additional keyword arguments for matplotlib figure (e.g., `{'figsize': (10, 8), 'dpi': 200}`) | `None` |
| `save` | `str`, `bool`, or `None` | File path to save figure, or `True` for default path | `None` |
| `save_kwargs` | `dict` | Additional arguments for saving (e.g., `{'dpi': 300, 'bbox_inches': 'tight'}`) | `None` |
| `fontsize` | `int` | Default font size for labels | `9` |
| `font_family` | `str` | Font family for text | `"Arial"` |
| `align_xaxis` | `bool` | Whether to align x-axes across all panels | `True` |
| `region_step` | `int` | Number of ticks to generate for x-axis when `align_xaxis=True` | `21` |
| `track_start_i` | `float` | X-axis offset for track-based plots | `0.0` |
| `verbose` | `bool` | Whether to show progress messages | `True` |
| `log` | `Log` | Logger instance for messages | `Log()` |

### Default Height Ratios

When `height_ratios` is not specified, default ratios are used based on panel type:

- `region`: 4.0 (larger for regional plots with scatter and gene track)
- `ld_block`: 4.0 (medium-large for LD block visualization)
- `track`: 2.0 (compact for genomic tracks)
- `arc`: 2.0 (compact for contact arcs)
- `chromatin`: 1.0 (very compact for chromatin states)
- `pipcs`: 3.0 (medium for PIP and Credible Sets plots)

Note: Panels with multiple axes (e.g., `region` and `ld_block`) automatically split their height ratio between their axes.

### Automatic Figure Sizing

The `plot_panels()` function automatically adjusts figure size and DPI based on:
- Number of panels
- Panel types (complex panels like `region` and `ld_block` get higher DPI)
- Total number of subplots

You can override this by providing `fig_kwargs` with explicit `figsize` and `dpi` values.

## Examples

### Example 1: Basic Multi-Panel Plot

Create a stacked figure with gene tracks, contact arcs, and a regional plot:

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("gwas_results.txt.gz")

# Define region
region = (1, 1000000, 2000000)

# Create panels
gene_panel = gl.Panel(
    "track",
    track_path="genes.gtf",
    region=region,
    color="#020080"
)

contact_panel = gl.Panel(
    "arc",
    bedpe_path="hic_contacts.bedpe.gz",
    region=region,
    color="#FF0000",
    alpha=0.3
)

regional_panel = mysumstats.Panel(
    "region",
    region=region,
    vcf_path="ld_ref.vcf.gz",
    build="38"
)

# Plot all panels
fig, axes = gl.plot_panels(
    [gene_panel, contact_panel, regional_panel],
    region=region,
    titles=["Genes", "Hi-C Contacts", "GWAS Regional Plot"],
    save="multi_panel_plot.png"
)
```

### Example 2: Panel with Chromatin States

Add chromatin state tracks to your visualization:

```python
import gwaslab as gl

# Create chromatin panel
chromatin_panel = gl.Panel(
    "chromatin",
    region_chromatin_files=["E098_15_coreMarks_mnemonics.bed.gz"],
    region_chromatin_labels=["E098 (Lung)"],
    region=(1, 1000000, 2000000)
)

# Combine with other panels
fig, axes = gl.plot_panels(
    [chromatin_panel, regional_panel],
    region=(1, 1000000, 2000000),
    titles=["Chromatin States", "Regional Plot"]
)
```

### Example 3: LD Block Panel

Visualize LD blocks alongside your regional plot:

```python
import gwaslab as gl

# Create LD block panel
ld_panel = mysumstats.Panel(
    "ld_block",
    region=(1, 1000000, 2000000),
    vcf_path="ld_ref.vcf.gz"
)

# Create regional panel
regional_panel = mysumstats.Panel(
    "region",
    region=(1, 1000000, 2000000),
    vcf_path="ld_ref.vcf.gz",
    build="38"
)

# Plot together
fig, axes = gl.plot_panels(
    [ld_panel, regional_panel],
    region=(1, 1000000, 2000000),
    titles=["LD Blocks", "Regional Plot"]
)
```

### Example 4: Custom Height Ratios and Styling

Customize panel heights and styling:

```python
import gwaslab as gl

# Create panels
panels = [
    gl.Panel("track", track_path="genes.gtf", region=region),
    gl.Panel("arc", bedpe_path="contacts.bedpe.gz", region=region),
    mysumstats.Panel("region", region=region, vcf_path="ld.vcf.gz")
]

# Plot with custom settings
fig, axes = gl.plot_panels(
    panels,
    region=region,
    height_ratios=[1.0, 1.0, 5.0],  # Make regional plot larger
    titles=["Genes", "Contacts", "Regional Plot"],
    title_pos="left",
    title_kwargs={"fontsize": 12, "fontweight": "bold"},
    fig_kwargs={"figsize": (12, 10), "dpi": 300},
    save="custom_panels.png"
)
```

### Example 5: PIP and Credible Sets Panel

Visualize fine-mapping results:

```python
import gwaslab as gl
import pandas as pd

# Load PIP/CS data (must have CHR, POS, PIP, CS columns)
pipcs_data = pd.read_csv("fine_mapping_results.txt", sep="\t")

# Create PIPCS panel
pipcs_panel = gl.Panel(
    "pipcs",
    pipcs_raw=pipcs_data,
    region=(1, 1000000, 2000000)
)

# Combine with regional plot
fig, axes = gl.plot_panels(
    [pipcs_panel, regional_panel],
    region=(1, 1000000, 2000000),
    titles=["Fine-mapping (PIP/CS)", "Regional Plot"]
)
```

## Panel Parameter Processing

Panel parameters are processed through GWASLab's parameter management system, which:

1. **Merges defaults**: Combines user-provided parameters with default values
2. **Filters allowed parameters**: Removes parameters not accepted by the plotting function
3. **Validates arguments**: Ensures required parameters are present

This means you can use the same parameters as the corresponding `plot_*` methods (e.g., `plot_track`, `plot_arc`, `plot_region`).

## Multi-Axis Panels

Some panel types require multiple axes:

- **`region`**: Uses 2 axes (main scatter plot + gene track)
- **`ld_block`**: Uses 2 axes (position bar + main LD block)

The `plot_panels()` function automatically handles these multi-axis panels and adjusts height ratios accordingly.

## X-Axis Alignment

By default, `plot_panels()` aligns x-axes across all panels so they share the same genomic coordinate system. This ensures:

- Consistent tick positions
- Shared x-axis limits
- Proper alignment of features across panels

You can disable this with `align_xaxis=False`, though this is generally not recommended.

## Tips and Best Practices

1. **Consistent regions**: All panels should use the same genomic region for proper alignment
2. **Parameter validation**: Panel parameters are validated when the Panel is created, not when plotted
3. **Data preservation**: Use `insumstats` instead of `sumstats` for regional panels to preserve original data
4. **File paths**: Ensure all file paths (GTF, BED, VCF, etc.) are accessible
5. **Memory management**: Large datasets may require significant memory for multi-panel plots
6. **Figure size**: Let the function auto-adjust figure size, or provide explicit `figsize` in `fig_kwargs`

## See Also

- [Regional Plot Documentation](RegionalPlot.md) - For details on regional plot parameters
- [Visualization Guide](Visualization.md) - For general visualization options
- [GTF Documentation](GTF.md) - For track file formats
- [VCF Documentation](VCF.md) - For LD reference formats
