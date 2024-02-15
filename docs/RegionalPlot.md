# Regional plots

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

!!! warning "Color issue"
    - gwaslab<=3.4.39 : the color assigned to each variant is actually the color for the lower LD r2 category. For example, variants with LD>0.8 will be colored with the color for 0.8>LD>0.6.
    - Solution: Update to new version (>=3.4.40) of gwaslab.

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

| Option                  | DataType     | Description                                                                                                                 | Default                                                                   |
|-------------------------|--------------|-----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `mode`                  | `r`          | specify regional plot mode                                                                                                  | -                                                                         |
| `region`                | `tuple`      | a three elements tuple (chr, start, end); for example, (7,156538803,157538803)                                              | -                                                                         |
| `vcf_path`              | `string`     | path to LD reference in VCF format: if None, LD information will not be plotted.                                            | `None`                                                                    |
| `region_ref`            | `boolean`    | the SNPID or rsID (if SNPID is not available) for reference variant; if None, lead variants will be selected                | `None`                                                                    |
| `region_ref2`           | `boolean`    | the SNPID or rsID for the second reference variant                                                                          | `None`                                                                    |
| `region_grid`           | `boolean`    | If True, plot the grid line                                                                                                 | `False`                                                                   |
| `region_grid_line`      | `dict`       | parameters for the grid line                                                                                                | `{"linewidth": 2,"linestyle":"--"}`                                       |
| `region_lead_grid`      | `string`     | If True, plot a line to show the reference variants                                                                         | -                                                                         |
| `region_lead_grid_line` | `string`     | parameters for the line to show the reference variants                                                                      | {"alpha":0.5,"linewidth" : 2,"linestyle":"--","color":"#FF0000"}          |
| `region_ld_threshold`   | `list`       | LD r2 categories                                                                                                            | `[0.2,0.4,0.6,0.8]`                                                       |
| `region_ld_colors`      | `list`       | LD r2 categories colors                                                                                                     | `["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"]` |
| `region_ld_colors1`     | `list`       | LD r2 categories colors for the first reference varaint (when region_ref2 is specified)                                     | `["#E4E4E4","#F8CFCF","#F5A2A5","#F17474","#EB4445","#E51819","#E51819"]` |
| `region_ld_colors2`     | `list`       | LD r2 categories colors for the second reference varaint (when region_ref2 is specified)                                    | `["#E4E4E4","#D8E2F2","#AFCBE3","#86B3D4","#5D98C4","#367EB7","#367EB7"]` |
| `region_hspace`         | `float`      | the space between the scatter plot and the gene track                                                                       | `0.02`                                                                    |
| `region_step`           | `int`        | number of X axis ticks                                                                                                      | `21`                                                                      |
| `region_recombination`  | `boolean`    |                                                                                                                             | `True`                                                                    |
| `tabix`                 | `string`     | path to tabix; if None, GWASLab will search in environmental path; Note: if tabix is available, the speed is much faster!!! | `None`                                                                    |
| `taf`                   | `list`       | a five-element list; number of gene track lanes, offset for gene track, font_ratio, exon_ratio, text_offset                 | `[4,0,0.95,1,1]`                                                          |
| `build`                 | `19` or `38` | reference genome build; `99` for unknown                                                                                    | `99`                                                                      |


!!! info "Calculation of LD r2"
    The calculation is based on [Rogers and Huff r implemented in scikit-alle](https://scikit-allel.readthedocs.io/en/stable/stats/ld.html). Variants in refernece vcf file should be biallelic format. Unphased data is acceptable. AF information is not needed. Variant ID is not required. Missing genotype is allowed.


## Examples

!!! example
    See [Regional plot](https://cloufield.github.io/gwaslab/visualization_regional/)


### FAQ

1. Why some genes are missing in the gene track?

We only included protein-coding genes in the reference GTF files for plotting the gene track.

2. Why some exons are missing in the gene track?

Sometimes the exon is too short to reach even 1 pixel in the plot. You can either increase the dpi or reduce the length of the region.

3. Why an error occurs even if both variants are in the reference VCF?

When the reference variant is mono-allelic in the reference VCF, LD can not be calculated.