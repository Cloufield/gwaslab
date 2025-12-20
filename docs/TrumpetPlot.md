# Trumpet plot

!!! info "Available since v3.4.20"

## Description

Scatter plot for visualization of the relationship between Minor Allele Frequency (MAF) and Effect size (BETA/ log(OR))

- X axis: MAF
- Y axis: Effect size
- additional curves: power

## .plot_trumpet()

```
.plot_trumpet()
```

### "q" mode for quantitative traits

| Option      | Type              | Description                                               | Default             |
|-------------|-------------------|-----------------------------------------------------------|---------------------|
| `sig_level` | `float`           | signiifcance level for power calculation                  | `5e-8`              |
| `n`         | `string` or `int` | column name for sample size or an integer for sample size | required            |
| `ts`        | `list`            | a list of power curve to draw                             | `[0.2,0.4,0.6,0.8]` |

### "b" mode for binary traits

| Option       | Type      | Description                                                                           | Default             |
|--------------|-----------|---------------------------------------------------------------------------------------|---------------------|
| `sig_level`  | `float`   | signiifcance level for power calculation                                              | `5e-8`              |
| `ncase`      | `int`     | number of cases                                                                       | required            |
| `ncontrol`   | `int`     | number of controls                                                                    | required            |
| `prevalence` | `float`   | disease prevalence in general population                                              | required            |
| `or_to_rr`   | `boolean` | if estimate RR using OR/beta and prevalence; for prevalence <10%, RR is similar to OR | "False"             |
| `ts`         | `list`    | a list of power curve to draw                                                         | `[0.2,0.4,0.6,0.8]` |

### Other options

| Option                | Type              | Description                                                                                     | Default                                |
|-----------------------|-------------------|-------------------------------------------------------------------------------------------------|----------------------------------------|
| `build`               | `string`          | `19` or `38`, which build to use for annotation                                                 | `"99"`                                 |
| `p_level`             | `float`           | upper limit of p values; variants with p values higher than this will be excluded from plotting | `5e-8`                                 |
| `anno`                | `string`          | which column in sumstats to use to annotate the variants                                        | `None`                                 |
| `anno_set`            | `list`            | list of SNPIDs/rsIDs to annotate (if None, auto-selects based on thresholds)                   | `None`                                 |
| `anno_alias`          | `dict`            | dictionary mapping SNPID/rsID to custom annotation text                                        | `None`                                 |
| `anno_style`          | `string`          | `expand`, `tight` or `right`                                                                    | `"expand"`                              |
| `anno_kwargs`         | `dict`            | arguments for annotation                                                                        | `None`                                 |
| `anno_source`         | `string`          | annotation source: `"ensembl"` or `"refseq"`                                                    | `"ensembl"`                             |
| `anno_x`              | `float`           | upper bound of abs(y) for annotation                                                            | `0.01`                                 |
| `anno_y`              | `float`           | lower bound of abs(y) for annotation                                                            | `1`                                    |
| `anno_d`              | `dict`            | dictionary to adjust annotation arm directions                                                  | `None`                                 |
| `arm_scale`           | `float`           | factor to adjust annotation arm height                                                          | `1`                                    |
| `repel_force`         | `float`           | determine the interval between each annotation                                                  | `0.01`                                 |
| `sort`                | `string`          | `beta` or `eaf`, which to use to determine the order for annotation                             | `"beta"`                               |
| `ylim`                | `tuple`           | y-axis limits                                                                                    | `None` (auto)                           |
| `xlim`                | `tuple`           | x-axis limits                                                                                    | `None` (auto)                           |
| `cmap`                | `string`          | matplotlib color map for power line                                                             | `"cool"`                               |
| `xscale`              | `string`          | `log` or `nonlog`                                                                               | `"log"`                                |
| `yscale_factor`       | `float`           | effect size will be multiplied by a factor                                                      | `1`                                    |
| `n_matrix`            | `int`             | The higher the value is, the smoother the power line will be                                    | `1000`                                 |
| `markercolor`         | `string`          | color for markers                                                                                | `"#597FBD"`                            |
| `hue`                 | `string`          | column name to use for color coding                                                              | `None`                                 |
| `size`                | `string` or `float`| column name or fixed size for markers                                                            | `None`                                 |
| `sizes`                | `tuple`           | (min, max) size range for markers                                                                | `None`                                 |
| `highlight`           | `list`            | list of SNPIDs/rsIDs to highlight                                                               | `None`                                 |
| `highlight_chrpos`    | `boolean`         | if True, highlight is interpreted as (CHR, POS) tuples                                           | `False`                                |
| `highlight_color`     | `string`          | color for highlighted variants                                                                   | `"#CB132D"`                            |
| `highlight_windowkb`  | `int`             | window size in kb for highlighting regions                                                       | `500`                                  |
| `pinpoint`            | `list`            | list of SNPIDs/rsIDs to pinpoint                                                                 | `None`                                 |
| `pinpoint_color`      | `string`          | color for pinpointed variants                                                                    | `"red"`                                |
| `ylabel`              | `string`          | Y axis label                                                                                     | `"Effect size"`                        |
| `xlabel`              | `string`          | X axis label                                                                                     | `"Minor allele frequency"`             |
| `xticks`              | `list`            | X axis tick positions                                                                            | `None`                                 |
| `xticklabels`         | `list`            | X axis tick labels                                                                               | `None`                                 |
| `yticks`              | `list`            | Y axis tick positions                                                                            | `None`                                 |
| `yticklabels`         | `list`            | Y axis tick labels                                                                               | `None`                                 |
| `title`               | `string`          | plot title                                                                                       | `None`                                 |
| `title_fontsize`      | `int`             | font size for title                                                                              | `15`                                   |
| `fontsize`            | `int`             | font size for labels                                                                             | `15`                                   |
| `font_family`         | `string`          | font family                                                                                      | `"Arial"`                              |
| `scatter_kwargs`      | `dict`            | additional arguments for scatter plot                                                            | `None`                                 |
| `fig_kwargs`          | `dict`            | additional arguments for figure creation                                                         | `None`                                 |
| `save`                | `string` or `bool`| if True, save to default path; if string, save to specified path                                | `False`                                |
| `save_kwargs`         | `dict`            | arguments for matplotlib savefig()                                                               | `None`                                 |
| `verbose`             | `boolean`         | print progress messages                                                                          | `True`                                 |

## example

!!! example "Trumpet plot for quantitative traits"
    ```python
    mysumstats.plot_trumpet(
        mode="q",
        ts=[0.2,0.4,0.6,0.8],
        anno="Gene",
        anno_style="expand",
        cmap="cool",
        sig_level=5e-8,
        build="19",
        anno_x=0.01,
        anno_y=1,
        p_level=5e-8,
        n_matrix=2000,
        fontsize=12,
        xscale="log",
        repel_force=0.15,
        yscale_factor=5.63,
        sort="eaf",
        ylim=(-5,4),
        save=True
    )
    ```
    ![image](https://github.com/Cloufield/gwaslab/assets/40289485/0b000467-4318-4045-b103-36b59aa3cd3d)


!!! example "Trumpet plot for binary traits"
    ```python
    mysumstats.get_lead(sig_level=5e-6, gls=True).plot_trumpet(
        mode="b",
        ncase=36614,
        ncontrol=155150,
        sig_level=5e-6,
        p_level=5e-6,
        ts=[0.2,0.4,0.6,0.8],
        anno="GENENAME",
        anno_style="right",
        cmap="cool",
        or_to_rr=True,
        build="19",
        anno_x=0.01,
        anno_y=0,
        n_matrix=2000,
        fontsize=12,
        xscale="log",
        repel_force=0.15,
        sort="eaf",
        ylim=(-5,4),
        save=True
    )
    ```
    ![image](https://github.com/Cloufield/gwaslab/assets/40289485/308f9e5b-386c-4f48-8ca4-dabe557ab472)

## Citation

!!! quote Trumpet plot
    This function was adapted from Corte, L., Liou, L., O'Reilly, P. F., & García-González, J. (2023). Trumpet plots: Visualizing The Relationship Between Allele Frequency And Effect Size In Genetic Association Studies. medRxiv, 2023-04. 
    R Shiny app : https://juditgg.shinyapps.io/shinytrumpets/

!!! Power calculation
    For power calculation: 

    - Sham, P. C., & Purcell, S. M. (2014). Statistical power and significance testing in large-scale genetic studies. Nature Reviews Genetics, 15(5), 335-346.
    - Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). Joint analysis is more efficient than replication-based analysis for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
