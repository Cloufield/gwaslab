# Trumpet plot

Available since v3.4.20

## Description

Scatter plot for visualization of the relationship between Minor Allele Frequency (MAF) and Effect size (BETA/ log(OR))

- X axis: MAF
- Y axis: Effect size
- additional curves : power

## Usage

### "q" mode for quantitative traits

|Option|Type|Description|Default|
|-|-|-|-|
|`sig_level`|`float`|signiifcance level for power calculation| `5e-8`| 
|`n`|`string` or `int`|column name for sample size or an integer for sample size| required |
|`ts`|`list`|a list of power curve to draw| `[0.2,0.4,0.6,0.8]`|

### "b" mode for binary traits

|Option|Type|Description|Default|
|-|-|-|-|
|`sig_level`|`float`|signiifcance level for power calculation| `5e-8`| 
|`scase`|`int`|number of cases|required|
|`scontrol`|`int`|number of controls|required|
|`prevalence`|`float`|disease prevalence in general population|required|
|`or_to_rr`|`boolean`|if estimate RR using OR/beta and prevalence; for prevalence <10%, RR is similar to OR|"False"|
|`ts`|`list`|a list of power curve to draw| `[0.2,0.4,0.6,0.8]`|

### other option

|Option|Type|Description|Default|
|-|-|-||
|`build`|`string`|`19` or `38`, which build to use for annotation|"99"|
|`p_level`|`float`|upper limit of p values; variants with p values higher than this will be excluded from plotting||
|`anno`|`string`|which column in sumstats to use to annotate the variants||
|`anno_style`|`string`|`expand`, `tight` or `right`||
|`anno_x`|`float`|upper bound of abs(y) for annotation||
|`anno_y`|`float`|lower bound of abs(y) for annotation||
|`repel_force`|`float`|determine the interval between each annotation||
|`sort`|`string`|`beta` or `maf`, which to use to determine the order for annotation||
|`ylim`|`tuple`|ylim||
|`cmap`|`string`|matplotlib color map for power line|"Reds"|
|`xscale`|`string`| `log` or `nonlog`|`log`|
|`yscale_factor`|`float`|in case effect size needs to be converted by a factor|1|
|`n_matrix`|`int`|the higher the value is, the soomther the power line will be||
|`markercolor`|`string`|color||
|`ylabel`|`string`|Y axis label|`Effect size`|
|`xlabel`|`string`|X axis label|`Minor allele frequency`|


## example

!!! example "Trumpet plot for quantatitive trait"
    ```
    a = mysumstats.plot_trumpet(mode="q",
                                ts=[0.2,0.4,0.6,0.8] ,
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
                                save=True)
    ```

!!! example "Trumpet plot for binary trait"
    ```
    mysumstats.get_lead(sig_level=5e-6,gls=True).plot_trumpet(mode="b",
                                scase=36614,
                                scontrol=155150,
                                sig_level=5e-6,
                                p_level=5e-6,
                                ts=[0.2,0.4,0.6,0.8] ,
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
                                save=True)
    ```

## Citation

!!! quote Trumpet plot
    This function was adapted from Corte, L., Liou, L., O'Reilly, P. F., & García-González, J. (2023). Trumpet plots: Visualizing The Relationship Between Allele Frequency And Effect Size In Genetic Association Studies. medRxiv, 2023-04. 
    R Shiny app : https://juditgg.shinyapps.io/shinytrumpets/

!!! Power calculation
    For power calculation: 

    - Sham, P. C., & Purcell, S. M. (2014). Statistical power and significance testing in large-scale genetic studies. Nature Reviews Genetics, 15(5), 335-346.
    - Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). Joint analysis is more efficient than replication-based analysis for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
