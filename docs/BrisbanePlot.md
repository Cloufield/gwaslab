# Brisbane plot

## Brisbane plot: GWAS signal density plot

GWASLab can create the [Brisbane plot](https://www.nature.com/articles/s41586-022-05275-y/figures/2) (GWAS signal density plot). Brisbane plot is a scatter plot that shows the signal density (number of variants within a specified window around each variant) for each variant, which is very useful for presenting the independent signals obtained from large-scale GWAS of complex traits. The signals are usually determined by other statistical methods such as conditional analysis.

**Key features:**

- **X-axis**: Genomic position (same as Manhattan plot)
- **Y-axis**: Signal density (number of neighboring variants within the window)
- **Purpose**: Visualize the density of independent signals across the genome

!!! note "Independent Signals"
    To create a meaningful Brisbane plot, you should load sumstats containing only independent signals (e.g., from conditional analysis). If you load the entire dataset, the plot will simply reflect the marker density for your sumstats. To investigate independent signals, please use other tools such as GCTA-COJO. GWASLab only calculates the density of all variants in the gl.Sumstats Object.

## .plot_mqq(mode="b")

```
mysumstats.plot_mqq(mode="b")
```

Creates a Brisbane plot showing signal density across the genome.

### Options

| Option              | DataType | Description                                                                                                 | Default |
|---------------------|----------|-------------------------------------------------------------------------------------------------------------|---------|
| `mode`              | `str`    | Plotting mode. Use `"b"` for Brisbane plot only        | `"b"`   |
| `bwindowsizekb`     | `int`    | Window size in kb (flanking region length on one side). Total window = 2 × bwindowsizekb                    | `100`   |
| `density_color`     | `bool`   | If True, color variants by density value. When True, `density_palette` overrides `colors`                  | `False` |
| `density_palette`    | `str`    | Color palette for density coloring (e.g., 'Reds', 'Blues', 'viridis')                                      | `'Reds'` |
| `density_range`     | `tuple`  | Color range for density plot. If None, auto-selected as (min(DENSITY), max(DENSITY))                       | `None`  |
| `density_threshold` | `int`    | Threshold for density coloring. Above threshold uses `density_palette`, below uses `density_tpalette`       | `5`     |
| `density_trange`    | `tuple`  | Threshold range for density plot (for variants below threshold)                                            | `(0,10)`|
| `density_tpalette`  | `str`    | Color palette for density below threshold                                                                   | `'Blues'`|
| `anno`              | `str`    | Annotation style. Options: `None`, `True` (chr:pos), `'GENENAME'`, or column name                          | `None`  |
| `anno_set`          | `list`   | List of variant IDs to annotate                                                                            | `None`  |
| `windowsizekb`      | `int`    | Window size in kb for determining lead variants for annotation                                              | `500`   |

### Example

!!! example "Basic Brisbane plot"
    ```
    mysumstats.plot_mqq(mode="b")
    ```

!!! example "Brisbane plot with custom window size"
    ```
    mysumstats.plot_mqq(mode="b", bwindowsizekb=250)
    ```

!!! example "Brisbane plot with density coloring"
    ```
    mysumstats.plot_mqq(
        mode="b",
        density_color=True,
        density_palette="Reds"
    )
    ```

!!! example "Manhattan-Brisbane layout"
    ```
    mysumstats.plot_mqq(mode="mb")
    ```

!!! example "Brisbane plot with annotations"
    ```
    mysumstats.plot_mqq(
        mode="b",
        anno="GENENAME",
        windowsizekb=500
    )
    ```

See more examples [here](https://cloufield.github.io/gwaslab/visualization_brisbane/)


## Calculate signal density

You can use `.get_density()` to calculate the density without creating a plot. This adds a `DENSITY` column to your sumstats.

```
mysumstats.get_density(windowsizekb=100)
```

### Parameters

| Option         | DataType | Description                                                                                                 | Default |
|----------------|----------|-------------------------------------------------------------------------------------------------------------|---------|
| `windowsizekb`  | `int`    | Window size in kb for calculation of signal density. Total window = 2 × windowsizekb (flanking on each side) | `100`   |
| `sig_list`     | `DataFrame`, optional | If provided, density is calculated based on significant variants in this list (for conditional analysis) | `None`  |

### How it works

The function calculates density by counting the number of variants within a window around each variant:
- For each variant at position P, it counts variants in the range [P - windowsizekb, P + windowsizekb]
- The density value represents how many neighboring variants are within this window
- Higher density values indicate regions with more clustered signals

### Output

The function adds a `DENSITY` column to `mysumstats.data` and prints summary statistics:
- Mean density
- Median density
- Standard deviation
- Maximum density and the variant with maximum density

!!! example "Calculate signal density"
    ```
    mysumstats.get_density(windowsizekb=100)
    
    # Output:
    # -Calculating DENSITY with windowsize of 100 kb
    # -Mean : 2.345 signals per 100 kb
    # -SD : 1.234
    # -Median : 2.0 signals per 100 kb
    # -Max : 15 signals per 100 kb at variant rs11555194
    
    mysumstats.data
    	SNPID	CHR	POS	P	STATUS	DENSITY
    0	rs2710888	1	959842	2.190000e-57	9999999	1
    1	rs3934834	1	1005806	2.440000e-29	9999999	1
    2	rs182532	1	1287040	1.250000e-18	9999999	1
    3	rs17160669	1	1305561	1.480000e-28	9999999	1
    4	rs9660106	1	1797947	1.860000e-12	9999999	0
    ...	...	...	...	...	...	...
    12106	rs9628283	22	50540766	5.130000e-15	9999999	1
    12107	rs28642259	22	50785718	1.140000e-13	9999999	1
    12108	rs11555194	22	50876662	2.000000e-15	9999999	2
    12109	rs762669	22	50943423	3.000000e-30	9999999	1
    12110	rs9628185	22	51109992	5.430000e-12	9999999	0
    ```

!!! example "Calculate density from significant variants list"
    ```
    # For conditional analysis: calculate density based on pre-identified significant variants
    significant_variants = mysumstats.data[mysumstats.data["P"] < 5e-8]
    mysumstats.get_density(windowsizekb=100, sig_list=significant_variants)
    ```

## Extract top variants by density

You can use `.get_top()` to extract variants with the highest density values within sliding windows. This is useful for identifying the most dense signal regions in your Brisbane plot.

```
mysumstats.get_top(by="DENSITY", windowsizekb=500)
```

### Parameters

| Option         | DataType | Description                                                                                                 | Default |
|----------------|----------|-------------------------------------------------------------------------------------------------------------|---------|
| `by`           | `str`    | Column name whose values are maximized to choose top variants. Default is `"DENSITY"`                      | `"DENSITY"` |
| `threshold`    | `float`, optional | If provided, only variants with `by` >= `threshold` are considered. If None, uses median of maximum values per chromosome | `None` |
| `windowsizekb` | `int`    | Sliding window size in kilobases used to determine locus boundaries                                        | `500`   |
| `bwindowsizekb` | `int`  | Window size in kb for calculating density (only used if DENSITY column doesn't exist and `by="DENSITY"`)    | `100`   |
| `anno`         | `bool`   | If True, annotate output with nearest gene names                                                           | `False` |
| `gls`          | `bool`   | If True, return a new Sumstats object instead of DataFrame                                                  | `False` |

### How it works

The function identifies top variants by:
1. **Calculating density** (if `by="DENSITY"` and DENSITY column doesn't exist)
2. **Filtering by threshold** (if provided, or using median of chromosome maxima)
3. **Sliding window approach**: Within each window on a chromosome, selects the variant with the highest value of the specified metric
4. **Non-overlapping windows**: Ensures selected variants are separated by at least `windowsizekb`

This is similar to `.get_lead()` but doesn't rely on P-values - it can use any metric column (DENSITY, MLOG10P, BETA, etc.).

### Return value

- **If `gls=False`**: Returns a pandas DataFrame containing the selected top variants
- **If `gls=True`**: Returns a new Sumstats object containing the selected top variants

!!! example "Extract top density variants"
    ```
    # Extract variants with highest density in 500kb windows
    top_density_variants = mysumstats.get_top(by="DENSITY", windowsizekb=500)
    
    # Output:
    # -Sliding window size for extracting top: 500 kb
    # -Using DENSITY threshold: 3.5
    # -Identified 25 top variants!
    ```

!!! example "Extract top density variants with custom threshold"
    ```
    # Only consider variants with density >= 5
    top_variants = mysumstats.get_top(
        by="DENSITY",
        threshold=5,
        windowsizekb=500
    )
    ```

!!! example "Extract top variants by other metrics"
    ```
    # Extract top variants by MLOG10P instead of DENSITY
    top_p_variants = mysumstats.get_top(by="MLOG10P", windowsizekb=500)
    
    # Extract top variants by BETA
    top_beta_variants = mysumstats.get_top(by="BETA", windowsizekb=500)
    ```

!!! example "Return as Sumstats object"
    ```
    # Get top density variants as a new Sumstats object
    top_sumstats = mysumstats.get_top(by="DENSITY", gls=True)
    top_sumstats.plot_mqq(mode="b")
    ```

!!! example "Extract top variants with gene annotation"
    ```
    # Extract top density variants and annotate with gene names
    top_variants = mysumstats.get_top(
        by="DENSITY",
        windowsizekb=500,
        anno=True
    )
    ```

## Reference
!!! quote "Citation for Brisbane plot"
    Yengo, L., Vedantam, S., Marouli, E., Sidorenko, J., Bartell, E., Sakaue, S., ... & Lee, J. Y. (2022). A saturated map of common genetic variants associated with human height. Nature, 1-16.