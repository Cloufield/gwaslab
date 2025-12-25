# Forest plot

!!! info "Available since v3.4.40"

GWASLab provides a function to create forest plots for meta-analysis results, displaying individual study effect estimates with confidence intervals and combined meta-analysis estimates.

## gl.plot_forest()

```python
gl.plot_forest(
    data,
    study_col,
    ...
)
```

## Options

| Option                    | DataType          | Description                                                                 | Default        |
|---------------------------|-------------------|-----------------------------------------------------------------------------|----------------|
| `data`                    | `pd.DataFrame` or `str` | DataFrame containing study data, or path to whitespace-separated file | required       |
| `study_col`               | `string`          | Column name containing study identifiers                                   | required       |
| `group_col`               | `string`, `bool`, or `None` | Column name for grouping studies. If False/None, all in "Group1"      | `None`         |
| `beta_col`                | `string`          | Column name for effect estimates (beta)                                    | `"beta"`       |
| `se_col`                  | `string`          | Column name for standard errors                                            | `"se"`         |
| `compact_factor`          | `float`           | Factor to adjust figure height (higher = more compact)                    | `1.0`          |
| `width_ratios`            | `list`            | Width ratios for subplot columns [group, plot, text]                       | `[2, 6, 2]`    |
| `sharex`                  | `string`          | Whether to share x-axis across rows                                        | `"col"`        |
| `meta`                    | `boolean`         | Whether to perform and display meta-analysis statistics                    | `True`         |
| `combine_effects_kwargs`  | `dict`            | Additional arguments passed to combine_effects function                    | `None`         |
| `fig_kwargs`              | `dict`            | Additional arguments passed to plt.subplots()                              | `None`         |
| `save`                    | `string`, `bool`, or `None` | Save path, True for default, or None/False to not save              | `None`         |
| `save_kwargs`             | `dict`            | Additional arguments for matplotlib savefig()                              | `None`         |
| `fontsize`                | `int`             | Font size for labels and text                                              | `12`           |
| `font_family`             | `string`          | Font family for text                                                       | `"Arial"`      |
| `colors`                  | `list`            | List of colors for alternating studies. If None, uses default colors       | `None`         |
| `verbose`                 | `boolean`         | Whether to print progress messages                                         | `True`         |

## Examples

!!! example "Basic forest plot"
    ```python
    import pandas as pd
    import gwaslab as gl
    
    # Create example data
    data = pd.DataFrame({
        'study': ['Study1', 'Study2', 'Study3', 'Study4'],
        'beta': [0.5, 0.7, 0.6, 0.65],
        'se': [0.2, 0.15, 0.18, 0.16]
    })
    
    # Create forest plot
    fig, axes = gl.plot_forest(data, study_col='study')
    ```

!!! example "Forest plot with grouping"
    ```python
    # Add grouping column
    data['group'] = ['Cohort1', 'Cohort1', 'Cohort2', 'Cohort2']
    
    # Create forest plot with groups
    fig, axes = gl.plot_forest(
        data, 
        study_col='study',
        group_col='group',
        save="forest_plot.png",
        save_kwargs={"dpi": 300}
    )
    ```

!!! example "Customized forest plot"
    ```python
    fig, axes = gl.plot_forest(
        data,
        study_col='study',
        beta_col='BETA',  # Custom column name
        se_col='SE',      # Custom column name
        compact_factor=1.2,
        width_ratios=[1.5, 6, 2],
        fontsize=14,
        font_family="Times New Roman",
        colors=["#597FBD", "#74BAD3"],
        meta=True,
        save="forest_custom.png"
    )
    ```

!!! example "Load from file"
    ```python
    # Load data from file
    fig, axes = gl.plot_forest(
        data="meta_analysis_results.txt",
        study_col="STUDY",
        beta_col="BETA",
        se_col="SE"
    )
    ```
