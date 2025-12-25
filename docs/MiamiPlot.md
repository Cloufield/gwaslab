# Miami plot

!!! info "Implemented since v3.3.4"

As a standalone function, GWASLab can plot Miami plot given a pair of sumstats files.

## gl.plot_miami2()

```python
gl.plot_miami2( 
          path1=None,
          path2=None,
          merged_sumstats=None,
          ...
          )
```

!!! info "Updated since v3.4.38"

GWASLab creates Miami plot by iteratively calling `.plot_mqq()`. The function name is `plot_miami2()` (not `plot_miami()`).

## Options

| Option            | DataType          | Description                                                                 | Default      |
|-------------------|-------------------|-----------------------------------------------------------------------------|--------------|
| `path1`           | `gl.Sumstats`     | First sumstats object for top plot                                         | `None`       |
| `path2`           | `gl.Sumstats`     | Second sumstats object for bottom plot                                     | `None`       |
| `merged_sumstats` | `pd.DataFrame`    | Merged sumstats DataFrame with suffixes (alternative to path1/path2)       | `None`       |
| `suffixes`        | `list`            | Suffixes for merged_sumstats columns (e.g., ["_1", "_2"])                  | `None`       |
| `mode`            | `string`          | Plot mode: "m" (Manhattan only) or "mqq" (Manhattan + QQ)                  | `"m"`        |
| `scaled`          | `boolean`         | Use MLOG10P for both plots                                                  | `False`      |
| `scaled1`         | `boolean`         | Use MLOG10P for top plot                                                    | `False`      |
| `scaled2`         | `boolean`         | Use MLOG10P for bottom plot                                                 | `False`      |
| `id0`             | `string`          | ID column for merged_sumstats                                              | `"TCHR+POS"` |
| `id1`             | `string`          | ID column for path1                                                        | `None`       |
| `id2`             | `string`          | ID column for path2                                                        | `None`       |
| `same_ylim`       | `boolean`         | Use same y-axis limits for both plots                                       | `False`      |
| `fig_kwargs`      | `dict`            | Arguments passed to plt.subplots()                                         | `None`       |
| `save`            | `string` or `bool`| Save path or True for default path                                         | `None`       |
| `save_kwargs`     | `dict`            | Arguments passed to matplotlib savefig()                                  | `None`       |
| `fontsize`        | `int`             | Font size for labels                                                       | `10`         |
| `font_family`     | `string`          | Font family                                                                | `"Arial"`    |
| `dpi`             | `int`             | Dots per inch for figure                                                   | `100`        |
| `verbose`         | `boolean`         | Print progress messages                                                    | `True`       |

`gl.plot_miami2()` supports most options from `plot_mqq()`.

- Add a suffix "1" to the options and it will be passed to the top Manhattan plot.
- Add a suffix "2" to the options and it will be passed to the bottom Manhattan plot.
- If no suffix is provided, it will be passed to both plots.

## Example

!!! example "Basic Miami plot"
    ```python
    import gwaslab as gl
    
    # Load two sumstats
    gl1 = gl.Sumstats("sumstats1.tsv.gz", fmt="auto")
    gl2 = gl.Sumstats("sumstats2.tsv.gz", fmt="auto")
    
    # Create Miami plot
    gl.plot_miami2(path1=gl1, path2=gl2)
    ```

!!! example "Customized Miami plot with annotations"
    ```python
    gl.plot_miami2(
        path1=gl1,
        path2=gl2,
        skip=2,
        cut1=20,
        cut2=15,
        id1="SNPID",
        id2="SNPID",
        anno1=True,
        anno2="GENENAME",
        additional_line1=[1e-14],
        anno_set1=["rs3798519"],
        pinpoint1=[["rs3798519","rs35560038"],["rs7933262","rs8098510"]],
        pinpoint_color1=["purple","black"],
        highlight1=["rs8098510"],
        highlight2=[["rs8098510","rs3798519"], ["rs1491850"]],
        highlight_color2=["red","yellow"],
        jagged=True,
        verbose1=False,
        verbose2=False,
        save="miami_plot.png",
        save_kwargs={"dpi": 300}
    )
    ```

!!! example "Using merged sumstats"
    ```python
    # Merge two sumstats with suffixes
    merged = gl.merge_sumstats([gl1, gl2], suffixes=["_1", "_2"])
    
    # Create Miami plot from merged data
    gl.plot_miami2(
        merged_sumstats=merged.data,
        suffixes=["_1", "_2"],
        id0="TCHR+POS"
    )
    ```

See also: [Miami plot gallery](https://cloufield.github.io/gwaslab/visualization_miami2/)
    