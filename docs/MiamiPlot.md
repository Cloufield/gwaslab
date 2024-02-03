# Miami plot

!!! info "Implemented since v3.3.4"

As a standalone function, GWASLab can plot Miami plot given a pair of sumstats files.

## gl.plot_miami()

```
gl.plot_miami( 
          path1,
          path2,
          ...
          )
```

!!! info "Updated since v3.4.38"

GWASLab creates Miami plot by iteratively calling `.plot_mqq()`.

## Options

| Option  | DataType      | Description        | Default |
|---------|---------------|--------------------|---------|
| `path1` | `gl.Sumstats` | gl.Sumstats Object | -       |
| `path2` | `gl.Sumstats` | gl.Sumstats Object | -       |

`gl.plot_miami()` now  supports most options in `plot_mqq()`.

Adding a suffix "1" to the options and it will be passed to the top Manhattan plot.
Adding a suffix "2" to the options and it will be passed to the bottom Manhattan plot.
If no suffix is provided, it will be passed to both plots.

## Example

!!! example

    ```
    gl.plot_miami2(path1= gl1,
                    path2= gl2,
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
                    verbose2=False
)
    ```

    See [Miami plot](https://cloufield.github.io/gwaslab/visualization_miami2/)
    