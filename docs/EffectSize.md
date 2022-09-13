# Scatter plot : effect size comparison

```python
gl.compare_effect( path1,
                   cols_name_list_1, effect_cols_list_1,
                   path2,
                   cols_name_list_2, effect_cols_list_2,
                   eaf=[],
                   maf_level=None,
                   label=["Sumstats_1","Sumstats_2","Both","None"],
                   snplist=None,
                   mode="beta",
                   anno=False,
                   null_beta=0,
                   is_q=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   legend_title=r'$ P < 5 x 10^{-8}$ in:',
                   legend_pos='upper left',
                   reg_box=dict(boxstyle='round', facecolor='white', alpha=1,edgecolor="grey"),
                   is_reg=True,
                   is_45_helper_line=True,
                   scatterargs={"s":20},
                   plt_args={"figsize":(8,8),"dpi":300},
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_args={"color":'black', "linestyle":'-',"lw":1},
                   fontargs={'family':'sans','fontname':'Arial','fontsize':12},
                   errargs={"ecolor":"#cccccc","elinewidth":1},
                   sep=["\t","\t"],
                   log = Log(),
                   verbose=False)
```

`path1` and `path2` : the paths to the sumstats.
`cols_name_list_1` and `cols_name_list_2` : list of column names for variants basic information, in the order of [snpid,p,ea,nea,chr,pos]
`effect_cols_list_1` and `effect_cols_list_1` : list of column names for effect size-related columns, in the order of [effect,se] or [OR,OR_95L,OR_95H]
`eaf` : optional, a list column names for effect allele frequency, in the order of [sumstats1_eaf, sumstats2_eaf]. It is needed when you need to filter by maf using `maf_level`.
`maf_level`: the maf filter for variants. Vairants with maf < maf_level will be removed from comparison.
`label` : a list of labels for the legend , in the order of ["Sumstats_1","Sumstats_2","Both","None"].
`snplist` : optional, specify the variants you want to compare. If None, gwaslab will automatically extract lead variants from both sumstats.
`anno` : if annotate the variants
`is_q` : if apply the heterogeneity tests by Cochran's Q test.
`q_level` : the significance threshold for Cochran's Q test.
`sig_level`: the significance level for auto-extracting lead variants.

Example:

```python
gl.compare()
```
