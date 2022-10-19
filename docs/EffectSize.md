# Scatter plot : effect size comparison

## gl.compare_effect()
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
                   verbose=False):
    

```

## Options
### path and column
- `path1` and `path2` : the paths to the sumstats.
- `cols_name_list_1` and `cols_name_list_2` : list of column names for variants basic information
- `effect_cols_list_1` and `effect_cols_list_1` : list of column names for effect size-related columns
- `mode` : use beta or OR 
-  examples:
    - [snpid,p,ea,nea]        ,[effect,se]
    - [snpid,p,ea,nea,chr,pos],[effect,se]
    - [snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]

### snplist
- `snplist` : optional, specify the variants you want to compare. If None, gwaslab will automatically extract lead variants from both sumstats.

### filter by maf: 
- `eaf` : optional, a list column names for effect allele frequency, in the order of [sumstats1_eaf, sumstats2_eaf]. It is needed when you need to filter by maf using `maf_level`.
- `maf_level`: the maf filter for variants. Vairants with maf < maf_level will be removed from comparison.

### label and annotation
- `label` : a list of labels for the legend , in the order of ["Sumstats_1","Sumstats_2","Both","None"].
- `anno` : if annotate the variants
- `sig_level`: the significance level for auto-extracting lead variants.
- `legend_title`: r'$ P < 5 x 10^{-8}$ in:',
- `legend_pos`: 'upper left'
- `xylabel_prefix`"Per-allele effect size in "
- `is_reg`
- `is_45_helper_line`

### heterogeneity test
- `is_q` : if apply the heterogeneity tests by Cochran's Q test.
- `q_level` : the significance threshold for Cochran's Q test.


## Example:
```
# Smoking behaviors : Smoking initiation (autosome, male)
!wget -O smoking_male.txt.gz http://jenger.riken.jp/16/

# Smoking behaviors : Smoking initiation (autosome, female)
!wget -O smoking_female.txt.gz http://jenger.riken.jp/18/
```

```python
# SNP	CHR	POS	A1	A2	A1Frq	Rsq	BETA	SE	P
# gwaslab will automatically extract significant variants from both sumstats. 
a = gl.compare_effect("smoking_female.txt.gz",
                      ["SNP","P","A1","A2","CHR","POS"],["BETA","SE"],
                      "smoking_male.txt.gz",
                      ["SNP","P","A1","A2","CHR","POS"],["BETA","SE"],
                      label=["Female","Male","Both","None"],
                      xylabel_prefix="Per-allele effect size for ",
                      sig_level=5e-6,
                      legend_title=r'$ P < 5 x 10^{-6}$ in:',
                      verbose=True)
```
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/196621266-4a7d2a83-68b0-4187-9d86-0e43258e268b.png">
