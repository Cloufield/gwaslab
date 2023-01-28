# Comparing effect sizes

## Scatter plot : effect size comparison
`gl.compare_effect()` will plot effect size comparison plot using two sets of sumstats. Alleles will be aligned to effect alleles in sumstats1.

```python
gl.compare_effect (path1,
                   cols_name_list_1, effect_cols_list_1,
                   path2,
                   cols_name_list_2, effect_cols_list_2,
                   eaf=[],
                   maf_level=None,
                   label=["Sumstats_1","Sumstats_2","Both","None"],
                   snplist=None,
                   mode="beta",
                   anno=False,
                   wc_correction=False, 
                   null_beta=0,
                   is_q=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   wc_sig_level=5e-8,
                   reg_box=dict(boxstyle='round', facecolor='white', alpha=1,edgecolor="grey"),
                   is_reg=True,
                   r2_se=False,
                   is_45_helper_line=True,
                   legend_title=r'$ P < 5 x 10^{-8}$ in:',
                   legend_pos='upper left',
                   scatterargs={"s":20},
                   plt_args={"figsize":(8,8),"dpi":300},
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_args={"color":'black', "linestyle":'-',"lw":1},
                   fontargs={'fontsize':12},
                   # 'family':'sans','fontname':'Arial',
                   errargs={"ecolor":"#cccccc","elinewidth":1},
                   sep=["\t","\t"],
                   log = Log(),
                   verbose=True)
    

```

## Options
### Path and column
- `path1` and `path2` : the paths to the sumstats.
- `cols_name_list_1` and `cols_name_list_2` : list of column names for variants basic information
- `effect_cols_list_1` and `effect_cols_list_1` : list of column names for effect size-related columns
- `mode` : use beta or OR 
-  examples:
    - [snpid,p,ea,nea]        ,[effect,se]
    - [snpid,p,ea,nea,chr,pos],[effect,se]
    - [snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]

### Snplist
- `snplist` : optional, specify the variants you want to compare. If None, gwaslab will automatically extract lead variants from both sumstats.

### Filter by maf: 
- `eaf` : optional, a list column names for effect allele frequency, in the order of [sumstats1_eaf, sumstats2_eaf]. It is needed when you need to filter by maf using `maf_level`.
- `maf_level`: the maf filter for variants. Vairants with maf < maf_level will be removed from comparison.

### Label and annotation
- `label` : a list of labels for the legend , in the order of ["Sumstats_1","Sumstats_2","Both","None"].
- `anno` : if annotate the variants
- `sig_level`: the significance level for auto-extracting lead variants.
- `legend_title`: r'$ P < 5 x 10^{-8}$ in:',
- `legend_pos`: legend position, default: 'upper left'
- `xylabel_prefix` : "Per-allele effect size in "
- `is_reg` : draw regression line or not
- `is_45_helper_line`: draw 45 degree line or not

### Heterogeneity test
- `is_q` : if apply the heterogeneity tests by Cochran's Q test.
- `q_level` : the significance threshold for Cochran's Q test.

### Winner's Curse correction

- `wc_correction`: `boolean` If True, wc correction will be performed before plotting. 

!!! info "Winner's Curse correction" 
    Reference: Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and confidence intervals for odds ratios in genome-wide association studies. Biostatistics, 9(4), 621-634.

### R2 SE

- `r2_se`: `boolean` If True, SE for r2 will be estimated using the jackknife method.

$$ s.e.(\hat{r^2}_{jack}) = \sqrt{ {{n-1}\over{n}} \sum_{i=1}^n(\hat{r^2_i} -\bar{r^2}_{jack} )^2 } $$

## Example:

!!! example "Male-specific and female-specific BMI Sumstats from JENGER"
    
    ```
    !wget -O bbj_bmi_male.txt.gz http://jenger.riken.jp/2analysisresult_qtl_download/
    !wget -O bbj_bmi_female.txt.gz http://jenger.riken.jp/4analysisresult_qtl_download/

    Headers of the files : # SNP	CHR	POS	REF	ALT	Frq	Rsq	BETA	SE	P
    ```
    
    ```python
    
    # gwaslab will automatically extract significant variants from both sumstats. 
    a = gl.compare_effect("bbj_bmi_female.txt.gz",
                          ["SNP","P","REF","ALT","CHR","POS"],["BETA","SE"],
                          "bbj_bmi_male.txt.gz",
                          ["SNP","P","REF","ALT","CHR","POS"],["BETA","SE"],
                          label=["Female","Male","Both","None"],
                          xylabel_prefix="Per-allele effect size for ",
                          sig_level=5e-6,
                          legend_title=r'$ P < 5 x 10^{-6}$ in:',
                          verbose=True)
    ```
    <img width="500" alt="image" src="https://user-images.githubusercontent.com/40289485/215021843-4572636d-b1a8-43f5-8f6e-070b09e5270f.png">
    
    ```python
    # wc_correction=True : perform winner's curse correction
    # r2_se=True : estimate the se for R2 using jackknife method
    
    a = gl.compare_effect("bbj_bmi_female.txt.gz",
                          ["SNP","P","REF","ALT","CHR","POS"],["BETA","SE"],
                          "bbj_bmi_male.txt.gz",
                          ["SNP","P","REF","ALT","CHR","POS"],["BETA","SE"],
                          label=["Female","Male","Both","None"],
                          xylabel_prefix="Per-allele effect size for ",
                          wc_correction=True,
                          r2_se=True,
                          sig_level=5e-6,
                          legend_title=r'$ P < 5 x 10^{-6}$ in:',
                          verbose=True
    )
    ```
    <img width="500" alt="image" src="https://user-images.githubusercontent.com/40289485/215021886-6ee4beb4-bf9f-42d1-a93b-2dca72841022.png">
    
    Reference: Akiyama, M., Okada, Y., Kanai, M., Takahashi, A., Momozawa, Y., Ikeda, M., ... & Kamatani, Y. (2017). Genome-wide association study identifies 112 new loci for body mass index in the Japanese population. Nature genetics, 49(10), 1458-1467.
