# Comparing effect sizes

## Scatter plot : effect size comparison 

!!! info "Available since v3.4.17"

`gl.compare_effect()` will plot effect size comparison plot using two sets of sumstats. Alleles will be aligned to effect alleles in sumstats1.


## gl.compare_effect()

```python
gl.compare_effect (path1,
                   path2,
                   cols_name_list_1=None, effect_cols_list_1=None,
                   cols_name_list_2=None, effect_cols_list_2=None,
                   eaf=[],
                   maf_level=None,
                   label=None,
                   snplist=None,
                   mode="beta",
                   anno=False,
                   anno_het=False,
                   anno_min=0,
                   anno_min1=0,
                   anno_min2=0,
                   anno_diff=0,
                   scaled=False,
                   scaled1=False,
                   scaled2=False,
                   wc_correction=False, 
                   null_beta=0,
                   is_q=False,
                   include_all=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   drop=False,
                   wc_sig_level=5e-8,
                   # reg
                   reg_box=None,
                   is_reg=True,
                   fdr=False,
                   allele_match=False,
                   r_se=False,
                   is_45_helper_line=True,
                   legend_mode="full",
                   legend_title=r'$ P < 5 x 10^{-8}$ in:',
                   legend_title2=r'Heterogeneity test:',
                   legend_pos='upper left',
                   scatterargs=None,
                   plt_args=None,
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_args=None,
                   fontargs=None,
                   errargs=None,
                   legend_args=None,
                   sep=["\t","\t"],
                   log = Log(),
                   save=False,
                   save_args=None,
                   verbose=False)
    

```

## Options

### Path and column

- `path1` and `path2` : the paths to the sumstats. Can also be `gl.Sumstats` Object or `pd.DataFrame` (from v3.4.17). If `gl.Sumstats` Objects are provided, there is no need to set cols_name_list and effect_cols_list.
- `cols_name_list_1` and `cols_name_list_2` : list of column names for variants basic information
- `effect_cols_list_1` and `effect_cols_list_2` : list of column names for effect size-related columns
- `mode` : use beta or OR 
-  `cols_name_list_x` and `effect_cols_list_x` examples:
    - `[snpid,p,ea,nea]`        ,`[effect,se]`
    - `[snpid,p,ea,nea,chr,pos]`,`[effect,se]`
    - `[snpid,p,ea,nea,chr,pos]`,`[OR,OR_l,OR_h]`

| Option               | Type                                      | Description                                                               | Default |
|----------------------|-------------------------------------------|---------------------------------------------------------------------------|---------|
| `path1`              | `string`,`gl.Sumstats`, or `pd.DataFrame` | path to the sumstats file, or gwaslab Sumstats object or pandas Dataframe | None    |
| `path2`              | `string`,`gl.Sumstats`, or `pd.DataFrame` | path to the sumstats file, or gwaslab Sumstats object or pandas Dataframe | None    |
| `cols_name_list_1`   | `list`                                    | "[snpid,p,ea,nea]" or "[snpid,p,ea,nea,chr,pos]"                          | None    |
| `cols_name_list_2`   | `list`                                    | "[snpid,p,ea,nea]" or "[snpid,p,ea,nea,chr,pos]"                          | None    |
| `effect_cols_list_1` | `list`                                    | "[effect,se]" or "[OR,OR_95l,OR_95h]"                                     | None    |
| `effect_cols_list_2` | `list`                                    | "[effect,se]" or "[OR,OR_95l,OR_95h]"                                     | None    |
| `mode`               | `beta` or `OR`                            | plot beta or OR                                                           | `beta`  |

!!! note 
    from v3.4.17, you need to specify the parameters using keywords instead of using them as positional arguments for `cols_name_list_1` and `cols_name_list_2`, `effect_cols_list_1` and `effect_cols_list_1` .

### Save figures
| Option      | Type     | Description                 | Default                                            |
|-------------|----------|-----------------------------|----------------------------------------------------|
| `save`      | `string` | path to the saved file      | `./Sumstats1_Sumstats2_effect_comparison_plot.png` |
| `save_args` | `dict`   | parametrs for plt.savefig() | `{"dpi":300,"facecolor":"white"}`                  |

### SNP List

| Option    | Type   | Description                                                                                                                       | Default |
|-----------|--------|-----------------------------------------------------------------------------------------------------------------------------------|---------|
| `snplist` | `list` | optional, specify the variants you want to compare. If None, GWASLab will automatically extract lead variants from both sumstats. | None    |

### Filter by MAF

| Option      | Type    | Description                                                                                                                                                                 | Default |
|-------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `eaf`       | `list`  | optional, a list column names for effect allele frequency, in the order of [sumstats1_eaf, sumstats2_eaf]. It is required when you need to filter by maf using `maf_level`. | `None`  |
| `maf_level` | `float` | the maf filter for variants. Variants with maf < maf_level will be removed from comparison.                                                                                 | `None`  |

### Label

| Option           | Type     | Description                                                                                                              | Default                                     |
|------------------|----------|--------------------------------------------------------------------------------------------------------------------------|---------------------------------------------|
| `label`          | `list`   | a list of labels for the legend , in the order of ["Sumstats_1","Sumstats_2","Both","None"]                              | `["Sumstats_1","Sumstats_2","Both","None"]` |
| `sig_level`      | `float`  | the significance level for auto-extracting lead variants. If `snplist` is provided, the auto-extraction will be skipped. | `5e-8`                                      |
| `legend_title`   | `string` | legend title                                                                                                             | `'$ P < 5 x 10^{-8}$ in:'`                  |
| `legend_pos`     | `string` | legend position                                                                                                          | `upper left`                                |
| `xylabel_prefix` | `string` | -                                                                                                                        | `"Per-allele effect size in "`              |

### Annotation

| Option              | Type      | Description                                                         | Default |
|---------------------|-----------|---------------------------------------------------------------------|---------|
| `is_reg`            | `boolean` | if true, draw regression line.                                      | `True`  |
| `is_45_helper_line` | `boolean` | if true, draw 45 degree line.                                       | `True`  |
| `anno`              | `boolean` | if true, annotate the variants with ID.                             | `False` |
| `anno_diff`         | `float`   | threshold of effect size difference for annotation.                 | `0`     |
| `anno_min1`         | `float`   | threshold of sumstats1 minimum absolute effect size for annotation. | `0`     |
| `anno_min2`         | `float`   | threshold of sumstats2 minimum absolute effect size for annotation. | `0`     |

### Heterogeneity test

| Option     | Type      | Description                                                    | Default |
|------------|-----------|----------------------------------------------------------------|---------|
| `is_q`     | `boolean` | if true, apply the heterogeneity tests by Cochran's Q test.    | `True`  |
| `q_level`  | `float`   | the significance threshold for Cochran's Q test (raw p value). | `0.05`  |
| `anno_het` |           | annotate only variants with Phet < `q_level`                   | `False` |

### R SE

| Option | Type      | Description                                                                                    | Default |
|--------|-----------|------------------------------------------------------------------------------------------------|---------|
| `r_se` | `boolean` | If True, SE for r will be estimated using the jackknife method. (Note: available from v3.4.17) | `False` |

$$ s.e.(\hat{r}_{jack}) = \sqrt{ {{n-1}\over{n}} \sum_{i=1}^n(\hat{r_i} -\bar{r}_{jack} )^2 } $$

## Examples


!!! example "Use gl.Sumstats object"
    ```
    gl1 = gl.Sumstats("bbj_bmi_female.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")
    gl2 = gl.Sumstats("bbj_bmi_male.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")

    # auto extract lead SNPs and compare the effect sizes
    a = gl.compare_effect(path1= gl1,
                          path2= gl2
    )

    ```
    <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/258d3316-7d73-408a-8e20-e560b831861c">

!!! example "Use pandas DataFrame"
    ```
    pd1 = pd.read_table("bbj_bmi_female.txt.gz",sep="\t")
    pd2 = pd.read_table("bbj_bmi_male.txt.gz",sep="\t")

    # cols_name_list should be SNPID, P, Effect Allele, Non-Effect allele, Chromosome and Position
    # effect_cols_list should be BETA,SE

    a = gl.compare_effect(path1 = pd1,
                        cols_name_list_1 = ["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_1= ["BETA","SE"],
                        path2 = pd2,
                        cols_name_list_2 = ["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_2= ["BETA","SE"]
                        )
    ```

!!! example "Heterogeneity test"
    ```
    pd1 = pd.read_table("bbj_bmi_female.txt.gz",sep="\t")
    pd2 = pd.read_table("bbj_bmi_male.txt.gz",sep="\t")

    # cols_name_list should be SNPID, P, Effect Allele, Non-Effect allele, Chromosome and Position
    # effect_cols_list should be BETA,SE
    # is_q : if true, apply the heterogeneity tests by Cochran's Q test.
    # q_level` : `float`, the significance threshold for Cochran's Q test (raw p value).

    a = gl.compare_effect(path1 = pd1,
                          cols_name_list_1 = ["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_1= ["BETA","SE"],
                          path2 = pd2,
                          cols_name_list_2 = ["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_2= ["BETA","SE"],
                          is_q=True,
                          q_level=0.05,
                          legend_mode="full"
                          )
    ```
     <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/b9c4511c-c1eb-4c36-8b72-788d2b21b790">
     
!!! example "Annotation"
    ```
    a = gl.compare_effect(path1="bbj_bmi_female.txt.gz",
                      cols_name_list_1=["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_1=["BETA","SE"],
                      path2="bbj_bmi_male.txt.gz",
                      cols_name_list_2=["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_2=["BETA","SE"],
                      label=["Female","Male","Both","None"],
                      xylabel_prefix="Per-allele effect size for ",
                      r_se=True, 
                      is_q=True, 
                      anno=True, 
                      anno_het=True, # only annotate variants with significant heterogeneity
                      anno_diff=0.015,    # only annotate variants with a difference in effect size > 0.015
                      sig_level=5e-8,
                      legend_title=r'$ P < 5 x 10^{-8}$ in:',
                      verbose=True)
                      
    ```
    <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/9e2ac224-63d6-44da-9dbc-1cd43596055f">
    
!!! example "Male-specific and female-specific BMI Sumstats from JENGER"
    
    ```
    !wget -O bbj_bmi_male.txt.gz http://jenger.riken.jp/2analysisresult_qtl_download/
    !wget -O bbj_bmi_female.txt.gz http://jenger.riken.jp/4analysisresult_qtl_download/

    Headers of the files : # SNP    CHR	POS	REF	ALT	Frq	Rsq	BETA	SE	P
    ```
    
    ```python
    # GWASLab will automatically extract significant variants from both sumstats. 
    a = gl.compare_effect(path1="bbj_bmi_female.txt.gz",
                      cols_name_list_1=["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_1=["BETA","SE"],
                      path2="bbj_bmi_male.txt.gz",
                      cols_name_list_2=["SNP","P","REF","ALT","CHR","POS"],effect_cols_list_2=["BETA","SE"],
                      label=["Female","Male","Both","None"],
                      xylabel_prefix="Per-allele effect size for ",
                      anno=True,
                      anno_het=True,
                      anno_diff=0.015,
                      sig_level=5e-8,
                      legend_title=r'$ P < 5 x 10^{-8}$ in:',
                      verbose=True,
                      save = "myplot.png",
                      save_args= {"dpi":300,"facecolor":"white"}
    )
    ```
    <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/ed462a3d-64cf-4f06-bdb5-13294b4182c8">
    
    Reference: Akiyama, M., Okada, Y., Kanai, M., Takahashi, A., Momozawa, Y., Ikeda, M., ... & Kamatani, Y. (2017). Genome-wide association study identifies 112 new loci for body mass index in the Japanese population. Nature genetics, 49(10), 1458-1467.
