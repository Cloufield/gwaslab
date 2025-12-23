# Comparing effect sizes

## Scatter plot : effect size comparison 

`gl.compare_effect()` will plot effect size comparison plot using two sets of sumstats. Alleles will be aligned to effect alleles in sumstats1.

!!! note "Important: Requires Sumstats Objects"
    From the latest version, `gl.compare_effect()` requires `gl.Sumstats` objects for both `path1` and `path2`. File paths or pandas DataFrames are no longer directly supported. Please load your sumstats using `gl.Sumstats()` first.


## gl.compare_effect()

```
gl.compare_effect(path1,
                   path2,
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
                   anno_kwargs=None,
                   wc_correction=False, 
                   null_beta=0,
                   is_q=False,
                   is_q_mc=False,
                   include_all=True,
                   q_level=0.05,
                   sig_level=5e-8,
                   get_lead_kwargs=None,
                   drop=False,
                   wc_sig_level=5e-8,
                   # reg
                   reg_box=None,
                   is_reg=True,
                   fdr=False,
                   reg_text="full",
                   allele_match=False,
                   r_se=False,
                   is_45_helper_line=True,
                   legend_mode="full",
                   legend_title=r'$\mathregular{ P < 5 x 10^{-8}}$ in:',
                   legend_title2=r'Heterogeneity test:',
                   legend_pos='upper left',
                   scatter_kwargs=None,
                   fig_kwargs=None,
                   xylabel_prefix="Per-allele effect size in ",
                   helper_line_kwargs=None,
                   adjust_text_kwargs=None,
                   adjust_text_kwargs_l=None,
                   adjust_text_kwargs_r=None,
                   font_kwargs=None,
                   build="19",
                   r_or_r2="r",
                   err_kwargs=None,
                   legend_kwargs=None,
                   clean_output=False,
                   log=Log(),
                   save=False,
                   save_kwargs=None,
                   verbose=False)
    

```

## Options

### Path and column

- `path1` and `path2` : **Required** `gl.Sumstats` objects. The function will automatically extract the necessary columns (SNPID, P/MLOG10P, EA, NEA, CHR, POS, BETA/OR, SE, EAF if available) from the Sumstats objects.
- `mode` : use beta or OR 

| Option               | Type                                      | Description                                                               | Default |
|----------------------|-------------------------------------------|---------------------------------------------------------------------------|---------|
| `path1`              | `gl.Sumstats`                             | GWASLab Sumstats object for the first dataset                             | Required |
| `path2`              | `gl.Sumstats`                             | GWASLab Sumstats object for the second dataset                            | Required |
| `mode`               | `beta` or `OR`                            | plot beta or OR                                                           | `beta`  |

!!! note 
    The function automatically detects whether sumstats use P or MLOG10P, and whether EAF columns are available. No need to specify column names manually.

### Save figures
| Option      | Type     | Description                 | Default                                            |
|-------------|----------|-----------------------------|----------------------------------------------------|
| `save`      | `string` or `bool` | path to the saved file or True to use default name | `False` |
| `save_kwargs` | `dict`   | parameters for plt.savefig() | `{"dpi":300,"facecolor":"white"}`                  |

### SNP List

| Option    | Type   | Description                                                                                                                       | Default |
|-----------|--------|-----------------------------------------------------------------------------------------------------------------------------------|---------|
| `snplist` | `list` | optional, specify the variants you want to compare. If None, GWASLab will automatically extract lead variants from both sumstats. | None    |

### Filter by MAF

| Option      | Type    | Description                                                                                                                                                                 | Default |
|-------------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `maf_level` | `float` | the maf filter for variants. Variants with maf < maf_level will be removed from comparison. Requires EAF columns in both sumstats.                                         | `None`  |

### Label and styling

| Option           | Type     | Description                                                                                                              | Default                                     |
|------------------|----------|--------------------------------------------------------------------------------------------------------------------------|---------------------------------------------|
| `label`          | `list`   | a list of labels for the legend , in the order of ["Sumstats_1","Sumstats_2","Both","None"]                              | `["Sumstats_1","Sumstats_2","Both","None"]` |
| `sig_level`      | `float`  | the significance level for auto-extracting lead variants. If `snplist` is provided, the auto-extraction will be skipped. | `5e-8`                                      |
| `legend_title`   | `string` | legend title (automatically adjusts if `sig_level` differs from 5e-8)                                                    | `r'$\mathregular{ P < 5 x 10^{-8}}$ in:'`   |
| `legend_title2`  | `string` | second legend title for heterogeneity test                                                                             | `r'Heterogeneity test:'`                    |
| `legend_pos`     | `string` | legend position                                                                                                          | `'upper left'`                              |
| `legend_mode`    | `string` | legend display mode: `"full"` for detailed legend with heterogeneity info                                                | `"full"`                                    |
| `legend_kwargs`  | `dict`   | additional parameters for legend                                                                                        | `None`                                      |
| `xylabel_prefix` | `string` | prefix for x and y axis labels                                                                                          | `"Per-allele effect size in "`              |
| `scatter_kwargs` | `dict`   | parameters for scatter plot (passed to plt.scatter)                                                                     | `None`                                      |
| `fig_kwargs`     | `dict`   | parameters for figure creation (passed to plt.subplots)                                                                 | `None`                                      |
| `font_kwargs`    | `dict`   | parameters for font styling                                                                                             | `None`                                      |
| `helper_line_kwargs` | `dict` | parameters for helper lines (45-degree line, x=0, y=0)                                                                  | `None`                                      |
| `err_kwargs`     | `dict`   | parameters for error bars                                                                                                | `None`                                      |

### Annotation

| Option              | Type      | Description                                                         | Default |
|---------------------|-----------|---------------------------------------------------------------------|---------|
| `is_reg`            | `boolean` | if true, draw regression line.                                      | `True`  |
| `is_45_helper_line` | `boolean` | if true, draw 45 degree line.                                       | `True`  |
| `anno`              | `boolean`, `str`, or `dict` | if `True`, annotate with SNPID; if `"GENENAME"`, annotate with gene names; if `dict`, annotate with custom labels (key=SNPID, value=label) | `False` |
| `anno_kwargs`       | `dict`    | parameters for text annotation (passed to plt.text)                | `None`  |
| `anno_diff`         | `float`   | threshold of effect size difference for annotation.                 | `0`     |
| `anno_min`          | `float`   | threshold of minimum absolute effect size for annotation (applies to both). | `0`     |
| `anno_min1`         | `float`   | threshold of sumstats1 minimum absolute effect size for annotation. | `0`     |
| `anno_min2`         | `float`   | threshold of sumstats2 minimum absolute effect size for annotation. | `0`     |
| `adjust_text_kwargs` | `dict`    | parameters for adjustText (applies to both left and right annotations) | `None`  |
| `adjust_text_kwargs_l` | `dict`  | parameters for adjustText for left-side annotations                 | `None`  |
| `adjust_text_kwargs_r` | `dict`  | parameters for adjustText for right-side annotations                | `None`  |

### Heterogeneity test

| Option     | Type      | Description                                                    | Default |
|------------|-----------|----------------------------------------------------------------|---------|
| `is_q`     | `boolean` | if true, apply the heterogeneity tests by Cochran's Q test.    | `False`  |
| `is_q_mc`  | `str` or `bool` | multiple correction method: `"fdr"` for FDR, `"bon"` for Bonferroni, `False`/`"non"` for no correction. If set to `"fdr"` or `"bon"`, `is_q` will be automatically set to `True`. | `False`  |
| `q_level`  | `float`   | the significance threshold for Cochran's Q test (raw p value or corrected if `is_q_mc` is set). | `0.05`  |
| `anno_het` | `boolean` | annotate only variants with Phet < `q_level`                   | `False` |

### Regression and correlation

| Option      | Type      | Description                                                                                    | Default |
|-------------|-----------|------------------------------------------------------------------------------------------------|---------|
| `reg_text`  | `str`     | regression text format: `"full"` for full equation, `"r"` for correlation only, `"r2"` for r² only | `"full"` |
| `r_se`      | `boolean` | If True, SE for r will be estimated using the jackknife method. (Note: available from v3.4.17) | `False` |
| `r_or_r2`   | `str`     | Display `"r"` or `"r2"` in regression text (when `reg_text` is not `"full"`)                  | `"r"`   |
| `null_beta` | `float`   | null hypothesis value for regression slope (for p-value calculation)                          | `0`     |

$$ s.e.(\hat{r}_{jack}) = \sqrt{ {{n-1}\over{n}} \sum_{i=1}^n(\hat{r_i} -\bar{r}_{jack} )^2 } $$

### Winner's Curse correction

| Option         | Type      | Description                                                                                    | Default |
|----------------|-----------|------------------------------------------------------------------------------------------------|---------|
| `wc_correction` | `bool` or `str` | Winner's Curse correction: `False` for no correction, `"all"` for all variants, `"sig"` for significant variants only, `"sumstats1"` for sumstats1 only, `"sumstats2"` for sumstats2 only | `False` |
| `wc_sig_level`  | `float`   | significance threshold for Winner's Curse correction (when `wc_correction` is `"sig"`, `"sumstats1"`, or `"sumstats2"`) | `5e-8`  |

### Additional options

| Option            | Type      | Description                                                                                    | Default |
|-------------------|-----------|------------------------------------------------------------------------------------------------|---------|
| `get_lead_kwargs` | `dict`    | additional parameters for lead SNP extraction (passed to `_get_sig`)                          | `None`  |
| `drop`            | `boolean` | if True, drop duplicates and NA values                                                       | `False` |
| `include_all`     | `boolean` | if True, include variants that are not significant in either dataset                          | `True`  |
| `fdr`             | `boolean` | if True, apply FDR correction to p-values (requires non-scaled p-values)                      | `False` |
| `allele_match`    | `boolean` | if True, only keep variants where effect alleles match after alignment                        | `False` |
| `build`           | `str`     | genome build version for gene annotation (when `anno="GENENAME"`)                             | `"19"`  |
| `clean_output`    | `boolean` | if True, return cleaned output with simplified column names                                   | `False` |
| `verbose`         | `boolean` | if True, print progress messages                                                               | `False` |
| `log`             | `Log`     | Log object for output messages                                                                 | `Log()` |

## Examples


!!! example "Use gl.Sumstats object (Recommended)"
    ```
    gl1 = gl.Sumstats("bbj_bmi_female.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")
    gl2 = gl.Sumstats("bbj_bmi_male.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")

    # auto extract lead SNPs and compare the effect sizes
    a = gl.compare_effect(path1=gl1,
                          path2=gl2
    )

    ```
    <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/258d3316-7d73-408a-8e20-e560b831861c">

!!! example "Heterogeneity test"
    ```
    gl1 = gl.Sumstats("bbj_bmi_female.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")
    gl2 = gl.Sumstats("bbj_bmi_male.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")

    # is_q : if true, apply the heterogeneity tests by Cochran's Q test.
    # q_level : float, the significance threshold for Cochran's Q test (raw p value or corrected if is_q_mc is set).
    # is_q_mc : multiple correction method ("fdr", "bon", or False/"non")

    a = gl.compare_effect(path1=gl1,
                          path2=gl2,
                          is_q=True,
                          is_q_mc="fdr",  # Use FDR correction for heterogeneity p-values
                          q_level=0.05,
                          legend_mode="full"
                          )
    ```
     <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/b9c4511c-c1eb-4c36-8b72-788d2b21b790">
     
!!! example "Annotation"
    ```
    gl1 = gl.Sumstats("bbj_bmi_female.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")
    gl2 = gl.Sumstats("bbj_bmi_male.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")

    a = gl.compare_effect(path1=gl1,
                          path2=gl2,
                          label=["Female","Male","Both","None"],
                          xylabel_prefix="Per-allele effect size for ",
                          r_se=True, 
                          is_q=True, 
                          anno=True,  # annotate with SNPID
                          anno_het=True, # only annotate variants with significant heterogeneity
                          anno_diff=0.015,    # only annotate variants with a difference in effect size > 0.015
                          sig_level=5e-8,
                          legend_title=r'$\mathregular{ P < 5 x 10^{-8}}$ in:',
                          verbose=True)
                      
    ```
    <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/9e2ac224-63d6-44da-9dbc-1cd43596055f">
    
!!! example "Male-specific and female-specific BMI Sumstats from JENGER"
    
    ```
    !wget -O bbj_bmi_male.txt.gz http://jenger.riken.jp/2analysisresult_qtl_download/
    !wget -O bbj_bmi_female.txt.gz http://jenger.riken.jp/4analysisresult_qtl_download/

    Headers of the files : # SNP    CHR	POS	REF	ALT	Frq	Rsq	BETA	SE	P
    ```
    
    ```
    # Load sumstats as Sumstats objects
    gl1 = gl.Sumstats("bbj_bmi_female.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")
    gl2 = gl.Sumstats("bbj_bmi_male.txt.gz",fmt="gwaslab",snpid="SNP",ea="REF",nea="ALT",sep="\t")
    
    # GWASLab will automatically extract significant variants from both sumstats. 
    a = gl.compare_effect(path1=gl1,
                          path2=gl2,
                          label=["Female","Male","Both","None"],
                          xylabel_prefix="Per-allele effect size for ",
                          anno=True,
                          anno_het=True,
                          anno_diff=0.015,
                          sig_level=5e-8,
                          legend_title=r'$\mathregular{ P < 5 x 10^{-8}}$ in:',
                          verbose=True,
                          save="myplot.png",
                          save_kwargs={"dpi":300,"facecolor":"white"}
    )
    ```
    <img width="500" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/ed462a3d-64cf-4f06-bdb5-13294b4182c8">
    
    Reference: Akiyama, M., Okada, Y., Kanai, M., Takahashi, A., Momozawa, Y., Ikeda, M., ... & Kamatani, Y. (2017). Genome-wide association study identifies 112 new loci for body mass index in the Japanese population. Nature genetics, 49(10), 1458-1467.
