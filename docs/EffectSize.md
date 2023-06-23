# Comparing effect sizes

## Scatter plot : effect size comparison 

## (WARNing testing!!! will be available from v3.4.17)

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
                   r_se=False,
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
- `path1` and `path2` : the paths to the sumstats. Can also be `gl.Sumstats` Object or `pd.DataFrame` (from v3.4.17). If `gl.Sumstats` Objects are provided, there is no need to set cols_name_list and effect_cols_list.
- `cols_name_list_1` and `cols_name_list_2` : list of column names for variants basic information
- `effect_cols_list_1` and `effect_cols_list_1` : list of column names for effect size-related columns
- `mode` : use beta or OR 
-  examples:
    - [snpid,p,ea,nea]        ,[effect,se]
    - [snpid,p,ea,nea,chr,pos],[effect,se]
    - [snpid,p,ea,nea,chr,pos],[OR,OR_l,OR_h]

!!! note 
    from v3.4.17, you need to specify the parameters using keywords instead of using them as positional arguments for `cols_name_list_1` and `cols_name_list_2`, `effect_cols_list_1` and `effect_cols_list_1` .

### Save
- `save`: path to the saved file.
- `saveargs`: parametrs for plt.savefig(). `{"dpi":300,"facecolor":"white"}`

### Snplist
- `snplist` : optional, specify the variants you want to compare. If None, GWASLab will automatically extract lead variants from both sumstats.

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
- `q_level` : the significance threshold for Cochran's Q test (raw p value).

### R SE
- `r_se`: `boolean` If True, SE for r will be estimated using the jackknife method. (Note: available from v3.4.17)

$$ s.e.(\hat{r}_{jack}) = \sqrt{ {{n-1}\over{n}} \sum_{i=1}^n(\hat{r_i} -\bar{r}_{jack} )^2 } $$

###

## Example:


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

!!! example "Use pandas DataFrame object"
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

!!! example "heterogeneity test"
    ```
    pd1 = pd.read_table("bbj_bmi_female.txt.gz",sep="\t")
    pd2 = pd.read_table("bbj_bmi_male.txt.gz",sep="\t")

    # cols_name_list should be SNPID, P, Effect Allele, Non-Effect allele, Chromosome and Position
    # effect_cols_list should be BETA,SE

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
                      saveargs= {"dpi":300,"facecolor":"white"}
    )
    ```
    <img width="500" alt="image" src="https://user-images.githubusercontent.com/40289485/215021843-4572636d-b1a8-43f5-8f6e-070b09e5270f.png">
    
    Reference: Akiyama, M., Okada, Y., Kanai, M., Takahashi, A., Momozawa, Y., Ikeda, M., ... & Kamatani, Y. (2017). Genome-wide association study identifies 112 new loci for body mass index in the Japanese population. Nature genetics, 49(10), 1458-1467.
