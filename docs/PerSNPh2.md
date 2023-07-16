# Per SNP Heritability

GWASLab provides a simple function to calculate the variance explained by each SNP.

!!! tip "Available since v3.4.20"

```
mysumstats.get_per_snp_r2()
```

It needs effect size `BETA` and effect allele frequency `EAF` for calculation.
If `N` is available, it will also calculate the F-statistics.

### Options

|`.get_per_snp_r2()` options|DataType|Description|Default|
|-|-|-|-|
|`mode`|`q` or `b`|`q`: quantitative trait; `b`: binary trait|`q`|

For quantitative traits (`mode="q"`), GWASLab will use `BETA`, `EAF` to calculate `SNPR2`.

|`.get_per_snp_r2()` options|DataType|Description|Default|
|-|-|-|-|
|`vary`|`float` or `se`| Var(Y); if `vary="se"`,Var(Y) will be estimated using `SE`|1|
|`k`|`int`|k for calculating F|1|

!!! info "For SNP `i`" 

    When `vary=1` and k=1:

    $$ h^2_{i} = 2 \times \beta_i^2 \times EAF \times (1 - EAF) $$

    $$ F = h^2_{i} \times (n -2) / (1 - h^2_{i}) $$

For binary traits  (`mode="b"`), `ncase`, `ncontrol` and `prevalence` are needed to estimate the variance of liability explained by variants:

|`.get_per_snp_r2()` options|DataType|Description|Default|
|-|-|-|-|
|`ncase`|`int`|number of cases|-|
|`ncontrol`|`int`|number of controls|-|
|`prevalence`|`float`|prevalence in general population|-|

!!! quote

    Equation 10 in Lee, S. H., Goddard, M. E., Wray, N. R., & Visscher, P. M. (2012). A better coefficient of determination for genetic profile analysis. Genetic epidemiology, 36(3), 214-224.
    
    - Implementation adopted from TwoSampleMR https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/add_rsq.r

!!! example
    ```
    SNPID	EAF	BETA	N	STATUS
    0	1:725932_G_A	0.9960	-0.0737	166718	9999999
    1	1:725933_A_G	0.0040	0.0737	166718	9999999
    2	1:737801_T_C	0.0051	0.0490	166718	9999999
    3	1:749963_T_TAA	0.8374	0.0213	166718	9999999
    4	1:751343_T_A	0.8593	0.0172	166718	9999999
    
    mysumstats.get_per_snp_r2()
    
    Mon Jul 17 01:39:02 2023 Start to calculate per-SNP heritability...
    Mon Jul 17 01:39:02 2023  -Calculating per-SNP rsq by 2 * (BETA**2) * AF * (1-AF) / Var(y)...
    Mon Jul 17 01:39:02 2023  -Var(y) is provided: 1...
    Mon Jul 17 01:39:02 2023  -Calculating F-statistic: F = [(N-k-1)/k] * (r2/1-r2)... where k = 1
    Mon Jul 17 01:39:02 2023  -For r2, SNPR2 is used.
    Mon Jul 17 01:39:02 2023 Finished calculating per-SNP heritibility!
    
    SNPID	EAF	BETA	N	STATUS	SNPR2	F
    0	1:725932_G_A	0.9960	-0.0737	166718	9999999	0.000043	7.215732
    1	1:725933_A_G	0.0040	0.0737	166718	9999999	0.000043	7.215732
    2	1:737801_T_C	0.0051	0.0490	166718	9999999	0.000024	4.062184
    3	1:749963_T_TAA	0.8374	0.0213	166718	9999999	0.000124	20.600305
    4	1:751343_T_A	0.8593	0.0172	166718	9999999	0.000072	11.927080
    ```
