# Per SNP Heritability

GWASLab provides a simple function to calculate the variance explained by each SNP.

```
mysumstats.get_per_snp_r2()
```

It needs effect size `BETA` and effect allele frequency `EAF` for calculation.
If `N` is available, it will also calculate the F-statistics.


!!! info "For SNP `i`" 

    $$ h^2_{i} = 2 \times \beta_i^2 \times EAF \times (1 - EAF) $$

    $$ F = h^2_{i} \times (n -2) / (1 - h^2_{i}) $$

!!! example
    ```
    SNPID	EAF	BETA	N	STATUS
    0	1:725932_G_A	0.9960	-0.0737	166718	9999999
    1	1:725933_A_G	0.0040	0.0737	166718	9999999
    2	1:737801_T_C	0.0051	0.0490	166718	9999999
    3	1:749963_T_TAA	0.8374	0.0213	166718	9999999
    4	1:751343_T_A	0.8593	0.0172	166718	9999999
    
    mysumstats.get_per_snp_r2()
    
    Tue Jan 31 22:44:34 2023 Start to calculate per-SNP heritibility...
    Tue Jan 31 22:44:34 2023  -Calculating per-SNP heritibility by 2 * (BETA**2) * AF * (1-AF)...
    Tue Jan 31 22:44:35 2023  -Calculating F-statistic...
    Tue Jan 31 22:44:35 2023 Finished calculating per-SNP heritibility!
    
    SNPID	EAF	BETA	N	STATUS	SNPH2	F
    0	1:725932_G_A	0.9960	-0.0737	166718	9999999	0.000043	7.215732
    1	1:725933_A_G	0.0040	0.0737	166718	9999999	0.000043	7.215732
    2	1:737801_T_C	0.0051	0.0490	166718	9999999	0.000024	4.062184
    3	1:749963_T_TAA	0.8374	0.0213	166718	9999999	0.000124	20.600305
    4	1:751343_T_A	0.8593	0.0172	166718	9999999	0.000072	11.927080
    ```
