# LDSC in gwaslab

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2025/12/25 22:12:30 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 22:12:30 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 22:12:30 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Loading and filter in only Hapmap3 SNPs

```python
t2d = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             beta="BETA",
             se="SE",
             p="P",
             direction="Dir",
             build="19",             
             n="N", verbose=False)
```

```python
t2d.filter_hapmap3(inplace=True)
```

**stdout:**
```
2025/12/25 22:12:53  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 22:12:53 Start to extract HapMap3 SNPs ...(v4.0.0)
2025/12/25 22:12:53  -Current Dataframe shape : 12557761 x 11 ; Memory usage: 937.22 MB
2025/12/25 22:12:53  -Loading Hapmap3 variants from built-in datasets...
2025/12/25 22:12:54  -Since rsID not in sumstats, CHR:POS( build 19) will be used for matching...
2025/12/25 22:13:14  -Checking if alleles are same...
2025/12/25 22:13:14  -Variants with macthed alleles: 1092430
2025/12/25 22:13:14  -Raw input contains 1092430 Hapmap3 variants based on CHR:POS...
2025/12/25 22:13:15  -Current Dataframe shape : 12557761 x 12 ; Memory usage: 1033.03 MB
2025/12/25 22:13:15 Finished extracting HapMap3 SNPs.
```

## Heritability estimation

```python
# available since v3.4.39
t2d.estimate_h2_by_ldsc(ref_ld_chr = "../../test/ref/eas_ldscores/", 
                         w_ld_chr = "../../test/ref/eas_ldscores/")
```

**stdout:**
```
2025/12/25 22:13:15  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 22:13:15 Start to extract HapMap3 SNPs ...(v4.0.0)
2025/12/25 22:13:15  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 101.03 MB
2025/12/25 22:13:15  -Loading Hapmap3 variants from built-in datasets...
2025/12/25 22:13:16  -rsID will be used for matching...
2025/12/25 22:13:18  -Raw input contains 1213752 Hapmap3 variants based on rsID...
2025/12/25 22:13:18  -Checking if alleles are same...
2025/12/25 22:13:19  -Filtered 0 Hapmap3 variants due to unmatech alleles...
2025/12/25 22:13:19 Finished extracting HapMap3 SNPs.
2025/12/25 22:13:19 Start to run LD score regression ...(v4.0.0)
2025/12/25 22:13:19  -Current Dataframe shape : 1213752 x 12 ; Memory usage: 111.90 MB
2025/12/25 22:13:19  -Run single variate LD score regression:
2025/12/25 22:13:19   -Adopted from LDSC source code: https://github.com/bulik/ldsc
2025/12/25 22:13:19   -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.
2025/12/25 22:13:19  -Arguments:
2025/12/25 22:13:19   -ref_ld_chr:../../test/ref/eas_ldscores/
2025/12/25 22:13:19   -w_ld_chr:../../test/ref/eas_ldscores/
2025/12/25 22:13:20  -LDSC log:
2025/12/25 22:13:20   -Reading reference panel LD Score from ../../test/ref/eas_ldscores/[1-22] ... (ldscore_fromlist)
2025/12/25 22:13:21   -Read reference panel LD Scores for 1208050 SNPs.
2025/12/25 22:13:21   -Removing partitioned LD Scores with zero variance.
2025/12/25 22:13:21   -Reading regression weight LD Score from ../../test/ref/eas_ldscores/[1-22] ... (ldscore_fromlist)
2025/12/25 22:13:22   -Read regression weight LD Scores for 1208050 SNPs.
2025/12/25 22:13:23   -After merging with reference panel LD, 1079390 SNPs remain.
2025/12/25 22:13:24   -After merging with regression SNP LD, 1079390 SNPs remain.
2025/12/25 22:13:24   -Using two-step estimator with cutoff at 30.
2025/12/25 22:13:25  -Results have been stored in .ldsc_h2
2025/12/25 22:13:25  -Current Dataframe shape : 1213752 x 13 ; Memory usage: 121.16 MB
2025/12/25 22:13:25 Finished running LD score regression.
```

```python
t2d.ldsc_h2
```

```
| h2_obs | h2_se | Lambda_gc | Mean_chi2 | Intercept Intercept_se | \ |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.10394433 | 0.00650644 | 1.32982693 | 1.49125406 | 1.09147712 |
|  | Ratio | Ratio_se Catagories |  |  |  |
| 0 | 0.18621142 | 0.02150169 | NA |  |  |
```

## Genetic correlation

```python
bmi_female = gl.Sumstats("../0_sample_data/bbj_bmi_female.txt.gz",fmt="auto",ea="REF",nea="ALT",rsid="SNP",n=70000, sep="\t",build="19",verbose=False)
bmi_male = gl.Sumstats("../0_sample_data/bbj_bmi_male.txt.gz",fmt="auto",ea="REF",nea="ALT",rsid="SNP",n=80000,sep="\t",build="19",verbose=False)
```

- other_traits : a list of gl.Sumstats object
- rg : alias for each trait including the main trait

```python
# available since v3.4.39
t2d.estimate_rg_by_ldsc(other_traits=[bmi_female,bmi_male], 
                               rg="T2D,BMI_female,BMI_male",
                               ref_ld_chr = "../../test/ref/eas_ldscores/", 
                               w_ld_chr = "../../test/ref/eas_ldscores/")
```

**stdout:**
```
2025/12/25 22:14:08  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 22:14:08 Start to extract HapMap3 SNPs ...(v4.0.0)
2025/12/25 22:14:08  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 101.03 MB
2025/12/25 22:14:08  -Loading Hapmap3 variants from built-in datasets...
2025/12/25 22:14:09  -rsID will be used for matching...
2025/12/25 22:14:11  -Raw input contains 1213752 Hapmap3 variants based on rsID...
2025/12/25 22:14:11  -Checking if alleles are same...
2025/12/25 22:14:13  -Filtered 0 Hapmap3 variants due to unmatech alleles...
2025/12/25 22:14:13 Finished extracting HapMap3 SNPs.
2025/12/25 22:14:13 Start to run LD score regression for genetic correlation ...(v4.0.0)
2025/12/25 22:14:13  -Run cross-trait LD score regression:
2025/12/25 22:14:13   -Adopted from LDSC source code: https://github.com/bulik/ldsc
2025/12/25 22:14:13   -Please cite LDSC: Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.
2025/12/25 22:14:13  -Arguments:
2025/12/25 22:14:13   -rg:T2D,BMI_female,BMI_male
2025/12/25 22:14:13   -ref_ld_chr:../../test/ref/eas_ldscores/
2025/12/25 22:14:13   -w_ld_chr:../../test/ref/eas_ldscores/
2025/12/25 22:14:13  -Processing sumstats with alias BMI_female (Study_1)
2025/12/25 22:14:14  -Processing sumstats with alias BMI_male (Study_1)
2025/12/25 22:14:15  -LDSC log:
2025/12/25 22:14:15   -Reading reference panel LD Score from ../../test/ref/eas_ldscores/[1-22] ... (ldscore_fromlist)
2025/12/25 22:14:16   -Read reference panel LD Scores for 1208050 SNPs.
2025/12/25 22:14:16   -Removing partitioned LD Scores with zero variance.
2025/12/25 22:14:16   -Reading regression weight LD Score from ../../test/ref/eas_ldscores/[1-22] ... (ldscore_fromlist)
2025/12/25 22:14:18   -Read regression weight LD Scores for 1208050 SNPs.
2025/12/25 22:14:18   -After merging with reference panel LD, 1079390 SNPs remain.
2025/12/25 22:14:19   -After merging with regression SNP LD, 1079390 SNPs remain.
2025/12/25 22:14:19   -Computing rg for phenotype 2/3
2025/12/25 22:14:19   -Read summary statistics for 5961600 SNPs.
2025/12/25 22:14:26   -After merging with summary statistics, 1021957 SNPs remain.
2025/12/25 22:14:27   -1021957 SNPs with valid alleles.
2025/12/25 22:14:28   -Heritability of phenotype 1
2025/12/25 22:14:28   -Total Observed scale h2: 0.09687973 (0.00885768)
                       -Lambda GC: 1.34889231
                       -Mean Chi^2: 1.50928684
                       -Intercept: 1.12024032 (0.02516315)
                       -Ratio: 0.23609547 (0.04940861)
2025/12/25 22:14:28   -Heritability of phenotype 2/3
2025/12/25 22:14:28   -Total Observed scale h2: 0.19319917 (0.01213263)
                       -Lambda GC: 1.2537503
                       -Mean Chi^2: 1.31732059
                       -Intercept: 1.0306093 (0.01022627)
                       -Ratio: 0.09646175 (0.03222695)
2025/12/25 22:14:28   -Genetic Covariance
2025/12/25 22:14:28   -Total Observed scale gencov: 0.04387074 (0.00729069)
                       -Mean z1*z2: 0.12628027
                       -Intercept: 0.01918354 (0.00945413)
2025/12/25 22:14:28   -Genetic Correlation
2025/12/25 22:14:28   -Genetic Correlation: 0.32066816 (0.06227918)
                       -Z-score: 5.14888223
                       -P: 2.62043342e-07

2025/12/25 22:14:28   -Computing rg for phenotype 3/3
2025/12/25 22:14:28   -Read summary statistics for 5961600 SNPs.
2025/12/25 22:14:34   -After merging with summary statistics, 1021957 SNPs remain.
2025/12/25 22:14:35   -1021957 SNPs with valid alleles.
2025/12/25 22:14:36   -Heritability of phenotype 3/3
2025/12/25 22:14:36   -Total Observed scale h2: 0.17595184 (0.01138382)
                       -Lambda GC: 1.25741577
                       -Mean Chi^2: 1.34453038
                       -Intercept: 1.04768706 (0.01071735)
                       -Ratio: 0.13841176 (0.03110712)
2025/12/25 22:14:36   -Genetic Covariance
2025/12/25 22:14:36   -Total Observed scale gencov: 0.02739806 (0.00705736)
                       -Mean z1*z2: 0.07617508
                       -Intercept: 0.00535343 (0.0118926)
2025/12/25 22:14:36   -Genetic Correlation
2025/12/25 22:14:36   -Genetic Correlation: 0.20984874 (0.06008801)
                       -Z-score: 3.49235659
                       -P: 0.00047878

2025/12/25 22:14:36 Summary of Genetic Correlation Results
 p1         p2       rg       se        z            p   h2_obs  h2_obs_se   h2_int  h2_int_se  gcov_int  gcov_int_se
T2D BMI_female 0.320668 0.062279 5.148882 2.620433e-07 0.193199   0.012133 1.030609   0.010226  0.019184     0.009454
T2D   BMI_male 0.209849 0.060088 3.492357 4.787786e-04 0.175952   0.011384 1.047687   0.010717  0.005353     0.011893

2025/12/25 22:14:36  -Results have been stored in .ldsc_rg
2025/12/25 22:14:36 Finished running LD score regression for genetic correlation.
```

```
| p1 | p2 | rg | se | z | p | h2_obs | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | T2D | BMI_female | 0.320668 | 0.062279 | 5.148882 | 2.620433e-07 | 0.193199 |
| 1 | T2D | BMI_male | 0.209849 | 0.060088 | 3.492357 | 4.787786e-04 | 0.175952 |
|  | h2_obs_se | h2_int | h2_int_se | gcov_int | gcov_int_se |  |  |
| 0 | 0.012133 | 1.030609 | 0.010226 | 0.019184 | 0.009454 |  |  |
| 1 | 0.011384 | 1.047687 | 0.010717 | 0.005353 | 0.011893 |  |  |
```

```python
t2d.ldsc_rg
```

```
| p1 | p2 | rg | se | z | p | h2_obs | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | T2D | BMI_female | 0.320668 | 0.062279 | 5.148882 | 2.620433e-07 | 0.193199 |
| 1 | T2D | BMI_male | 0.209849 | 0.060088 | 3.492357 | 4.787786e-04 | 0.175952 |
|  | h2_obs_se | h2_int | h2_int_se | gcov_int | gcov_int_se |  |  |
| 0 | 0.012133 | 1.030609 | 0.010226 | 0.019184 | 0.009454 |  |  |
| 1 | 0.011384 | 1.047687 | 0.010717 | 0.005353 | 0.011893 |  |  |
```

### visualization using plot_rg

```python
fig, ax, log, df = gl.plot_rg(t2d.ldsc_rg, p="p")
```

**stdout:**
```
2025/12/25 22:15:04 Start to create ldsc genetic correlation heatmap...
2025/12/25 22:15:04 Configured plot style for plot_rg:None
2025/12/25 22:15:04 Raw dataset records: 2
2025/12/25 22:15:04  -Raw dataset non-NA records: 2
2025/12/25 22:15:04 Filling diagnal line and duplicated pair for plotting...
2025/12/25 22:15:04 Valid unique trait pairs: 2
2025/12/25 22:15:04  -Valid unique trait1: 1
2025/12/25 22:15:04  -Valid unique trait2: 2
2025/12/25 22:15:04  -Significant correlations with P < 0.05: 2
2025/12/25 22:15:04  -Significant correlations after Bonferroni correction: 2
2025/12/25 22:15:04  -Significant correlations with FDR <0.05: 2
2025/12/25 22:15:04 Plotting heatmap...
2025/12/25 22:15:04 Full cell : fdr-corrected P == 0.05
2025/12/25 22:15:04 P value annotation text (Order: Bon -> FDR -> Pnom): 
2025/12/25 22:15:04  -* : non-corrected P < 0.05 
2025/12/25 22:15:04  -** : fdr-corrected P < 0.05 
2025/12/25 22:15:04  -*** : bon-corrected P < 0.05 
2025/12/25 22:15:04 Finished creating ldsc genetic correlation heatmap!
```

![Output image](images/notebooks/ldsc_in_gwaslab_img_0.png)

## Partitioned h2

```python
t2d.estimate_partitioned_h2_by_ldsc(       ref_ld_chr = "/home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.", 
                               w_ld_chr = "/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.",
                               frqfile_chr= "/home/yunye/tools/ldsc/eas_frq/1000G.EAS.QC.",
                               overlap_annot = True, 
                               print_coefficients = True, 
                               print_delete_vals=True)
```

**stdout:**
```
2025/12/25 22:15:07  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 22:15:07 Start to extract HapMap3 SNPs ...(v4.0.0)
2025/12/25 22:15:07  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 101.03 MB
2025/12/25 22:15:07  -Loading Hapmap3 variants from built-in datasets...
2025/12/25 22:15:08  -rsID will be used for matching...
2025/12/25 22:15:11  -Raw input contains 1213752 Hapmap3 variants based on rsID...
2025/12/25 22:15:11  -Checking if alleles are same...
2025/12/25 22:15:12  -Filtered 0 Hapmap3 variants due to unmatech alleles...
2025/12/25 22:15:12 Finished extracting HapMap3 SNPs.
2025/12/25 22:15:12 Start to run LD score regression ...(v4.0.0)
2025/12/25 22:15:12  -Run partitioned LD score regression:
2025/12/25 22:15:12   -Adopted from LDSC source code: https://github.com/bulik/ldsc
2025/12/25 22:15:12   -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.
2025/12/25 22:15:12  -Arguments:
2025/12/25 22:15:12   -ref_ld_chr:/home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.
2025/12/25 22:15:12   -w_ld_chr:/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.
2025/12/25 22:15:12   -frqfile_chr:/home/yunye/tools/ldsc/eas_frq/1000G.EAS.QC.
2025/12/25 22:15:12   -overlap_annot:True
2025/12/25 22:15:12   -print_coefficients:True
2025/12/25 22:15:12   -print_delete_vals:True
2025/12/25 22:15:12  -LDSC log:
2025/12/25 22:15:13   -Reading reference panel LD Score from /home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.[1-22] ... (ldscore_fromlist)
```

```python
t2d.estimate_h2_cts_by_ldsc(ref_ld_chr = "/home/yunye/tools/ldsc/eas_baseline/baseline1_2/baseline.", 
                            ref_ld_chr_cts = "/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Multi_tissue_gene_expr.EAS.ldcts.new",
                            w_ld_chr = "/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.")
```

```python
# available since v3.4.40
t2d.estimate_partitioned_h2_by_ldsc(       ref_ld_chr = "/home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.", 
                               w_ld_chr = "/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.",
                               frqfile_chr= "/home/yunye/tools/ldsc/eas_frq/1000G.EAS.QC.",
                               overlap_annot = True, 
                               print_coefficients = True, 
                               print_delete_vals=True)
```

**stdout:**
```
2024/12/20 12:55:59 Start to extract HapMap3 SNPs...v3.5.4
2024/12/20 12:55:59  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 118.32 MB
2024/12/20 12:55:59  -Loading Hapmap3 variants from built-in datasets...
2024/12/20 12:56:00  -rsID will be used for matching...
2024/12/20 12:56:00  -Raw input contains 1092430 Hapmap3 variants based on rsID...
2024/12/20 12:56:00 Start to run LD score regression...v3.5.4
2024/12/20 12:56:00  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 118.32 MB
2024/12/20 12:56:00  -Run partitioned LD score regression:
2024/12/20 12:56:00   -Adopted from LDSC source code: https://github.com/bulik/ldsc
2024/12/20 12:56:00   -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.
2024/12/20 12:56:00  -Arguments:
2024/12/20 12:56:00   -ref_ld_chr:/home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.
2024/12/20 12:56:00   -w_ld_chr:/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.
2024/12/20 12:56:00   -frqfile_chr:/home/yunye/tools/ldsc/eas_frq/1000G.EAS.QC.
2024/12/20 12:56:00   -overlap_annot:True
2024/12/20 12:56:00   -print_coefficients:True
2024/12/20 12:56:00   -print_delete_vals:True
2024/12/20 12:56:00  -LDSC log:
2024/12/20 12:56:00   -Reading reference panel LD Score from /home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.[1-22] ... (ldscore_fromlist)
2024/12/20 12:56:11   -Read reference panel LD Scores for 1071448 SNPs.
2024/12/20 12:56:11   -Removing partitioned LD Scores with zero variance.
2024/12/20 12:56:12   -Reading regression weight LD Score from /home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.[1-22] ... (ldscore_fromlist)
2024/12/20 12:56:13   -Read regression weight LD Scores for 1071448 SNPs.
2024/12/20 12:56:13   -After merging with reference panel LD, 1061153 SNPs remain.
2024/12/20 12:56:15   -After merging with regression SNP LD, 1061153 SNPs remain.
2024/12/20 12:56:15   -Removed 49 SNPs with chi^2 > 191.764 (1061104 SNPs remain)
2024/12/20 12:56:22   -Printing block jackknife delete values to ldsc.delete.
2024/12/20 12:56:22   -Printing partitioned block jackknife delete values to ldsc.part_delete.
2024/12/20 12:56:22   -Reading annot matrix from /home/yunye/tools/ldsc/eas_baseline/baselineLD2_2/baselineLD.[1-22] ... (annot)
2024/12/20 12:57:14  -Results have been stored in .ldsc_h2
2024/12/20 12:57:14 Finished running LD score regression.
```

```python
t2d.ldsc_partitioned_h2_summary
```

```
| h2_obs | h2_se | Lambda_gc | Mean_chi2 | Intercept Intercept_se | \ |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.11773265 | 0.00727336 | 1.33210306 | 1.47668034 | 1.07649298 |
|  | Ratio | Ratio_se |  |  |  |
| 0 | 0.16047018 | 0.02577482 |  |  |  |
```

```python
t2d.ldsc_partitioned_h2_results
```

```
| Category | Prop._SNPs | Prop._h2 | \ |
| --- | --- | --- | --- |
| 0 | baseL2_0 | 1.000000 | 1.000000 |
| 1 | Coding_UCSC.bedL2_0 | 0.014196 | 0.068254 |
| 2 | Coding_UCSC.bed.flanking.500L2_0 | 0.049260 | 0.072198 |
| 3 | Conserved_LindbladToh.bedL2_0 | 0.024551 | 0.088877 |
| 4 | Conserved_LindbladToh.bed.flanking.500L2_0 | 0.305295 | 0.649833 |
| ... | ... | ... | ... |
| 92 | 3.149492e-07 | 1.103395e-07 | 2.854366 |
| 93 | 5.790740e-09 | 6.996917e-08 | 0.082761 |
| 94 | 4.797773e-09 | 5.983189e-09 | 0.801875 |
| 95 | 2.645270e-07 | 3.249973e-07 | 0.813936 |
| 96 -1.734589e-08 | 6.804492e-07 | -0.025492 |  |

*[97 rows x 10 columns]*
```

## Cell type specific 

```python
# available since v3.4.40
t2d.estimate_h2_cts_by_ldsc(ref_ld_chr = "/home/yunye/tools/ldsc/eas_baseline/baseline1_2/baseline.", 
                            ref_ld_chr_cts = "/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Multi_tissue_gene_expr.EAS.ldcts.new",
                            w_ld_chr = "/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.")
```

**stdout:**
```
2024/12/20 12:57:14 Start to extract HapMap3 SNPs...v3.5.4
2024/12/20 12:57:14  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 118.32 MB
2024/12/20 12:57:14  -Loading Hapmap3 variants from built-in datasets...
2024/12/20 12:57:15  -rsID will be used for matching...
2024/12/20 12:57:15  -Raw input contains 1092430 Hapmap3 variants based on rsID...
2024/12/20 12:57:15 Start to run LD score regression...v3.5.4
2024/12/20 12:57:15  -Current Dataframe shape : 1092430 x 12 ; Memory usage: 118.32 MB
2024/12/20 12:57:15  -Run cell type specific LD score regression:
2024/12/20 12:57:15   -Adopted from LDSC source code: https://github.com/bulik/ldsc
2024/12/20 12:57:15   -Please cite LDSC: Finucane, H. K., Reshef, Y. A., Anttila, V., Slowikowski, K., Gusev, A., Byrnes, A., ... & Price, A. L. (2018). Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature genetics, 50(4), 621-629.
2024/12/20 12:57:15  -Arguments:
2024/12/20 12:57:15   -ref_ld_chr:/home/yunye/tools/ldsc/eas_baseline/baseline1_2/baseline.
2024/12/20 12:57:15   -ref_ld_chr_cts:/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Multi_tissue_gene_expr.EAS.ldcts.new
2024/12/20 12:57:15   -w_ld_chr:/home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.
2024/12/20 12:57:15  -LDSC log:
2024/12/20 12:57:15   -Reading reference panel LD Score from /home/yunye/tools/ldsc/eas_baseline/baseline1_2/baseline.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:22   -Read reference panel LD Scores for 1071448 SNPs.
2024/12/20 12:57:22   -Removing partitioned LD Scores with zero variance.
2024/12/20 12:57:22   -Reading regression weight LD Score from /home/yunye/tools/ldsc/eas_weights/weights.EAS.hm3_noMHC.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:23   -Read regression weight LD Scores for 1071448 SNPs.
2024/12/20 12:57:24   -After merging with reference panel LD, 1061153 SNPs remain.
2024/12/20 12:57:25   -After merging with regression SNP LD, 1061153 SNPs remain.
2024/12/20 12:57:25   -Removed 49 SNPs with chi^2 > 191.764 (1061104 SNPs remain)
2024/12/20 12:57:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.1.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:27   -Performing regression.
2024/12/20 12:57:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.2.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:32   -Performing regression.
2024/12/20 12:57:36   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.3.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:38   -Performing regression.
2024/12/20 12:57:40   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.4.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:43   -Performing regression.
2024/12/20 12:57:45   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.5.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:48   -Performing regression.
2024/12/20 12:57:50   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.6.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:53   -Performing regression.
2024/12/20 12:57:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.7.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:57:58   -Performing regression.
2024/12/20 12:58:00   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.8.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:02   -Performing regression.
2024/12/20 12:58:04   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.9.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:06   -Performing regression.
2024/12/20 12:58:08   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.10.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:10   -Performing regression.
2024/12/20 12:58:12   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.11.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:14   -Performing regression.
2024/12/20 12:58:16   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.12.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:19   -Performing regression.
2024/12/20 12:58:21   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.13.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:23   -Performing regression.
2024/12/20 12:58:26   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.14.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:29   -Performing regression.
2024/12/20 12:58:31   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.15.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:34   -Performing regression.
2024/12/20 12:58:36   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.16.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:38   -Performing regression.
2024/12/20 12:58:39   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.17.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:41   -Performing regression.
2024/12/20 12:58:44   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.18.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:46   -Performing regression.
2024/12/20 12:58:49   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.19.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:51   -Performing regression.
2024/12/20 12:58:54   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.20.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:58:56   -Performing regression.
```

**stdout:**
```
2024/12/20 12:58:59   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.21.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:01   -Performing regression.
2024/12/20 12:59:04   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.22.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:06   -Performing regression.
2024/12/20 12:59:08   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.23.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:10   -Performing regression.
2024/12/20 12:59:12   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.24.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:14   -Performing regression.
2024/12/20 12:59:16   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.25.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:18   -Performing regression.
2024/12/20 12:59:20   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.26.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:22   -Performing regression.
2024/12/20 12:59:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.27.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:27   -Performing regression.
2024/12/20 12:59:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.28.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:32   -Performing regression.
2024/12/20 12:59:35   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.29.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:37   -Performing regression.
2024/12/20 12:59:39   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.30.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:41   -Performing regression.
2024/12/20 12:59:43   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.31.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:45   -Performing regression.
2024/12/20 12:59:47   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.32.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:49   -Performing regression.
2024/12/20 12:59:51   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.33.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:54   -Performing regression.
2024/12/20 12:59:56   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.34.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 12:59:59   -Performing regression.
2024/12/20 13:00:01   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.35.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:03   -Performing regression.
2024/12/20 13:00:05   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.36.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:07   -Performing regression.
2024/12/20 13:00:10   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.37.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:12   -Performing regression.
2024/12/20 13:00:14   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.38.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:16   -Performing regression.
2024/12/20 13:00:17   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.39.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:20   -Performing regression.
2024/12/20 13:00:21   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.40.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:23   -Performing regression.
2024/12/20 13:00:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.41.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:27   -Performing regression.
2024/12/20 13:00:29   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.42.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:32   -Performing regression.
2024/12/20 13:00:34   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.43.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:36   -Performing regression.
2024/12/20 13:00:39   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.44.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:41   -Performing regression.
2024/12/20 13:00:43   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.45.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:45   -Performing regression.
2024/12/20 13:00:47   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.46.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:49   -Performing regression.
2024/12/20 13:00:51   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.47.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
```

**stdout:**
```
2024/12/20 13:00:53   -Performing regression.
2024/12/20 13:00:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.48.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:00:58   -Performing regression.
2024/12/20 13:01:00   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.49.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:03   -Performing regression.
2024/12/20 13:01:05   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.50.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:08   -Performing regression.
2024/12/20 13:01:09   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.51.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:11   -Performing regression.
2024/12/20 13:01:14   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.52.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:16   -Performing regression.
2024/12/20 13:01:19   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.53.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/GTEx.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:22   -Performing regression.
2024/12/20 13:01:24   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.1.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:27   -Performing regression.
2024/12/20 13:01:29   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.2.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:32   -Performing regression.
2024/12/20 13:01:34   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.3.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:37   -Performing regression.
2024/12/20 13:01:40   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.4.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:42   -Performing regression.
2024/12/20 13:01:45   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.5.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:47   -Performing regression.
2024/12/20 13:01:49   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.6.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:51   -Performing regression.
2024/12/20 13:01:53   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.7.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:55   -Performing regression.
2024/12/20 13:01:57   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.8.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:01:59   -Performing regression.
2024/12/20 13:02:02   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.9.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:04   -Performing regression.
2024/12/20 13:02:07   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.10.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:09   -Performing regression.
2024/12/20 13:02:11   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.11.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:14   -Performing regression.
2024/12/20 13:02:17   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.12.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:19   -Performing regression.
2024/12/20 13:02:20   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.13.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:23   -Performing regression.
2024/12/20 13:02:24   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.14.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:26   -Performing regression.
2024/12/20 13:02:28   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.15.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:31   -Performing regression.
2024/12/20 13:02:33   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.16.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:36   -Performing regression.
2024/12/20 13:02:39   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.17.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:41   -Performing regression.
2024/12/20 13:02:44   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.18.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:46   -Performing regression.
2024/12/20 13:02:49   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.19.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:51   -Performing regression.
2024/12/20 13:02:52   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.20.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:02:55   -Performing regression.
2024/12/20 13:02:56   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.21.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
```

**stdout:**
```
2024/12/20 13:02:58   -Performing regression.
2024/12/20 13:03:01   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.22.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:04   -Performing regression.
2024/12/20 13:03:07   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.23.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:09   -Performing regression.
2024/12/20 13:03:12   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.24.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:15   -Performing regression.
2024/12/20 13:03:18   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.25.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:20   -Performing regression.
2024/12/20 13:03:22   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.26.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:24   -Performing regression.
2024/12/20 13:03:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.27.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:28   -Performing regression.
2024/12/20 13:03:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.28.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:33   -Performing regression.
2024/12/20 13:03:37   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.29.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:39   -Performing regression.
2024/12/20 13:03:43   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.30.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:45   -Performing regression.
2024/12/20 13:03:47   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.31.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:49   -Performing regression.
2024/12/20 13:03:51   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.32.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:54   -Performing regression.
2024/12/20 13:03:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.33.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:03:57   -Performing regression.
2024/12/20 13:03:59   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.34.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:02   -Performing regression.
2024/12/20 13:04:04   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.35.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:06   -Performing regression.
2024/12/20 13:04:09   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.36.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:12   -Performing regression.
2024/12/20 13:04:15   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.37.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:17   -Performing regression.
2024/12/20 13:04:21   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.38.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:23   -Performing regression.
2024/12/20 13:04:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.39.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:27   -Performing regression.
2024/12/20 13:04:28   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.40.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:31   -Performing regression.
2024/12/20 13:04:32   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.41.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:35   -Performing regression.
2024/12/20 13:04:36   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.42.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:39   -Performing regression.
2024/12/20 13:04:40   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.43.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:43   -Performing regression.
2024/12/20 13:04:46   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.44.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:48   -Performing regression.
2024/12/20 13:04:51   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.45.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:53   -Performing regression.
2024/12/20 13:04:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.46.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:04:57   -Performing regression.
2024/12/20 13:05:03   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.47.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:05   -Performing regression.
2024/12/20 13:05:08   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.48.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
```

**stdout:**
```
2024/12/20 13:05:11   -Performing regression.
2024/12/20 13:05:13   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.49.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:15   -Performing regression.
2024/12/20 13:05:18   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.50.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:20   -Performing regression.
2024/12/20 13:05:22   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.51.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:24   -Performing regression.
2024/12/20 13:05:26   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.52.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:28   -Performing regression.
2024/12/20 13:05:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.53.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:32   -Performing regression.
2024/12/20 13:05:35   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.54.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:37   -Performing regression.
2024/12/20 13:05:40   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.55.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:42   -Performing regression.
2024/12/20 13:05:45   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.56.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:47   -Performing regression.
2024/12/20 13:05:49   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.57.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:51   -Performing regression.
2024/12/20 13:05:54   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.58.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:05:57   -Performing regression.
2024/12/20 13:05:58   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.59.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:01   -Performing regression.
2024/12/20 13:06:02   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.60.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:04   -Performing regression.
2024/12/20 13:06:06   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.61.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:08   -Performing regression.
2024/12/20 13:06:10   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.62.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:12   -Performing regression.
2024/12/20 13:06:15   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.63.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:17   -Performing regression.
2024/12/20 13:06:19   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.64.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:21   -Performing regression.
2024/12/20 13:06:24   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.65.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:26   -Performing regression.
2024/12/20 13:06:28   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.66.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:30   -Performing regression.
2024/12/20 13:06:31   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.67.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:33   -Performing regression.
2024/12/20 13:06:36   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.68.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:39   -Performing regression.
2024/12/20 13:06:41   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.69.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:43   -Performing regression.
2024/12/20 13:06:46   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.70.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:48   -Performing regression.
2024/12/20 13:06:51   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.71.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:53   -Performing regression.
2024/12/20 13:06:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.72.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:06:57   -Performing regression.
2024/12/20 13:06:59   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.73.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:01   -Performing regression.
2024/12/20 13:07:02   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.74.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:05   -Performing regression.
2024/12/20 13:07:07   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.75.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
```

**stdout:**
```
2024/12/20 13:07:10   -Performing regression.
2024/12/20 13:07:12   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.76.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:15   -Performing regression.
2024/12/20 13:07:17   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.77.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:19   -Performing regression.
2024/12/20 13:07:21   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.78.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:24   -Performing regression.
2024/12/20 13:07:26   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.79.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:29   -Performing regression.
2024/12/20 13:07:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.80.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:33   -Performing regression.
2024/12/20 13:07:34   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.81.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:36   -Performing regression.
2024/12/20 13:07:38   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.82.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:40   -Performing regression.
2024/12/20 13:07:43   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.83.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:45   -Performing regression.
2024/12/20 13:07:48   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.84.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:50   -Performing regression.
2024/12/20 13:07:53   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.85.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:55   -Performing regression.
2024/12/20 13:07:57   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.86.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:07:59   -Performing regression.
2024/12/20 13:08:00   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.87.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:03   -Performing regression.
2024/12/20 13:08:04   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.88.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:06   -Performing regression.
2024/12/20 13:08:08   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.89.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:10   -Performing regression.
2024/12/20 13:08:13   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.90.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:15   -Performing regression.
2024/12/20 13:08:16   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.91.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:19   -Performing regression.
2024/12/20 13:08:20   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.92.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:23   -Performing regression.
2024/12/20 13:08:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.93.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:27   -Performing regression.
2024/12/20 13:08:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.94.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:32   -Performing regression.
2024/12/20 13:08:34   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.95.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:36   -Performing regression.
2024/12/20 13:08:38   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.96.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:40   -Performing regression.
2024/12/20 13:08:41   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.97.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:44   -Performing regression.
2024/12/20 13:08:45   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.98.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:47   -Performing regression.
2024/12/20 13:08:50   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.99.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:52   -Performing regression.
2024/12/20 13:08:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.100.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:08:57   -Performing regression.
2024/12/20 13:09:00   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.101.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:02   -Performing regression.
2024/12/20 13:09:04   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.102.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
```

**stdout:**
```
2024/12/20 13:09:06   -Performing regression.
2024/12/20 13:09:08   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.103.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:10   -Performing regression.
2024/12/20 13:09:11   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.104.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:13   -Performing regression.
2024/12/20 13:09:16   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.105.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:18   -Performing regression.
2024/12/20 13:09:22   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.106.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:24   -Performing regression.
2024/12/20 13:09:27   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.107.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:29   -Performing regression.
2024/12/20 13:09:32   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.108.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:34   -Performing regression.
2024/12/20 13:09:36   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.109.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:38   -Performing regression.
2024/12/20 13:09:39   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.110.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:41   -Performing regression.
2024/12/20 13:09:44   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.111.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:47   -Performing regression.
2024/12/20 13:09:49   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.112.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:52   -Performing regression.
2024/12/20 13:09:55   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.113.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:09:57   -Performing regression.
2024/12/20 13:09:59   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.114.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:01   -Performing regression.
2024/12/20 13:10:04   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.115.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:06   -Performing regression.
2024/12/20 13:10:07   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.116.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:10   -Performing regression.
2024/12/20 13:10:11   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.117.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:14   -Performing regression.
2024/12/20 13:10:16   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.118.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:19   -Performing regression.
2024/12/20 13:10:21   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.119.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:24   -Performing regression.
2024/12/20 13:10:25   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.120.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:27   -Performing regression.
2024/12/20 13:10:30   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.121.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:32   -Performing regression.
2024/12/20 13:10:34   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.122.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:37   -Performing regression.
2024/12/20 13:10:38   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.123.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:40   -Performing regression.
2024/12/20 13:10:42   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.124.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:44   -Performing regression.
2024/12/20 13:10:47   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.125.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:49   -Performing regression.
2024/12/20 13:10:52   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.126.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:54   -Performing regression.
2024/12/20 13:10:57   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.127.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:10:59   -Performing regression.
2024/12/20 13:11:02   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.128.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:04   -Performing regression.
2024/12/20 13:11:07   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.129.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
```

**stdout:**
```
2024/12/20 13:11:10   -Performing regression.
2024/12/20 13:11:11   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.130.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:13   -Performing regression.
2024/12/20 13:11:15   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.131.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:17   -Performing regression.
2024/12/20 13:11:20   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.132.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:22   -Performing regression.
2024/12/20 13:11:24   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.133.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:26   -Performing regression.
2024/12/20 13:11:28   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.134.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:31   -Performing regression.
2024/12/20 13:11:33   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.135.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:36   -Performing regression.
2024/12/20 13:11:38   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.136.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:40   -Performing regression.
2024/12/20 13:11:42   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.137.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:44   -Performing regression.
2024/12/20 13:11:46   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.138.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:48   -Performing regression.
2024/12/20 13:11:51   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.139.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:53   -Performing regression.
2024/12/20 13:11:56   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.140.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:11:58   -Performing regression.
2024/12/20 13:12:01   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.141.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:03   -Performing regression.
2024/12/20 13:12:06   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.142.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:08   -Performing regression.
2024/12/20 13:12:11   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.143.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:13   -Performing regression.
2024/12/20 13:12:15   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.144.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:17   -Performing regression.
2024/12/20 13:12:18   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.145.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:21   -Performing regression.
2024/12/20 13:12:23   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.146.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:25   -Performing regression.
2024/12/20 13:12:27   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.147.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:29   -Performing regression.
2024/12/20 13:12:32   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.148.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:34   -Performing regression.
2024/12/20 13:12:37   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.149.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:39   -Performing regression.
2024/12/20 13:12:41   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.150.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:44   -Performing regression.
2024/12/20 13:12:45   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.151.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:47   -Performing regression.
2024/12/20 13:12:49   -Reading cts reference panel LD Score from /home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.152.,/home/yunye/tools/ldsc/Multi_tissue_gene_expr_EAS_1000G_v3_ldscores/Franke.EAS.control.[1-22] ... (ldscore_fromlist)
2024/12/20 13:12:51   -Performing regression.
2024/12/20 13:12:54  -Results have been stored in .ldsc_partitioned_h2
2024/12/20 13:12:54 Finished running LD score regression.
```

```python
t2d.ldsc_h2_cts
```

```
| Name | Coefficient | \ |
| --- | --- | --- |
| 40 | Pancreas | 9.460678e-09 |
| 186 | A03.556.875.Upper.Gastrointestinal.Tract | 9.761378e-09 |
| 149 | A03.556.124.526.767.Rectum | 9.773394e-09 |
| 70 | A03.556.249.249.209.Cecum | 9.289977e-09 |
| 152 | A03.556.875.875.Stomach | 8.648790e-09 |
| ... | ... | ... |
| 91 | 2.852965e-09 | 0.984537 |
| 134 | 2.841565e-09 | 0.985709 |
| 18 | 2.177073e-09 | 0.993719 |
| 135 | 2.721569e-09 | 0.997003 |
| 55 | 2.520982e-09 | 0.998254 |

*[205 rows x 4 columns]*
```

```python
extra_javascript:
  - javascripts/mathjax.js
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
  - https://unpkg.com/tablesort@5.3.0/dist/tablesort.min.js
  - javascripts/tablesort.js
  - javascripts/hide-sidebar-tutorial.js

extra_css:
  - stylesheets/extra.css

plugins:
  - search
```

```python
     line_spans: __span
```

```python
      anchor_linenums: true
  - pymdownx.inlinehilite
```

```python
  - pymdownx.superfences
```

```python
markdown_extensions:
  - toc:
      toc_depth: 3
  - admonition
  - pymdownx.details
  - pymdownx.arithmatex:
      generic: true
  - pymdownx.superfences
  - pymdownx.highlight:
      anchor_linenums: true
      line_spans: __span
      pygments_lang_class: true
  - pymdownx.inlinehilite
  - pymdownx.snippets
  - tables
```

```python
extra_javascript:
  - javascripts/mathjax.js
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
  - https://unpkg.com/tablesort@5.3.0/dist/tablesort.min.js
  - javascripts/tablesort.js
  - javascripts/hide-sidebar-tutorial.js

extra_css:
  - stylesheets/extra.css
  - pymdownx.superfences
```

```python
markdown_extensions:
  - toc:
      toc_depth: 3
  - admonition
  - pymdownx.details
  - pymdownx.arithmatex:
      generic: true
  - pymdownx.highlight:
      anchor_linenums: true
  - pymdownx.superfences 
  - pymdownx.snippets
  - tables
```

```python
extra_javascript:
  - javascripts/mathjax.js
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
  - https://unpkg.com/tablesort@5.3.0/dist/tablesort.min.js
  - javascripts/tablesort.js
  - javascripts/hide-sidebar-tutorial.js
```

```python


extra_css:
  - stylesheets/extra.css
```

```python
extra_javascript:
  - https://polyfill.io/v3/polyfill.min.js?features=es6
  - https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js
  - https://unpkg.com/tablesort@5.3.0/dist/tablesort.min.js
  - javascripts/tablesort.js
  - javascripts/hide-sidebar-tutorial.js

extra_css:
  - stylesheets/extra.css
```
