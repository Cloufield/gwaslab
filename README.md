![gwaslab_logo](https://cloufield.github.io/gwaslab/images/index_logo.jpg)
A collection of handy python scripts for GWAS. 
Just want to make life eaiser and save myself from repetitive work.

For usage, please check GWASLab document at https://cloufield.github.io/gwaslab/ .

## What you can do with gwaslab:
### Standardization, Normalization & Harmonization

- CHR and POS notation standardization

- Variant POS and allele normalization

- Genome build : Liftover

- Reference allele alignment using a reference genome sequence

- rsID assignment based on CHR, POS, REF and ALT

- CHR POS assignment based on rsID using a reference VCF

- Palindromic SNPs and indels strand inference using a reference VCF

- Check allele frequency discrepancy using a reference VCF

### Quality control, Value conversion & Filtering

- Statistics sanity check

- Equivalent statistics conversion
  
  - BETA/SE , OR/OR_95L/OR_95U
  
  - P, Z, CHISQ, MLOG10

- Extract/exclude hapmap3 variants 

- Extract/exclude MHC variants

- Filtering values.

### Visualization

- Mqq plot : Manhattan plot and QQ plot side by side

- Heatmap : ldsc-rg genetic correlation matrix

- Scatter Plot : variant effect size comparison with sumstats

- Scatter Plot : allele frequency comparison 

### Other Utilities

- Read ldsc h2 or rg outputs directly as DataFrames

- Extract lead SNPs given a window size

- Logging : keep a complete record of manipulations from raw data to munged data

- Formating GWAS sumstats in certain formats
  
  - LDSC / MAGMA / METAL / MR-MEGA / FUMA
...




## Requirements:
- Python >= 3
- pyVCF >= 0.6.8
- Biopython >= 1.79
- liftover >= 1.1.13
- pandas >= 1.2.4
- numpy >= 1.21.2
- matplotlib>3.5
- seaborn >= 0.11.1
- scipy >= 1.6.2
- adjustText

## Install:
```
pip install gwaslab
```
Current version: 3.0.0

# Usage:

Input: pandas dataframe or directly import form text files.

### Create Manhattan plot and QQ plot with just one line
```
import gwaslab as gl

## load gwaslab Sumstats object
AF = gl.Sumstats("./AF_bbj.txt.gz",
                   snpid="SNP",
                   eaf="FREQ1",
                   chrom="CHR",
                   pos="POS",
                   ea="A1",
                   nea="A2",
                   n=12121,
                   p="PVALUE",
                   beta="EFFECT1",
                   se="STDERR")

## creat qqplot and manhattan plot with just one line
myplot = AF.plot_mqq(
        snpid="MARKERNAME",
        mode="mqq",
        stratified=True,
        eaf="EAF",
        anno=True,
        cut=20,
        highlight=["rs7434417","rs12044963"], #the lead SNPs to highlight
        highlight_color="#33FFA0", 
        maf_bin_colors = ["#f0ad4e","#5cb85c", "#7878BA","#000042"])
```
Or you can plot it separately.

### Calculate genomic inflation factor
```
AF.get_gc()
```

### Extract top snps given a sliding window size

```
AF.get_lead()
```
Ref:
Zhou, Wei, and Global Biobank Meta-analysis Initiative. "Global Biobank Meta-analysis Initiative: Powering genetic discovery across human diseases." medRxiv (2021).

### Converting observed scale heritability to liability scale heritability
```
gl.h2_obs_to_liab(h2_obs, P, K)

gl.h2_obs_to_liab(h2_obs, P, K, se_obs=None)
```
Ref: 
Equation 23
Lee, Sang Hong, et al. "Estimating missing heritability for disease from genome-wide association studies." The American Journal of Human Genetics 88.3 (2011): 294-305.


### Read ldsc results in to pandas DataFrame

Directly read ldsc -h2 or -rg into pandas dataframe...

```
pathlist=["./test.results.log","./test2.results.log"]

ldsc_h2 = gl.read_ldsc(pathlist, mode="h2")
ldsc_rg = gl.read_ldsc(pathlist, mode="rg")

ldsc_h2
Filename	h2_obs	h2_se	Lambda_gc	Mean_chi2	Intercept	Intercept_se	Ratio	Ratio_se
test.results.log	42.9954	8.657	1.2899	1.3226	0.0098	0.0098	0.6538	0.0304
test2.results.log	NA	NA	1.2899	1.3226	0.0098	0.0098	Ratio < 0	NA

ldsc_rg
p1	p2	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se
./test.results.log	./test.results.log	0.2317	0.0897	2.5824	0.0098	0.3305	0.0571	0.9612	0.009	-0.0001	0.0062
./test.results.log	./test2.results.log	0.2317	0.0897	2.5824	0.0098	0.3305	0.0571	0.9612	0.009	-0.0001	0.0062

```

### Compare effect sizes of selected variants from two sumstats
```
gl.compare_effect()
```

### preformat your sumstats for a qc workflow



--------------------------
# Log
- 1.0.0 implemented Sumstats object

- 0.0.5 - 0.0.6
- added  compare_effect, read_ldsc 

- 0.0.4  
  -  added mqqplot feature
  -  fixed gtesig algorithm
  -  recreated mplot and qqplot


For more information: 
https://gwaslab.com/
