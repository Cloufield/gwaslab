# gwaslab
A collection of handy python scripts for GWAS. 

Just want to make lif eaiser and save myself from repetitive work.

## What you can do with gwaslab:
1. [Side-by-side Manhattan and QQ plot](#create-manhattan-plot-and-qq-plot-with-just-one-line)
2. [Manhattan plot](#manhattan-plot)
3. [QQ plot](#qq-plot)
4. [Calculate lamda GC](#calculate-genomic-inflation-factor)
5. [Select top SNPs based on a given window size.]
6. Convert beta/se <-> OR/95%L_U/95%L_L
7. Select hapmap3 SNPs from sumstats and convert to ldsc format
8. [Convert Observed scale heritability to liability scale heritability](#converting-observed-scale-heritability-to-liability-scale-heritability)
9. [read ldsc log and extract numeric results directly into a pandas dataframe.](#read-ldsc-results-in-to-pandas-dataframe)
10. compare the effect size of select variants / or automatically detected lead variants from two sumstats.

![屏幕截图 2022-03-28 235029](https://user-images.githubusercontent.com/40289485/160425251-63549fcd-73b8-4328-b00d-d7b83e26d4e9.jpg)


## Requirements:
1. Python>3  
2. "scipy"  
3. "numpy"  
4. "pandas"  
5. "matplotlib"  
6. "seaborn" 
7. "adjustText"

## Install:
```
pip install gwaslab
```
Current version: 1.0.4

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
                   se="STDERR",
                   other=["WALDCHISQ"])

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
