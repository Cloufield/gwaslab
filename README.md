# gwaslab
A collection of handy python scripts for GWAS. 
This package is based on matplotlib and seaborn.
Just want to save myself from repetitive work.

## What you can do with gwaslab:
1. [Side-by-side Manhattan and QQ plot](#create-manhattan-plot-and-qq-plot-with-just-one-line)
2. [Manhattan plot](#manhattan-plot)
3. [QQ plot](#qq-plot)
4. [Calculate lamda GC](#calculate-genomic-inflation-factor)
5. [Select top SNPs based on a given window size.]
6. Convert beta/se <-> OR/95%L_U/95%L_L
7. Select hapmap3 SNPs from sumstats
8. [Convert Observed scale heritability to liability scale heritability](#converting-observed-scale-heritability-to-liability-scale-heritability)

![manhattan_qq_plot](https://user-images.githubusercontent.com/40289485/154832769-eddaf72e-9664-4f33-86e9-199e8fe92e56.png)

## Requirements:
1. Python>3  2. "scipy"  3. "numpy"  4. "pandas"  5. "matplotlib"  6. "seaborn"

## Install:
```
pip install gwaslab
```
Current version: 0.0.4

# Usage:

Input: pandas dataframe

### Create Manhattan plot and QQ plot with just one line
```
import gwaslab as gl

## creat qqplot and manhattan plot with just one line
## pass a dataframe in, and specify the column name for chromosome, base pair position, and also the p values.
gl.mqqplot(sumstats,"CHR","POS","PVALUE")

## adjust the plot, select top snps and add annotation sutomatically.
gl.mqqplot(sumstats,"CHR","POS","PVALUE",cut=20,cutfactor=10,anno=True,verbose=True,save=True,title="gwaslab")

## all options
gl.mqqplot(insumstats,
          chrom,
          pos,
          p,
          scaled=False,
          cut=0,
          cutfactor=10,
          cut_line_color="#ebebeb",
          windowsizekb=500,
          anno=None,
          sig_level=5e-8,
          sig_line_color="grey",
          suggestive_sig_level=5e-6,
          title =None,
          mtitle=None,
          qtitle=None,
          figsize =(15,5),
          fontsize = 10,
          colors = ["#000042", "#7878BA"],
          verbose=True,
          repel_force=0.03,
          gc=True,
          save=None,
          saveargs={"dpi":300,"facecolor":"white"}
          )
```
Or you can plot it separately.
### Manhattan plot
```
gl.mplot()
```
### QQ plot
```
gl.qqplot()
```

### Calculate genomic inflation factor
```
gc(insumstats{"PVALUE"},mode="p",level=0.5)
gc(insumstats["Z"],mode="z",level=0.5)
gc(insumstats["chi2"],mode="chi2",level=0.5)
```

### Extract top snps given a sliding window size
```
gl.getsig(insumstats,id,chrom,pos,p)

gl.getsig(insumstats,id,chrom,pos,p,windowsizekb=500,verbose=True,sig_level=5e-8)
```

### Converting observed scale heritability to liability scale heritability
```
gl.h2_obs_to_liab(h2_obs, P, K)

gl.h2_obs_to_liab(h2_obs, P, K, se_obs=None)
```
Ref: 

--------------------------
# Log
- 0.0.4  
  -  added mqqplot feature
  -  fixed gtesig algorithm
  -  recreated mplot and qqplot

# Next 
- beta to OR
- OR to beta 
- (Possibly intergrate ldsc's munge.py)

For more information: 
https://gwaslab.com/
