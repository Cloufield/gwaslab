# gwaslab
A collection of handy python scripts for GWAS.

1. Manhattan plot
2. QQ plot
3. Calculate lamda GC
4. Select top SNPs based on a given window size.
5. Convert beta/se <-> OR/95%L_U/95%L_L
6. Select hapmap3 SNPs from sumstats

Requirements:
1. Python>3
2. scipy
3. numpy
4. pandas
5. matplotlib
6. seaborn

Install:
```
pip install gwaslab
```
--

Usage:

```
import gwaslab as gl

gl.mplot(df,chr,pos,p,cut,cutfactor)

gl.qqplot()

gl.gc()

gl.getsig()
```
