# gwaslab
Current version: 0.0.4
A collection of handy python scripts for GWAS. This package is based on matplotlib and seaborn.

What you can do with gwaslab:
1. Manhattan plot
2. QQ plot
3. Calculate lamda GC
4. Select top SNPs based on a given window size.
5. Convert beta/se <-> OR/95%L_U/95%L_L
6. Select hapmap3 SNPs from sumstats
7. Convert Observed scale heritability to liability scale heritability 

![manhattan_qq_plot](https://user-images.githubusercontent.com/40289485/154832769-eddaf72e-9664-4f33-86e9-199e8fe92e56.png)

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

gl.mqqplot(sumstats,"CHR","POS","PVALUE",cut=20,cutfactor=10,anno=True,verbose=True,save=True,title="gwaslab")

gl.mplot()

gl.qqplot()

gl.gc()

gl.getsig()
```

--------------------------
log
- 0.0.4  
  -  added mqqplot feature
  -  fixed gtesig algorithm
  -  recreated mplot and qqplot

For more information: 
https://gwaslab.com/
