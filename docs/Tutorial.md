# Quick Start

Using a jupyter notebook, we first import gwaslab package:

```python
import gwaslab as gl
```

The sample sumstats we use in this study: 

```bash
!wget -O t2d_bbj.txt.gz http://jenger.riken.jp/14/
```

Let's import this raw sumstats into the gwaslab Sumstats Object by specifying the necessary columns, and all data are imported as strings.
Note: you can either specify eaf (effect allele frequency) or neaf(non-effect allele frequency), if neaf is specified, it will be converted to eaf when loading sumstats.
```python
mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             neaf="Frq",
             beta="BETA",
             se="SE",
             p="P",
             direction="Dir",
             n="N",
             build="19")
```
See details in [SumstatsObject](SumstatsObject.md).

Maybe the first thing you want to check is the manhattan plot, you can do this with one line of code, gwaslab will perform a minimum QC for just the plotting.

```python
mysumstats.plot_mqq()
```

![tutorial_mqq](images/tutorial_mqq.jpg)

# Sanity check

Looks good, but we need to perform QC to make sure there are no unexpected errors, let's check the statitsics first.

```
mysumstats.check_sanity()
```

## Filtering

There are more than 10 million variants in the original sumstats and it will take long to process the entrie dataset. 

So, let's just filter-in (include) variants with P<0.00005 and filter-out (exclude) variants on ChrX.  This could also be used for filtering other columns like INFO,N and so forth if you need. 

```python
mysumstats.filter_out(gt={"P":0.005},eq={"CHR":"X"})
```

See details in [QC&Filtering](QC&Filtering.md).

## Standardization & normalization
It is also needed to check ID,CHR,POS and alleles:

simply run:

```python
mysumstats.basic_check()
```

`.basic_check()` is a wrapper of all the following basic functions, you can use these separately.

```python
mysumstats.fix_ID()
mysumstats.fix_chr()
mysumstats.fix_pos()
mysumstats.fix_allele()
mysumstats.check_sanity()
mysumstats.normalize_allele()
```
See details in [Standardization](Standardization.md).
## 



## Extract lead variants

Let's extract the lead variants in each significant loci to check our data. 

The significant loci are detected based on a sliding window (default window size: 500kb)

```python
mysumstats.get_lead()
```

See details in [ExtractLead](ExtractLead.md).


## Customized manhattan plot

GWASlab can plot more complicated manhattan plot: (not finished yet)

```python
mysumstats.plot_mqq(snpid="SNPID",mode="mqq",
                  cut=20,skip=3, eaf="EAF",
                  anno=True,anno_set=["9:22132729_A_G","6:20688121_T_A","9:22132729_A_G","15:62394264_G_C"] ,
                  pinpoint=["9:22132729_A_G","5:176513896_C_A"], 
                  highlight=["7:127253550_C_T","19:46166604_C_T"],
                  highlight_windowkb =1000,
                  stratified=True,
                  marker_size=(5,10),
                  figargs={"figsize":(15,5),"dpi":300})
```

![mqq_customized](images/mqq_customized.jpg)

See details in [Visualization](Visualization.md).

## 

## Harmonise the sumstats

Coming soon.

-------------

## Liftover

```python
mysumstats.liftover(n_cores=1,from_build="19", to_build="38")
```

Gwaslab only liftover CHR and POS, and when lifted, the last two digits status code will be rolled back to 99. Since for difference reference genome, the reference allele or strand might be reverse, so it is need to align and check agin. 

See details in [Harmonization](Harmonization.md).
