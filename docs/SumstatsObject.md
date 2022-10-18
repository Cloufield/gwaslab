# Sumstats Object in gwaslab

In gwaslab, sumstats were stored in a `Sumstats Object`，which is built on `pandas Dataframe`. All other function are designed as methods of this Sumstats Object. 

To load any sumstats into the object, simply specify the column name and load the raw GWAS summary statsitics from a pandas dataframe or specifying file path. All raw data will be loaded as "string" datatype. 

```python
mysumstats = gl.Sumstats(
             sumstats,
             snpid=None,
             rsid=None,
             chrom=None,
             pos=None,
             ea=None,
             nea=None,
             eaf=None,
             n=None,
             beta=None,
             se=None,
             chisq=None,
             z=None,
             p=None,
             mlog10p=None,
             info=None,
             OR=None,
             OR_95L=None,
             OR_95U=None,
             status=None,
             other=[],
             direction=None,
             verbose=True,
             build="00",
             **args
)
```

`sumstats`: either a file path or a pandas DataFrame

Currently, gwaslab supports the following columns:
- `snpid `: variant ID column name, preferably in chr:pos:ea:nea format.
- `rsid `: dbSNP rsID column name

The minimum required columns are just either `rsid `or `snpid`. 
All other columns are optional.
- `chrom `: chromosome column name
- `pos`: basepair position column name
- `ea`: effect allele column name 
- `nea`: non-effect allele column name
- `eaf`: effect allele frequency
- `n`: sample size column name or just input a single  `integer` 
- `beta`: effect size beta column name
- `se`: standard error column name
- `chisq`: chi square column name
- `z`: z score column name
- `p`: p value column name
- `mlog10p`: -log10(P) column name
- `info`: imputation info or rsq column name
- `OR`: odds ratio column name
- `OR_95L`:odds ratio lower 95% ci column name
- `OR_95U`:odds ratio upper 95% ci column name
- `direction`: direction column name in METAL format (e.g. "++--+?+")
- `other`: a list  of other column names you want to keep with the core columns, probably some annotations.
- `status`: gwaslab 5-digit vairants status code. For details, please check status code page.
- `verbose`: if true: output log 
- `build`:  `str `genome build ("19","38")
- `**arg `: additional parameters for pl.read_table function. 

After loading, the raw data columns will be renamed to new columns without ambiguity and the dataframe is store in .data :

```python
mysumstats.data
```

You can simply save the processed data using pandas saving functions, for example
```
mysumstats.data.to_csv("./mysumstats.csv")
```  
or convert the sumstats to other sumstats:
please check [https://cloufield.github.io/gwaslab/Format/](https://cloufield.github.io/gwaslab/Format/)


All manipulation conducted to the sumstats will be logged for reproducibility and traceability. The log is stored in a gl.Log object . You can check it by` .log.show() `and save it using `.log.save()`

```
mysumstats.log.show()

mysumstats.log.save()
```

You can check the meta information of this sumstats by:

```python
mysumstats.summary()
```

Other functions of gwaslab is implemented as the methods of Sumstats Object.

```python
mysumstats.basic_check()

mysumstats.plot_mqq()
```
