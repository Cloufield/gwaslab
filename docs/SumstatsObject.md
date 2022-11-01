# Sumstats Object in gwaslab

In gwaslab, sumstats were stored in a `Sumstats Object`，which is built on `pandas Dataframe`. All other function are designed as methods of this Sumstats Object. 

To load any sumstats into the object, simply specify the column name and load the raw GWAS summary statsitics from a pandas dataframe or specifying file path. All raw data will be loaded as "string" datatype. 

```python
mysumstats = gl.Sumstats(
             sumstats,
             fmt=None,
             snpid=None,
             rsid=None,
             chrom=None,
             pos=None,
             ea=None,
             nea=None,
             ref=None,
             alt=None,
             eaf=None,
             neaf=None,
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

`sumstats`: either a file path `string` or a pandas `DataFrame`

Currently, gwaslab supports the following columns:

* `snpid `: variant ID column name, preferably in chr:pos:ea:nea format.
* `rsid `: dbSNP rsID column name

The minimum required columns are just either `rsid `or `snpid`. 
All other columns are optional.

* `fmt`: input sumstats format : For formats supported by gwaslab, please check [https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook)
* `chrom `: chromosome column name
* `pos`: basepair position column name
* `ea`: effect allele column name. 
* `nea`: non-effect allele column name.
* `ref`: reference allele column name.
* `alt`: alternative allele column name , when `ea`,`ref` and `alt` are specified, `nea` will be inferred.
* `eaf`: effect allele frequency
* `neaf`: non-effect allele frequency, neaf will be converted to eaf (eaf = 1 - neaf) while loading.
* `n`: sample size column name or just input a single  `integer` 
* `beta`: effect size beta column name
* `se`: standard error column name
* `chisq`: chi square column name
* `z`: z score column name
* `p`: p value column name
* `mlog10p`: -log10(P) column name
* `info`: imputation info or rsq column name
* `OR`: odds ratio column name
* `OR_95L`:odds ratio lower 95% ci column name
* `OR_95U`:odds ratio upper 95% ci column name
* `direction`: direction column name. gwaslab uses METAL format (e.g. "++--+?+")
* `other`: a list  of other column names you want to keep with the core columns, probably some annotations.
* `status`: gwaslab 7-digit vairants status code. For details, please check status code page.
* `verbose`: if true: output log 
* `build`:  `str `genome build ("19","38")
* `**arg `: additional parameters for [pd.read_table()](https://pandas.pydata.org/docs/reference/api/pandas.read_table.html) function. Some common options include : `sep`,`nrows`, `skiprows` and `na_values`.

## Loading sumstats
you can load the path by specifying the columns like:
```
mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
             snpid="SNPID",
             chrom="CHR",
             pos="POS",
             ea="Allele2",
             nea="Allele1",
             eaf="AF_Allele2",
             beta="BETA",
             se="SE",
             p="p.value",
             n="N")
```
or just specify the format:

```
mysumstats = gl.Sumstats("t2d_bbj.txt.gz", format="saige")
```
gwaslab uses manully curated format conversion dictionary in [https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook)

Supported formats:

* `ssf`: GWAS-SSF
* `gwascatalog` : GWAS Catalog format
* `pgscatalog` : PGS Catalog format
* `plink`: PLINK output format
* `plink2`: PLINK2 output format
* `saige`: SAIGE output format
* `regenie`: output format
* `fastgwa`: output format
* `metal`: output format
* `mrmega`: output format
* `fuma`: input format
* `ldsc`: input format
* `locuszoom`: input format
* `vcf`: gwas-vcf format

### Loading sumstats from chromosome-separated files
gwaslab support loading sumstats from chromosome-separated files (file names need to be in the same pattern.). Just use @ to replace chromosome number. 

For example:
```
mysumstats = gl.Sumstats("t2d.chr@.txt.gz",fmt="metal")
```


## Check and save sumstats
After loading, the raw data columns will be renamed to new columns without ambiguity and the dataframe is store in .data :
```python
mysumstats.data
```

You can simply save the processed data using pandas saving functions, for example:
```
mysumstats.data.to_csv("./mysumstats.csv")
```  
or convert the sumstats to other sumstats using gwaslab `to_format()` function (**stringly recommended**):

```
mysumstats.to_format("./mysumstats", fmt="ldsc",hapmap3=True, exclude_hla=True, build="19")
```
Please check [https://cloufield.github.io/gwaslab/Format/](https://cloufield.github.io/gwaslab/Format/)


## Logging
All manipulation conducted to the sumstats will be logged for reproducibility and traceability. The log is stored in a gl.Log object . You can check it by` .log.show() `and save it using `.log.save()`

```
mysumstats.log.show()

mysumstats.log.save()
```

## Sumstats summary
You can check the meta information of this sumstats by:

```python
mysumstats.summary()
```

## Other functions
Other functions of gwaslab is implemented as the methods of Sumstats Object.

```python
mysumstats.basic_check()

mysumstats.plot_mqq()

...
```
