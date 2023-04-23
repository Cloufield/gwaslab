# Sumstats Object in GWASLab

In GWASLab, sumstats were stored in a `Sumstats Object`，which is built on `pandas Dataframe`. All other function are designed as methods of this Sumstats Object. 

To load any sumstats into the object, simply specify the column name and load the raw GWAS summary statsitics from a pandas dataframe or specifying file path. All raw data will be loaded as `"string"` datatype. 

## Usage

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
             build="99",
             **args
)
```

## Parameters

`sumstats`: either a file path `string` or a pandas `DataFrame`

Currently, GWASLab supports the following columns:

* `snpid `: variant ID column name, preferably in chr:pos:ea:nea format.
* `rsid `: dbSNP rsID column name

The minimum required columns are just either `rsid `or `snpid`. 
All other columns and options are optional.

* `fmt`: `string`, input sumstats format. For formats supported by GWASLab, please check [https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook)
* `chrom `: `string`, chromosome column name
* `pos`: `string`, basepair position column name
* `ea`: `string`, effect allele column name. 
* `nea`: `string`, non-effect allele column name.
* `ref`: `string`, reference allele column name.
* `alt`: `string`, alternative allele column name , when `ea`,`ref` and `alt` are specified, `nea` will be inferred.
* `eaf`: `string`, effect allele frequency
* `neaf`: `string`, non-effect allele frequency. NEAF will be converted to EAF (EAF = 1 - NEAF) while loading.
* `n` `string` or `integer`, sample size column name or just input a single `integer` as sample size for all variants.
* `beta`: `string`, effect size beta column name
* `se`: `string`, standard error column name
* `chisq`: `string`, chi square column name
* `z`: `string`, z score column name
* `p`: `string`, p value column name
* `mlog10p`: `string`, -log10(P) column name
* `info`: `string`, imputation info or rsq column name
* `OR`: `string`, odds ratio column name
* `OR_95L`: `string`, odds ratio lower 95% ci column name
* `OR_95U`:`string`, odds ratio upper 95% ci column name
* `direction`: `string`, direction column name. GWASLab uses METAL format (e.g. "++--+?+")
* `other`: `list`, a list  of other column names you want to keep with the core columns (probably some annotations).
* `status`: `string`, status code column name. GWASLab uses a 7-digit vairant status code. For details, please check status code page.
* `verbose`: `boolean`, if True, print log. 
* `build`:  `string `, genome build. `19` for hg19, `38` for hg38 and `99` for unknown.
* `**arg `: additional parameters for [pd.read_table()](https://pandas.pydata.org/docs/reference/api/pandas.read_table.html) function. Some common options include : `sep`,`nrows`, `skiprows` and `na_values`.

## Loading sumstats

You can load the sumstats by specifying the columns like:

!!! example "Load sumstats by manually specifying columns"
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

or just specify the sumstats format:

!!! example "Load sumstats by simply specifying the format"
    ```
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz", fmt="saige")
    ```
    
!!! note formatbook
    GWASLab uses manually curated format conversion dictionary in [https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook)

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
    * `bolt-lmm`: output format
    
### Loading sumstats from chromosome-separated files
GWASLab support loading sumstats from chromosome-separated files (file names need to be in the same pattern.). Just use @ to replace chromosome number. 

!!! example
    ```
    mysumstats = gl.Sumstats("t2d.chr@.txt.gz",fmt="metal")
    ```


## Check and save sumstats
After loading, the raw data columns will be renamed to new columns without ambiguity and the dataframe is store in `.data` :

!!! example
    ```python
    mysumstats.data
    ```

You can simply save the processed data using pandas saving functions, for example:

!!! example
    ```
    mysumstats.data.to_csv("./mysumstats.csv")
    ```  
    
or convert the sumstats to other sumstats using GWASLab `to_format()` function (**strongly recommended**):

!!! example 
    ```
    mysumstats.to_format("./mysumstats", fmt="ldsc",hapmap3=True, exclude_hla=True, build="19")
    ```
    
Please check [GWASLab - Format](https://cloufield.github.io/gwaslab/Format/) for more details.

## Saving half-finished Sumstats Object

If the pipeline is very long, and you need to temporarily save the Sumstats Object, you can use the `.dump_pickle()` method to temporarily save the Sumstats Object.

Please check [GWASLab - Pickle](https://cloufield.github.io/gwaslab/Pickle/) for more details.

## Logging

All manipulation conducted to the sumstats will be logged for reproducibility and traceability. 

The log is stored in a `gl.Log()` object . You can check it by `.log.show() `and save it using `.log.save()`

!!! example 
    ```
    mysumstats.log.show()
    
    mysumstats.log.save()
    ```

## Sumstats summary

You can check the meta information of sumstats by:

!!! example 
    ```python
    mysumstats.summary()
    ```

## Other functions

Other functions of GWASLab are implemented as the methods of Sumstats Object.

!!! example
    ```python
    mysumstats.basic_check()
    
    mysumstats.plot_mqq()
    
    ...
    ```
