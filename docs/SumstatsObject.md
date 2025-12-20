# Sumstats Object in GWASLab

In GWASLab, sumstats were stored in a `Sumstats Object`, which is built on `pandas Dataframe`. All other functions are designed as methods of this Sumstats Object. 

To load any sumstats into the object, simply specify the column name and load the raw GWAS summary statistics from a pandas DataFrame or specify a file path. All raw data will be loaded as `"string"` datatype. 

## gl.Sumstats()

```py
mysumstats = gl.Sumstats(
             sumstats,
             fmt=None,
             snpid=None,
             rsid=None,
             chrom=None,
             pos=None,
             ea=None,
             nea=None,
             ...
)
```

## Options

`sumstats`: either a file path `string` or a pandas `DataFrame`

Currently, GWASLab supports the following columns:

| Option  | DataType | Description                                                  | Header in GWASLab |
|---------|----------|--------------------------------------------------------------|-------------------|
| `snpid` | `string` | variant ID column name, preferably in chr:pos:ea:nea format. | `SNPID`           |
| `rsid`  | `string` | dbSNP rsID column name                                       | `rsID`            |

The minimum required columns are just either `rsid `or `snpid`. 
All other columns and options are optional.

| Option      | DataType              | Description                                                                                                                                                                                         | Header in GWASLab                |
|-------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------|
| `fmt`       | `string`              | input sumstats format. For formats supported by GWASLab, please check [https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook)                                            | -                                |
| `chrom`     | `string`              | chromosome column name                                                                                                                                                                              | `CHR`                            |
| `pos`       | `string`              | basepair position column name                                                                                                                                                                       | `POS`                            |
| `ea`        | `string`              | effect allele column name; BETA, OR, HR, EAF are in reference to EA...                                                                                                                              | `EA`                             |
| `nea`       | `string`              | non-effect allele column name                                                                                                                                                                       | `NEA`                            |
| `ref`       | `string`              | reference allele column name; the allele on reference genome                                                                                                                                        | `REF`                            |
| `alt`       | `string`              | alternative allele column name; the allele that is not on reference genome; when `ea`,`ref` and `alt` are specified, `nea` will be inferred.                                                        | `ALT`                            |
| `eaf`       | `string`              | effect allele frequency                                                                                                                                                                             | `EAF`                            |
| `neaf`      | `string`              | non-effect allele frequency. NEAF will be converted to EAF (EAF = 1 - NEAF) while loading.                                                                                                          | `EAF`                            |
| `n`         | `string` or `integer` | sample size column name or just input a single `integer` as sample size for all variants                                                                                                            | `N`                              |
| `beta`      | `string`              | effect size beta column name; in reference to EA                                                                                                                                                    | `BETA`                           |
| `se`        | `string`              | standard error column name                                                                                                                                                                          | `SE`                             |
| `chisq`     | `string`              | chi square column name                                                                                                                                                                              | `CHISQ`                          |
| `z`         | `string`              | z score column name                                                                                                                                                                                 | `Z`                              |
| `p`         | `string`              | p value column name                                                                                                                                                                                 | `P`                              |
| `mlog10p`   | `string`              | -log10(P) column name                                                                                                                                                                               | `MLOG10P`                        |
| `info`      | `string`              | imputation info or rsq column name                                                                                                                                                                  | `INFO`                           |
| `OR`        | `string`              | odds ratio column name; in reference to EA                                                                                                                                                          | `OR`                             |
| `OR_95L`    | `string`              | odds ratio lower 95% CI column name                                                                                                                                                                 | `OR_95L`                         |
| `OR_95U`    | `string`              | odds ratio upper 95% CI column name                                                                                                                                                                 | `OR_95U`                         |
| `direction` | `string`              | direction column name. GWASLab uses METAL format (e.g. `++--+?+0`)                                                                                                                                  | `DIRECTION`                      |
| `other`     | `list`                | a list  of other column names you want to keep with the core columns (probably some annotations).                                                                                                   | -                                |
| `ncontrol`  | `string` or `integer` | sample size column name for controls or just input a single `integer` as sample size for all variants                                                                                               | `N_CONTROL`                      |
| `ncase`     | `string` or `integer` | sample size column name for cases or just input a single `integer` as sample size for all variants                                                                                                  | `N_CASE`                         |
| `HR`        | `string`              | hazrad ratio column name                                                                                                                                                                            | `HR`                             |
| `HR_95U`    | `string`              | hazrad ratio upper 95% CI column name                                                                                                                                                               | `HR_95U`                         |
| `HR_95L`    | `string`              | hazrad ratio lower 95% CI column name                                                                                                                                                               | `HR_95L`                         |
| `beta_95U`  | `string`              | beta upper 95% CI column name                                                                                                                                                                       | `BETA_95U`                       |
| `beta_95L`  | `string`              | beta lower 95% CI column name                                                                                                                                                                       | `BETA_95L`                       |
| `i2`        | `string`              | I2 column name                                                                                                                                                                                      | `I2_HET`                         |
| `phet`      | `string`              | heterogeneity test P value                                                                                                                                                                          | `P_HET`                          |
| `dof`       | `string` or `integer` | degree of freedom                                                                                                                                                                                   | `DOF`                            |
| `snpr2`     | `string`              | column name for proportion of phenotypic variance explained by each variant                                                                                                                         | `SNPR2`                          |
| `maf`       | `string`              | minor allele frequency column header                                                                                                                                                                | `MAF`                            |
| `f`         | `string`              | F statistics column header                                                                                                                                                                          | `F`                              |
| `status`    | `string`              | status code column name. GWASLab uses a 7-digit variant status code. For details, please check status code page.                                                                                    | `STATUS`                         |
| `verbose`   | `boolean`             | if True, print log.                                                                                                                                                                                 | -                                |
| `build`     | `string`              | genome build. `19` for hg19, `38` for hg38 and `99` for unknown.                                                                                                                                    | The first two digits of `STATUS` |
| `**arg `    | `string`              | additional parameters for [pd.read_table()](https://pandas.pydata.org/docs/reference/api/pandas.read_table.html) function. Some common options include : `sep`,`nrows`, `skiprows` and `na_values`. | -                                |

## Loading sumstats

### Load by specifying columns

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

### Load by specifying formats

GWASLab supports common sumstats formats and you can load sumstats by specifying `fmt`.  
    
GWASLab uses a manually curated format conversion dictionary in [https://github.com/Cloufield/formatbook](https:/github.com/Cloufield/formatbook). Currently, it supports the following formats:

| `Keyword`     | `Description`        |
|---------------|----------------------|
| `ssf`         | GWAS-SSF             |
| `gwascatalog` | GWAS Catalog format  |
| `pgscatalog`  | PGS Catalog format   |
| `plink`       | PLINK output format  |
| `plink2`      | PLINK2 output format |
| `saige`       | SAIGE output format  |
| `regenie`     | output format        |
| `fastgwa`     | output format        |
| `metal`       | output format        |
| `mrmega`      | output format        |
| `fuma`        | input format         |
| `ldsc`        | input format         |
| `locuszoom`   | input format         |
| `vcf`         | gwas-vcf format      |
| `bolt_lmm`    | output format        |

!!! info "Update Formatbook using `gl.update_formaybook()`"
    ```
    gl.update_formatbook()
    Mon Jul 17 17:38:11 2023 Updating formatbook from: https://raw.github.com/Cloufield/formatbook/main/formatbook.json
    Mon Jul 17 17:38:12 2023 Overwrite formatbook to :  /home/yunye/gwaslab/gwaslab/src/gwaslab/data/formatbook.json
    Mon Jul 17 17:38:12 2023 Available formats: auto,bolt_lmm,fastgwa,gwascatalog,gwascatalog_hm,gwaslab,ldsc,metal,mrmega,mtag,pgscatalog,pgscatalog_hm,pheweb,plink,plink2,regenie,saige,ssf,template,vcf
    Mon Jul 17 17:38:12 2023 Formatbook has been updated!
    ```

!!! example "Load sumstats by simply specifying the format"
    ```
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz", fmt="saige")
    ```

### Load sumstats with `auto` mode

GWASLab also provides an auto mode (`fmt="auto"`; available since v3.4.21) which assumes A1 or alternative allele (ALT) is the effect allele (EA) and Frq refers to the allele frequency of effect allele (EAF). Common headers will be detected. You can find the conversion table [here](https://github.com/Cloufield/formatbook/blob/main/formats/auto.json)

!!! example "Load sumstats with auto mode"
    ```
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz", fmt="auto")
    ```
    
### Load sumstats from chromosome-separated files
GWASLab supports loading sumstats from chromosome-separated files (file names need to be in the same pattern.). Just use @ to replace the chromosome numbers. 

!!! example
    ```
    mysumstats = gl.Sumstats("t2d.chr@.txt.gz",fmt="metal")
    ```


## Check and save sumstats
After loading, the raw data columns will be renamed to new columns without ambiguity and the DataFrame is stored in `.data` :

!!! example
    ```python
    mysumstats.data
    ```

You can simply save the processed data using pandas saving functions, for example:

!!! example
    ```
    mysumstats.data.to_csv("./mysumstats.csv")
    ```  
    
or convert the sumstats to other sumstats using GWASLab `to_format()` function (recommended):

!!! example 
    ```
    mysumstats.to_format("./mysumstats", fmt="ldsc",hapmap3=True, exclude_hla=True, build="19")
    ```
    
Please check [GWASLab - Format](https://cloufield.github.io/gwaslab/Format/) for more details.

## Saving half-finished Sumstats Object

If the pipeline is very long, and you need to temporarily save the Sumstats Object, you can use the `.to_pickle()` method (recommended) or the `gl.dump_pickle()` function to temporarily save the Sumstats Object.

```python
# Recommended method (object method)
mysumstats.to_pickle("./mysumstats.pickle", overwrite=True)

# Alternative method (standalone function)
gl.dump_pickle(mysumstats, "./mysumstats.pickle", overwrite=True)
```

Please check [GWASLab - Pickle](https://cloufield.github.io/gwaslab/Pickle/) for more details.

## Logging

All manipulation conducted to the sumstats will be logged for reproducibility and traceability. 

The log is stored in a `gl.Log()` object. You can check it by `.log.show() `and save it using `.log.save()`

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
