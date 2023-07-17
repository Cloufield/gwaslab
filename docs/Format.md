# Output sumstats in certain formats

GWASLab provides a flexible formatting and saving function.

## Usage

```
.to_format(
          path="./sumstats",
          fmt="ldsc",   
          ...
          )
```

|`.to_format()` options|DataType|Description|Default|
|-|-|-|-|
|`path`|`string`|the path for the output file; only prefix is needed.|`"./sumstats"`|
|`fmt`|`string`|output format for sumstats. Currently support `plink` ,`plink2`, `ldsc`, `saige`, `fastgwa`, `regenie` and so forth. For details , please check [ https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook).|`"gwaslab"`|
|`cols`|`list`|list of additional columns to linclude in the output|`None`|
|`extract`|`list`|a list of variant SNPIDs to include.|`None`|
|`exclude`|`list`|a list of variant SNPIDs to exclude.|`None`|
|`id_use`|`SNPID` or `rsID`| specify which ID to use when excluding or extracting variants.|`rsID`|
|`hapmap3`|`boolean`|If True, only output Hapmap3 SNPs.|`False`|
|`exclude_hla`|`boolean`|if True, exclude variants in the MHC region from the output.|`False`|
|`hla_range`|`tuple`|a tuple of 2 numbers (MBp) indicating the start and the end position of the HLA region. |`(25,34)`|
|`build`|`string`| reference genome build.|`"19"`|
|`xymt_number`|`boolean`|if True, output sex chromosomes as X/Y/MT |"False"|
|`xymt`|`list`| 3-element list of sex chromosome notations to indicate how to convert integers to sex chromosome|`["X","Y","MT"]`|
|`chr_prefix`|`string`|Add a prefix to chromosomes. For example, 6 -> Chr6.|`""`|
|`bgzip`|`boolean`|If True, bgzip the output file. Only works for bed and vcf format.|-|
|`tabix`|`boolean`|If True, use tabix to index the bgzipped output file. Only works for bed and vcf format.|-|
|`md5sum`|`boolean`|If True, calculate and outpit the file MD5 hashes|`False`|
|`to_csvargs`|`dict`|extra parameters for pd.to_csv()|`None`|
|`float_formats`|`dict`|a dictionary to specify the float format for each column.|`None`|
|`verbose`|`boolean`|If True, print logs.|`True`|
|`output_log`|`boolean`|If True, save the log to a file.|`True`|
|`ssfmeta`|`boolean`|If True, output a gwas-ssf-style meta file.|`False`|

## Format dictionary

Using `float_formats`, you can specify the formats for numbers.

##info "Default formats for floating-point numbers"
    ```
    {'EAF': '{:.4g}', 'BETA': '{:.4f}', 'Z': '{:.4f}','CHISQ': '{:.4f}','SE': '{:.4f}','OR': '{:.4f}','OR_95U': '{:.4f}','OR_95L': '{:.4f}','INFO': '{:.4f}','P': '{:.4e}','MLOG10P': '{:.4f}','DAF': '{:.4f}'}
    ```

## Example 1: common tabular format

GWASLab support commonly used tabular formats, which are listed in a companion repository `formatbook`.

!!! quote "formatbook"
    for more details, please check [formatbook](https://github.com/Cloufield/formatbook)
    ssf: GWAS-SSF format
    gwascatalog : GWAS Catalog format
    pgscatalog : PGS Catalog format
    plink: PLINK output format
    plink2: PLINK2 output format
    saige: SAIGE output format
    regenie: output format
    fastgwa: output format
    metal: output format
    mrmega: output format
    fuma: input format
    ldsc: input format
    locuszoom: input format
    vcf: gwas-vcf format
    bolt_lmm : output format


!!! example "output metal format"
    ```
    import gwaslab as gl
    # load your raw sumstats
    mysumstats = gl.Sumstats(...)
    # basic QC
    mysumstats.basic_check()
    # output metal format
    mysumstats.to_format("./test",fmt="metal")
    
    Tue Sep 13 18:00:41 2022 Start to format the output sumstats in:  metal  format
    Tue Sep 13 18:00:41 2022  -Formatting statistics ...
    Tue Sep 13 18:00:41 2022  - Float statistics formats:
    Tue Sep 13 18:00:41 2022   - Columns: ['EAF', 'BETA', 'SE', 'P']
    Tue Sep 13 18:00:41 2022   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
    Tue Sep 13 18:00:41 2022  - Start outputting sumstats in metal format...
    Tue Sep 13 18:00:41 2022  -metal format will be loaded...
    Tue Sep 13 18:00:41 2022  -metal format meta info:
    Tue Sep 13 18:00:41 2022   - format_name  :  metal
    Tue Sep 13 18:00:41 2022   - format_source  :  https://genome.sph.umich.edu/wiki/METAL_Documentation
    Tue Sep 13 18:00:41 2022   - format_version  :  20220726
    Tue Sep 13 18:00:41 2022  -gwaslab to metal format dictionary:
    Tue Sep 13 18:00:41 2022   - gwaslab keys: ['SNPID', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'DIRECTION']
    Tue Sep 13 18:00:41 2022   - metal values: ['MarkerName', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr', 'P-value', 'Direction']
    Tue Sep 13 18:00:41 2022  -Output columns: Index(['MarkerName', 'Allele1', 'Allele2', 'Freq1', 'Effect', 'StdErr',
           'P-value'],
          dtype='object')
    Tue Sep 13 18:00:41 2022  -Output path: ./test.metal.tsv.gz
    Tue Sep 13 18:00:41 2022  -Saving log file: ./test.metal.log
    Tue Sep 13 18:00:41 2022 Finished outputting successfully!
    ```

## Example 2 : LDSC
!!! example "LDSC format, extract hapmap3 SNPs and exclude SNPs in HLA region"
    ```
    ## format the sumstats to ldsc format
    ## extract only hapmap3 SNPs
    ## exclude SNPs in HLA region
    mysumstats.to_format("./test",fmt="ldsc", hapmap3=True, exclude_hla=False, build="19")
    ```

## Example 3 : BED
!!! example bed-like format
    ```
    # output 1-based bed-like files for vep
    mysumstats.to_format("./test",fmt="vep",xymt_number=True,chr_prefix="Chr")
    
    # output 0-based bed-like file, and then bgzip and index the file.
    mysumstats.to_format("./test",fmt="bed",bgzip=True,tabix=True)
    ```

## Example 4: VCF 
!!! example vcf format
    ```
    # output vcf file, and then bgzip and index the file.
    mysumstats.to_format("./test",fmt="vcf",bgzip=True,tabix=True)
    ```

## Example 5: GWAS-ssf
!!! example GWAS-ssf format
    ```
    # output  GWAS-ssf format
    mysumstats.to_format("./test",fmt="ssf")
    ```

