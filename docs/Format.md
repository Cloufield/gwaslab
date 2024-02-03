# Output sumstats in certain formats

GWASLab provides a flexible formatting and saving function.

## .to_format()

```
.to_format(
          path="./sumstats",
          fmt="ldsc",   
          ...
          )
```

## Options

| `.to_format()` options | DataType          | Description                                                                                                                                                                                                                          | Default          |
|------------------------|-------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
| `path`                 | `string`          | the path for the output file; only prefix is needed.                                                                                                                                                                                 | `"./sumstats"`   |
| `fmt`                  | `string`          | output format for sumstats. Currently support `plink` ,`plink2`, `ldsc`, `saige`, `fastgwa`, `regenie` and so forth. For details , please check [ https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook). | `"gwaslab"`      |
| `cols`                 | `list`            | list of additional columns to linclude in the output                                                                                                                                                                                 | `None`           |
| `extract`              | `list`            | a list of variant SNPIDs to include.                                                                                                                                                                                                 | `None`           |
| `exclude`              | `list`            | a list of variant SNPIDs to exclude.                                                                                                                                                                                                 | `None`           |
| `id_use`               | `SNPID` or `rsID` | specify which ID to use when excluding or extracting variants.                                                                                                                                                                       | `rsID`           |
| `hapmap3`              | `boolean`         | If True, only output Hapmap3 SNPs.                                                                                                                                                                                                   | `False`          |
| `exclude_hla`          | `boolean`         | if True, exclude variants in the MHC region from the output.                                                                                                                                                                         | `False`          |
| `hla_range`            | `tuple`           | a tuple of 2 numbers (MBp) indicating the start and the end position of the HLA region.                                                                                                                                              | `(25,34)`        |
| `build`                | `string`          | reference genome build.                                                                                                                                                                                                              | `"19"`           |
| `xymt_number`          | `boolean`         | if True, output sex chromosomes as X/Y/MT                                                                                                                                                                                            | "False"          |
| `xymt`                 | `list`            | 3-element list of sex chromosome notations to indicate how to convert integers to sex chromosome                                                                                                                                     | `["X","Y","MT"]` |
| `chr_prefix`           | `string`          | Add a prefix to chromosomes. For example, 6 -> Chr6.                                                                                                                                                                                 | `""`             |
| `bgzip`                | `boolean`         | If True, bgzip the output file. Only works for bed and vcf format.                                                                                                                                                                   | -                |
| `tabix`                | `boolean`         | If True, use tabix to index the bgzipped output file. Only works for bed and vcf format.                                                                                                                                             | -                |
| `tabix_indexargs`      | `dict`            | extra parameters for pysam.tabix_index()                                                                                                                                                                                             | `{}`             |
| `md5sum`               | `boolean`         | If True, calculate and output the file MD5 hashes                                                                                                                                                                                    | `False`          |
| `to_csvargs`           | `dict`            | extra parameters for pd.to_csv()                                                                                                                                                                                                     | `None`           |
| `float_formats`        | `dict`            | a dictionary to specify the float format for each column.                                                                                                                                                                            | `None`           |
| `verbose`              | `boolean`         | If True, print logs.                                                                                                                                                                                                                 | `True`           |
| `output_log`           | `boolean`         | If True, save the log to a file.                                                                                                                                                                                                     | `True`           |
| `ssfmeta`              | `boolean`         | If True, output a gwas-ssf-style meta file.                                                                                                                                                                                          | `False`          |

## Format dictionary

Using `float_formats`, you can specify the formats for numbers.

!!! info "Default formats for floating-point numbers"
    ```
    {'EAF': '{:.4g}', 'BETA': '{:.4f}', 'Z': '{:.4f}','CHISQ': '{:.4f}','SE': '{:.4f}','OR': '{:.4f}','OR_95U': '{:.4f}','OR_95L': '{:.4f}','INFO': '{:.4f}','P': '{:.4e}','MLOG10P': '{:.4f}','DAF': '{:.4f}'}
    ```

## Examples

GWASLab supports commonly used tabular formats, which are listed in a companion repository `formatbook`.

!!! quote "formatbook"
    For more details, please check [formatbook](https://github.com/Cloufield/formatbook)

!!! example
    See [Output sumstats](https://cloufield.github.io/gwaslab/format_load_save/)
