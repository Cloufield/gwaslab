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
| `fmt`                  | `string`          | output format for sumstats. Currently support `plink`, `plink2`, `ldsc`, `saige`, `fastgwa`, `regenie` and so forth. For details, please check [https://github.com/Cloufield/formatbook](https://github.com/Cloufield/formatbook). | `"gwaslab"`      |
| `tab_fmt`              | `string`          | tabular format type when `fmt` is not 'vcf', 'bed', or 'annovar'. Options: `tsv`, `csv`, `parquet`                                                                                                                                   | `"tsv"`          |
| `cols`                 | `list`            | list of additional columns to include in the output                                                                                                                                                                                 | `None`           |
| `extract`              | `list`            | a list of variant SNPIDs to include.                                                                                                                                                                                                 | `None`           |
| `exclude`              | `list`            | a list of variant SNPIDs to exclude.                                                                                                                                                                                                 | `None`           |
| `id_use`               | `SNPID` or `rsID` | specify which ID to use when excluding or extracting variants.                                                                                                                                                                       | `rsID`           |
| `hapmap3`              | `boolean`         | If True, only output Hapmap3 SNPs.                                                                                                                                                                                                   | `False`          |
| `exclude_hla`          | `boolean`         | If True, exclude variants in the MHC region from the output.                                                                                                                                                                         | `False`          |
| `hla_range`            | `tuple`           | a tuple of 2 numbers (Mbp) indicating the start and the end position of the HLA region.                                                                                                                                              | `(25,34)`        |
| `build`                | `string`          | reference genome build.                                                                                                                                                                                                              | `None`           |
| `n`                    | `float`           | sample size to add as 'N' column.                                                                                                                                                                                                    | `None`           |
| `no_status`            | `boolean`         | If True, exclude 'STATUS' column from output.                                                                                                                                                                                       | `False`          |
| `xymt_number`          | `boolean`         | If True, output sex chromosomes as numeric codes (23, 24, 25 for X, Y, MT).                                                                                                                                                         | `False`          |
| `xymt`                 | `list`            | 3-element list of sex chromosome notations to indicate how to convert integers to sex chromosome                                                                                                                                     | `["X","Y","MT"]` |
| `chr_prefix`           | `string`          | Add a prefix to chromosomes. For example, 6 -> chr6.                                                                                                                                                                                 | `""`             |
| `gzip`                 | `boolean`         | If True, gzip compress the output file.                                                                                                                                                                                              | `True`           |
| `bgzip`                | `boolean`         | If True, bgzip the output file. Only works for bed and vcf format.                                                                                                                                                                   | `False`          |
| `tabix`                | `boolean`         | If True, use tabix to index the bgzipped output file. Only works for bed and vcf format. Requires `bgzip=True`.                                                                                                                       | `False`          |
| `tabix_indexargs`      | `dict`            | extra parameters for pysam.tabix_index()                                                                                                                                                                                             | `{}`             |
| `md5sum`               | `boolean`         | If True, calculate and output the file MD5 hashes                                                                                                                                                                                    | `False`          |
| `to_csvargs`           | `dict`            | extra parameters for pd.to_csv()                                                                                                                                                                                                     | `None`           |
| `to_tabular_kwargs`    | `dict`            | extra parameters for tabular format output (tsv, csv, parquet)                                                                                                                                                                        | `None`           |
| `float_formats`        | `dict`            | a dictionary to specify the float format for each column.                                                                                                                                                                            | `None`           |
| `validate`             | `boolean`         | If True, use gwas-ssf CLI tool for validation (only for SSF format).                                                                                                                                                                | `False`          |
| `verbose`              | `boolean`         | If True, print logs.                                                                                                                                                                                                                 | `True`           |
| `output_log`           | `boolean`         | If True, save the log to a file.                                                                                                                                                                                                     | `True`           |
| `ssfmeta`              | `boolean`         | If True, output a gwas-ssf-style meta file.                                                                                                                                                                                          | `False`          |

## Format dictionary

Using `float_formats`, you can specify the formats for numbers.

!!! info "Default formats for floating-point numbers"
    ```
    {'EAF': '{:.4g}', 'BETA': '{:.4f}', 'Z': '{:.4f}','CHISQ': '{:.4f}','SE': '{:.4f}','OR': '{:.4f}','OR_95U': '{:.4f}','OR_95L': '{:.4f}','INFO': '{:.4f}','P': '{:.4e}','MLOG10P': '{:.4f}','DAF': '{:.4f}'}
    ```

## Output File Naming

The output filename is automatically constructed based on the format and tabular format:
- Pattern: `{path}.{fmt}.{tab_fmt}[.gz]`
- Examples:
  - `path="./sumstats"`, `fmt="gwaslab"`, `tab_fmt="tsv"`, `gzip=True` → `sumstats.gwaslab.tsv.gz`
  - `path="./output"`, `fmt="ldsc"`, `tab_fmt="csv"`, `gzip=False` → `output.ldsc.csv`
  - `path="./data"`, `fmt="vcf"`, `bgzip=True` → `data.vcf.bcf.gz` (if bgzip) or `data.vcf.gz` (if gzip)

## Examples

GWASLab supports commonly used tabular formats, which are listed in a companion repository `formatbook`.

!!! quote "formatbook"
    For more details, please check [formatbook](https://github.com/Cloufield/formatbook)

!!! example "Basic format conversion"
    ```python
    # Convert to LDSC format
    mysumstats.to_format(path="./output", fmt="ldsc")
    # Output: output.ldsc.tsv.gz
    
    # Convert to PLINK format with CSV
    mysumstats.to_format(path="./output", fmt="plink", tab_fmt="csv", gzip=False)
    # Output: output.plink.csv
    
    # Convert to VCF with bgzip and tabix index
    mysumstats.to_format(path="./output", fmt="vcf", bgzip=True, tabix=True)
    # Output: output.vcf.bcf.gz and output.vcf.bcf.gz.tbi
    ```

!!! example "Filtering and formatting"
    ```python
    # Extract HapMap3 SNPs only
    mysumstats.to_format(path="./hapmap3", fmt="ldsc", hapmap3=True)
    
    # Exclude HLA region
    mysumstats.to_format(path="./no_hla", fmt="gwaslab", exclude_hla=True, hla_range=(25, 34))
    
    # Extract specific variants
    mysumstats.to_format(path="./subset", fmt="gwaslab", extract=["rs123", "rs456"], id_use="rsID")
    
    # Add chromosome prefix and N column
    mysumstats.to_format(path="./formatted", fmt="gwaslab", chr_prefix="chr", n=10000)
    ```

!!! example "Advanced options"
    ```python
    # Custom float formatting
    float_formats = {'P': '{:.2e}', 'BETA': '{:.6f}', 'SE': '{:.6f}'}
    mysumstats.to_format(path="./formatted", fmt="gwaslab", float_formats=float_formats)
    
    # Output as Parquet format
    mysumstats.to_format(path="./parquet", fmt="gwaslab", tab_fmt="parquet", gzip=False)
    
    # Generate MD5 checksum
    mysumstats.to_format(path="./checksummed", fmt="gwaslab", md5sum=True)
    ```

!!! example "CLI usage"
    ```bash
    # Basic format conversion
    gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt ldsc
    
    # With filtering options
    gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab \
      --hapmap3 --exclude-hla --n 10000
    
    # Output as CSV without compression
    gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab \
      --tab-fmt csv --no-gzip
    ```

See also: [Output sumstats](https://cloufield.github.io/gwaslab/format_load_save/)
