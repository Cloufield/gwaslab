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
    ```
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
    ```
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
    ```
    # Custom float formatting
    float_formats = {'P': '{:.2e}', 'BETA': '{:.6f}', 'SE': '{:.6f}'}
    mysumstats.to_format(path="./formatted", fmt="gwaslab", float_formats=float_formats)
    
    # Output as Parquet format
    mysumstats.to_format(path="./parquet", fmt="gwaslab", tab_fmt="parquet", gzip=False)
    
    # Generate MD5 checksum
    mysumstats.to_format(path="./checksummed", fmt="gwaslab", md5sum=True)
    ```

!!! example "CLI usage"
    ```
    # Basic format conversion
    gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt ldsc
    
    # With filtering options
    gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab \
      --hapmap3 --exclude-hla --n 10000
    
    # Output as CSV without compression
    gwaslab --input sumstats.tsv --fmt gwaslab --out output --to-fmt gwaslab \
      --tab-fmt csv --no-gzip
    ```

## Special Format Specifications

GWASLab supports several specialized formats for variant annotation tools. These formats have specific coordinate conventions that are automatically handled.

### BED Format (0-based)

The BED (Browser Extensible Data) format is used by the UCSC Genome Browser and other tools. GWASLab outputs BED format with **0-based, half-open coordinates**.

**Output columns:** `CHR`, `START`, `END`, `NEA/EA`, `STRAND`, `SNPID`

**Coordinate conventions:**
- **SNPs**: `START = POS - 1`, `END = POS - 1 + len(NEA)` (0-based)
- **Insertions**: `START = POS`, `END = POS` (0-based)
- **Deletions**: `START = POS`, `END = POS + len(NEA) - 1` (0-based)

**Example:**
```python
mysumstats.to_format(path="./output", fmt="bed", bgzip=True, tabix=True)
# Output: output.bed.gz (bgzipped and tabix-indexed)
```

**Source:** [UCSC BED Format Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

---

### VEP Format (1-based)

The VEP (Variant Effect Predictor) format is used by Ensembl's VEP tool for variant annotation. GWASLab outputs VEP format with **1-based coordinates**.

**Output columns:** `CHR`, `START`, `END`, `NEA/EA`, `STRAND`, `SNPID`

**Coordinate conventions:**
- **SNPs**: `START = END = POS + (len(NEA) - 1)` (1-based)
- **Insertions**: `START = POS + 1`, `END = POS` (VEP convention: `START > END` indicates insertion)
- **Deletions**: `START = POS + 1`, `END = POS + (len(NEA) - 1)` (1-based)

!!! note "VEP Insertion Convention"
    VEP format uses `START > END` for insertions as a special convention to indicate an insertion between positions. This is intentional and correct according to VEP specifications.

**Example:**
```python
mysumstats.to_format(path="./output", fmt="vep", bgzip=True, tabix=True)
# Output: output.vep.gz (bgzipped and tabix-indexed)
```

**Source:** [Ensembl VEP Format Documentation](https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html)

---

### Annovar Format (1-based)

The Annovar format is used by the ANNOVAR tool for functional annotation of genetic variants. GWASLab outputs Annovar format with **1-based coordinates**.

**Output columns:** `CHR`, `START`, `END`, `NEA_out`, `EA_out`, `SNPID`

**Format specification:**
According to the [ANNOVAR documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/input/), the input format requires:
- **First 5 columns**: Chromosome, Start, End, Reference Allele, Alternative Allele
- **Coordinate system**: 1-based (by default)
- **Insertions**: Use `-` for reference allele
- **Deletions**: Use `-` for alternative allele

**Coordinate conventions:**
- **SNPs**: `START = POS`, `END = POS - 1 + len(NEA)` 
  - Example: SNP A/G at POS=100 → `START=100, END=100` (since len(NEA)=1)
  - Matches ANNOVAR example: `1 948921 948921 T C`
- **Insertions**: `START = END = POS` (matches ANNOVAR specification)
  - Example: Insertion A/AT at POS=200 → `START=200, END=200, REF=-, ALT=TC`
  - Matches ANNOVAR example: `1 11403596 11403596 - AT` (START=END=POS)
- **Deletions**: `START = POS`, `END = POS - 1 + len(NEA)`
  - Example: Deletion AT/A at POS=300 → `START=300, END=302, REF=TC, ALT=-`
  - Matches ANNOVAR example: `1 13211293 13211294 TC -` (2-bp deletion)

**Examples from ANNOVAR documentation:**
```
1 948921 948921 T C comments: rs15842, a SNP in 5' UTR
1 13211293 13211294 TC - comments: rs59770105, a 2-bp deletion
1 11403596 11403596 - AT comments: rs35561142, a 2-bp insertion
```

**Example:**
```python
mysumstats.to_format(path="./output", fmt="annovar", bgzip=True, tabix=True)
# Output: output.annovar.gz (bgzipped and tabix-indexed)
```

**Source:** [ANNOVAR Input Format Documentation](https://annovar.openbioinformatics.org/en/latest/user-guide/input/)

---

### Format Comparison

| Format | Coordinate System | Insertion Convention | Deletion Convention |
|--------|------------------|---------------------|-------------------|
| **BED** | 0-based, half-open | `START = END = POS` | `START = POS`, `END = POS + len(NEA) - 1` |
| **VEP** | 1-based | `START = POS + 1`, `END = POS` (START > END) | `START = POS + 1`, `END = POS + len(NEA) - 1` |
| **Annovar** | 1-based | `START = END = POS` | `START = POS`, `END = POS - 1 + len(NEA)` |

!!! tip "Automatic Coordinate Conversion"
    GWASLab automatically handles coordinate conversion for all variant types (SNPs, insertions, deletions) based on the selected format. You don't need to manually adjust coordinates.

---

### SSF Format (GWAS-SSF v0.1)

The SSF (Summary Statistics Format) is a standardized format for GWAS summary statistics proposed to improve interoperability and reproducibility. GWASLab supports outputting and validating SSF format files.

**Format Specification:**
- **Source**: [GWAS-SSF v0.1](https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1.full)
- **Separator**: Tab (`\t`)
- **Missing values**: `#NA`
- **File extension**: `.tsv` or `.tsv.gz`

**Required columns:**
- `chromosome` - Chromosome number
- `base_pair_location` - Base pair position
- `effect_allele` - Effect allele
- `other_allele` - Non-effect allele
- `standard_error` - Standard error
- `effect_allele_frequency` - Effect allele frequency
- `p_value` - P-value

**Optional columns:**
- `beta`, `odds_ratio`, or `hazard_ratio` - At least one effect measure required
- `neg_log_10_p_value` - Alternative to p_value
- `rsid` - dbSNP rsID
- `variant_id` - Variant identifier
- `info` - Imputation quality score
- `ref_allele` - Reference allele
- `n` - Sample size
- `ci_upper`, `ci_lower` - Confidence interval bounds

**Column order:** SSF format has a strict column order requirement. GWASLab automatically orders columns according to the SSF specification.

**Example:**
```python
# Output in SSF format
mysumstats.to_format(path="./output", fmt="ssf")

# Output SSF format with metadata file
mysumstats.to_format(path="./output", fmt="ssf", ssfmeta=True)

# Output and validate SSF format
mysumstats.to_format(path="./output", fmt="ssf", validate=True)
```

**SSF Metadata:**
When `ssfmeta=True`, GWASLab generates a YAML metadata file alongside the SSF file containing:
- Study information
- Sample characteristics
- File checksums (MD5)
- Format version information

---

### SSF Validation

GWASLab provides built-in validation for SSF format files to ensure compliance with the GWAS-SSF specification. The built-in validator replicates the behavior of `gwas-sumstats-validator` (also known as `gwas-ssf` CLI tool) and performs the same validation checks.

**Validation checks:**
1. **File extension**: Must be `.tsv` or `.tsv.gz`
2. **Required columns**: All 7 required columns must be present
3. **Effect field**: At least one of `beta`, `odds_ratio`, or `hazard_ratio` must be present
4. **P-value field**: Either `p_value` or `neg_log_10_p_value` must be present
5. **Column order**: Columns must follow the SSF specification order
6. **Chromosome coverage**: Validates presence of all autosomes (1-22)
7. **Data validation**: Checks data types, ranges, and consistency
8. **Minimum rows**: Requires at least 100,000 variants (configurable)

**Usage:**
```python
# Validate SSF file during output
mysumstats.to_format(path="./output", fmt="ssf", validate=True)

# Validation uses built-in validator by default
# If gwas-ssf CLI is available, it will be used instead
```

**Validation methods:**
- **Primary**: Uses `gwas-ssf` CLI tool (gwas-sumstats-validator) if available (external dependency)
- **Fallback**: Uses built-in GWASLab validator that replicates the same validation logic as `gwas-sumstats-validator` (no external dependencies)

**Validation output:**
- Success: `✓ SSF validation successful`
- Failure: Lists specific validation errors and issues found

**Example validation output:**
```
✓ SSF validation successful: File passes all validation checks
```

or

```
✗ SSF validation failed: Missing required columns: ['standard_error']
  - Missing required columns: ['standard_error']
```

See also: [Output sumstats](https://cloufield.github.io/gwaslab/format_load_save/)
