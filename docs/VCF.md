# VCF I/O in GWASLab

GWASLab provides functions for working with VCF (Variant Call Format) and BCF files, including automatic chromosome detection, LD matrix calculation, and format conversion.

## Supported File Formats

GWASLab supports:
- VCF files (`.vcf`)
- Compressed VCF files (`.vcf.gz`)
- BCF files (`.bcf`)
- Indexed VCF/BCF files (with `.tbi` or `.csi` indices)

## VCF File Detection and Chromosome Handling

### `is_vcf_file()`

Check if a file is a VCF or BCF file by examining headers.

```
from gwaslab.io.io_vcf import is_vcf_file

if is_vcf_file("variants.vcf.gz"):
    print("This is a VCF/BCF file")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `path` | `str` | Path to the file to check | Required |

**Returns:**
- `bool` - True if the file is a VCF/BCF file, False otherwise

---

### `auto_check_vcf_chr_dict()`

Automatically determine chromosome naming convention used in VCF/BCF files.

```
from gwaslab.io.io_vcf import auto_check_vcf_chr_dict
from gwaslab.info.g_Log import Log

log = Log()
vcf_chr_dict = auto_check_vcf_chr_dict(
    vcf_path="variants.vcf.gz",
    vcf_chr_dict=None,
    verbose=True,
    log=log
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `vcf_path` | `str` or `None` | Path to the VCF/BCF file to check. If None, returns vcf_chr_dict | Required |
| `vcf_chr_dict` | `dict` or `None` | Optional pre-defined chromosome dictionary. If None, the function will attempt to determine the appropriate dictionary | `None` |
| `verbose` | `bool` | If True, print detailed progress messages | Required |
| `log` | `gwaslab.g_Log.Log` | Logging object for recording process information | Required |

**Returns:**
- `dict` - A chromosome dictionary mapping that matches the chromosome naming convention used in the VCF/BCF file

**Detection order:**
1. RefSeq IDs (for hg19 or hg38 builds)
2. Chromosome prefixes (e.g., "chr1", "Chr1", "CHR1")
3. Standard numeric chromosomes (e.g., 1, 2, 3, ...)

!!! note "Automatic Detection"
    The function automatically detects the chromosome naming convention and filters the dictionary to only include contigs present in the VCF file.

---

### `check_vcf_chr_prefix()`

Check for chromosome prefix in VCF/BCF file headers.

```
from gwaslab.io.io_vcf import check_vcf_chr_prefix
from gwaslab.info.g_Log import Log

log = Log()
prefix = check_vcf_chr_prefix("variants.vcf.gz", log, verbose=True)
if prefix:
    print(f"Chromosome prefix detected: {prefix}")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `vcf_bcf_path` | `str` | Path to VCF/BCF file | Required |
| `log` | `gwaslab.g_Log.Log` | Logging object | Required |
| `verbose` | `bool` | Whether to log detailed messages | Required |

**Returns:**
- `str` or `None` - Detected chromosome prefix (e.g., "chr", "Chr", "CHR") if found, otherwise None

---

### `check_vcf_chr_NC()`

Check for RefSeq chromosome IDs in VCF/BCF file headers.

```
from gwaslab.io.io_vcf import check_vcf_chr_NC
from gwaslab.info.g_Log import Log

log = Log()
chr_dict = check_vcf_chr_NC("variants.vcf.gz", log, verbose=True)
if chr_dict:
    print("RefSeq IDs detected")
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `vcf_bcf_path` | `str` | Path to VCF/BCF file | Required |
| `log` | `gwaslab.g_Log.Log` | Logging object | Required |
| `verbose` | `bool` | Whether to log detailed messages | Required |

**Returns:**
- `dict` or `None` - Chromosome mapping dictionary for detected build (hg19/hg38) if found, otherwise None

---

## LD Matrix Calculation from VCF

### `_get_ld_matrix_from_vcf()`

Calculate full LD matrix from VCF file and return both the LD matrix and corresponding sumstats.

```
from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
from gwaslab.info.g_Log import Log
import pandas as pd

# Prepare sumstats DataFrame
sumstats = pd.DataFrame({
    "CHR": [1, 1, 1],
    "POS": [1000, 2000, 3000],
    "NEA": ["A", "C", "G"],
    "EA": ["T", "G", "A"]
})

log = Log()
region = (1, 1000, 5000)  # (chromosome, start, end)

matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
    sumstats=sumstats,
    vcf_path="reference.vcf.gz",
    region=region,
    log=log,
    verbose=True,
    pos="POS",
    nea="NEA",
    ea="EA",
    vcf_chr_dict=None,
    tabix=None,
    export_path=None
)
```

| Parameter | DataType | Description | Default |
|-----------|----------|-------------|---------|
| `sumstats` | `pandas.DataFrame` | Summary statistics dataframe | Required |
| `vcf_path` | `str` | Path to the VCF file | Required |
| `region` | `tuple` | Region specification (chromosome, start, end) | Required |
| `log` | `Log` | Logging object | `Log()` |
| `verbose` | `bool` | Verbose flag | `True` |
| `pos` | `str` | Position column name | `"POS"` |
| `nea` | `str` | Non-effect allele column name | `"NEA"` |
| `ea` | `str` | Effect allele column name | `"EA"` |
| `vcf_chr_dict` | `dict` | Chromosome dictionary for VCF | `None` |
| `tabix` | `bool` | Whether to use tabix indexing | `None` (auto-detect) |
| `export_path` | `str` | Path to export results. If provided, sumstats and LD matrix will be saved | `None` |

**Returns:**
- `tuple` - (matched_sumstats, ld_matrix) where:
  - `matched_sumstats`: Subset of sumstats with valid variants, ordered to match ld_matrix
  - `ld_matrix`: Full LD matrix as numpy array (r² values)

**Matching algorithm:**
The function matches variants between sumstats and VCF using:
1. Position matching
2. Allele matching (handles strand flips)
3. Handles multiple variants at the same position

**LD calculation:**
- Uses `allel` library's `rogers_huff_r_between()` function
- Calculates pairwise LD (r²) for all matched variants
- Returns r² values (squared correlation coefficients)

!!! note "Export functionality"
    If `export_path` is provided, the function automatically saves:
    - Matched sumstats: `{prefix}_sumstats.tsv.gz`
    - LD matrix: `{prefix}_ldmatrix.npy`
    The prefix is automatically constructed from the region (e.g., `chr1_1000_5000`).

---

## Output to VCF Format

GWASLab can output summary statistics in GWAS-VCF format using the `.to_format()` method. See the [Format documentation](Format.md) for details.

```
# Output sumstats as VCF
mysumstats.to_format(path="./output", fmt="vcf", bgzip=True, tabix=True)
```

**GWAS-VCF format features:**
- Standard VCF format with GWAS-specific fields
- Automatic header generation with metadata
- Support for bgzip compression and tabix indexing
- Allele frequency information in INFO field
- Summary statistics in FORMAT column

---

## Examples

!!! example "Check VCF file and detect chromosome format"
    ```
    from gwaslab.io.io_vcf import is_vcf_file, auto_check_vcf_chr_dict
    from gwaslab.info.g_Log import Log
    
    vcf_path = "reference.vcf.gz"
    
    if is_vcf_file(vcf_path):
        log = Log()
        vcf_chr_dict = auto_check_vcf_chr_dict(
            vcf_path=vcf_path,
            vcf_chr_dict=None,
            verbose=True,
            log=log
        )
        print(f"Chromosome dictionary: {vcf_chr_dict}")
    ```

!!! example "Calculate LD matrix from VCF"
    ```
    from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
    from gwaslab.info.g_Log import Log
    import pandas as pd
    import numpy as np
    
    # Prepare sumstats
    sumstats = pd.DataFrame({
        "CHR": [1, 1, 1, 1],
        "POS": [1000, 2000, 3000, 4000],
        "NEA": ["A", "C", "G", "T"],
        "EA": ["T", "G", "A", "C"],
        "P": [1e-5, 1e-6, 1e-7, 1e-8]
    })
    
    log = Log()
    region = (1, 500, 5000)  # chr1:500-5000
    
    matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
        sumstats=sumstats,
        vcf_path="reference.vcf.gz",
        region=region,
        log=log,
        verbose=True
    )
    
    print(f"Matched variants: {len(matched_sumstats)}")
    print(f"LD matrix shape: {ld_matrix.shape}")
    print(f"LD between first two variants: {ld_matrix[0, 1]}")
    ```

!!! example "Export LD results"
    ```
    from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
    from gwaslab.info.g_Log import Log
    import pandas as pd
    
    sumstats = pd.DataFrame({
        "CHR": [1, 1, 1],
        "POS": [1000, 2000, 3000],
        "NEA": ["A", "C", "G"],
        "EA": ["T", "G", "A"]
    })
    
    log = Log()
    region = (1, 1000, 5000)
    
    # Export results to directory
    matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
        sumstats=sumstats,
        vcf_path="reference.vcf.gz",
        region=region,
        log=log,
        verbose=True,
        export_path="./ld_results"
    )
    
    # Files saved:
    # - ./ld_results/chr1_1000_5000_sumstats.tsv.gz
    # - ./ld_results/chr1_1000_5000_ldmatrix.npy
    ```

!!! example "Output sumstats as VCF"
    ```
    # Output in GWAS-VCF format
    mysumstats.to_format(
        path="./output",
        fmt="vcf",
        bgzip=True,
        tabix=True
    )
    
    # Output files:
    # - output.vcf.gz (bgzipped)
    # - output.vcf.gz.tbi (tabix index)
    ```

!!! example "Workflow: VCF to LD matrix for regional plot"
    ```
    from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf, auto_check_vcf_chr_dict
    from gwaslab.info.g_Log import Log
    import pandas as pd
    
    # Load your sumstats
    sumstats = pd.read_csv("sumstats.tsv.gz", sep="\t")
    
    # Define region of interest
    region = (6, 25000000, 35000000)  # HLA region
    
    # Filter sumstats to region
    region_sumstats = sumstats[
        (sumstats["CHR"] == region[0]) &
        (sumstats["POS"] >= region[1]) &
        (sumstats["POS"] <= region[2])
    ]
    
    log = Log()
    
    # Auto-detect chromosome format
    vcf_chr_dict = auto_check_vcf_chr_dict(
        vcf_path="reference.vcf.gz",
        vcf_chr_dict=None,
        verbose=True,
        log=log
    )
    
    # Calculate LD matrix
    matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
        sumstats=region_sumstats,
        vcf_path="reference.vcf.gz",
        region=region,
        log=log,
        verbose=True,
        vcf_chr_dict=vcf_chr_dict
    )
    
    # Use for regional plot
    # matched_sumstats and ld_matrix can be used with visualization functions
    ```

---

## Performance Tips

1. **Use tabix indexing**: For large VCF files, ensure they are indexed with tabix (`.tbi` or `.csi` files) for fast region-based queries.

2. **Filter by region**: When calculating LD matrices, filter your sumstats to the region of interest before calling `_get_ld_matrix_from_vcf()` to reduce computation time.

3. **Auto-detect chromosome format**: Use `auto_check_vcf_chr_dict()` to automatically handle different chromosome naming conventions without manual configuration.

4. **Export for reuse**: If you need to reuse LD matrices, use the `export_path` parameter to save results, avoiding recalculation.

5. **Memory considerations**: LD matrix calculation can be memory-intensive for large numbers of variants. Consider processing in smaller regions if memory is limited.

---

## VCF Format Specification

VCF (Variant Call Format) is a text file format for storing genetic variation data. The format consists of:

1. **Header lines**: Starting with `##`, containing metadata
2. **Column header**: Starting with `#CHROM`, defining column names
3. **Data lines**: Tab-separated variant information

**Standard VCF columns:**
- `#CHROM`: Chromosome
- `POS`: Position (1-based)
- `ID`: Variant identifier
- `REF`: Reference allele
- `ALT`: Alternate allele(s)
- `QUAL`: Quality score
- `FILTER`: Filter status
- `INFO`: Additional information
- `FORMAT`: Format of genotype data
- Sample columns: Genotype data for each sample

**GWAS-VCF format:**
GWASLab outputs summary statistics in GWAS-VCF format, which includes:
- Standard VCF columns
- Allele frequency in INFO field (`AF=...`)
- Summary statistics (BETA, SE, P, etc.) in FORMAT column
- Study metadata in header

For more details, see the [VCF format specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf) and [GWAS-VCF format](https://github.com/MRCIEU/gwas-vcf-spec).

