# Standardization and Normalization

After loading raw sumstats into a GWASLab Sumstats object, the first step is typically to standardize variant-related notations, normalize indels, and check for unexpected errors in the statistics. When checking is finished, the status code will be automatically updated to reflect the quality and processing state of each variant.

## Overview

Standardization ensures that your sumstats data follows consistent formats and conventions:

- **IDs**: **SNPID** and **rsID** validation and formatting
- **Chromosomes**: Standardized notation (integers for autosomes, special handling for X, Y, MT)
- **Positions**: Validated base-pair positions within expected ranges
- **Alleles**: Standardized to ATCG notation, normalized for indels
- **Coordinates**: Sorted by genomic position
- **Columns**: Reordered to GWASLab default order

## Methods Summary

### Standardization Methods

| Method | Key Options | Description |
|--------|------------|-------------|
| `.fix_id()` | `fixchrpos=False`, `fixid=False`, `fixsep=False`, `overwrite=False`, `forcefixid=False` | Check and fix **rsID** or **SNPID** (chr:pos:ref:alt), or use **SNPID** to fix **CHR** and **POS** |
| `.fix_chr()` | `remove=False` | Standardize chromosome notation (species-aware) |
| `.fix_pos()` | `remove=False`, `limit=250000000`, `lower_limit=0`, `upper_limit=None` | Standardize base-pair position notation and filter out bad values |
| `.fix_allele()` | `remove=False` | Standardize base notations to ATCG |
| `.normalize_allele()` | `threads=1`, `mode="s"`, `chunk=3000000` | Normalize indels (left-alignment and parsimony principle) |
| `.sort_coordinate()` | | Sort variants by genomic coordinates (**CHR**, **POS**) |
| `.sort_column()` | `order=None` | Sort columns to GWASLab default order |
| `.basic_check()` | `remove=False`, `remove_dup=False`, `normalize=True`, `threads=1`, ... | All-in-one function to perform the complete standardization pipeline |

## IDs - SNPID and rsID

GWASLab requires at least one ID column for sumstats, either in the form of **SNPID** or **rsID** (or both). GWASLab will automatically check if **SNPID** is mixed in **rsID**.

### ID Formats

- **SNPID**: User-provided IDs, or in `**CHR**:**POS**:**REF**:**ALT**` format. Delimiters can be `":"`, `"_"`, or `"-"`.

  - Examples: `"1:12345:A:G"`, `"1_12345_A_G"`, `"1-12345-A-G"`
- **rsID**: dbSNP **rsID**s in the format `rs[numbers]`

  - Examples: `"rs123456"`, `"rs789012"`

### `.fix_id()`

Check or fix **SNPID** and **rsID**. This method can:

- Extract **CHR** and **POS** information from `**CHR**:**POS**:**REF**:**ALT**` formatted **SNPID**s
- Reconstruct **SNPID** from **CHR**, **POS**, **EA**, and **NEA**
- Standardize **SNPID** delimiters
- Validate **rsID** format

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `fixchrpos` | `bool` | `False` | If `True`, extract **CHR** and **POS** from **SNPID** (**CHR**:**POS**:**NEA**:**EA**) to fill **CHR** and **POS** columns |
| `fixid` | `bool` | `False` | If `True`, use **CHR**/**POS**/**NEA**/**EA** to reconstruct **SNPID**. For variants not aligned with reference genome, only **CHR**/**POS** will be used |
| `forcefixid` | `bool` | `False` | If `True`, use **CHR**/**POS**/**NEA**/**EA** to reconstruct **SNPID** without checking if the variant is aligned |
| `fixsep` | `bool` | `False` | If `True`, fix **SNPID** delimiter (e.g., `1:123_A_C` to `1:123:A:C`) |
| `overwrite` | `bool` | `False` | If `True`, overwrite existing data |
| `fixprefix` | `bool` | `False` | If `True`, remove 'chr' prefix in **SNPID** |
| `fixeanea` | `bool` | `False` | If `True`, extract **EA** and **NEA** from **SNPID** |
| `reversea` | `bool` | `False` | If `True`, reverse alleles in **SNPID** |

**Notes:**

- **SNPID** will be fixed to `**CHR**:**POS**:**NEA**:**EA**` format only when variants are already aligned with the reference genome
- Otherwise, a temporary **SNPID** in the format `**CHR**:**POS**` will be assigned
- The method automatically validates **rsID** format (must start with "rs" followed by numbers)

!!! example "Extract CHR and POS from SNPID"
    ```python
    # Extract chromosome and position from SNPID
    mysumstats.fix_id(fixchrpos=True)
    ```

!!! example "Reconstruct SNPID from coordinates and alleles"
    ```python
    # Reconstruct SNPID from CHR, POS, EA, NEA
    mysumstats.fix_id(fixid=True)
    
    # Force reconstruction without alignment check
    mysumstats.fix_id(forcefixid=True)
    ```

!!! example "Standardize SNPID format"
    ```python
    # Fix delimiter separators
    mysumstats.fix_id(fixsep=True)
    
    # Remove 'chr' prefix
    mysumstats.fix_id(fixprefix=True)
    ```

## CHR - Chromosomes

### `.fix_chr()`

Standardize chromosome notation. This method:
- Removes prefixes like "chr" from chromosome labels
- Converts string datatype to integers for autosomes
- Maps sex chromosomes and mitochondrial DNA to standardized integers
- Automatically uses species-specific chromosome mappings from the Sumstats object
- Removes variants with unrecognized chromosome notations (if `remove=True`)
- Removes variants with CHR < minimum autosome number

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `remove` | `bool` | `False` | If `True`, remove variants with invalid or unrecognized chromosome labels |
| `add_prefix` | `str` | `""` | Prefix to prepend to chromosome labels (rarely used) |

**Species-Aware Behavior:**
Chromosome mappings are automatically derived from the Sumstats object's `chromosomes` attribute (initialized based on the `species` parameter when creating the Sumstats object).

**For Human (default):**
- Integers `1-22` for autosomes
- `23` for X chromosome
- `24` for Y chromosome  
- `25` for MT (mitochondrial DNA)

**For Other Species:**
The function automatically uses the correct chromosomes based on the species:

```python
# Human (default)
mysumstats = gl.Sumstats("data.txt", species="homo sapiens")
mysumstats.fix_chr()  # Uses X, Y, MT -> 23, 24, 25

# Mouse
mysumstats = gl.Sumstats("data.txt", species="mouse")
mysumstats.fix_chr()  # Uses X, Y, MT -> 23, 24, 25 (19 autosomes)

# Chicken (ZW sex determination system)
mysumstats = gl.Sumstats("data.txt", species="chicken")
mysumstats.fix_chr()  # Uses Z, W, MT -> 23, 24, 25 (28 autosomes)

# Zebrafish (no sex chromosomes)
mysumstats = gl.Sumstats("data.txt", species="zebrafish")
mysumstats.fix_chr()  # Uses MT -> 25 (25 autosomes, no sex chromosomes)
```

**Supported Species:**
- Human (homo sapiens)
- Mouse (mus musculus)
- Rat (rattus norvegicus)
- Chicken (gallus gallus) - uses Z, W instead of X, Y
- Zebrafish (danio rerio) - no sex chromosomes
- Fruit fly (drosophila melanogaster)
- Pig (sus scrofa)
- Cattle (bos taurus)
- Dog (canis lupus familiaris)
- Horse (equus caballus)
- Rice (oryza sativa)
- Arabidopsis (arabidopsis thaliana)

!!! example "Standard chromosome fixing"
    ```python
    # Use species-specific chromosome mappings (automatic)
    mysumstats.fix_chr()
    
    # Remove invalid chromosomes
    mysumstats.fix_chr(remove=True)
    ```

!!! tip "Species-Specific Chromosomes"
    Chromosome information is managed by the `Chromosomes` class, which is automatically initialized when creating a Sumstats object. The `chromosomes` attribute contains all chromosome identifiers for the species.
    
    ```python
    # Access chromosome information
    mysumstats.chromosomes.chromosomes  # All chromosomes
    mysumstats.chromosomes.autosomes   # Autosomes only
    mysumstats.chromosomes.sex_chromosomes  # Sex chromosomes
    mysumstats.chromosomes.mitochondrial    # Mitochondrial chromosome
    ```

## POS - Base-pair Positions

### `.fix_pos()`

Check and fix values in **POS**. This method:

- Validates that **POS** values are positive integers
- Converts base-pair positions to integers (handles string formats with thousands separators)
- Converts invalid **POS** values to NA
- Removes **POS** outliers outside the specified range
- Updates status codes to reflect position validation state

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `remove` | `bool` | `False` | If `True`, remove variants with invalid or out-of-range positions |
| `limit` | `int` | `250000000` | Default upper limit for position validation (longest human chromosome is ~250 Mb) |
| `lower_limit` | `int` | `0` | Minimum acceptable genomic position |
| `upper_limit` | `int` | `None` | Maximum acceptable genomic position. If `None`, uses `limit` value |

**Notes:**

- Handles string-formatted positions with thousands separators (e.g., `"1,234,567"` → `1234567`)
- The default limit of 250,000,000 bp covers the longest human chromosome (chromosome 1)
- For other species, adjust `limit` or `upper_limit` as needed
- Invalid positions are flagged in the STATUS column

!!! example "Standard position fixing"
    ```python
    # Use default limits (0 to 250,000,000)
    mysumstats.fix_pos()
    
    # Remove invalid positions
    mysumstats.fix_pos(remove=True)
    
    # Custom limits for a different species
    mysumstats.fix_pos(lower_limit=1, upper_limit=500000000)
    ```

## Alleles

### Allele Notation Standardization

### `.fix_allele()`

Standardize allele representations to ATCG notation. This method:
- Validates that alleles contain only valid nucleotide characters (`A`, `T`, `C`, `G`)
- Converts lowercase letters to UPPERCASE
- Uses categorical data types for efficient storage
- Classifies variants as SNPs, indels, normalized, or not normalized
- Updates status codes to reflect allele validation state

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `remove` | `bool` | `False` | If `True`, remove variants with invalid allele representations |

**Notes:**

- Currently supports SNPs and INDELs only
- Copy number variants (CNV) like `<CN0>` won't be recognized
- Variants with invalid alleles (containing characters other than A, T, C, G) are flagged
- Status codes are updated to indicate allele quality

!!! example "Standardize alleles"
    ```python
    # Standardize allele notation
    mysumstats.fix_allele()
    
    # Remove variants with invalid alleles
    mysumstats.fix_allele(remove=True)
    ```

### Variant Normalization

### `.normalize_allele()`

Normalize indels according to the left-alignment and parsimony principle. This follows VCF normalization standards, removing common suffixes and prefixes from both alleles and adjusting positions accordingly.

**Example:** `chr1:123456:ATG:AT` will be normalized to `chr1:123457:TG:T`

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `threads` | `int` | `1` | Number of threads to use for parallel processing |
| `mode` | `str` | `"s"` | Normalization mode (`"s"` for standard, `"v"` for variant) |
| `chunk` | `int` | `3000000` | Size of chunks for parallel processing |

**How it works:**
- Only processes variants that need normalization (indels, status code digit 5 = 4 or 5)
- Left-aligns variants by removing common sequence context
- Adjusts positions to reflect the normalized representation
- Updates status codes to indicate normalization state

!!! warning "Reference Genome Limitation"
    Currently, normalization is implemented without checking a reference genome. This means it cannot normalize variants like `chr1:123456:G:-` where the missing allele information needs to be obtained from a reference genome.

!!! quote "Variant Normalization"
    For details on variant normalization principles, see: [Variant Normalization](https://genome.sph.umich.edu/wiki/Variant_Normalization)

!!! example "Normalize indels"
    ```python
    # Normalize with single thread
    mysumstats.normalize_allele(threads=1)
    
    # Normalize with multiple threads (faster for large datasets)
    mysumstats.normalize_allele(threads=4)
    ```
    
    **Before normalization:**
    
    <img width="345" alt="image" src="https://user-images.githubusercontent.com/40289485/212257048-1a5517a4-dd4d-4210-9e19-1dcf7747b7f5.png">
    
    **After normalization:**
    
    <img width="345" alt="image" src="https://user-images.githubusercontent.com/40289485/212256576-7808a1ec-5cf2-42ec-819e-a6f1c9c200bb.png">
    

## Coordinate and Column Sorting

### `.sort_coordinate()`

Sort variants by genomic coordinates (CHR, then POS). This is essential for many downstream analyses and ensures consistent ordering.

**Requirements:**
- CHR and POS must be fixed beforehand (use `fix_chr()` and `fix_pos()` first)

!!! example "Sort by coordinates"
    ```python
    # Fix coordinates first
    mysumstats.fix_chr()
    mysumstats.fix_pos()
    
    # Then sort
    mysumstats.sort_coordinate()
    ```

### `.sort_column()`

Sort columns to GWASLab default order. This ensures consistent column ordering across different datasets.

**Default Column Order:**
```python
"SNPID", "rsID", "CHR", "POS", "EA", "NEA", "EAF", "MAF", 
"BETA", "SE", "BETA_95L", "BETA_95U", "Z", "CHISQ", "P", "MLOG10P", 
"OR", "OR_95L", "OR_95U", "HR", "HR_95L", "HR_95U", 
"INFO", "N", "N_CASE", "N_CONTROL", "DIRECTION", 
"I2", "P_HET", "DOF", "SNPR2", "STATUS"
```
Additional columns are appended after the default columns.

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `order` | `list` | `None` | Custom column order. If `None`, uses GWASLab default order |

!!! example "Sort columns"
    ```python
    # Use default GWASLab column order
    mysumstats.sort_column()
    
    # Custom column order
    mysumstats.sort_column(order=["CHR", "POS", "SNPID", "P", "BETA"])
    ```

## All-in-One: `.basic_check()`

The `basic_check()` method performs a comprehensive standardization and quality control pipeline in a single call. It's the recommended starting point for processing new sumstats.

### Basic Check Workflow

The `basic_check()` process follows a sequential workflow to standardize and validate summary statistics:

```
Load Sumstats
     │
     ▼
Fix ID (fix_id)
     │
     ├─→ remove_dup=True ─→ Remove Duplicates ─┐
     │                                          │
     └─→ remove_dup=False ──────────────────────┘
                         │
                         ▼
                    Fix CHR (fix_chr)
                         │
                         ├─→ remove=True ─→ Remove bad CHR ─┐
                         │                                  │
                         └─→ remove=False ──────────────────┘
                                         │
                                         ▼
                                    Fix POS (fix_pos)
                                         │
                                         ├─→ remove=True ─→ Remove bad POS ─┐
                                         │                                  │
                                         └─→ remove=False ──────────────────┘
                                                         │
                                                         ▼
                                                    Fix Allele (fix_allele)
                                                         │
                                                         ├─→ remove=True ─→ Remove bad alleles ─┐
                                                         │                                      │
                                                         └─→ remove=False ──────────────────────┘
                                                                         │
                                                                         ▼
                                                                    Normalize Allele
                                                                    (if normalize=True)
                                                                         │
                                                                         ▼
                                                                    Check Sanity
                                                                         │
                                                                         ▼
                                                                    Check Consistency
                                                                         │
                                                                         ▼
                                                                    Sort Coordinate
                                                                         │
                                                                         ▼
                                                                    Sort Column
                                                                         │
                                                                         ▼
                                                              Standardized Sumstats
```

### Pipeline Order

The `basic_check()` method executes the following steps in order:

1. `fix_id()` - Fix SNPID and rsID
2. `remove_dup()` - Remove duplicates (if `remove_dup=True`)
3. `fix_chr()` - Standardize chromosome notation
4. `fix_pos()` - Standardize base-pair positions
5. `fix_allele()` - Standardize allele notation
6. `normalize_allele()` - Normalize indels (if `normalize=True`)
7. `check_sanity()` - Sanity check for statistics
8. `check_data_consistency()` - Check data consistency
9. `sort_coordinate()` - Sort by genomic coordinates
10. `sort_column()` - Sort columns to default order

**Parameters:**

| Parameter | DataType | Default | Description |
|-----------|------|---------|-------------|
| `remove` | `bool` | `False` | If `True`, remove bad quality variants detected in `fix_chr()`, `fix_pos()`, and `fix_allele()` |
| `remove_dup` | `bool` | `False` | If `True`, remove duplicated or multi-allelic variants using `remove_dup()` |
| `normalize` | `bool` | `True` | If `True`, perform indel normalization |
| `threads` | `int` | `1` | Number of threads to use for parallel processing (affects `normalize_allele()`) |
| `fix_id_kwargs` | `dict` | `{}` | Keyword arguments passed to `fix_id()` |
| `remove_dup_kwargs` | `dict` | `{}` | Keyword arguments passed to `remove_dup()` |
| `fix_chr_kwargs` | `dict` | `{}` | Keyword arguments passed to `fix_chr()` |
| `fix_pos_kwargs` | `dict` | `{}` | Keyword arguments passed to `fix_pos()` |
| `fix_allele_kwargs` | `dict` | `{}` | Keyword arguments passed to `fix_allele()` |
| `sanity_check_stats_kwargs` | `dict` | `{}` | Keyword arguments passed to `check_sanity()` |
| `consistency_check_kwargs` | `dict` | `{}` | Keyword arguments passed to `check_data_consistency()` |
| `normalize_allele_kwargs` | `dict` | `{}` | Keyword arguments passed to `normalize_allele()` |
| `verbose` | `bool` | `True` | If `True`, print log messages |

**Notes:**

- By default, `basic_check()` does not remove any variants (`remove=False`, `remove_dup=False`)
- All operations update the STATUS column to reflect data quality
- For details on `remove_dup()` and `check_sanity()`, see [QC and Filtering](https://cloufield.github.io/gwaslab/QC%26Filtering/)

!!! example "Basic QC workflow"
    ```python
    import gwaslab as gl
    
    # Load sumstats
    mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
    
    # Run comprehensive QC (doesn't remove variants by default)
    mysumstats.basic_check()
    
    # Run QC and remove bad variants
    mysumstats.basic_check(remove=True, remove_dup=True)
    
    # Customize individual steps
    mysumstats.basic_check(
        remove=True,
        remove_dup=True,
        normalize=True,
        threads=4,
        fix_chr_kwargs={"remove": True},
        fix_pos_kwargs={"limit": 300000000},
        sanity_check_stats_kwargs={"p": (1e-300, 1)}
    )
    ```

## Quick Reference: Standardization Workflow

### Recommended Pipeline

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")

# Option 1: Use basic_check() (recommended)
mysumstats.basic_check(remove=True, remove_dup=True, normalize=True)

# Option 2: Step-by-step (for fine control)
mysumstats.fix_id(fixchrpos=True)
mysumstats.fix_chr(remove=True)
mysumstats.fix_pos(remove=True)
mysumstats.fix_allele(remove=True)
mysumstats.normalize_allele(threads=4)
mysumstats.check_sanity()
mysumstats.check_data_consistency()
mysumstats.sort_coordinate()
mysumstats.sort_column()

# Check results
mysumstats.summary()
```

## Status Code Updates

All standardization methods automatically update the STATUS column to reflect:
- **Digit 1-2**: Genome build information
- **Digit 3**: rsID and SNPID validation status
- **Digit 4**: CHR and POS validation status
- **Digit 5**: EA and NEA standardization status
- **Digit 6**: REF-NEA alignment status
- **Digit 7**: Palindromic SNP/indel status

For details on status codes, see [Status Code](https://cloufield.github.io/gwaslab/StatusCode/) documentation.

