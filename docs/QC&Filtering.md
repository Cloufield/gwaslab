# QC and Filtering

GWASLab provides comprehensive quality control (QC) and filtering functions for summary statistics. These tools help you identify and remove problematic variants, check data consistency, and filter variants based on various criteria.

## Overview

GWASLab's QC and filtering capabilities include:

- **Quality checks**: Sanity checks for statistics, data consistency validation
- **Duplicate removal**: Remove duplicated and multiallelic variants
- **Variant filtering**: Filter by value, region, variant type, and more
- **All-in-one QC**: `basic_check()` performs multiple QC steps in one call

## Methods Summary

### Quality Control Methods

| Method | Options | Description |
|--------|---------|-------------|
| `.basic_check()` | `remove=False`, `remove_dup=False`, `normalize=True`, `threads=1`, ... | All-in-one QC function that performs multiple checks and fixes |
| `.check_sanity()` | `n`, `ncase`, `ncontrol`, `beta`, `se`, `eaf`, `p`, `mlog10p`, ... | Sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N... |
| `.check_data_consistency()` | | Check if `BETA/SE-derived P/MLOG10P = original P/MLOG10P`, `N = N_CASE + N_CONTROL`... |
| `.remove_dup()` | `mode="md"`, `keep='first'`, `keep_col="P"`, `remove=False` | Remove duplicated, multiallelic or NA variants |

### Filtering Methods

| Method | Options | Description |
|--------|---------|-------------|
| `.filter_value()` | `expr`, `inplace=False` | Filter variants based on expression (pandas query syntax) |
| `.filter_flanking_by_id()` | `snpid`, `windowsizekb=500`, `inplace=False` | Filter variants in flanking regions around specified SNPID/rsID |
| `.filter_flanking_by_chrpos()` | `chrpos`, `windowsizekb=500`, `inplace=False` | Filter variants in flanking regions around specified CHR:POS coordinates |
| `.filter_region_in()` | `path`, `high_ld=False`, `build="19"`, `inplace=False` | Filter in variants in regions defined by BED file |
| `.filter_region_out()` | `path`, `high_ld=False`, `build="19"`, `inplace=False` | Filter out variants in regions defined by BED file |
| `.filter_palindromic()` | `mode="out"`, `inplace=False` | Filter palindromic variants (keep or remove) |
| `.filter_snp()` | `mode="in"`, `inplace=False` | Filter SNPs (keep or remove) |
| `.filter_indel()` | `mode="out"`, `inplace=False` | Filter indels (keep or remove) |
| `.filter_hapmap3()` | `inplace=False` | Filter to HapMap3 variants only |
| `.exclude_hla()` | `inplace=False` | Exclude variants in HLA region (chr6:25-35 Mb) |

## Statistics Sanity Check

`.check_sanity()` performs basic sanity checks on statistics to identify `extreme values` or `values out of expected ranges`.

Comparison is performed with `float_tolerance = 1e-7` for any float type statistics. For example, `eaf=(0, 1)` will be converted to `eaf=(-1e-7, 1 + 1e-7)` to account for floating-point precision.

**Default Parameter Ranges:**

| Parameter | Type | Default Range | Description |
|-----------|------|--------------|-------------|
| `float_tolerance` | `float` | `1e-7` | Tolerance for floating-point comparisons |
| `n` | `tuple` | `(0, 2**31-1)` | Sample size: 0 < N < 2³¹-1 |
| `ncase` | `tuple` | `(0, 2**31-1)` | Number of cases: 0 < N_CASE < 2³¹-1 |
| `ncontrol` | `tuple` | `(0, 2**31-1)` | Number of controls: 0 < N_CONTROL < 2³¹-1 |
| `mac` | `tuple` | `(0, 2**31-1)` | Minor allele count: MAC ≥ 0 |
| `eaf` | `tuple` | `(0, 1)` | Effect allele frequency: 0 < EAF < 1 |
| `maf` | `tuple` | `(0, 0.5)` | Minor allele frequency: 0 < MAF < 0.5 |
| `chisq` | `tuple` | `(0, float("Inf"))` | Chi-square statistic: CHISQ > 0 |
| `z` | `tuple` | `(-9999, 9999)` | Z-score: -9999 < Z < 9999 |
| `t` | `tuple` | `(-99999, 99999)` | T-statistic: -99999 < T < 99999 |
| `f` | `tuple` | `(0, float("Inf"))` | F-statistic: F > 0 |
| `p` | `tuple` | `(0, 1)` | P-value: 0 < P < 1 (P=0 will cause a warning) |
| `mlog10p` | `tuple` | `(0, 99999)` | Negative log10 p-value: 0 < MLOG10P < 99999 |
| `beta` | `tuple` | `(-100, 100)` | Effect size: -100 < BETA < 100 |
| `se` | `tuple` | `(0, float("Inf"))` | Standard error: SE > 0 |
| `OR` | `tuple` | `(0, 100)` | Odds ratio: 0 < OR < 100 |
| `OR_95L` | `tuple` | `(0, float("Inf"))` | OR 95% CI lower bound: OR_95L > 0 |
| `OR_95U` | `tuple` | `(0, float("Inf"))` | OR 95% CI upper bound: OR_95U > 0 |
| `HR` | `tuple` | `(0, 100)` | Hazard ratio: 0 < HR < 100 |
| `HR_95L` | `tuple` | `(0, float("Inf"))` | HR 95% CI lower bound: HR_95L > 0 |
| `HR_95U` | `tuple` | `(0, float("Inf"))` | HR 95% CI upper bound: HR_95U > 0 |
| `info` | `tuple` | `(0, 2)` | Imputation info score: 0 < INFO < 2 |

**Usage:**

```python
# Use default ranges
mysumstats.check_sanity()

# Customize ranges
mysumstats.check_sanity(
    p=(1e-300, 1),           # Allow very small p-values
    beta=(-50, 50),          # Tighter beta range
    info=(0.8, 2),           # Require INFO > 0.8
    float_tolerance=1e-6     # Stricter tolerance
)
```

**Notes:**
- Columns not present in the sumstats are automatically skipped
- Invalid values are flagged but not automatically removed (use with `basic_check(remove=True)` to remove)
- The `direction` column is checked to ensure it only contains `"+"`, `"-"`, `"0"`, or `"?"`

## Remove duplicated or multiallelic variants

After standardizing and normalizing the sumstats, you can also remove duplicated or multiallelic variants using:

```
.remove_dup(mode="md")
```

- `mode=d` , remove duplicate variants.
    - remove duplicate SNPs based on  1. SNPID, 
    - remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
    - remove duplicate SNPs based on  3. rsID
- `mode=s` ,remove duplicate variants.
    - remove duplicate SNPs based on  1. SNPID
- `mode=c` ,remove duplicate variants.
    - remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
- `mode=r` ,remove duplicate variants.
    - remove duplicate SNPs based on  3. rsID
- `mode=m`, remove multiallelic variants.
    - remove multiallelic SNPs based on  4. CHR, POS
- `remove=True` : remove NAs 
- `keep_col` : use which column to sort the values (`keep_ascend=True`: ascending order)
- `keep`: keep 'first' or 'last'.

!!! example
    ```python
    sumstats.remove_dup(mode="md",keep='first',keep_col="P",remove=False)
        
    Fri Jan 13 17:34:38 2023 Start to sort the sumstats using P...
    Fri Jan 13 17:34:38 2023 Start to remove duplicated variants based on snpid...
    Fri Jan 13 17:34:38 2023  -Current Dataframe shape : 9  x  11
    Fri Jan 13 17:34:38 2023  -Which variant to keep:  first
    Fri Jan 13 17:34:38 2023  -Removed  1  based on SNPID...
    Fri Jan 13 17:34:38 2023 Start to remove duplicated variants based on rsID...
    Fri Jan 13 17:34:38 2023  -Removed  1  based on rsID...
    Fri Jan 13 17:34:38 2023 Start to remove duplicated variants based on CHR,POS,EA and NEA...
    Fri Jan 13 17:34:38 2023  -Current Dataframe shape : 7  x  11
    Fri Jan 13 17:34:38 2023  -Which variant to keep:  first
    Fri Jan 13 17:34:38 2023  -Removed  1  based on CHR,POS,EA and NEA...
    Fri Jan 13 17:34:38 2023 Start to remove multiallelic variants based on chr:pos...
    Fri Jan 13 17:34:38 2023  -Which variant to keep:  first
    Fri Jan 13 17:34:38 2023  -Removed  0  multiallelic variants...
    Fri Jan 13 17:34:38 2023  -Removed  3  variants in total.
    Fri Jan 13 17:34:38 2023  -Sort the coordinates...
    Fri Jan 13 17:34:38 2023 Finished removing successfully!
    ```
    
    This will remove duplicated and multiallelic variants and keep the one with the lowest P.
    
    Before:
    
    <img width="525" alt="image" src="https://user-images.githubusercontent.com/40289485/212273929-330531bc-ed85-4e65-8eeb-0263b9250204.png">
    
    After
    
    <img width="525" alt="image" src="https://user-images.githubusercontent.com/40289485/212274043-fe37a99e-1fed-4340-944a-e731126e51f3.png">

## Check_data consistency

```
.check_data_consistency()
```

GWASLab checks if `BETA/SE-derived P/MLOG10P = original P/MLOG10P` or `N = N_CASE + N_CONTROL`.

## Filtering Variants by Expression: `filter_value()`

Filter variants using pandas query syntax. This is the most flexible filtering method and allows complex conditional filtering.

```python
mysumstats.filter_value(expr, inplace=False)
```

**Parameters:**

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `expr` | `str` | Query string using pandas syntax. For example: `'BETA > 0 & N > 10000'` | Required |
| `inplace` | `bool` | If `False`, return a new Sumstats object. If `True`, filter in place. | `False` |

**Features:**
- Supports method chaining (returns new object by default)
- Uses pandas `DataFrame.query()` syntax
- Can combine multiple conditions with `&` (and), `|` (or), `~` (not)

!!! quote "pandas.DataFrame.query()"
    For detailed query syntax, see: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html

!!! example "Filter by p-value and effect size"
    ```python
    # Keep only significant variants with positive effect
    significant = mysumstats.filter_value('P < 5e-8 & BETA > 0')
    
    # Filter and plot in one chain
    mysumstats.filter_value('BETA < 0 & CHR == 1').plot_mqq()
    ```

!!! example "Complex filtering"
    ```python
    # Multiple conditions
    filtered = mysumstats.filter_value(
        'P < 1e-5 & INFO > 0.8 & N > 5000 & EAF > 0.01 & EAF < 0.99'
    )
    
    # Filter by chromosome
    chr1_variants = mysumstats.filter_value('CHR == 1')
    
    # Filter by range
    beta_range = mysumstats.filter_value('-1 < BETA < 1')
    ```
 
## Filtering Variants in Flanking Regions

Extract variants within a specified window around target variants. Useful for regional association plots and fine-mapping.

!!! info "Available since v3.4.37"

### `.filter_flanking_by_id()`

Filter variants in flanking regions around specified SNPID or rsID.

```python
mysumstats.filter_flanking_by_id(
    snpid=["rs123", "rs456", "1:12345:A:G"],
    windowsizekb=500,
    inplace=False
)
```

**Parameters:**

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `snpid` | `list` | List of SNPID or rsID. Examples: `["rs123"]`, `["1:12345:A:G"]`, `["rs123", "rs456"]` | Required |
| `windowsizekb` | `int` | Flanking window size in kilobases (kb) on each side | `500` |
| `inplace` | `bool` | If `False`, return a new Sumstats object. If `True`, filter in place. | `False` |

### `.filter_flanking_by_chrpos()`

Filter variants in flanking regions around specified chromosome and position coordinates.

```python
mysumstats.filter_flanking_by_chrpos(
    chrpos=[(1, 12345), (2, 67891)],
    windowsizekb=500,
    inplace=False
)
```

**Parameters:**

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `chrpos` | `list` | List of (CHR, POS) tuples. Example: `[(1, 12345), (2, 67891)]` | Required |
| `windowsizekb` | `int` | Flanking window size in kilobases (kb) on each side | `500` |
| `inplace` | `bool` | If `False`, return a new Sumstats object. If `True`, filter in place. | `False` |

!!! example "Extract variants around lead SNPs"
    ```python
    # Extract variants ±500kb around a lead SNP
    regional = mysumstats.filter_flanking_by_id(
        snpid=["rs123456"],
        windowsizekb=500
    )
    
    # Extract variants around multiple positions
    multiple_regions = mysumstats.filter_flanking_by_chrpos(
        chrpos=[(1, 12345678), (2, 87654321)],
        windowsizekb=200  # ±200kb window
    )
    ```

## Filtering Variants in Regions Using BED Files

Filter variants based on genomic regions defined in BED (Browser Extensible Data) format files or built-in regions.

### `.filter_region_in()`

Keep only variants within specified regions.

```python
mysumstats.filter_region_in(
    path="regions.bed",
    high_ld=False,
    build="19",
    inplace=False
)
```

### `.filter_region_out()`

Remove variants within specified regions.

```python
mysumstats.filter_region_out(
    path="regions.bed",
    high_ld=False,
    build="19",
    inplace=False
)
```

**Parameters:**

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `path` | `str` | Path to BED file. BED format: `chr\tstart\tend` (0-based, half-open) | `None` |
| `high_ld` | `bool` | If `True`, use built-in high-LD regions (commonly excluded regions) | `False` |
| `build` | `str` | Genome build for high-LD regions (`"19"` for hg19, `"38"` for hg38) | `"19"` |
| `inplace` | `bool` | If `False`, return a new Sumstats object. If `True`, filter in place. | `False` |

**BED File Format:**
```
chr1    1000000    2000000
chr2    5000000    6000000
chr3    10000000   11000000
```

!!! example "Filter by BED file"
    ```python
    # Keep variants in regions defined by BED file
    in_regions = mysumstats.filter_region_in(path="my_regions.bed")
    
    # Remove variants in high-LD regions (built-in)
    no_high_ld = mysumstats.filter_region_out(high_ld=True, build="38")
    
    # Remove variants in custom regions
    filtered = mysumstats.filter_region_out(path="exclude_regions.bed")
    ```

!!! example "Common use cases"
    ```python
    # Exclude high-LD regions (often done before LDSC)
    mysumstats.filter_region_out(high_ld=True, build="38")
    
    # Keep only variants in gene regions
    gene_regions = mysumstats.filter_region_in(path="genes.bed")
    ```

## Quick Reference: Common QC Workflows

### Standard QC Pipeline

```python
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")

# Comprehensive QC
mysumstats.basic_check(remove=True, remove_dup=True, normalize=True)

# Additional filtering
mysumstats.filter_palindromic(mode="out")  # Remove palindromic variants
mysumstats.exclude_hla()                    # Exclude HLA region
mysumstats.filter_value('INFO > 0.8 & EAF > 0.01 & EAF < 0.99')  # Quality filters

# Check results
mysumstats.summary()
```

### Pre-LDSC Filtering

```python
# Standard filtering for LDSC
mysumstats.filter_hapmap3()                 # HapMap3 variants only
mysumstats.filter_palindromic(mode="out")   # Remove palindromic
mysumstats.exclude_hla()                    # Exclude HLA
mysumstats.filter_region_out(high_ld=True, build="38")  # Exclude high-LD regions
mysumstats.filter_value('INFO > 0.9 & MAF > 0.01')      # Quality filters
```

### Regional Analysis

```python
# Extract variants around lead SNPs for regional plots
regional = mysumstats.filter_flanking_by_id(
    snpid=["rs123456", "rs789012"],
    windowsizekb=500
)

# Plot regional association
regional.plot_regional()
```

## Notes

- **Inplace vs New Object**: Most filtering methods return a new Sumstats object by default (`inplace=False`), allowing method chaining. Set `inplace=True` to modify the current object.
- **Method Chaining**: You can chain multiple filtering operations: `mysumstats.filter_value('P < 1e-5').filter_palindromic(mode="out").plot_mqq()`
- **Performance**: For large datasets, consider using `inplace=True` to avoid creating multiple copies of the data.
- **Status Codes**: Filtering operations may update the STATUS column. Check status codes after filtering to understand variant quality.
- **Logging**: All QC and filtering operations are logged. Use `mysumstats.log.show()` to review the operations performed.


## All-in-One QC: `basic_check()`

The `basic_check()` method is a comprehensive QC function that performs multiple quality control steps in a single call. It's recommended for initial QC of your sumstats.

```python
mysumstats.basic_check(
    remove=False,           # Remove bad quality variants
    remove_dup=False,       # Remove duplicated variants
    normalize=True,         # Normalize indels
    threads=1,              # Number of threads
    verbose=True
)
```

**What it does:**
- Fixes SNPID and rsID (`fix_id`)
- Standardizes chromosome notation (`fix_chr`)
- Standardizes position notation (`fix_pos`)
- Standardizes allele notation (`fix_allele`)
- Performs sanity checks on statistics (`check_sanity`)
- Checks data consistency (`check_data_consistency`)
- Normalizes indels (`normalize_allele`)
- Removes duplicates (if `remove_dup=True`)
- Sorts coordinates and columns

**Parameters:**
- `remove` (bool): Remove bad quality variants detected during fixing
- `remove_dup` (bool): Remove duplicated or multiallelic variants
- `normalize` (bool): Perform indel normalization
- `threads` (int): Number of threads for parallel processing
- `fix_id_kwargs`, `fix_chr_kwargs`, etc.: Keyword arguments passed to individual functions

!!! example "Basic QC workflow"
    ```python
    import gwaslab as gl
    
    # Load sumstats
    mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")
    
    # Run comprehensive QC
    mysumstats.basic_check(remove=True, remove_dup=True)
    
    # Check results
    mysumstats.summary()
    ```

## Additional Filtering Methods

### Filter by Variant Type

#### `.filter_palindromic()`

Filter palindromic variants (e.g., A/T, C/G SNPs that are ambiguous on the reverse strand).

```python
# Remove palindromic variants (recommended for most analyses)
mysumstats.filter_palindromic(mode="out", inplace=False)

# Keep only palindromic variants
mysumstats.filter_palindromic(mode="in", inplace=False)
```

**Parameters:**
- `mode` (str): `"in"` to keep palindromic variants, `"out"` to remove them (default: `"out"`)
- `inplace` (bool): If `False`, return a new Sumstats object. If `True`, filter in place.

#### `.filter_snp()` and `.filter_indel()`

Filter variants by type (SNPs vs indels).

```python
# Keep only SNPs (remove indels)
mysumstats.filter_snp(mode="in", inplace=False)

# Remove indels
mysumstats.filter_indel(mode="out", inplace=False)
```

**Parameters:**
- `mode` (str): `"in"` to keep the variant type, `"out"` to remove it
- `inplace` (bool): If `False`, return a new Sumstats object. If `True`, filter in place.

### Filter by Reference Panel

#### `.filter_hapmap3()`

Filter to variants present in HapMap3 reference panel (commonly used for LDSC and other analyses).

```python
mysumstats.filter_hapmap3(inplace=False)
```

**Parameters:**
- `inplace` (bool): If `False`, return a new Sumstats object. If `True`, filter in place.

### Exclude Specific Regions

#### `.exclude_hla()`

Exclude variants in the HLA (Human Leukocyte Antigen) region on chromosome 6 (chr6:25-35 Mb). This region has complex LD structure and is often excluded from certain analyses.

```python
mysumstats.exclude_hla(inplace=False)
```

**Parameters:**
- `inplace` (bool): If `False`, return a new Sumstats object. If `True`, filter in place.

!!! example "Common filtering workflow"
    ```python
    # Remove palindromic variants
    mysumstats.filter_palindromic(mode="out")
    
    # Remove indels
    mysumstats.filter_indel(mode="out")
    
    # Exclude HLA region
    mysumstats.exclude_hla()
    
    # Filter to HapMap3 variants (for LDSC)
    mysumstats.filter_hapmap3()
    ```

## Filtering with Thresholds: `filter_in()` and `filter_out()`

These methods provide a convenient way to filter variants based on threshold values.

```python
# Filter in variants (keep variants meeting criteria)
mysumstats.filter_in(
    gt={"INFO": 0.9, "N": 10000},      # Keep variants with INFO > 0.9 and N > 10000
    lt={"P": 1e-5},                     # Keep variants with P < 1e-5
    eq={"CHR": 1},                      # Keep variants on chromosome 1
    inplace=False
)

# Filter out variants (remove variants meeting criteria)
mysumstats.filter_out(
    gt={"P": 0.05},                     # Remove variants with P > 0.05
    lt={"INFO": 0.3},                   # Remove variants with INFO < 0.3
    eq={"CHR": "X"},                    # Remove variants on chromosome X
    inplace=False
)
```

**Parameters:**
- `gt` (dict): Dictionary of `{column: threshold}` for "greater than" filtering
- `lt` (dict): Dictionary of `{column: threshold}` for "less than" filtering
- `eq` (dict): Dictionary of `{column: value}` for "equal to" filtering
- `inplace` (bool): If `False`, return a new Sumstats object. If `True`, filter in place.

!!! note "Note on deprecated methods"
    The `filter_in()` and `filter_out()` methods are still available and functional, but `filter_value()` with pandas query syntax is generally more flexible and recommended for complex filtering.
