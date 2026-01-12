# Meta-Analysis

GWASLab provides comprehensive meta-analysis functionality for combining summary statistics from multiple GWAS studies. The implementation supports both fixed-effects and random-effects (DerSimonian-Laird) models.

## Overview

Meta-analysis in GWASLab can be performed using:

1. **`SumstatsPair`**: For combining two studies
2. **`SumstatsMulti`**: For combining three or more studies

Both classes provide a unified interface through the `.run_meta_analysis()` method.

## Quick Start

### Two Studies (SumstatsPair)

```python
import gwaslab as gl

# Load two sumstats
sumstats1 = gl.Sumstats("study1.txt.gz", verbose=False)
sumstats2 = gl.Sumstats("study2.txt.gz", verbose=False)

# Create SumstatsPair object
pair = gl.SumstatsPair(sumstats1, sumstats2, keep_all_variants=False)

# Run fixed-effects meta-analysis
result = pair.run_meta_analysis(random_effects=False)

# Run random-effects meta-analysis
result = pair.run_meta_analysis(random_effects=True)
```

### Multiple Studies (SumstatsMulti)

```python
import gwaslab as gl

# Load multiple sumstats
sumstats_list = [
    gl.Sumstats("study1.txt.gz", verbose=False),
    gl.Sumstats("study2.txt.gz", verbose=False),
    gl.Sumstats("study3.txt.gz", verbose=False)
]

# Create SumstatsMulti object
multi = gl.SumstatsMulti(sumstats_list, engine="pandas", verbose=False)

# Run meta-analysis
result = multi.run_meta_analysis(random_effects=False)
```

## SumstatsPair

### Creating a SumstatsPair Object

```python
pair = gl.SumstatsPair(
    sumstatsObject1,      # First Sumstats object
    sumstatsObject2,      # Second Sumstats object
    study=None,           # Optional study name
    suffixes=("_1", "_2"), # Column suffixes for each study
    keep_all_variants=True, # If True (default), keep variants present in only one study
    verbose=True
)
```

**Parameters:**

- `sumstatsObject1`, `sumstatsObject2`: `Sumstats` objects to combine
- `keep_all_variants`: 
  - `True` (default): Keep all variants from both studies (outer merge)
  - `False`: Only keep variants present in both studies (inner merge)
  - `True`: Keep all variants from both studies (outer merge)
- `suffixes`: Tuple of suffixes to append to column names (default: `("_1", "_2")`)

### Running Meta-Analysis

```python
result = pair.run_meta_analysis(
    random_effects=False,  # Use random-effects model
    verbose=True
)
```

**Returns:** A `Sumstats` object containing meta-analysis results

## SumstatsMulti

### Creating a SumstatsMulti Object

```python
multi = gl.SumstatsMulti(
    sumstatsObjects,      # List of Sumstats objects
    engine="pandas",      # "pandas" or "polars"
    merge_mode="outer",   # "inner" or "outer" (ignored if keep_all_variants is set)
    keep_all_variants=False,  # If True, keep variants present in only some studies
    verbose=True
)
```

**Parameters:**

- `sumstatsObjects`: List of `Sumstats` objects (minimum 2)
- `engine`: 
  - `"pandas"`: Uses pandas DataFrame (default)
  - `"polars"`: Uses polars DataFrame (faster for large datasets)
- `keep_all_variants`: 
  - `True` (default): Keep all variants from all studies (outer merge)
  - `False`: Only keep variants present in all studies (inner merge)
- `merge_mode`: Legacy parameter, overridden by `keep_all_variants` if provided

### Running Meta-Analysis

```python
result = multi.run_meta_analysis(
    random_effects=False,  # Use random-effects model
    verbose=True
)
```

**Returns:** 
- `Sumstats` object (pandas engine)
- `pl.DataFrame` (polars engine)

## Meta-Analysis Models

### Fixed-Effects Model

The fixed-effects model assumes all studies estimate the same underlying effect size. The weighted average is calculated as:

$$\beta_{meta} = \frac{\sum w_i \beta_i}{\sum w_i}$$

where $w_i = 1/SE_i^2$ is the inverse-variance weight.

**Standard error:**
$$SE_{meta} = \sqrt{\frac{1}{\sum w_i}}$$

### Random-Effects Model (DerSimonian-Laird)

The random-effects model accounts for between-study heterogeneity using the DerSimonian-Laird estimator.

**Between-study variance (tau²):**
$$\tau^2 = \max\left(0, \frac{Q - df}{C}\right)$$

where:
- $Q = \sum w_i \beta_i^2 - \frac{(\sum w_i \beta_i)^2}{\sum w_i}$ (Cochran's Q statistic)
- $df = k - 1$ (degrees of freedom, where $k$ is the number of studies)
- $C = \sum w_i - \frac{\sum w_i^2}{\sum w_i}$

**Random-effects weights:**
$$w_i^{RE} = \frac{1}{SE_i^2 + \tau^2}$$

**Random-effects estimate:**
$$\beta_{RE} = \frac{\sum w_i^{RE} \beta_i}{\sum w_i^{RE}}$$

## Output Columns

### Fixed-Effects Columns

| Column | Description |
|--------|-------------|
| `BETA` | Meta-analysis effect size (weighted average) |
| `SE` | Standard error of meta-analysis estimate |
| `Z` | Z-score (BETA / SE) |
| `P` | P-value (two-tailed) |
| `N` | Total sample size (sum across studies) |
| `EAF` | Weighted average effect allele frequency |
| `DOF` | Degrees of freedom (number of studies - 1) |
| `Q` | Cochran's Q statistic for heterogeneity |
| `P_HET` | P-value for heterogeneity test (chi-square) |
| `I2` | I² statistic (percentage of variance due to heterogeneity) |
| `DIRECTION` | Direction of effect in each study (e.g., "++" for both positive) |

### Random-Effects Columns (when `random_effects=True`)

Additional columns:

| Column | Description |
|--------|-------------|
| `BETA_RANDOM` | Random-effects meta-analysis effect size |
| `SE_RANDOM` | Standard error of random-effects estimate |
| `Z_RANDOM` | Z-score for random-effects model |
| `P_RANDOM` | P-value for random-effects model |

## Heterogeneity Statistics

### Cochran's Q

$$Q = \sum w_i \beta_i^2 - \frac{(\sum w_i \beta_i)^2}{\sum w_i}$$

Q follows a chi-square distribution with $df = k - 1$ degrees of freedom under the null hypothesis of no heterogeneity.

### I² Statistic

$$I^2 = \max\left(0, \frac{Q - df}{Q}\right) \times 100\%$$

I² represents the percentage of total variation across studies due to heterogeneity rather than chance:
- 0-25%: Low heterogeneity
- 25-50%: Moderate heterogeneity
- 50-75%: High heterogeneity
- 75-100%: Very high heterogeneity

## Handling Missing Data

GWASLab automatically handles missing data:

1. **Deduplication**: Variants are deduplicated by `CHR:POS:EA:NEA` before analysis
2. **Valid Study Mask**: Only studies with valid data are included:
   - `BETA` not null
   - `SE` not null and > 0
   - `N` not null
   - `EAF` not null
3. **NA Handling**: Variants with missing data in some studies are handled gracefully:
   - Fixed-effects: Only studies with valid data contribute to the meta-analysis
   - DOF is adjusted based on the number of studies with valid data

## Edge Cases

### Single Study

When only one study has valid data for a variant:
- `DOF = 0` (no heterogeneity test possible)
- `I2 = 0` (no heterogeneity)
- `P_HET = NaN` (cannot compute)
- Meta-analysis result equals the single study result

### Zero or Negative Q

When $Q \leq df$:
- `I2 = 0` (no heterogeneity detected)
- `tau² = 0` (random-effects equals fixed-effects)

### Division by Zero Guards

GWASLab includes guards to prevent division by zero:
- I² calculation: Returns 0 when $Q \leq 0$
- tau² calculation: Returns 0 when $C \leq 0$ or $Q \leq df$

## Examples

### Example 1: Basic Two-Study Meta-Analysis

```python
import gwaslab as gl

# Load studies
study1 = gl.Sumstats("study1.txt.gz", verbose=False)
study2 = gl.Sumstats("study2.txt.gz", verbose=False)

# Create pair
pair = gl.SumstatsPair(study1, study2, keep_all_variants=False)

# Run fixed-effects meta-analysis
result_fe = pair.run_meta_analysis(random_effects=False)

# Check results
print(result_fe.data[["SNPID", "BETA", "SE", "P", "I2"]].head())

# Run random-effects meta-analysis
result_re = pair.run_meta_analysis(random_effects=True)

# Compare fixed vs random effects
print(result_re.data[["SNPID", "BETA", "BETA_RANDOM", "SE", "SE_RANDOM"]].head())
```

### Example 2: Multiple Studies with Outer Merge

```python
import gwaslab as gl

# Load multiple studies
studies = [
    gl.Sumstats("study1.txt.gz", verbose=False),
    gl.Sumstats("study2.txt.gz", verbose=False),
    gl.Sumstats("study3.txt.gz", verbose=False)
]

# Create multi-object
multi = gl.SumstatsMulti(studies, engine="pandas", verbose=False)

# Run meta-analysis (keeps all variants)
result = multi.run_meta_analysis(random_effects=True)

# Filter for significant variants
significant = result.data[result.data["P_RANDOM"] < 5e-8]
print(f"Found {len(significant)} significant variants")
```

### Example 3: Checking Heterogeneity

```python
import gwaslab as gl

# Run meta-analysis
pair = gl.SumstatsPair(study1, study2)
result = pair.run_meta_analysis(random_effects=True)

# Check heterogeneity
heterogeneous = result.data[
    (result.data["I2"] > 50) & 
    (result.data["P_HET"] < 0.05)
]

print(f"Found {len(heterogeneous)} variants with significant heterogeneity")
print(heterogeneous[["SNPID", "I2", "P_HET", "BETA", "BETA_RANDOM"]].head())
```

## MR-MEGA (Meta-Regression for Genome-wide Association Studies)

MR-MEGA is a meta-regression method that accounts for genetic ancestry differences between cohorts using Multidimensional Scaling (MDS). Unlike standard fixed-effects or random-effects meta-analysis, MR-MEGA explicitly models ancestry-related heterogeneity by incorporating principal components derived from EAF (Effect Allele Frequency) distances between cohorts.

### Overview

MR-MEGA performs meta-regression by:

1. **Capturing genetic diversity**: Uses MDS to extract principal components from EAF correlation distances between cohorts
2. **Meta-regression model**: Fits a weighted least squares regression model:
   $$Y = \beta_0 + \beta_1 \cdot PC_1 + \beta_2 \cdot PC_2 + \ldots + \varepsilon$$
   where $Y$ is the effect size from each cohort, and $PC_1, PC_2, \ldots$ are MDS principal components
3. **Statistical tests**: Provides association, ancestry heterogeneity, and residual heterogeneity tests

**Reference**: Mägi, R., Horikoshi, M., Sofer, T., et al. (2017). Trans-ethnic meta-regression of genome-wide association studies accounting for ancestry increases power for discovery and improves fine-mapping resolution. *Human molecular genetics*, 26(18), 3639-3650.

### Quick Start

```python
import gwaslab as gl
from gwaslab.extension.mrmega import meta_regress_mrmega

# Load multiple sumstats
sumstats_list = [
    gl.Sumstats("study1.txt.gz", verbose=False),
    gl.Sumstats("study2.txt.gz", verbose=False),
    gl.Sumstats("study3.txt.gz", verbose=False),
    gl.Sumstats("study4.txt.gz", verbose=False),
    gl.Sumstats("study5.txt.gz", verbose=False)
]

# Create SumstatsMulti object
multi = gl.SumstatsMulti(sumstats_list, verbose=False)

# Perform MR-MEGA meta-regression
result = meta_regress_mrmega(
    multi,
    num_pcs=2,                    # Number of principal components
    use_genomic_control=True,     # Apply genomic control correction
    use_gco=True,                 # Apply second genomic control correction
    min_maf=0.01                  # Minimum MAF for marker selection
)

# Check results
print(result.data[["SNPID", "BETA", "SE", "P", "P_HET"]].head())
```

### Parameters

```python
result = meta_regress_mrmega(
    sumstats_multi,               # SumstatsMulti object
    num_pcs=2,                    # Number of principal components (must satisfy: num_pcs < ncohorts - 2)
    use_genomic_control=True,     # Apply genomic control correction to standard errors
    use_gco=True,                 # Apply second genomic control correction to output
    min_maf=0.01,                 # Minimum minor allele frequency for marker selection
    verbose=True                  # Print progress messages
)
```

**Parameters:**

- `sumstats_multi`: `SumstatsMulti` object containing summary statistics from multiple studies (minimum 3 studies recommended)
- `num_pcs`: Number of principal components to use. Must satisfy: `num_pcs < ncohorts - 2` (strict inequality)
- `use_genomic_control`: If `True`, applies genomic control correction to standard errors before regression
- `use_gco`: If `True`, applies a second genomic control correction to association chi-squares
- `min_maf`: Minimum minor allele frequency threshold for marker selection (default: 0.01 = 1%)


### Output Columns

| Column | Description |
|--------|-------------|
| `BETA` | Main association effect (intercept $\beta_0$) |
| `SE` | Standard error of main association effect |
| `Z` | Z-score (BETA / SE) |
| `P` | P-value for association test |
| `Q` | Residual heterogeneity chi-square statistic |
| `DOF` | Degrees of freedom for residual heterogeneity test |
| `P_HET` | P-value for residual heterogeneity test |
| `I2` | I² statistic (percentage of variance due to residual heterogeneity) |
| `EAF` | Weighted average effect allele frequency |
| `N` | Total sample size |
| `Ncohort` | Number of cohorts with valid data |

Additional MR-MEGA specific columns:
- `beta_0`, `se_0`: Intercept (main effect)
- `beta_1`, `se_1`, ..., `beta_{num_pcs}`, `se_{num_pcs}`: PC coefficients
- `chisq_association`, `ndf_association`, `P-value_association`: Association test statistics
- `chisq_ancestry_het`, `ndf_ancestry_het`, `P-value_ancestry_het`: Ancestry heterogeneity test statistics
- `chisq_residual_het`, `ndf_residual_het`, `P-value_residual_het`: Residual heterogeneity test statistics
- `lnBF`: Natural logarithm of Bayes factor

### When to Use MR-MEGA

MR-MEGA is particularly useful when:

1. **Trans-ethnic meta-analysis**: Combining studies from different ancestral populations
2. **Ancestry-related heterogeneity**: When genetic ancestry differences may confound associations
3. **Improved fine-mapping**: MR-MEGA can improve fine-mapping resolution by accounting for ancestry
4. **Multiple cohorts**: Works best with 5+ cohorts (minimum 3 required)

### Examples

#### Example 1: Basic MR-MEGA Analysis

```python
import gwaslab as gl
from gwaslab.extension.mrmega import meta_regress_mrmega

# Load multiple studies
studies = [
    gl.Sumstats("study1.txt.gz", verbose=False),
    gl.Sumstats("study2.txt.gz", verbose=False),
    gl.Sumstats("study3.txt.gz", verbose=False),
    gl.Sumstats("study4.txt.gz", verbose=False),
    gl.Sumstats("study5.txt.gz", verbose=False)
]

# Create SumstatsMulti object
multi = gl.SumstatsMulti(studies, verbose=False)

# Run MR-MEGA with 2 principal components
result = meta_regress_mrmega(multi, num_pcs=2)

# Check significant associations
significant = result.data[result.data["P"] < 5e-8]
print(f"Found {len(significant)} genome-wide significant variants")
```

#### Example 2: Comparing MR-MEGA with Standard Meta-Analysis

```python
import gwaslab as gl
from gwaslab.extension.mrmega import meta_regress_mrmega

# Create SumstatsMulti object
multi = gl.SumstatsMulti(studies, verbose=False)

# Standard fixed-effects meta-analysis
result_fe = multi.run_meta_analysis(random_effects=False)

# MR-MEGA meta-regression
result_mrmega = meta_regress_mrmega(multi, num_pcs=2)

# Compare results
comparison = result_fe.data.merge(
    result_mrmega.data[["SNPID", "BETA", "SE", "P"]],
    on="SNPID",
    suffixes=("_FE", "_MRMEGA")
)

print(comparison[["SNPID", "BETA_FE", "BETA_MRMEGA", "P_FE", "P_MRMEGA"]].head())
```

#### Example 3: Checking Ancestry Heterogeneity

```python
import gwaslab as gl
from gwaslab.extension.mrmega import meta_regress_mrmega

# Run MR-MEGA
result = meta_regress_mrmega(multi, num_pcs=2)

# Check for ancestry-related heterogeneity
ancestry_het = result.data[
    (result.data["P-value_ancestry_het"] < 0.05) &
    (result.data["P"] < 5e-8)
]

print(f"Found {len(ancestry_het)} significant variants with ancestry heterogeneity")
print(ancestry_het[["SNPID", "P-value_ancestry_het", "BETA", "P"]].head())
```

### Constraints and Limitations

1. **Minimum number of cohorts**: Requires at least 3 cohorts (recommended: 5+)
2. **PC constraint**: `num_pcs` must satisfy: `num_pcs < ncohorts - 2` (strict inequality)
3. **Marker selection**: Requires sufficient markers with valid EAF across all cohorts
4. **Computational cost**: More computationally intensive than standard meta-analysis due to MDS and regression steps