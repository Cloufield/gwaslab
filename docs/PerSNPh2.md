# Per SNP Heritability

GWASLab provides a function to calculate the variance explained by each SNP (per-SNP heritability or R²), which is useful for Mendelian randomization, polygenic risk scores, and understanding the contribution of individual variants to trait variance.

!!! tip "Available since v3.4.20"

## .get_per_snp_r2()

Calculates the proportion of variance explained by each SNP and optionally computes F-statistics for instrument strength assessment.

```
mysumstats.get_per_snp_r2(**kwargs)
```

### Required Columns

- **For quantitative traits** (`mode="q"`): `BETA` and `EAF` (effect allele frequency)
- **For binary traits** (`mode="b"`): `BETA`, `EAF`, plus `ncase`, `ncontrol`, and `prevalence` parameters
- **For F-statistics**: `N` (sample size) column

### Parameters

#### General Parameters

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `mode` | `str` | Trait type: `"q"` for quantitative, `"b"` for binary | `"q"` |
| `beta` | `str` | Column name for effect size (beta coefficient) | `"BETA"` |
| `af` | `str` | Column name for effect allele frequency | `"EAF"` |
| `n` | `str` | Column name for sample size | `"N"` |
| `se` | `str` | Column name for standard error (used when `vary="se"`) | `"SE"` |
| `adjuested` | `bool` | If `True`, calculate adjusted R² | `False` |
| `verbose` | `bool` | Print progress messages | `True` |

#### Parameters for Quantitative Traits (`mode="q"`)

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `vary` | `float` or `"se"` | Variance of the phenotype Y. If `"se"`, Var(Y) is estimated from SE, N, and MAF | `1` |
| `k` | `int` or `"all"` | Number of parameters for F-statistic calculation. Use `"all"` to set k = number of SNPs | `1` |

#### Parameters for Binary Traits (`mode="b"`)

| Parameter | DataType | Description | Default |
|-----------|------|-------------|---------|
| `ncase` | `int` | Number of cases in the study | Required |
| `ncontrol` | `int` | Number of controls in the study | Required |
| `prevalence` | `float` | Disease prevalence in the general population (0-1) | Required |

### Output Columns

The function adds the following columns to the sumstats:

- **`SNPR2`**: Proportion of variance explained by each SNP (R²)
- **`ADJUESTED_SNPR2`**: Adjusted R² (only if `adjuested=True`)
- **`F`**: F-statistic for instrument strength (only if **N** column is available)
- **`_VAR(BETAX)`**: Variance of **BETA** × X (intermediate calculation for quantitative traits)
- **`_SIGMA2`**: Error variance (intermediate calculation when `vary="se"`)
- **`_POPAF`**: Population allele frequency (intermediate calculation for binary traits)
- **`_VG`**: Genetic variance (intermediate calculation for binary traits)

### Quantitative Traits (`mode="q"`)

For quantitative traits, the variance explained by SNP $i$ is calculated as:

$$ R^2_i = \frac{Var(\beta_i \times X_i)}{Var(Y)} $$

Where:

- $Var(\beta_i \times X_i) = 2 \times \beta_i^2 \times **EAF**_i \times (1 - **EAF**_i)$
- $Var(Y)$ depends on the `vary` parameter:

  - If `vary` is a number: $Var(Y) = \text{vary}$
  - If `vary="se"`: $Var(Y) = Var(\beta_i \times X_i) + **SE**_i^2 \times 2 \times **N** \times **EAF**_i \times (1 - **EAF**_i)$

**F-statistic** (when `N` is available):

$$ F_i = \frac{R^2_i \times (N - k - 1)}{(1 - R^2_i) \times k} $$

Where $k$ is the number of parameters (default: 1).

!!! info "When `vary=1` and `k=1`"

    The simplified formulas are:
    
    $$ R^2_i = 2 \times \beta_i^2 \times EAF_i \times (1 - EAF_i) $$
    
    $$ F_i = \frac{R^2_i \times (N - 2)}{1 - R^2_i} $$

!!! quote "Reference for quantitative traits"

    Shim, H., Chasman, D. I., Smith, J. D., Mora, S., Ridker, P. M., Nickerson, D. A., ... & Stephens, M. (2015). A multivariate genome-wide association analysis of 10 LDL subfractions, and their response to statin treatment, in 1868 Caucasians. *PloS one*, 10(4), e0120758.

### Binary Traits (`mode="b"`)

For binary traits, the function estimates the variance of liability explained by each variant using a liability threshold model.

The calculation involves:

1. Estimating population allele frequency from sample allele frequency, case-control ratio, odds ratio, and population prevalence
2. Calculating genetic variance: $V_G = \beta^2 \times p \times (1-p)$ where $p$ is the population allele frequency
3. Calculating R²: $R^2 = \frac{V_G}{V_G + V_E}$ where $V_E = \frac{\pi^2}{3} \approx 3.29$ (error variance for logistic regression)

!!! quote "Reference for binary traits"

    Equation 10 in Lee, S. H., Goddard, M. E., Wray, N. R., & Visscher, P. M. (2012). A better coefficient of determination for genetic profile analysis. *Genetic epidemiology*, 36(3), 214-224.
    
    Implementation adopted from [TwoSampleMR](https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/add_rsq.r).

### Examples

**Basic usage for quantitative trait:**

```
import gwaslab as gl

# Load sumstats
mysumstats = gl.Sumstats("t2d_bbj.txt.gz", 
                         snpid="SNP",
                         chrom="CHR",
                         pos="POS",
                         ea="ALT",
                         nea="REF",
                         neaf="Frq",
                         beta="BETA",
                         se="SE",
                         n=170000)

# Calculate per-SNP R² with default settings (vary=1)
mysumstats.get_per_snp_r2()

# View results
print(mysumstats.data[["SNP", "BETA", "EAF", "SNPR2", "F"]].head())
```

**Quantitative trait with estimated Var(Y):**

```
# Estimate Var(Y) from SE, N, and MAF
mysumstats.get_per_snp_r2(vary="se")
```

**Quantitative trait with custom Var(Y):**

```
# Provide known variance of the phenotype
mysumstats.get_per_snp_r2(vary=2.5)
```

**Binary trait:**

```
# For binary traits, provide case/control counts and prevalence
mysumstats.get_per_snp_r2(mode="b",
                          ncase=8500,
                          ncontrol=161500,
                          prevalence=0.10)
```

**With adjusted R²:**

```
# Calculate adjusted R²
mysumstats.get_per_snp_r2(adjuested=True)
```

**Custom k for F-statistic:**

```
# Use k=10 for F-statistic calculation
mysumstats.get_per_snp_r2(k=10)

# Or use all SNPs as k
mysumstats.get_per_snp_r2(k="all")
```

### Use Cases

- **Mendelian Randomization**: F-statistics > 10 indicate strong instruments
- **Polygenic Risk Scores**: R² values help assess individual variant contributions
- **Power Analysis**: Understanding variance explained by each variant
- **Quality Control**: Identifying variants with unexpectedly high or low R² values