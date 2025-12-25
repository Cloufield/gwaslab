# Heritability conversion 

GWASLab provides functions to convert heritability estimates between observed-scale and liability-scale, which is essential for comparing heritability estimates across studies with different case-control ratios or for binary traits.

## When to use liability-scale conversion

For binary (case-control) traits, heritability estimates depend on the case-control ratio in the sample. To compare heritability estimates across studies or to estimate population-level heritability, conversion to the liability scale is necessary. The liability scale assumes an underlying continuous liability distribution with a threshold model for disease.

!!! quote
    Conversion formula (Equation 23 from Lee et al. 2011):
    $$
    h^2_{liability-scale} = h^2_{observed-scale} \times \frac{K^2(1-K)^2}{P(1-P) \times Z^2}
    $$
    
    Where:

    - $K$ : Population disease prevalence (proportion in the general population).
    - $P$ : Sample disease prevalence (proportion of cases in the study sample).
    - $Z$ : The height of the standard normal probability density function at threshold $T$. Calculated as `scipy.stats.norm.pdf(T, loc=0, scale=1)`.
    - $T$ : The threshold on the liability scale. Calculated as `scipy.stats.norm.isf(K)` (inverse survival function).
    
    Reference: Lee, S. H., Wray, N. R., Goddard, M. E., & Visscher, P. M. (2011). Estimating missing heritability for disease from genome-wide association studies. *The American Journal of Human Genetics*, 88(3), 294-305. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/

## gl.h2_obs_to_liab()

Converts heritability on the observed scale in an ascertained sample to heritability on the liability scale in the population.

```python
gl.h2_obs_to_liab(h2_obs, P, K, se_obs=None)
```

### Parameters

- `h2_obs` : float
    - Heritability on the observed scale in an ascertained sample.
- `P` : float in (0,1)
    - Prevalence of the phenotype in the sample (proportion of cases).
- `K` : float in (0,1)
    - Prevalence of the phenotype in the population (proportion in general population).
- `se_obs` : float, optional
    - Standard error of `h2_obs`. If provided, the function will also return the standard error on the liability scale.

### Returns

- If `se_obs` is `None`: Returns a single float value `h2_liab` (heritability on liability scale).
- If `se_obs` is provided: Returns a tuple `(h2_liab, se_liab)` where:
    - `h2_liab` : float - Heritability on the liability scale.
    - `se_liab` : float - Standard error of heritability on the liability scale.

### Examples

**Basic conversion without standard error:**

```python
import gwaslab as gl

# Example: Type 2 diabetes
# Observed-scale heritability from LDSC
h2_obs = 0.15

# Sample prevalence (50% cases in the study)
P = 0.5

# Population prevalence (estimated 10% in general population)
K = 0.10

# Convert to liability scale
h2_liab = gl.h2_obs_to_liab(h2_obs, P, K)
print(f"Liability-scale heritability: {h2_liab:.4f}")
```

**Conversion with standard error:**

```python
# With standard error
h2_obs = 0.15
se_obs = 0.02
P = 0.5
K = 0.10

h2_liab, se_liab = gl.h2_obs_to_liab(h2_obs, P, K, se_obs=se_obs)
print(f"Liability-scale heritability: {h2_liab:.4f} ± {se_liab:.4f}")
```

**Handling missing prevalence values:**

If both `P` and `K` are `NaN`, the function returns `h2_obs` unchanged (no conversion applied).

!!! quote
    Implementation note: Codes were adopted from [LDSC](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldscore/regressions.py). 
    
    Reference: Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., Ripke, S., Yang, J., Patterson, N., ... & Neale, B. M. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature genetics*, 47(3), 291-295.

## gl.h2_se_to_p()

Converts a heritability estimate and its standard error to a p-value using a one-sided test.

```python
gl.h2_se_to_p(h2, se)
```

### Parameters

- `h2` : float
    - Heritability estimate (can be on observed or liability scale).
- `se` : float
    - Standard error of the heritability estimate.

### Returns

- `p_value` : float
    - One-sided p-value testing whether the heritability is significantly greater than zero.

### Example

```python
import gwaslab as gl

# Heritability estimate and standard error
h2 = 0.15
se = 0.02

# Calculate p-value
p_value = gl.h2_se_to_p(h2, se)
print(f"P-value: {p_value:.4e}")
```

!!! note
    This function uses a one-sided test (survival function) to test if heritability is significantly greater than zero. The test statistic is $z = |h^2| / SE(h^2)$.
