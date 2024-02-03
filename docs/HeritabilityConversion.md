# Heritability conversion 

GWASLab can convert Observed-scale heritability to Liability-scale heritability. 

!!! quote
    Conversion formula (Equation 23 from Lee. 2011):
    $$
    h^2_{liability-scale} = h^2_{observed-scale} * {{K(1-K)}\over{Z^2}} *  {{K(1-K)}\over{P(1-P)}}
    $$
    
    - $K$ : Population disease prevalence.
    - $P$ : Sample disease prevalence.
    - $Z$ : The height of the standard normal probability density function at threshold T. `scipy.stats.norm.pdf(T, loc=0, scale=1)`.
    - $T$ : The threshold. `scipy.stats.norm.ppf(1 - K, loc=0, scale=1)` or `scipy.stats.norm.isf(K)`.
    
    Reference: Estimating Missing Heritability for Disease from Genome-wide Association Studies https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/

## gl.h2_obs_to_liab()

```
gl.h2_obs_to_liab(h2_obs, P, K, se_obs=None)
```

## Parameters

- `h2_obs` : float. Heritability on the observed scale in an ascertained sample. 
- `P` : float in (0,1). Prevalence of the phenotype in the sample. 
- `K` : float in (0,1) . Prevalence of the phenotype in the population. 
- `se_obs` : float. se of h2_obs.

!!! quote
    Codes were adopted from [LDSC](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldscore/regressions.py). 
    
    Reference : Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., Ripke, S., Yang, J., Patterson, N., ... & Neale, B. M. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature genetics, 47(3), 291-295.
    
