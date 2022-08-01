# Heritabilty conversion (Observed-scale -> Liability-scale)

```
gl.h2_obs_to_liab(h2_obs, P, K, se_obs=None)
```

`h2_obs` : float. Heritability on the observed scale in an ascertained sample. `P` : float in (0,1). Prevalence of the phenotype in the sample. `K` : float in (0,1) . Prevalence of the phenotype in the population. `se_obs` : float. se of h2_obs.

Adopted from LDSC.

Reference: Estimating Missing Heritability for Disease from Genome-wide Association Studies https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/
