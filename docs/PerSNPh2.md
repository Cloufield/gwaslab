# Per SNP Heritability

GWASLab provides a simple function to calculate the variance explained by each SNP.

!!! tip "Available since v3.4.20"


## .get_per_snp_r2()

```
mysumstats.get_per_snp_r2()
```

It needs effect size `BETA` and effect allele frequency `EAF` for calculation.
If `N` is available, it will also calculate the F-statistics.

### Options

| `.get_per_snp_r2()` options | DataType   | Description                                | Default |
|-----------------------------|------------|--------------------------------------------|---------|
| `mode`                      | `q` or `b` | `q`: quantitative trait; `b`: binary trait | `q`     |

For quantitative traits (`mode="q"`), GWASLab will use `BETA`, `EAF` to calculate `SNPR2`.

| `.get_per_snp_r2()` options | DataType        | Description                                                | Default |
|-----------------------------|-----------------|------------------------------------------------------------|---------|
| `vary`                      | `float` or `se` | Var(Y); if `vary="se"`,Var(Y) will be estimated using `SE` | 1       |
| `k`                         | `int`           | k for calculating F                                        | 1       |

!!! info "For SNP `i`" 

    When `vary=1` and k=1:

    $$ h^2_{i} = 2 \times \beta_i^2 \times EAF \times (1 - EAF) $$

    $$ F = h^2_{i} \times (n -2) / (1 - h^2_{i}) $$

For binary traits  (`mode="b"`), `ncase`, `ncontrol` and `prevalence` are needed to estimate the variance of liability explained by variants:

| `.get_per_snp_r2()` options | DataType | Description                      | Default |
|-----------------------------|----------|----------------------------------|---------|
| `ncase`                     | `int`    | number of cases                  | -       |
| `ncontrol`                  | `int`    | number of controls               | -       |
| `prevalence`                | `float`  | prevalence in general population | -       |

!!! quote

    Equation 10 in Lee, S. H., Goddard, M. E., Wray, N. R., & Visscher, P. M. (2012). A better coefficient of determination for genetic profile analysis. Genetic epidemiology, 36(3), 214-224.
    
    - Implementation adopted from TwoSampleMR https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/add_rsq.r


## Examples

!!! example

    ```    
    mysumstats.get_per_snp_r2()
    ```

    See [Data conversion](https://cloufield.github.io/gwaslab/utility_data_conversion/#calculate-per-snp-r2)