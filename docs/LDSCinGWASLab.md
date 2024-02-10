# LDSC in GWASLab (BETA version)

!!! info "Available since v3.4.39"

LD score regression has become one of the most common methods to evaluate the inflation caused by confounding factors and evaluate the genetic correlation across traits in GWAS. 

The original [LDSC software](https://github.com/bulik/ldsc) was implemented in Python2 and was only available for the command line interface. 

GWASLab integrates the core functions of LDSC into the gl.Sumstats object, which makes the LD score regression much more convenient to conduct.

!!! info "The difference between original LDSC and LDSC in GWASLab"
    - GWASLab will automatically extract Hapmap3 SNPs based on CHR:POS and EA and NEA if rsID not avaiable in sumstats
    - Codes have been adjusted to be compatible with Python3. (`map`, `xrange` and so forth)
    - Sumstats were supplied by GWASLab instead of reading from files.
    - Log system has been replaced by GWASLab.Log
    - Fixed minor errors

!!! warning "LICENSE change"
    Since LDSC was integrated into GWASLab, LICENSE for GWASLab has also been changed from [MIT](https://github.com/Cloufield/gwaslab/blob/main/LICENSE_before_v3.4.39) to [GPL-3.0](https://github.com/Cloufield/gwaslab/?tab=GPL-3.0-1-ov-file#readme) license to be compatible with LDSC's LICENSE.

## Single variate LD score regression

!!! quote "Single variate LD score regression"
    Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.

```
.estimate_h2_by_ldsc()
```

| `.estimate_h2_by_ldsc()` options | DataType | Description                              | Default |
|----------------------------------|----------|------------------------------------------|---------|
| `ref_ld_chr`                     | `string` | ld score reference file directory        | -       |
| `w_ld_chr`                       | `string` | ld score weight reference file directory | -       |

Results (a pd.DataFrame) will be stored in `.ldsc_h2`. 

!!! example

    ```
    mysumstats.basic_check()
    mysumstats.estimate_h2_by_ldsc(ref_ld_chr = "/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr = "/home/yunye/tools/ldsc/ldscores/eas_ldscores/")
    mysumstats.ldsc_h2
    ```

For more examples, see [LDSC in gwaslab](https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/)

## Cross-trait LD score regression

!!! quote "Cross-trait LD score regression"
    Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.

```
.estimate_rg_by_ldsc()
```

| `.estimate_rg_by_ldsc()` options | DataType | Description                                    | Default |
|----------------------------------|----------|------------------------------------------------|---------|
| `other_traits`                   | `list`   | a list of gl.Sumstats objects for other traits | -       |
| `rg`                             | `string` | alias for each traits separated by commas      | -       |
| `ref_ld_chr`                     | `string` | ld score reference file directory              | -       |
| `w_ld_chr`                       | `string` | ld score weight reference file directory       | -       |
| `samp_prev`                      | `string` | prevalences in samples separated by commas     | -       |
| `pop_prev`                       | `string` | prevalences in population separated by commas  | -       |

Results (a pd.DataFrame) will be stored in `.ldsc_rg`. 

!!! example
    ```
    mysumstats.estimate_rg_by_ldsc(other_traits=[bmi_female,bmi_male], 
                               rg="T2D,BMI_female,BMI_male",
                               ref_ld_chr = "/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                               w_ld_chr = "/home/yunye/tools/ldsc/ldscores/eas_ldscores/")
    mysumstats.ldsc_rg
    ```

For more examples, see [LDSC in gwaslab](https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/)