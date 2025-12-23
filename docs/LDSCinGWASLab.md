# LDSC in GWASLab (BETA version)

!!! info "Available since v3.4.39"

LD score regression has become one of the most common methods to evaluate the inflation caused by confounding factors and evaluate the genetic correlation across traits in GWAS. 

The original [LDSC software](https://github.com/bulik/ldsc) was implemented in Python2 and was only available for the command line interface. 

GWASLab integrates the core functions of LDSC into the gl.Sumstats object, which makes the LD score regression much more convenient to conduct.

!!! info "The difference between original LDSC and LDSC in GWASLab"
    - GWASLab will automatically extract Hapmap3 SNPs based on CHR:POS and EA and NEA if rsID not available in sumstats
    - Codes have been adjusted to be compatible with Python3. (`map`, `xrange` and so forth)
    - Sumstats were supplied by GWASLab instead of reading from files.
    - Log system has been replaced by GWASLab.Log
    - Integrated munging workflow based on original LDSC munge_sumstats.py
    - Automatic column name handling (works with both GWASLab and LDSC standard column names)
    - Fixed minor errors

!!! warning "LICENSE change"
    Since LDSC was integrated into GWASLab, LICENSE for GWASLab has also been changed from [MIT](https://github.com/Cloufield/gwaslab/blob/main/LICENSE_before_v3.4.39) to [GPL-3.0](https://github.com/Cloufield/gwaslab/?tab=GPL-3.0-1-ov-file#readme) license to be compatible with LDSC's LICENSE.

## Munging (Filtering and Harmonization)

!!! info "Munging workflow"
    GWASLab implements the LDSC munging workflow based on the original [munge_sumstats.py](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py). Munging applies standard filtering and harmonization procedures to prepare summary statistics for LDSC analysis.

The munging process includes:

1. **Column name mapping**: Automatically maps GWASLab column names to LDSC standard format:
   - `EA` → `A1` (effect allele)
   - `NEA` → `A2` (non-effect allele)
   - `EAF` → `FRQ` (frequency)
   - `rsID` → `SNP` (variant ID)

2. **P-value filtering**: Removes SNPs with P-values outside (0, 1] with warnings

3. **INFO score filtering**: Filters SNPs with INFO < threshold (default 0.9), warns if INFO outside [0, 1.5]

4. **MAF filtering**: Converts EAF to MAF (minor allele frequency) and filters by MAF threshold (default 0.01), warns if frequency outside [0, 1]

5. **Allele filtering**: Keeps only strand-unambiguous SNPs (A/T, C/G, A/C, A/G, T/C, T/G)

6. **Palindromic SNP removal**: Optionally removes palindromic SNPs (default: True)

7. **Sample size filtering**: Filters by N using 90th percentile / 1.5 threshold (LDSC default) or user-specified value

8. **P to Z conversion**: Creates Z-scores from P-values (prefers BETA/SE if available for more accurate conversion)

9. **Duplicate removal**: Removes duplicate SNPs based on SNP ID

10. **Optional exclusions**: Can exclude HLA region and sex chromosomes

Munging can be enabled for `estimate_h2_by_ldsc()` using the `munge=True` parameter. All LDSC estimation functions are compatible with pre-munged data (they automatically handle both munged and non-munged column formats).

| Munging parameter | DataType | Description                              | Default |
|-------------------|----------|------------------------------------------|---------|
| `munge`           | `bool`   | If `True`, apply munging procedures (available for `estimate_h2_by_ldsc()`) | `False` |
| `munge_kwargs`    | `dict`   | Additional munging parameters            | `None`  |

!!! note "Munging compatibility"
    All LDSC estimation functions (`estimate_h2_by_ldsc`, `estimate_rg_by_ldsc`, `estimate_h2_cts_by_ldsc`, `estimate_partitioned_h2_by_ldsc`) can work with pre-munged data. The functions automatically detect and handle both GWASLab column names (EA/NEA/rsID) and LDSC standard names (A1/A2/SNP).

### Munging parameters

When `munge=True`, you can customize munging behavior using `munge_kwargs`:

| Parameter         | DataType | Description                              | Default |
|-------------------|----------|------------------------------------------|---------|
| `info`            | `float`  | Minimum INFO score threshold             | `0.9`   |
| `maf`             | `float`  | Minimum minor allele frequency threshold | `0.01`  |
| `n`               | `float` or `None` | Minimum sample size. If `None`, uses 90th percentile / 1.5 | `None` |
| `nopalindromic`   | `bool`   | If `True`, remove palindromic SNPs        | `True`  |
| `exclude_hla`     | `bool`   | If `True`, exclude HLA region             | `True`  |
| `exclude_sexchr`  | `bool`   | If `True`, exclude sex chromosomes        | `True`  |

## Single variate LD score regression

!!! quote "Single variate LD score regression"
    Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.

```
.estimate_h2_by_ldsc(build=None, verbose=True, match_allele=True, how="right", **kwargs)
```

| `.estimate_h2_by_ldsc()` options | DataType | Description                              | Default |
|----------------------------------|----------|------------------------------------------|---------|
| `ref_ld_chr`                     | `string` | **Required**. LD score reference file directory (e.g., `"/path/to/ldscores/"`) | -       |
| `w_ld_chr`                       | `string` | LD score weight reference file directory. Often the same as `ref_ld_chr` | -       |
| `build`                          | `str`    | Genome build version (e.g., `"19"`, `"38"`). If `None`, uses the build from sumstats metadata | `None`  |
| `verbose`                        | `bool`   | If `True`, print detailed progress messages | `True`  |
| `match_allele`                   | `bool`   | If `True`, match alleles with reference panel | `True`  |
| `how`                            | `str`    | Merge strategy for allele matching: `"left"`, `"right"`, `"inner"`, `"outer"` | `"right"` |
| `munge`                          | `bool`   | If `True`, apply standard munging procedures (filtering, harmonization, QC) | `False` |
| `munge_kwargs`                   | `dict`   | Additional parameters for munging (e.g., `info=0.9`, `maf=0.01`) | `None`  |
| `samp_prev`                      | `float`  | Sample prevalence (case proportion) for case-control studies | Auto from metadata |
| `pop_prev`                       | `float`  | Population prevalence for case-control studies | Auto from metadata |
| `print_coefficients`            | `str`    | Print coefficient results. Set to `"ldsc"` to enable | `"ldsc"` |

Results (a pd.DataFrame) will be stored in `.ldsc_h2`. If `print_coefficients` is enabled, coefficient results will be stored in `.ldsc_h2_results`. 

!!! example "Basic heritability estimation"

    ```
    mysumstats.basic_check()
    mysumstats.estimate_h2_by_ldsc(ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/")
    mysumstats.ldsc_h2
    ```

!!! example "With munging and additional options"

    ```
    mysumstats.estimate_h2_by_ldsc(ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/",
                                   munge=True,
                                   munge_kwargs={"info": 0.9, "maf": 0.01, "exclude_hla": True},
                                   print_coefficients="ldsc")
    # Access results
    mysumstats.ldsc_h2  # Summary results
    mysumstats.ldsc_h2_results  # Coefficient results (if print_coefficients is enabled)
    ```

!!! example "Custom munging parameters"

    ```
    # Customize munging with stricter filters
    mysumstats.estimate_h2_by_ldsc(ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/",
                                   munge=True,
                                   munge_kwargs={
                                       "info": 0.95,        # Stricter INFO filter
                                       "maf": 0.05,         # Higher MAF threshold
                                       "n": 5000,           # Minimum sample size
                                       "nopalindromic": True,
                                       "exclude_hla": True,
                                       "exclude_sexchr": True
                                   })
    mysumstats.ldsc_h2
    ```

!!! example "Case-control study with prevalence"

    ```
    mysumstats.estimate_h2_by_ldsc(ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/",
                                   samp_prev=0.5,  # 50% cases in sample
                                   pop_prev=0.1)  # 10% prevalence in population
    mysumstats.ldsc_h2
    ```

For more examples, see [LDSC in gwaslab](https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/)

## Cross-trait LD score regression

!!! quote "Cross-trait LD score regression"
    Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.

```
.estimate_rg_by_ldsc(build=None, verbose=True, match_allele=True, how="right", get_hm3=True, **kwargs)
```

| `.estimate_rg_by_ldsc()` options | DataType | Description                                    | Default |
|----------------------------------|----------|------------------------------------------------|---------|
| `other_traits`                   | `list`   | **Required**. A list of `gl.Sumstats` objects for other traits to compare | -       |
| `ref_ld_chr`                     | `string` | **Required**. LD score reference file directory (e.g., `"/path/to/ldscores/"`) | -       |
| `w_ld_chr`                       | `string` | LD score weight reference file directory. Often the same as `ref_ld_chr` | -       |
| `rg`                             | `string` | Alias for each trait separated by commas (e.g., `"T2D,BMI_female,BMI_male"`). If not provided, uses study names from metadata | Auto |
| `build`                          | `str`    | Genome build version (e.g., `"19"`, `"38"`). If `None`, uses the build from sumstats metadata | `None`  |
| `verbose`                        | `bool`   | If `True`, print detailed progress messages | `True`  |
| `match_allele`                   | `bool`   | If `True`, match alleles with reference panel | `True`  |
| `how`                            | `str`    | Merge strategy for allele matching: `"left"`, `"right"`, `"inner"`, `"outer"` | `"right"` |
| `get_hm3`                        | `bool`   | If `True`, filter to HapMap3 SNPs before analysis | `True`  |
| `samp_prev`                      | `string` | Sample prevalences separated by commas (e.g., `"0.5,0.3,0.4"`). Auto-detected from metadata if available | Auto |
| `pop_prev`                       | `string` | Population prevalences separated by commas (e.g., `"0.1,0.2,0.15"`). Auto-detected from metadata if available | Auto |

Results (a pd.DataFrame) will be stored in `.ldsc_rg`. 

!!! example "Basic genetic correlation"

    ```
    # Load other traits as Sumstats objects
    bmi_female = gl.Sumstats("bmi_female.txt.gz", fmt="gwaslab")
    bmi_male = gl.Sumstats("bmi_male.txt.gz", fmt="gwaslab")
    
    mysumstats.estimate_rg_by_ldsc(other_traits=[bmi_female, bmi_male], 
                                   ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/")
    mysumstats.ldsc_rg
    ```

!!! example "With custom trait aliases"

    ```
    mysumstats.estimate_rg_by_ldsc(other_traits=[bmi_female, bmi_male], 
                                   rg="T2D,BMI_female,BMI_male",  # Custom aliases
                                   ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/")
    mysumstats.ldsc_rg
    ```

!!! example "Case-control studies with prevalence"

    ```
    mysumstats.estimate_rg_by_ldsc(other_traits=[bmi_female, bmi_male], 
                                   ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/",
                                   samp_prev="0.5,0.3,0.4",  # Sample prevalences for each trait
                                   pop_prev="0.1,0.2,0.15")  # Population prevalences for each trait
    mysumstats.ldsc_rg
    ```

!!! example "Without HapMap3 filtering"

    ```
    mysumstats.estimate_rg_by_ldsc(other_traits=[bmi_female, bmi_male], 
                                   ref_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/", 
                                   w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/",
                                   get_hm3=False)  # Use all SNPs, not just HapMap3
    mysumstats.ldsc_rg
    ```

For more examples, see [LDSC in gwaslab](https://cloufield.github.io/gwaslab/ldsc_in_gwaslab/)

## Cell type specific heritability

!!! quote "Cell type specific heritability"
    Finucane, H. K., Reshef, Y. A., Anttila, V., Slowikowski, K., Gusev, A., Byrnes, A., ... & Price, A. L. (2018). Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature genetics, 50(4), 621-629.

```
.estimate_h2_cts_by_ldsc(build=None, verbose=True, match_allele=True, how="right", **kwargs)
```

| `.estimate_h2_cts_by_ldsc()` options | DataType | Description                              | Default |
|--------------------------------------|----------|------------------------------------------|---------|
| `ref_ld_chr_cts`                     | `string` | **Required**. LD score reference file directory for cell type specific analysis (e.g., `"/path/to/baseline/baseline."`) | -       |
| `build`                              | `str`    | Genome build version (e.g., `"19"`, `"38"`). If `None`, uses the build from sumstats metadata | `None`  |
| `verbose`                            | `bool`   | If `True`, print detailed progress messages | `True`  |
| `match_allele`                       | `bool`   | If `True`, match alleles with reference panel | `True`  |
| `how`                                | `str`    | Merge strategy for allele matching: `"left"`, `"right"`, `"inner"`, `"outer"` | `"right"` |
| `print_all_cts`                      | `bool`   | If `True`, print all cell type specific results | `False` |

Results (a pd.DataFrame) will be stored in `.ldsc_h2_cts`.

!!! example "Cell type specific heritability"

    ```
    mysumstats.estimate_h2_cts_by_ldsc(ref_ld_chr_cts="/home/yunye/tools/ldsc/eas_baseline/baseline1_2/baseline.")
    mysumstats.ldsc_h2_cts
    ```

## Partitioned heritability

!!! quote "Partitioned heritability"
    Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.

```
.estimate_partitioned_h2_by_ldsc(build=None, verbose=True, match_allele=True, how="right", **kwargs)
```

| `.estimate_partitioned_h2_by_ldsc()` options | DataType | Description                              | Default |
|----------------------------------------------|----------|------------------------------------------|---------|
| `ref_ld_chr`                                 | `string` | **Required**. LD score reference file directory with annotations (e.g., `"/path/to/annotations/"`) | -       |
| `w_ld_chr`                                   | `string` | LD score weight reference file directory. Often the same as `ref_ld_chr` | -       |
| `build`                                      | `str`    | Genome build version (e.g., `"19"`, `"38"`). If `None`, uses the build from sumstats metadata | `None`  |
| `verbose`                                    | `bool`   | If `True`, print detailed progress messages | `True`  |
| `match_allele`                               | `bool`   | If `True`, match alleles with reference panel | `True`  |
| `how`                                        | `str`    | Merge strategy for allele matching: `"left"`, `"right"`, `"inner"`, `"outer"` | `"right"` |
| `samp_prev`                                  | `float`  | Sample prevalence (case proportion) for case-control studies | Auto from metadata |
| `pop_prev`                                   | `float`  | Population prevalence for case-control studies | Auto from metadata |

Results will be stored in `.ldsc_partitioned_h2_summary` and `.ldsc_partitioned_h2_results`.

!!! example "Partitioned heritability"

    ```
    mysumstats.estimate_partitioned_h2_by_ldsc(ref_ld_chr="/home/yunye/tools/ldsc/annotations/", 
                                                w_ld_chr="/home/yunye/tools/ldsc/ldscores/eas_ldscores/")
    mysumstats.ldsc_partitioned_h2_summary  # Summary results
    mysumstats.ldsc_partitioned_h2_results   # Detailed results by annotation
    ```