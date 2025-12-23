# Statistics conversion

GWASLab can convert equivalent statistics, including:

| Target stats               | Original stats              | Implementation                                                                                                                                                                                                                                     |
|----------------------------|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| MLOG10P                    | P                           | `sumstats["MLOG10P"] = -np.log10(sumstats["P"])`                                                                                                                                                                                                   |
| P                          | MLOG10P                     | `sumstats["P"] = np.power(10,-sumstats["MLOG10P"])`                                                                                                                                                                                                |
| P                          | Z                           | `sumstats["P"] = ss.chi2.sf(sumstats["Z"]**2, 1)` (equivalent to two-sided normal test)                                                                                                                                                            |
| P                          | CHISQ                       | `sumstats["P"] = ss.chi2.sf(sumstats["CHISQ"], 1)`                                                                                                                                                                                                 |
| OR<br />OR_95L<br />OR_95U | BETA<br />SE                | `sumstats["OR"]   = np.exp(sumstats["BETA"])`, <br /> `sumstats["OR_95L"] = np.exp(sumstats["BETA"]-ss.norm.ppf(0.975)*sumstats["SE"])`, <br /> `sumstats["OR_95U"] = np.exp(sumstats["BETA"]+ss.norm.ppf(0.975)*sumstats["SE"])`                  |
| BETA <br /> SE             | OR <br />OR_95L<br />OR_95U | `sumstats["BETA"]  = np.log(sumstats["OR"])  `, <br /> `sumstats["SE"]=(np.log(sumstats["OR"]) - np.log(sumstats["OR_95L"]))/ss.norm.ppf(0.975)`, <br /> `sumstats["SE"]=(np.log(sumstats["OR_95U"]) - np.log(sumstats["OR"]))/ss.norm.ppf(0.975)` |
| Z                          | BETA/SE                     | `sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]`                                                                                                                                                                                                  |
| CHISQ                      | P                           | `sumstats["CHISQ"] = ss.chi2.isf(sumstats["P"], 1)`                                                                                                                                                                                                |
| CHISQ                      | Z                           | `sumstats["CHISQ"] = (sumstats["Z"])**2`                                                                                                                                                                                                           |
| MAF                        | EAF                         | ` sumstats["MAF"] =  sumstats["EAF"].apply(lambda x: min(x,1-x) if pd.notnull(x) else np.nan)`                                                                                                                                                     |


!!! info "Extreme P values"
    For extreme P values (P < 1e-308), set `extreme=True` to overcome float64 precision limitations. MLOG10P will be calculated using Z-scores (or BETA/SE) or CHISQ using the method described [here](https://stackoverflow.com/questions/46416027/how-to-compute-p-values-from-z-scores-in-r-when-the-z-score-is-large-pvalue-muc/46416222#46416222):
    
    ```
    mysumstats.fill_data(to_fill=["MLOG10P"], extreme=True)
    ```
    
    <img width="446" alt="image" src="https://github.com/Cloufield/gwaslab/assets/40289485/7c1c96e3-f6e0-4232-a13f-af7fc551cc24">
    
    When `extreme=True`:
    - Z-scores (or BETA and SE) or CHISQ will be used to calculate MLOG10P directly, bypassing P-value calculation
    - Two additional columns `P_MANTISSA` and `P_EXPONENT` will be added to represent P values in scientific notation
    - This allows handling of extremely small P-values that exceed standard floating-point precision
    
    **Formulas for Extreme MLOG10P Calculation:**
    
    **From Z-scores (or BETA/SE):**
    ```
    log_pvalue = log(2) + norm.logsf(|Z|)  # two-sided test
    log10_pvalue = log_pvalue / log(10)
    MLOG10P = -log10_pvalue
    P_MANTISSA = 10^(log10_pvalue mod 1)
    P_EXPONENT = floor(log10_pvalue)
    ```
    
    Where:
    - `norm.logsf(|Z|)` is the log survival function of the standard normal distribution
    - `log(2)` accounts for the two-sided test
    - `P_MANTISSA` and `P_EXPONENT` represent P in scientific notation: P = P_MANTISSA × 10^P_EXPONENT
    
    **From CHISQ (with degrees of freedom):**
    ```
    log_pvalue = chi2.logsf(CHISQ, df)
    log10_pvalue = log_pvalue / log(10)
    MLOG10P = -log10_pvalue
    P_MANTISSA = 10^(log10_pvalue mod 1)
    P_EXPONENT = floor(log10_pvalue)
    ```
    
    Where:
    - `chi2.logsf(CHISQ, df)` is the log survival function of the chi-square distribution
    - `df` is the degrees of freedom (typically 1 for GWAS)
    

!!! note
    The conversion is implemented using scipy and numpy.
    
    - ss : `import scipy.stats as ss`
    - np : `import numpy as np`

See examples [here.](https://cloufield.github.io/gwaslab/utility_data_conversion/)

## fill_data()

```
mysumstats.fill_data( 
    to_fill=None,
    df=None,
    overwrite=False,
    verbose=True,
    only_sig=False,
    sig_level=5e-8,
    extreme=False
)
```

## Options

| Option      | DataType          | Description                                                                                                 | Default   |
|-------------|-------------------|-------------------------------------------------------------------------------------------------------------|-----------|
| `to_fill`   | `str` or `list`   | Column name(s) to fill. Valid values: `"OR"`, `"OR_95L"`, `"OR_95U"`, `"BETA"`, `"SE"`, `"P"`, `"Z"`, `"CHISQ"`, `"MLOG10P"`, `"MAF"`, `"SIG"`. Note: `"SIG"` creates a `"SIGNIFICANT"` column (boolean) based on P or MLOG10P threshold | `None`    |
| `df`        | `str`             | Column name containing degrees of freedom for chi-square tests (only used when filling CHISQ)              | `None`    |
| `overwrite` | `boolean`         | If True, overwrite existing values in target columns                                                       | `False`   |
| `verbose`   | `boolean`         | If True, display progress messages                                                                          | `True`    |
| `only_sig`  | `boolean`         | If True, only fill data for significant variants (P < sig_level)                                           | `False`   |
| `sig_level` | `float`           | Significance threshold for P-value filtering (used when `only_sig=True` or when filling `SIG` column)        | `5e-8`    |
| `extreme`   | `boolean`         | If True, use extreme value calculations for MLOG10P (helpful when P < 1e-300)                              | `False`   |

## Conversion Priority

GWASLab uses the following priority order when multiple source columns are available:

- **For P**: MLOG10P → Z → CHISQ
- **For MLOG10P**: P → Z → CHISQ (or BETA/SE if `extreme=True`)
- **For BETA/SE**: OR/OR_95L/OR_95U
- **For OR/OR_95L/OR_95U**: BETA/SE
- **For Z**: BETA/SE
- **For CHISQ**: Z → P
- **For MAF**: EAF (MAF = min(EAF, 1-EAF))
- **For SIG**: P or MLOG10P (creates SIGNIFICANT column: True if P < sig_level or MLOG10P > -log10(sig_level))

!!! note "Iterative Filling Process"
    The function performs iterative filling in multiple rounds:
    - **Round 1**: Attempts to fill all requested columns using available source data
    - **Subsequent rounds**: Newly filled columns may enable additional conversions
    - The process continues until all columns are filled or no further progress can be made
    - This allows complex conversions like: BETA/SE → Z → P → MLOG10P

!!! note "Column Handling"
    - If a target column already exists and `overwrite=False`, it will be skipped
    - When `extreme=True`, MLOG10P is calculated from Z-scores (or BETA/SE) to handle P-values below float64 precision limits
    - Intermediate columns created during conversion (but not requested) are automatically removed

## Examples

!!! example "Basic conversion"
    ```
    # Raw data with BETA, SE, P
    # SNPID	CHR	POS	EA	NEA	EAF	BETA	SE	P	STATUS
    # 1:725932_G_A	1	725932	G	A	0.9960	-0.0737	0.1394	0.5970	9999999
    
    # Fill missing statistics
    # GWASLab will automatically search for equivalent statistics
    mysumstats.fill_data(to_fill=["MLOG10P", "Z", "OR", "OR_95L", "OR_95U"])
    
    # Output:
    # Start filling data using existing columns...
    #  -Raw input columns:  ['SNPID', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'STATUS']
    #  -Overwrite mode:  False
    #  -Skipping columns:  []
    # Filling columns:  ['MLOG10P', 'OR', 'OR_95L', 'OR_95U', 'Z']
    #   - Filling OR using BETA column...
    #   - Filling OR_95L/OR_95U using BETA/SE columns...
    #   - Filling MLOG10P using P column...
    #   - Filling Z using BETA/SE columns...
    # Finished filling data using existing columns.
    ```

!!! example "Fill only significant variants"
    ```
    # Fill statistics only for variants with P < 5e-8
    mysumstats.fill_data(
        to_fill=["MLOG10P", "Z"],
        only_sig=True,
        sig_level=5e-8
    )
    ```

!!! example "Overwrite existing columns"
    ```
    # Overwrite existing MLOG10P values
    mysumstats.fill_data(
        to_fill=["MLOG10P"],
        overwrite=True
    )
    ```

!!! example "Fill single column"
    ```
    # Fill a single column (can pass string instead of list)
    mysumstats.fill_data(to_fill="Z")
    ```

!!! example "Fill MAF from EAF"
    ```
    # Calculate MAF from EAF
    mysumstats.fill_data(to_fill=["MAF"])
    # MAF = min(EAF, 1-EAF)
    ```

!!! example "Extreme P values"
    ```
    # Handle extreme P values using Z-scores
    mysumstats.fill_data(
        to_fill=["MLOG10P"],
        extreme=True
    )
    # Creates P_MANTISSA and P_EXPONENT columns
    # MLOG10P is calculated from Z-scores (or BETA/SE) to avoid precision loss
    ```

!!! example "Fill SIGNIFICANT column"
    ```
    # Create a SIGNIFICANT column based on P-value threshold
    mysumstats.fill_data(
        to_fill=["SIG"],
        sig_level=5e-8
    )
    # Creates SIGNIFICANT column: True if P < 5e-8 or MLOG10P > 7.3
    ```