# Statistics conversion

GWASLab can convert equvalent statistics, including:

|Target stats|Original stats|Implementation|
|-|-|-|
|MLOG10P|P|`sumstats["MLOG10P"] = -np.log10(sumstats["P"])` |
|P|MLOG10P|`sumstats["P"] = np.power(10,-sumstats["MLOG10P"])`|
|P|Z|`sumstats["P"] = ss.norm.sf(np.abs(sumstats["Z"])) * 2`|
|P|CHISQ|`sumstats["P"] = ss.chi2.sf(sumstats["CHISQ"], 1)`|
|OR<br />OR_95L<br />OR_95U|BETA<br />SE|`sumstats["OR"]   = np.exp(sumstats["BETA"])`, <br /> `sumstats["OR_95L"] = np.exp(sumstats["BETA"]-ss.norm.ppf(0.975)*sumstats["SE"])`, <br /> `sumstats["OR_95U"] = np.exp(sumstats["BETA"]+ss.norm.ppf(0.975)*sumstats["SE"])`|
| BETA <br /> SE|OR <br />OR_95L<br />OR_95U|`sumstats["BETA"]  = np.log(sumstats["OR"])  `, <br /> `sumstats["SE"]=(np.log(sumstats["OR"]) - np.log(sumstats["OR_95L"]))/ss.norm.ppf(0.975)`, <br /> `sumstats["SE"]=(np.log(sumstats["OR_95U"]) - np.log(sumstats["OR"]))/ss.norm.ppf(0.975)`|
|Z|BETA/SE|`sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]`|
|CHISQ|P|`sumstats["CHISQ"] = ss.chi2.isf(sumstats["P"], 1)`|
|CHISQ|Z|`sumstats["CHISQ"] = (sumstats["Z"])**2`|
|MAF|EAF|` sumstats["MAF"] =  sumstats["EAF"].apply(lambda x: min(x,1-x) if pd.notnull(x) else np.nan)`|


!!! info "Extreme P values"
    For extreme P, `extreme=True` can be added to overcome the limitation of extreme P values (P<1e-308):
    
    ```mysumstats.fill_data(to_fill=["MLOG10P"], extreme=True)```
    
    Z socres (or BETA and SE) will be used to calculate MLOG10P, two additional columns `P_MANTISSA` and `P_EXPONENT` will be added to present p values. 


!!! note
    The conversion is implemented using scipy and numpy.
    
    - ss : `import scipy.stats as ss`
    - np : `import numpy as np`

See examples [here.](https://cloufield.github.io/gwaslab/utility_data_conversion/)

## fill_data()

```
mysumstats.fill_data( 
    to_fill=[],
    df=None,
    overwrite=False,
    only_sig=False
    )
```

## Options

- `to_fill`: the columns to fill. ["OR","OR_95L","OR_95U","BETA","SE","P","MLOG10P","Z","CHISQ"]
- `df` : columns name for degree of freedom
- `overwrite`: if overwrite when the specified column existed
- `only_sig` : fill the data only for significant variants

## Priority

- For P : using MLOG10P, Z, CHISQ 
- For MLOG10P : using P, MLOG10P, Z, CHISQ 
- For BETA/SE : using OR/OR_95L/OR_95U
- For OR/OR_95L/OR_95U : using BETA/SE
- For Z : using BETA/SE
- For CHISQ : using  Z, P

## Example

!!! example
    ```
    # raw data
    #SNPID	CHR	POS	EA	NEA	EAF	BETA	SE	P	STATUS
    #1:725932_G_A	1	725932	G	A	0.9960	-0.0737	0.1394	0.5970	9999999
    #1:725933_A_G	1	725933	G	A	0.0040	0.0737	0.1394	0.5973	9999999
    #1:737801_T_C	1	737801	C	T	0.0051	0.0490	0.1231	0.6908	9999999
    
    # let's fill "MLOG10P","Z","OR","OR_95L","OR_95U"
    # gwaslab will automatically search for equivalent statistics
    
    mysumstats.fill_data(to_fill=["MLOG10P","Z","OR","OR_95L","OR_95U"])
    
    Wed Oct 19 10:13:30 2022 Start filling data using existing columns...
    Wed Oct 19 10:13:30 2022  -Raw input columns:  ['SNPID', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'STATUS']
    Wed Oct 19 10:13:30 2022  -Overwrite mode:  False
    Wed Oct 19 10:13:30 2022   - Skipping columns:  []
    Wed Oct 19 10:13:30 2022 Filling columns:  ['MLOG10P', 'OR', 'OR_95L', 'OR_95U']
    Wed Oct 19 10:13:30 2022   - Filling OR using BETA column...
    Wed Oct 19 10:13:31 2022   - Filling OR_95L/OR_95U using BETA/SE columns...
    Wed Oct 19 10:13:32 2022   - Filling MLOG10P using P column...
    Wed Oct 19 10:13:38 2022 Finished filling data using existing columns.
    ```

For more examples:
[Utility_data_conversion](https://github.com/Cloufield/gwaslab/blob/main/examples/Utility_data_conversion.ipynb)
