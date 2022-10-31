# Statistics conversion

gwaslab can convert equvalent statistics, including:
- P => MLOG10P
- MLOG10P => P
- Z => P
- CHISQ => P
- BETA/SE => OR/OR_95L/OR_95U
- OR/OR_95L/OR_95U => BETA/SE
- BETA/SE => Z
- P => CHISQ
- Z => CHISQ

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
More examples:
[Utility_data_conversion](https://github.com/Cloufield/gwaslab/blob/main/examples/Utility_data_conversion.ipynb)
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


