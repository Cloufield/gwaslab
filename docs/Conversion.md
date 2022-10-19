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

## Example
```
mysumstats.fill_data(to_fill=["MLOG10P","Z","OR","OR_95L","OR_95U"])

Wed Oct 19 00:10:53 2022 Start filling data using existing columns...
Wed Oct 19 00:10:53 2022  -Raw input columns:  ['CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'Z', 'P', 'MLOG10P', 'INFO', 'SNPID', 'STATUS']
Wed Oct 19 00:10:53 2022  -Overwrite mode:  False
Wed Oct 19 00:10:53 2022   - Skipping columns:  ['MLOG10P']
Wed Oct 19 00:10:53 2022 Filling columns:  ['Z', 'OR', 'OR_95L', 'OR_95U']
Wed Oct 19 00:10:53 2022   - Filling OR using BETA column...
Wed Oct 19 00:10:53 2022   - Filling OR_95L/OR_95U using BETA/SE columns...
Wed Oct 19 00:10:55 2022   - Filling Z using BETA/SE column...
Wed Oct 19 00:10:58 2022 Finished filling data using existing columns.
```
