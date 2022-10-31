# QC and filtering

## Methods Summary

| Sumstats Methods  | Options                  | Description                                                             |
| ----------------- | ------------------------ | ----------------------------------------------------------------------- |
| `.check_sanity()` |                          | sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N... |
| `.remove_dup()`  |  `mode="md"`, <br/>` keep='first'`, <br/>`keep_col="P"`, <br/>`remove=False` | remove duplicated, multiallelic or NA variants |
| `.filter_in()`    | lt, gt, eq, inplace      |                                                                         |
| `.filter_out()`   |                          |                                                                         |
| `.filter_region_in()`   | path , inplace=True ,high_ld=False, build="19"                         |                                                                         |
| `.filter_region_out()`   |  path , inplace=True ,high_ld=False, build="19"                       |                                                                         |


## Statistics Sanity Check

`.check_sanity()`: Basic sanity check will. be performed on statistics to check if there are any `extreme values` or `values out of expected range`.

`BETA/SE`: float, -10<BETA<10, -10<log(OR)<10

`OR/OR_95L/OR_95U`: float, 0<OR<10, OR_95L>0, OR_95U>0

`EAF`: 0<= EAF <=1, if EAF of >95% of valid variants is less than 0.5, a warning will be sent.

`P` : float, 0<P<5e-300

`MLOG10` : float, MLOG10>0

`Z` : float

`CHISQ` : float , CHISQ>0

`N` : interger, N>0

`Direction` : string, only contains "+","-" ,"0"or "?"

```python
sumstats.check_sanity()
```

## Remove duplicated or multiallelic variants

after standardize the sumstats, you can also remove duplicated or multiallelic variants using :

- `mode=d` ,remove duplicate.
    - remove duplicate SNPs based on  1. SNPID, 
    - remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
    - remove duplicate SNPs based on  3. rsID
- 'mode=m', remove multiallelic variants.
     - remove multiallelic SNPs based on  4. CHR, POS
- `remove=True` : remove NAs 
- `keep_col` : use which column to sort the values (keep_ascend=True: ascending order)
- `keep`: keep 'first' or 'last'.

```python
sumstats.remove_dup(mode="md",keep='first',keep_col="P",remove=False)
```
 -> remove duplicated and multiallelic variants and keep the one with lowest P.


## FIltering

```python
.filter_in(gt={},lt={},eq={},inplace=True)

.filter_out(gt={},lt={},eq={},inplace=True)
```

`gt`: greater than

`lt`: less than

`eq`: equal to

`inplace`: True or False. If False, return a dataframe. If true, the Sumstats object will be filtered.



## Filtering region
```
.filter_region_in(path="./abc.bed")
.filter_region_out(high_ld=True)
```

