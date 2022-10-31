# QC and filtering

## Methods Summary

| Sumstats Methods  | Options                  | Description                                                             |
| ----------------- | ------------------------ | ----------------------------------------------------------------------- |
| `.check_sanity()` |  `n=(0,float("Inf"))`, <br/>`eaf=(0,1)`, <br/>`mac=(5,float("Inf"))`, <br/>`chisq=(0,float("Inf"))`, <br/>`p=(5e-300,1)`, <br/>`mlog10p=(0,float("Inf"))`, <br/>`beta=(-10,10)`, <br/>`z=(-37.5,37.5)`, <br/>`se=(0,float("Inf"))`, <br/>`OR=(-10,10)` , <br/>`OR_95L=(0,float("Inf"))`, <br/>`OR_95U=(0,float("Inf"))`, <br/>`info=(0,float("Inf"))`   | sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N... |
| `.remove_dup()`   |  `mode="md"`, <br/>` keep='first'`, <br/>`keep_col="P"`, <br/>`remove=False` | remove duplicated, multiallelic or NA variants |
| `.filter_in()`    |  lt, gt, eq, inplace     |    filter in variants base on given threshold                                                                      |
| `.filter_out()`   |  lt, gt, eq, inplace     |       filter out variants base on given threshold                                                                      |
| `.filter_region_in()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                         |      filter in variants in the specified region define by a bed file                                                                   |
| `.filter_region_out()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                        |      filter out variants in the specified region define by a bed file                                                                  |


## Statistics Sanity Check

`.check_sanity()`: Basic sanity check will. be performed on statistics to check if there are any `extreme values` or `values out of expected range`.

- `n=(0,float("Inf"))` : interger, N>0
- `eaf=(0,1)` : float ,0<= EAF <=1
- `mac=(5,float("Inf"))` : mac>=5
- `chisq=(0,float("Inf"))` : float , CHISQ>0
- `p=(5e-300,1)` : float, 5e-300<P<=1
- `mlog10p=(0,float("Inf"))` : float, MLOG10>0
- `beta=(-10,10)` : float, -10<BETA<10
- `z=(-37.5,37.5)`: float, -37.5<z<37.5
- `se=(0,float("Inf"))` : float, SE>0
- `OR=(-10,10)` : float, -10<log(OR)<10
- `OR_95L=(0,float("Inf"))` :float, OR_95L>0
- `OR_95U=(0,float("Inf"))` :float, OR_95U>0
- `info=(0,float("Inf"))` : float, INFO>0 
- `Direction` : string, only contains "+","-" ,"0"or "?"

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

