# QC and filtering
GWASLab provides all-in-one functions and customizable functions for sumstats QC and filtering.

## 1. Methods Summary

| Sumstats Methods               | Options                                                                     | Description                                                                           |
|--------------------------------|-----------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| `.check_sanity()`              | `n`.`ncase`,`ncontrol`,`beta`,`se`,`eaf` ...                                | sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N...               |
| `.check_data_consistency()`    |                                                                             | check if `BETA/SE-drived P/MLOG10P = original P/MLOG10P`, `N = N_CASE + N_CONTROL`... |
| `.remove_dup()`                | `mode="md"`, <br/>` keep='first'`, <br/>`keep_col="P"`, <br/>`remove=False` | remove duplicated, multiallelic or NA variants                                        |
| `.filter_value()`              | `expr` , <br/> `inplace=False`                                              | filter in variants base on expr                                                       |
| `.filter_flanking_by_id()`     | `snpid` , <br/> `inplace=False`                                             | filter in variants in the specified flacking regions (SNPID/rsID)                     |
| `.filter_flanking_by_chrpos()` | `chrpos` , <br/> `inplace=False`                                            | filter in variants in the specified flacking regions (CHR POS tuples)                 |
| `.filter_region_in()`          | `path` , <br/> `inplace=False` , <br/>`high_ld=False`, <br/> `build="19"`   | filter in variants in the specified region define by a bed file                       |
| `.filter_region_out()`         | `path` , <br/> `inplace=False` , <br/>`high_ld=False`, <br/> `build="19"`   | filter out variants in the specified region define by a bed file                      |

## 2. Statistics Sanity Check

Note: Default parameters have been updated since v3.3.36

`.check_sanity()`: Basic sanity check will. be performed on statistics to check if there are any `extreme values` or `values out of expected ranges`.

Comparison will be performed with `float_tolerence = 1e-7` for any float type statistics. For example, `eaf=(0, 1)` will be converted to `eaf=(-1e-7, 1 + 1e-7)`.

| Parameters                | Type       | Range                                    |
|---------------------------|------------|------------------------------------------|
| `float_tolerence`         | `float`    | tolerence for comparison                 |
| `n=(0,2**31-1))`          | `interger` | 0<N< $2^{31}-1$                          |
| `ncase=(0,2**31-1)`       | `interger` | 0<N< $2^{31}-1$                          |
| `ncontrol=(0,2**31-1)`    | `interger` | 0<N< $2^{31}-1$                          |
| `mac=(0,2**31-1)`         | `interger` | MAC>=0                                   |
| `eaf=(0,1)`               | `float`    | 0<EAF<1                                  |
| `chisq=(0,float("Inf"))`  | `float`    | CHISQ>0                                  |
| `p=(0,1)`                 | `float`    | 0<P<1   (Any P=0 will cause a warning)   |
| `mlog10p=(0,9999)`        | `float`    | 0<MLOG10P<9999                           |
| `beta=(-100,100)`         | `float`    | -10<BETA<10                              |
| `z=(-9999,9999)`          | `float`    | -9999<z<9999                             |
| `se=(0,float("Inf"))`     | `float`    | SE>0                                     |
| `OR=(-100,100)`           | `float`    | -100<log(OR)<100                         |
| `OR_95L=(0,float("Inf"))` | `float`    | OR_95L>0                                 |
| `OR_95U=(0,float("Inf"))` | `float`    | OR_95U>0                                 |
| `HR=(-100,100)`           | `float`    | -100<log(HR)<100                         |
| `HR_95L=(0,float("Inf"))` | `float`    | HR_95L>0                                 |
| `HR_95U=(0,float("Inf"))` | `float`    | HR_95U>0                                 |
| `info=(0,1)`              | `float`    | 0<INFO<1                                 |
| `direction`               | `string`   | only contains `"+"`,`"-"` ,`"0"`or `"?"` |

## 3. Remove duplicated or multiallelic variants

After standardizing and normalizing the sumstats, you can also remove duplicated or multiallelic variants using

```
.remove_dup(mode="md")
```

- `mode=d` , remove duplicate variants.
    - remove duplicate SNPs based on  1. SNPID, 
    - remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
    - remove duplicate SNPs based on  3. rsID
- `mode=s` ,remove duplicate variants.
    - remove duplicate SNPs based on  1. SNPID
- `mode=c` ,remove duplicate variants.
    - remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
- `mode=r` ,remove duplicate variants.
    - remove duplicate SNPs based on  3. rsID
- `mode=m`, remove multiallelic variants.
    - remove multiallelic SNPs based on  4. CHR, POS
- `remove=True` : remove NAs 
- `keep_col` : use which column to sort the values (`keep_ascend=True`: ascending order)
- `keep`: keep 'first' or 'last'.

!!! example
    ```python
    sumstats.remove_dup(mode="md",keep='first',keep_col="P",remove=False)
        
    Fri Jan 13 17:34:38 2023 Start to sort the sumstats using P...
    Fri Jan 13 17:34:38 2023 Start to remove duplicated variants based on snpid...
    Fri Jan 13 17:34:38 2023  -Current Dataframe shape : 9  x  11
    Fri Jan 13 17:34:38 2023  -Which variant to keep:  first
    Fri Jan 13 17:34:38 2023  -Removed  1  based on SNPID...
    Fri Jan 13 17:34:38 2023 Start to remove duplicated variants based on rsID...
    Fri Jan 13 17:34:38 2023  -Removed  1  based on rsID...
    Fri Jan 13 17:34:38 2023 Start to remove duplicated variants based on CHR,POS,EA and NEA...
    Fri Jan 13 17:34:38 2023  -Current Dataframe shape : 7  x  11
    Fri Jan 13 17:34:38 2023  -Which variant to keep:  first
    Fri Jan 13 17:34:38 2023  -Removed  1  based on CHR,POS,EA and NEA...
    Fri Jan 13 17:34:38 2023 Start to remove multiallelic variants based on chr:pos...
    Fri Jan 13 17:34:38 2023  -Which variant to keep:  first
    Fri Jan 13 17:34:38 2023  -Removed  0  multiallelic variants...
    Fri Jan 13 17:34:38 2023  -Removed  3  variants in total.
    Fri Jan 13 17:34:38 2023  -Sort the coordinates...
    Fri Jan 13 17:34:38 2023 Finished removing successfully!
    ```
    
    This will remove duplicated and multiallelic variants and keep the one with the lowest P.
    
    Before:
    
    <img width="525" alt="image" src="https://user-images.githubusercontent.com/40289485/212273929-330531bc-ed85-4e65-8eeb-0263b9250204.png">
    
    After
    
    <img width="525" alt="image" src="https://user-images.githubusercontent.com/40289485/212274043-fe37a99e-1fed-4340-944a-e731126e51f3.png">


## 4. Filtering by condition

Filter the sumstats by `expr` (a wrapper of `pandas.DataFrame.query`), and return a new Sumstats Object by default. This allows method chaining. For example, you can filter certain variants first and then create a Mahanttan plot like `mysumstats.filter_value('BETA<0 & CHR==1').plot_mqq()`.

```
.filter_value(expr,inplace=False)
```

| Options   | DataType  | Description                                                                                             | Default |
|-----------|-----------|---------------------------------------------------------------------------------------------------------|---------|
| `expr`    | `string`  | the query string used fot filtering. For example: '1>BETA>0 & N>10000'                                  |         |
| `inplace` | `boolean` | if False, return a new Sumstats object. If true, the current Sumstats object will be filtered in place. | `False` |

!!! quote pd.DataFrame.query()
    Please check https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html

!!! example
    ```
    mysumstats.filter_value('BETA<0 & CHR==1').plot_mqq()
    Wed Dec  7 01:32:06 2022 Start filtering values by condition: BETA<0 & CHR==1
    Wed Dec  7 01:32:06 2022  -Removing 12075769 variants not meeting the conditions: BETA<0 & CHR==1
    Wed Dec  7 01:32:06 2022 Finished filtering values.
    Wed Dec  7 01:32:06 2022 Start to plot manhattan/qq plot with the following basic settings:
    Wed Dec  7 01:32:06 2022  -Genome-wide significance level is set to 5e-08 ...
    Wed Dec  7 01:32:06 2022  -Raw input contains 481992 variants...
    Wed Dec  7 01:32:06 2022  -Plot layout mode is : mqq
    Wed Dec  7 01:32:06 2022 Finished loading specified columns from the sumstats.
    Wed Dec  7 01:32:06 2022 Start conversion and sanity check:
    Wed Dec  7 01:32:06 2022  -Removed 0 variants with nan in CHR or POS column ...
    Wed Dec  7 01:32:06 2022  -Removed 0 variants with nan in P column ...
    Wed Dec  7 01:32:06 2022  -P values are being converted to -log10(P)...
    Wed Dec  7 01:32:06 2022  -Sanity check after conversion: 0 variants with P value outside of (0,1] will be removed...
    Wed Dec  7 01:32:06 2022  -Sanity check: 0 na/inf/-inf variants will be removed...
    Wed Dec  7 01:32:07 2022  -Maximum -log10(P) values is 10.598771832501887 .
    Wed Dec  7 01:32:07 2022 Finished data conversion and sanity check.
    Wed Dec  7 01:32:07 2022 Start to create manhattan plot with 481992 variants:
    Wed Dec  7 01:32:09 2022  -Found 2 significant variants with a sliding window size of 500 kb...
    Wed Dec  7 01:32:09 2022 Finished creating Manhattan plot successfully!
    Wed Dec  7 01:32:09 2022  -Skip annotating
    Wed Dec  7 01:32:09 2022 Start to create QQ plot with 481992 variants:
    Wed Dec  7 01:32:09 2022  -Calculating GC using P : 1.2372836632762023
    Wed Dec  7 01:32:09 2022 Finished creating QQ plot successfully!
    ```
    ![image](https://user-images.githubusercontent.com/40289485/211584317-6c1583b5-53e4-4aae-9141-af5781e2439b.png)
 
## 5. Flanking regions

Available since v3.3.37

```
.filter_flanking_by_chrpos(snpid)

.filter_flanking_by_id(chrpos)
```

Extract variants in the specified flacking regions.

| Options        | DataType  | Description                                                                                             | Default |
|----------------|-----------|---------------------------------------------------------------------------------------------------------|---------|
| `snpid`        | `list`    | a list of reference SNPID or rsID. `["rs123", "1:123:A:G"]`                                             |         |
| `snpid`        | `list`    | a list of reference CHR, POS tuples `[(1, 12345), (2, 67891)]`                                          |         |
| `windonsizekb` | `int`     | flanking window size in kb `["rs123", "1:123:A:G"]`                                                     | `500`   |
| `inplace`      | `boolean` | if False, return a new Sumstats object. If true, the current Sumstats object will be filtered in place. | `False` |

## 6. Filtering regions

```
.filter_region_in()
.filter_region_out()
```

Filter variants in  predifined regions or regions defined in bed files.

| Options   | DataType  | Description                                                                                             | Default |
|-----------|-----------|---------------------------------------------------------------------------------------------------------|---------|
| `path`    | `string`  | path to the bed files                                                                                   | `None`  |
| `high_ld` | `boolean` | if True, filter high ld regions using built-in data                                                     | `False` |
| `inplace` | `boolean` | if False, return a new Sumstats object. If true, the current Sumstats object will be filtered in place. | `False` |

!!! example
    ```
    mysumstats.filter_region_in(high_ld=True,inplace=False).data

    mysumstats.filter_region_out(high_ld=True,inplace=False).data
    ```


## 7. filter_in & filter_out (deprecated)

```python
.filter_in(gt={},lt={},eq={},inplace=False)

.filter_out(gt={},lt={},eq={},inplace=False)
```
- `gt`: greater than
- `lt`: less than
- `eq`: equal to
- `inplace`: True or False. If False, return a dataframe. If true, the Sumstats object will be filtered.
