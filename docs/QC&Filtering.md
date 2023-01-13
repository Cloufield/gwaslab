# QC and filtering
gwaslab provides all-in-one function and customizable functions for sumstats QC and filtering.

See examples [here.](https://cloufield.github.io/gwaslab/quality_control_and_filtering/)

## Methods Summary

| Sumstats Methods  | Options                  | Description                                                             |
| ----------------- | ------------------------ | ----------------------------------------------------------------------- |
| `.check_sanity()` |  `n=(0,float("Inf"))`, <br/>`eaf=(0,1)`, <br/>`mac=(5,float("Inf"))`, <br/>`chisq=(0,float("Inf"))`, <br/>`p=(5e-300,1)`, <br/>`mlog10p=(0,float("Inf"))`, <br/>`beta=(-10,10)`, <br/>`z=(-37.5,37.5)`, <br/>`se=(0,float("Inf"))`, <br/>`OR=(-10,10)` , <br/>`OR_95L=(0,float("Inf"))`, <br/>`OR_95U=(0,float("Inf"))`, <br/>`info=(0,float("Inf"))`   | sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N... |
| `.remove_dup()`   |  `mode="md"`, <br/>` keep='first'`, <br/>`keep_col="P"`, <br/>`remove=False` | remove duplicated, multiallelic or NA variants |
| `.filter_value()`    |  expr     |    filter in variants base on expr                                                                    |
| `.filter_in()`    |  `lt`, `gt`, `eq`, `inplace`     |    filter in variants base on given threshold                                                                      |
| `.filter_out()`   |  `lt`, `gt`, `eq`, `inplace`     |       filter out variants base on given threshold                                                                      |
| `.filter_region_in()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                         |      filter in variants in the specified region define by a bed file                                                                   |
| `.filter_region_out()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                        |      filter out variants in the specified region define by a bed file                                                                  |


## Statistics Sanity Check

`.check_sanity()`: Basic sanity check will. be performed on statistics to check if there are any `extreme values` or `values out of expected range`.
|Parameters|Type|Range|
|-|-|-|
|`n=(0,2**64-1))` | `interger`| 0<N<2**64-1|
|`eaf=(0,1)` | `float` | 0<=EAF<=1|
|`mac=(5,float("Inf"))`| `float`| mac>=5|
|`chisq=(0,float("Inf"))` | `float` | CHISQ>0|
|`p=(5e-300,1)` | `float`| 5e-300<P<=1|
|`mlog10p=(0,float("Inf"))` | `float`| MLOG10>0|
|`beta=(-10,10)` | `float`| -10<BETA<10|
|`z=(-37.5,37.5)`| `float`| -37.5<z<37.5|
|`se=(0,float("Inf"))` | `float`| SE>0|
|`OR=(-10,10)` | `float`| -10<log(OR)<10|
|`OR_95L=(0,float("Inf"))` |`float`| OR_95L>0|
|`OR_95U=(0,float("Inf"))` |`float`| OR_95U>0|
|`info=(0,float("Inf"))` | `float`| INFO>0 |
|`direction` | `string`| only contains "+","-" ,"0"or "?"|

!!! example
    ```python
    sumstats.check_sanity()
    
    Wed Dec  7 01:31:10 2022 Start sanity check for statistics ...
    Wed Dec  7 01:31:10 2022  -Current Dataframe shape : 12557761  x  12
    Wed Dec  7 01:31:18 2022  -Checking if  0 <=N<= inf  ...
    Wed Dec  7 01:31:26 2022  -Removed 0 variants with bad N.
    Wed Dec  7 01:31:26 2022  -Checking if  0 <=EAF<= 1  ...
    Wed Dec  7 01:31:30 2022  -Removed 0 variants with bad EAF.
    Wed Dec  7 01:31:30 2022  -Checking if  5 <=MAC<= inf  ...
    Wed Dec  7 01:31:33 2022  -Removed 0 variants with bad MAC.
    Wed Dec  7 01:31:33 2022  -Checking if  5e-300 <= P <= 1  ...
    Wed Dec  7 01:31:35 2022  -Removed 0 variants with bad P.
    Wed Dec  7 01:31:35 2022  -Checking if  -10 <BETA)< 10  ...
    Wed Dec  7 01:31:36 2022  -Removed 0 variants with bad BETA.
    Wed Dec  7 01:31:36 2022  -Checking if  0 <SE< inf  ...
    Wed Dec  7 01:31:38 2022  -Removed 0 variants with bad SE.
    Wed Dec  7 01:31:38 2022  -Checking STATUS...
    Wed Dec  7 01:31:40 2022  -Coverting STAUTUS to interger.
    Wed Dec  7 01:31:42 2022  -Removed 0 variants with bad statistics in total.
    Wed Dec  7 01:31:42 2022 Finished sanity check successfully!
    ```

## Remove duplicated or multiallelic variants

After standardizing the sumstats, you can also remove duplicated or multiallelic variants using:

- `mode=d` ,remove duplicate.
    - remove duplicate SNPs based on  1. SNPID, 
    - remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
    - remove duplicate SNPs based on  3. rsID
- `mode=m`, remove multiallelic variants.
    - remove multiallelic SNPs based on  4. CHR, POS
- `remove=True` : remove NAs 
- `keep_col` : use which column to sort the values (`keep_ascend=True`: ascending order)
- `keep`: keep 'first' or 'last'.

!!! example
    ```python
    sumstats.remove_dup(mode="md",keep='first',keep_col="P",remove=False)
    ```

-> remove duplicated and multiallelic variants and keep the one with lowest P.


## Filtering by condition

Filtering sumstats by `expr` (a wrapper of `pd.query`), return a new Sumstats Object by default. 

```
.filter_value(expr,inplace=False)
```

- `filter_value`, is a wrapper of pd.DataFrame query. This can be used for value filtering.
- `expr`: the conditions used fot filtering. (for example: '1>BETA>0 & N>10000')
- `inplace`: True or False. If False, return a new  Sumstats Object (so chain calling is possible). If true, the current Sumstats Object will be filtered in place.

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
 

## Filtering regions
```
.filter_region_in(path="./abc.bed")
.filter_region_out(high_ld=True)
```
filter variants in certain regions defined in bed files.

- `path` : path to the bed files
- `high_ld` : high-ld region (built-in data, no additional bed needed)

!!! example
    ```
    mysumstats.filter_region_in(high_ld=True,inplace=False).data
    mysumstats.filter_region_out(high_ld=True,inplace=False).data
    ```
   
## filter_in & filter_out (deprecated)

```python
.filter_in(gt={},lt={},eq={},inplace=False)

.filter_out(gt={},lt={},eq={},inplace=False)
```
- `gt`: greater than
- `lt`: less than
- `eq`: equal to
- `inplace`: True or False. If False, return a dataframe. If true, the Sumstats object will be filtered.
