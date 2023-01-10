# Standardization and Normalization
After loading raw sumstats into gwaslab Sumstats Object, the first thing we probably want to do is to standardize the variant-related notations, normalize the indels and check if there are any unexpected errors in the statistics. When checking is finished, the status code will be automatically changed.

See examples [here](https://cloufield.github.io/gwaslab/standardization_workflow/).

## Methods Summary

| Sumstats Methods      | Options                                                      | Description                                                                    |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| `.fix_id()`           | `fixchrpos=False`, <br/>`fixid=False`, <br/>`fixsep=False`,<br/>`overwrite=False`,<br/>`forcefixid=False` | check and  fix rsID or SNPID(chr:pos:ref:alt), or use snpid to fix CHR and POS |
| `.fix_CHR()`          | `remove=False`                                               | standardize chromsome notation                                                 |
| `.fix_POS()`          | `remove=False` , `limit=250000000`          | standardize basepair posituion notation and filter out bad values              |
| `.fix_allele()`       | `remove=False`, `x="X"`, `y="Y"`, `mt="MT"`         | standardize base notation to ATCG                                              |
| `.normalize_allele()` | `n_cores=1`                                                  | normalize indels (only support ATA:AA -> AT:A but not -:T)                     |
| `.sort_coordinate()`  |                                                              | sort the variant coordinates                                                   |
| `.sort_column()`  |                                                              | sort the column order to gwaslab default                                                   |
| `.basic_check()`  |                                                              | all-in-one function to perform the defaul pipeline                                             |

## 1. IDs

Gwaslab requires at least one ID columns for sumstats, either in the form of SNPID or rsID, (or both). Gwaslab will automatically check if SNPID is mixed in rsID.

- `SNPID` : it could be user provided IDs, or in `CHR:POS:REF:ALT` format, delimiter can be":","_" or "-"
- `rsID`: dbSNP rsIDs

Gwaslab will check if the IDs is valid SNPID orrsID.  It can also extract CHR and POS information from the CHR:POS:REF:ALT formatted IDs using `.fix_id()` method.

SNPID will be fixed by `CHR:POS:NEA:EA`  only when the variants is already aligned with reference genome. Otherwise, a temporary SNPID in the format of `CHR:POS` will be given.

`.fix_id()` : check or fix SNPID and rsID.   

|`.fix_id()` options|Type|Description|
|-|-|-|
|`fixchrpos`|`Boolean`|If True, extract CHR and POS from SNPID (CHR:POS:NEA:EA) to fill CHR and POS columns (deaful:`False`)|
|`fixid`|`Boolean`|If True, use CHR/POS/NEA/EA to reconstruct the SNPID. For variant that are not aligned with reference genome. Only CHR/POS will be used. (deaful:`False`)|
|`forcefixid`|`Boolean`|If True, use CHR/POS/NEA/EA to reconstruct the SNPID without checking if the variant is aligned. (deaful:`False`)|
|`fixsep`|`Boolean`|If True, fix SNPID delimiter (For example: 1:123_A_C to 1:123:A:C) (deaful:`False`)|
|`overwrite`|`Boolean`|If True, overwrite existing data. (deaful:`False`)|

!!! example
    ```python
    sumstats.fix_id(fixchrpos=False,
                    fixid=False,
                    fixsep=False,
                    forcefixid=False,
                    overwrite=False)
    ```

## 2. CHR

`.fix_CHR()`

CHR will be standardized to integers. For human chromosomes, CHR will be converted to `1-22` for autosomes, `23` for `X`, `24` for `Y`,and `25` for `MT`.

-  autosomes: 1-22
- `x="X"` : 23
- `y="Y"` : 24
- `mt="MT"` : 25

Leading "chr" or leading 0s will be stripped. Before conversion, CHR will be converted to uppercase.

!!! example
    ```python
    sumstats.fix_CHR(remove=False,x="X",y="Y",mt="MT")
    ```

## 3. POS

`.fix_POS()`

Values in POS must be positive integer numbers. Basepair position will be force converted to integers. Invalid pos will be converted to NA. 

`limit=250000000` : After conversion, gwaslab will also sanity check if POS is in the range of 1 to 250,000,000. (the longest chromosome, CHR1, is around 250,000,000bp long)  

!!! example
    ```python
    sumstats.fix_POS(remove=False, limit=250000000)
    ```

## 4. Allele

### 4.1 Allele notation Standardization

`.fix_allele()`

- Currently, gwaslab only support processing SNPs and INDELs. 
- All alleles will be checked if containing letters other than `A`,`T`,`C`,`G`.
- Copy number variant (CNV) like `<CN0>` won't be recognized.
- Lower cases will converted to UPPERCASES.

!!! example
    ```python
    sumstats.fix_allele(remove=False)
    ```

### 4.2 Variant Normalization

`.normalize_allele()`

Alleles will be normalized accroding to left alignment and parsimony principal. For example, chr1:123456:ATG:AT will be normalized to chr1:123455:TG:T.

`n_cores` : threads to use for normalization.

!!! note
    Currently, the normalizeation is implemented without checking reference, which means it can not normalize variants like chr1:123456:G:- if the missing information need to be obtained from a reference genome.  
    
!!! quote
    For details on variant normalization, please check: [https://genome.sph.umich.edu/wiki/Variant_Normalization](https://genome.sph.umich.edu/wiki/Variant_Normalization) )

!!! example
    ```python
    sumstats.normalize_allele(n_cores=1)
    ```

## 5. Coordinate sorting

Sort genomic coordinates， 1-25 (23:X, 24:Y ,25:MT)

!!! example
    ```python
    sumstats.sort_coordinate()
    ```

## 6. Column sorting

The default column order is `"SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "Z",  "CHISQ", "P", "MLOG10P", "OR", "OR_95L", "OR_95U", "INFO", "N","DIRECTION","STATUS" and other additional columns`.

`order`:`list` , column order.

!!! example
    ```python
    sumstats.sort_column()
    ```

## Example
Please check [https://cloufield.github.io/gwaslab/standardization_workflow/](https://cloufield.github.io/gwaslab/standardization_workflow/)
