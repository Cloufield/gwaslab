# Standardization and Normalization
After loading raw sumstats into gwaslab.Sumstats Object, the first thing we probably want to do is to standardize the variant-related notations, normalize the indels and check if there are any unexpected errors in the statistics. When checking is finished, the status code will be automatically changed.

See examples [here](https://cloufield.github.io/gwaslab/standardization_workflow/).

## Methods Summary

| Sumstats Methods      | Options                                                      | Description                                                                    |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| `.fix_id()`           | `fixchrpos=False`, <br/>`fixid=False`, <br/>`fixsep=False`,<br/>`overwrite=False`,<br/>`forcefixid=False` | check and  fix rsID or SNPID(chr:pos:ref:alt), or use snpid to fix CHR and POS |
| `.fix_chr()`          | `remove=False`, `x=("X",23)`, `y=("Y",24)`, `mt=("MT",25)`, `chrom_list = gl.get_chr_list()` | standardize chromsome notation                                                 |
| `.fix_pos()`          | `remove=False` , `limit=250000000`          | standardize basepair posituion notation and filter out bad values              |
| `.fix_allele()`       | `remove=False`        | standardize base notation to ATCG                                              |
| `.normalize_allele()` | `n_cores=1`                                                  | normalize indels (only support ATA:AA -> AT:A but not -:T)                     |
| `.sort_coordinate()`  |                                                              | sort the variant coordinates                                                   |
| `.sort_column()`  |                                                              | sort the column order to GWASLab default                                                   |
| `.basic_check()`  |                                                              | all-in-one function to perform the defaul pipeline                                             |

## 1. IDs

GWASLab requires at least one ID column for sumstats, either in the form of SNPID or rsID, (or both). GWASLab will automatically check if SNPID is mixed in rsID.

- `SNPID` : it could be user-provided IDs, or in `CHR:POS:REF:ALT` format, delimiters can be `":"`,`"_"` or `"-"`.
- `rsID`: dbSNP rsIDs

GWASLab will check if the IDs is valid SNPID or rsID.  It can also extract CHR and POS information from the `CHR:POS:REF:ALT` formatted IDs using `.fix_id()` method.

SNPID will be fixed by `CHR:POS:NEA:EA`  only when the variants are already aligned with the reference genome. Otherwise, a temporary SNPID in the format of `CHR:POS` will be assigned.

`.fix_id()`: check or fix SNPID and rsID.   

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

```
.fix_chr(
    x=("X",23), 
    y=("Y",24), 
    mt=("MT",25),
    chrom_list=gl.get_chr_list()
)
```

CHR will be standardized to integers. 

- `"x","y","mt"`: `tuple` , how to convert sex chromosomes. For examplem `x = ("X",23)`
- `chrom_list` : `list` , vaild chromosome notations for filtering. Datatype should be `string` or `object`.

(by default) For human chromosomes, CHR will be converted to integers for autosomes, `23` for `X`, `24` for `Y`, and `25` for `MT`. 

- `x=("X",23)` 
- `y=("Y",24)`  
- `mt=("MT",25)`

For other species, you can change it like:

- `x=("X",26)` 
- `y=("Y",27)`  
- `mt=("MT",28)`
- `chrom_list=gl.get_chr_list(n=28)`

!!! example
    ```python
    sumstats.fix_chr(x=("X",23))
    ```

!!! tip gl.get_chr_list()
    Get a chromosome list for n autosomes and 'X', 'Y', 'M', 'MT. ("string" data type)
    ```
    gl.get_chr_list(n=10)
    ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'X', 'Y', 'M', 'MT']
    ```

## 3. POS

`.fix_pos()`

Values in POS must be positive integer numbers. Basepair position will be force converted to integers. Invalid pos will be converted to NA. 

`limit=250000000` : After conversion, GWASLab will also sanity check if POS is in the range of 1 to 250,000,000. (the longest chromosome, CHR1, is around 250,000,000bp long)  

!!! example
    ```python
    sumstats.fix_pos(remove=False, limit=250000000)
    ```

## 4. Allele

### 4.1 Allele notation Standardization

`.fix_allele()`

- Currently, GWASLab only supports processing SNPs and INDELs. 
- All alleles will be checked if containing letters other than `A`,`T`,`C`,`G`.
- Copy number variant (CNV) like `<CN0>` won't be recognized.
- Lower cases will be converted to UPPERCASES.

!!! example
    ```python
    sumstats.fix_allele(remove=False)
    ```

### 4.2 Variant Normalization

`.normalize_allele()`

Alleles will be normalized according to the left alignment and parsimony principle. For example, chr1:123456:ATG:AT will be normalized to chr1:123457:TG:T.

`n_cores`: threads to use for normalization.

!!! note
    Currently, the normalization is implemented without checking reference, which means it can not normalize variants like chr1:123456:G:- if the missing information needs to be obtained from a reference genome.  
    
!!! quote
    For details on variant normalization, please check: [https://genome.sph.umich.edu/wiki/Variant_Normalization](https://genome.sph.umich.edu/wiki/Variant_Normalization) )

!!! example

    ```python
    sumstats.normalize_allele(n_cores=1)
    ```
    
    Before:
    
    <img width="345" alt="image" src="https://user-images.githubusercontent.com/40289485/212257048-1a5517a4-dd4d-4210-9e19-1dcf7747b7f5.png">

    After:
    
    <img width="345" alt="image" src="https://user-images.githubusercontent.com/40289485/212256576-7808a1ec-5cf2-42ec-819e-a6f1c9c200bb.png">
    

## 5. Coordinate sorting

Sort genomic coordinates. Make sure CHR and POS are fixed beforehand.

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

## All-in-one `.basic_check()`

```
.basic_check()
```
Fix data without dropping any variant.

!!! note "`basic_check()` pipeline"
    
    - fix_id()
    - fix_chr()
    - fix_pos()
    - fix_allele()
    - normalize_allele()
    - sort_coordinate()
    - sort_column()

## Example
Please check [https://cloufield.github.io/gwaslab/standardization_workflow/](https://cloufield.github.io/gwaslab/standardization_workflow/)
