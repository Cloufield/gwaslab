# Standardization and Normalization
After loading raw sumstats into gwaslab.Sumstats Object, the first thing we probably want to do is to standardize the variant-related notations, normalize the indels and check if there are any unexpected errors in the statistics. When checking is finished, the status code will be automatically changed.

## 1. Methods Summary

| Sumstats Methods      | Options                                                      | Description                                                                    |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| `.fix_id()`           | `fixchrpos=False`, <br/>`fixid=False`, <br/>`fixsep=False`,<br/>`overwrite=False`,<br/>`forcefixid=False` | check and  fix rsID or SNPID(chr:pos:ref:alt), or use snpid to fix CHR and POS |
| `.fix_chr()`          | `remove=False`, `x=("X",23)`, `y=("Y",24)`, `mt=("MT",25)`, `chrom_list = gl.get_chr_list()` | standardize chromsome notation                                                 |
| `.fix_pos()`          | `remove=False` , `limit=250000000`          | standardize basepair posituion notation and filter out bad values              |
| `.fix_allele()`       | `remove=False`        | standardize base notations to ATCG                                              |
| `.normalize_allele()` | `n_cores=1`                                                  | normalize indels (only support ATA:AA -> AT:A but not -:T)                     |
| `.sort_coordinate()`  |                                                              | sort the variant coordinates                                                   |
| `.sort_column()`  |     `order`                                                         | sort the column order to GWASLab default                                                   |
| `.basic_check()`  |                                                              | all-in-one function to perform the default pipeline                                             |

## 2. IDs - SNPID and rsID

GWASLab requires at least one ID column for sumstats, either in the form of SNPID or rsID, (or both). GWASLab will automatically check if SNPID is mixed in rsID.

- `SNPID` : it could be user-provided IDs, or in `CHR:POS:REF:ALT` format, delimiters can be `":"`,`"_"` or `"-"`.
- `rsID`: dbSNP rsIDs

GWASLab will check if the IDs are valid SNPID or rsID.  It can also extract CHR and POS information from the `CHR:POS:REF:ALT` formatted IDs using `.fix_id()` method.

SNPID will be fixed by `CHR:POS:NEA:EA`  only when the variants are already aligned with the reference genome. Otherwise, a temporary SNPID in the format of `CHR:POS` will be assigned.

`.fix_id()`: check or fix SNPID and rsID.   

|`.fix_id()` options|Type|Description|
|-|-|-|
|`fixchrpos`|`Boolean`|If True, extract CHR and POS from SNPID (CHR:POS:NEA:EA) to fill CHR and POS columns (default:`False`)|
|`fixid`|`Boolean`|If True, use CHR/POS/NEA/EA to reconstruct the SNPID. For variant that are not aligned with reference genome. Only CHR/POS will be used. (default:`False`)|
|`forcefixid`|`Boolean`|If True, use CHR/POS/NEA/EA to reconstruct the SNPID without checking if the variant is aligned. (default:`False`)|
|`fixsep`|`Boolean`|If True, fix SNPID delimiter (For example: 1:123_A_C to 1:123:A:C) (default:`False`)|
|`overwrite`|`Boolean`|If True, overwrite existing data. (default:`False`)|

!!! example
    ```python
    sumstats.fix_id(fixchrpos=False,
                    fixid=False,
                    fixsep=False,
                    forcefixid=False,
                    overwrite=False)
    ```

## 3. CHR - Chromosomes

```
.fix_chr()
```

CHR will be standardized. 

- Prefixes like "chr" will be removed.
- String datatype will be converted to integers.
- Sex chromosomes will be converted to integers.
- Variants with unrecognized chromosome notations will be removed. (notations not in `chrom_list`)
- Variants with CHR<=0 will be removed. 

Options:

|`.fix_chr()` options|DataType|Description|Default|
|-|-|-|-|
|`x`|`tuple`|how to convert X chromosomes. For examplem `x = ("X",23)`|`("X",23)`|
|`y`|`tuple`|how to convert Y chromosomes.|`("Y",24)`|
|`mt`|`tuple`|how to convert Mitochondrial DNA.|`("MT",25)`|
|`chrom_list`|`list`|vaild chromosome notations for filtering. Datatype should be `string`|`gl.get_chr_list()`|

(by default) For human chromosomes, CHR will be converted to integers for autosomes, `23` for `X`, `24` for `Y`, and `25` for `MT`. 

For other species, you can change it like:
- `x=("X",26)` 
- `y=("Y",27)`  
- `mt=("MT",28)`
- `chrom_list=gl.get_chr_list(n=28)`

!!! example
    ```python
    sumstats.fix_chr()
    ```

!!! tip gl.get_chr_list()
    Get a chromosome list for n autosomes plus 'X', 'Y', 'M', 'MT. (`string`` data type)
    
    ```python
    gl.get_chr_list()
    ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','X','Y','M','MT']

    gl.get_chr_list(n=10)
    ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'X', 'Y', 'M', 'MT']
    ```

## 4. POS - Base-pair positions

```
.fix_pos()
```

Check and fix values in POS.

- Values in POS must be positive integer numbers. 
- Basepair position will be force converted to integers. 
- Invalid POS values will be converted to NA. 
- POS outliers will be removed.

|`.fix_pos()` options|DataType|Description|Default|
|-|-|-|-|
|`limit`|`integer`|After conversion, GWASLab will also perform a sanity check for POS (if is in the range of 1 to 250,000,000). (The longest chromosome for human, namely chromosome 1, is less than 250,000,000bp long) |`250000000`|

!!! example
    ```python
    sumstats.fix_pos()
    ```

## 5. Allele

### 5.1 Allele notation Standardization

`.fix_allele()`

- Currently, GWASLab only supports processing SNPs and INDELs. 
- All alleles will be checked if containing letters other than `A`,`T`,`C`,`G`.
- Copy number variants (CNV) like `<CN0>` won't be recognized.
- Lower cases will be converted to UPPERCASES.

!!! example
    ```python
    sumstats.fix_allele()
    ```

### 5.2 Variant Normalization

`.normalize_allele()`

Alleles will be normalized according to the left alignment and parsimony principle. For example, chr1:123456:ATG:AT will be normalized to chr1:123457:TG:T.

|`.normalize_allele()` options|DataType|Description|Default|
|-|-|-|-|
|`n_cores`|`integer`|threads to use for normalization.|`1`| 

!!! warning
    Currently, the normalization is implemented without checking reference, which means it can not normalize variants like chr1:123456:G:- if the missing allele information needs to be obtained from a reference genome.  
    
!!! quote "Variant Normalization"
    For details on variant normalization, please check: [https://genome.sph.umich.edu/wiki/Variant_Normalization](https://genome.sph.umich.edu/wiki/Variant_Normalization)

!!! example

    ```python
    sumstats.normalize_allele(n_cores=1)
    ```
    
    Before:
    
    <img width="345" alt="image" src="https://user-images.githubusercontent.com/40289485/212257048-1a5517a4-dd4d-4210-9e19-1dcf7747b7f5.png">

    After:
    
    <img width="345" alt="image" src="https://user-images.githubusercontent.com/40289485/212256576-7808a1ec-5cf2-42ec-819e-a6f1c9c200bb.png">
    

## 6. Genome coordinate sorting

Sort genomic coordinates. Make sure CHR and POS are fixed beforehand.

!!! example
    ```python
    sumstats.sort_coordinate()
    ```

## 7. Column sorting

```
.sort_column()
```

|`.sort_column()` options|DataType|Description|Default|
|-|-|-|-|
|`order`|`list`|sort the columns|The default column order is `"SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "MAF", "BETA", "SE","BETA_95L","BETA_95U", "Z", "CHISQ", "P", "MLOG10P", "OR", "OR_95L", "OR_95U","HR", "HR_95L", "HR_95U","INFO", "N","N_CASE","N_CONTROL","DIRECTION","I2","P_HET","DOF","SNPR2","STATUS"` and other additional columns.| 

!!! example
    ```python
    sumstats.sort_column()
    ```

## 8. All-in-one `.basic_check()`

```
.basic_check()
```
Fix data without dropping any variant.

!!! note "`basic_check()` pipeline"
    
    - fix_id()
    - remove_dup()
    - fix_chr()
    - fix_pos()
    - fix_allele()
    - normalize_allele()
    - check_sanity()  
    - sort_coordinate()
    - sort_column()

For `remove_dup()` and `check_sanity()`, please check [QC and filtering](https://cloufield.github.io/gwaslab/QC%26Filtering/)

|`.basic_check()` options|DataType|Description|Default|
|-|-|-|-|
|`remove`|`boolean`|if True, remove duplicate and multi-allelic variants|`False`| 
|`removedup_args`|`dict`|options for remove_dup()| `{}`| 
|`fixid_args`|`dict`|options for `.fix_id()`|`{}`| 
|`fixchr_agrs`|`dict`|options for `.fix_chr()`|`{}`| 
|`fixpos_args`|`dict`|options for `.fix_pos()`|`{}`| 
|`fixallele_args`|`dict`|options for `.fix_allele()`|`{}`| 
|`sanitycheckstats_args`|`dict`|options for `.check_sanity()`|`{}`|
|`normalizeallele_args`|`dict`|options for `.normalize_allele()` |`{}`| 
|`verbose`|`boolean`|if True, print log|`True`| 

