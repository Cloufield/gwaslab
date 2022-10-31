# Standardization and normalization

```python
import gwaslab as gl
sumstats = gl.Sumstats(...)
```

After loading raw sumstats into gwaslab Sumstats Object, the first thing we probably want to do is to standardize the variant-related notations and check if there are any unexpected errors in the statistics. When checking is finished, the status code will be automatically changed.

## Methods Summary

| Sumstats Methods      | Options                                                      | Description                                                                    |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| `.fix_ID()`           | `fixchrpos=False`, <br/>`fixid=False`, <br/>`fixsep=False`,<br/>`overwrite=False`,<br/>`forcefixid=False` | check and  fix rsID or SNPID(chr:pos:ref:alt), or use snpid to fix CHR and POS |
| `.fix_CHR()`          | `remove=False`                                               | standardize chromsome notation                                                 |
| `.fix_POS()`          | `remove=False`                                               | standardize basepair posituion notation and filter out bad values              |
| `.fix_allele()`       | `remove=False`                                               | standardize base notation to ATCG                                              |
| `.normalize_allele()` | `n_cores=1`                                                  | normalize indels (only support ATA:AA -> AT:A but not -:T)                     |
| `.sort_coordinate()`  |                                                              | sort the variant coordinates                                                   |

## 1. IDs

Gwaslab requires at least one ID columns for sumstats, either in the form of SNPID or rsID, (or both). Gwaslab will automatically check if SNPID is mixed in rsID.

 `SNPID` : it could be user provided IDs, or in `CHR:POS:REF:ALT` format, delimiter can be":","_" or "-"

 `rsID`: dbSNP rsIDs

gwaslab checks if the IDs you provided is valid SNPID or rsID.  It can also extract CHR and POS information from the CHR:POS:REF:ALT formatted IDs using `.fix_ID()` method.

SNPID will be fixed by `CHR:POS:NEA:EA`  only when the variants is already aligned with reference genome. Otherwise, a temporary SNPID in the format of `CHR:POS` will be given.

`.fix_ID()` : check or fix SNPID and rsID.   

```python
sumstats.fixID(fixchrpos=False,
               fixid=False,
               fixsep=False,
               overwrite=False)
```
`.fix_ID()` options:

- `fixchrpos` : `Boolean`, extract CHR and POS from SNPID (CHR:POS:NEA:EA) to fill CHR and POS columns (deaful:`False`)
- `fixid` : `Boolean` ,  use CHR/POS/NEA/EA to reconstruct the SNPID. For variant that are not aligned with reference genome. Only CHR/POS will be used. (deaful:`False`)
- `forcefixid`: `Boolean` ,  use CHR/POS/NEA/EA to reconstruct the SNPID without checking if the variant is aligned. (deaful:`False`)
- `fixsep` : `Boolean`, fix SNPID delimiter (For example: 1:123_A_C to 1:123:A:C) (deaful:`False`)
- `overwrite` : `Boolean`, whether overwrite existing data. (deaful:`False`)

## 2. CHR

`.fix_CHR()`

CHR will be standardized to `1-22,X,Y,MT`

Leading "chr" and leading 0s will be stripped.

```python
sumstats.fix_CHR(remove=False)
```

## 3. POS

`.fix_POS()`

Values in POS must be positive integer numbers.

Basepair position will be force converted to integers.

Invalid pos will be converted to NA. 

(not implemented yet) After conversion, gwaslab will also sanity check if POS is in the range of 1 to 300,000,000. (the longest chromosome, CHR1, is around 250,000,000bp long)  

```python
sumstats.fix_POS(remove=False)
```

## 4. Allele

### 4.1 Standardization

`.fix_allele()`

Currently, gwaslab only support processing SNPs and INDELs.

All alleles will be checked if containing letters other than `ATCG`.

Copy number variant (CNV) like `<CN0>` won't be recognized.

Lower cases will converted to UPPERCASES.

```python
sumstats.fix_allele(remove=False)
```

### 4.2 Normalization

`.normalize_allele()`

Alleles will be normalized accroding to left alignment and parsimony principal. (For details: add link here )

For example, chr1:123456:ATG:AT will be normalized to chr1:123455:TG:T.

Note: Currently, the normalizeation is implemented without checking reference, which means it can not normalize variants like chr1:123456:G:- if the missing information need to be obtained from a reference genome. 

```python
sumstats.normalize_allele(n_cores=1)
```

## 5. Remove duplicated or multiallelic variants

- remove duplicate SNPs based on  1. SNPID, 
- remove duplicate SNPs based on  2. CHR, POS, EA, and NEA
- remove duplicate SNPs based on  3. rsID for non-NA variants
- remove multiallelic SNPs based on  4. CHR, POS

```
sumstats.remove_dup()
```

## 5. Coordinate sorting

Sort genomic coordinates， 1-22 X Y MT

```python
sumstats.sort_coordinate()
```

## 6. Column sorting

The default column order is "SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "Z",  "CHISQ", "P", "MLOG10P", "OR", "OR_SE", "OR_95L", "OR_95U", "INFO", "N","DIRECTION","STATUS" and other additional columns.

```
sumstats.sort_columns()
```
