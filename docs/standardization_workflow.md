# Standardization Workflow

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2025/12/25 21:54:45 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:54:45 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:54:45 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Load sample data

```python
mysumstats = gl.Sumstats("../0_sample_data/toy_data/dirty_sumstats.tsv",fmt="gwaslab",other=["NOTE"])
```

**stdout:**
```
2025/12/25 21:54:45 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:54:45 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:54:45 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/12/25 21:54:45 Start to load format from formatbook....
2025/12/25 21:54:45  -gwaslab format meta info:
2025/12/25 21:54:45   - format_name  : gwaslab
2025/12/25 21:54:45   - format_source  : https://cloufield.github.io/gwaslab/
2025/12/25 21:54:45   - format_version  : 20231220_v4
2025/12/25 21:54:45 Start to initialize gl.Sumstats from file :../0_sample_data/toy_data/dirty_sumstats.tsv
2025/12/25 21:54:45  -Reading columns          : NOTE,NEA,POS,N_CASE,N_CONTROL,DIRECTION,MLOG10P,Z,OR_95U,OR,P,SE,SNPID,EAF,EA,BETA,N,CHISQ,CHR,OR_95L
2025/12/25 21:54:45  -Renaming columns to      : NOTE,NEA,POS,N_CASE,N_CONTROL,DIRECTION,MLOG10P,Z,OR_95U,OR,P,SE,SNPID,EAF,EA,BETA,N,CHISQ,CHR,OR_95L
2025/12/25 21:54:45  -Current Dataframe shape : 63  x  20
2025/12/25 21:54:45  -Initiating a status column: STATUS ...
2025/12/25 21:54:45 #WARNING! Version of genomic coordinates is unknown...
2025/12/25 21:54:45 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:54:45  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,CHISQ,Z,P,MLOG10P,DIRECTION,N,N_CASE,N_CONTROL,NOTE
2025/12/25 21:54:45 Finished reordering the columns.
2025/12/25 21:54:45  -Trying to convert datatype for CHR: string -> Int64...Failed...
2025/12/25 21:54:45  -Trying to convert datatype for POS: object -> int64...Failed...
2025/12/25 21:54:45  -Trying to convert datatype for N: float64 -> Int64...Failed...
2025/12/25 21:54:45  -Column  : SNPID  CHR    POS    EA       NEA      STATUS EAF     BETA    SE      OR      OR_95U  OR_95L  CHISQ   Z       P       MLOG10P DIRECTION N       N_CASE N_CONTROL NOTE  
2025/12/25 21:54:45  -DType   : object string object category category int64  float64 float64 float64 float64 float64 float64 float64 float64 float64 float64 object    float64 int64  int64     object
2025/12/25 21:54:45  -Verified: T      F      F      T        T        T      T       T       T       T       T       T       T       T       T       T       T         F       T      T         NA    
2025/12/25 21:54:45 #WARNING! Columns with possibly incompatible dtypes: CHR,POS,N
2025/12/25 21:54:45 #WARNING! Consider using Sumstats.fix_chr() to fix CHR dtype
2025/12/25 21:54:45 #WARNING! Consider using Sumstats.fix_pos() to fix POS dtype
2025/12/25 21:54:45  -Current Dataframe memory usage: 0.01 MB
2025/12/25 21:54:45 Finished loading data successfully!
```

Dirty sumstats with issues specified in NOTE column

```python
mysumstats.data
```

```
| SNPID CHR POS | EA | NEA | STATUS | EAF | BETA | SE | OR | ... | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1_G_A | 1 | 1 | A | G | 9999999 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:1_A_G | 1 | 1 | A | G | 9999999 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:1_A_G | 1 | 1 | A | G | 9999999 | 0.996 | 0.0603 | 0.0103 |
| 3 | 1:2 | 1 | 2 | T | C | 9999999 | 0.996 | 0.0603 | 0.0103 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9999999 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 58 | 160000.0 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |  |
| 59 | 160000.0 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |  |
| 60 | 160000.0 | 120000 | 40000 | MLOG10P missing |  |  |  |  |  |
| 61 | 160000.0 | 120000 | 40000 | Clean sumstats |  |  |  |  |  |
| 62 | 160000.0 | 120000 | 40000 | Clean sumstats |  |  |  |  |  |

*[63 rows x 21 columns]*
```

## All in one function

```python
mysumstats.basic_check(remove=True,remove_dup=True)
```

**stdout:**
```
2025/12/25 21:54:45 Start to check SNPID/rsID ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 63 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Checking SNPID data type...
2025/12/25 21:54:46  -Converted datatype for SNPID: object -> string
2025/12/25 21:54:46  -Checking NA strings :na,NA,Na,Nan,NaN,<NA>,null,NULL,#N/A,#VALUE!,N/A,n/a,missing,
2025/12/25 21:54:46  -Checking if SNPID contains NA strings...
2025/12/25 21:54:46  -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)
2025/12/25 21:54:46 Finished checking SNPID/rsID.
2025/12/25 21:54:46 Start to fix chromosome notation (CHR) ...(v4.0.0)
2025/12/25 21:54:46  -Checking CHR data type...
2025/12/25 21:54:46  -Variants with standardized chromosome notation: 56
2025/12/25 21:54:46  -Variants with fixable chromosome notations: 4
2025/12/25 21:54:46  -Variants with NA chromosome notations: 1
2025/12/25 21:54:46  -Variants with invalid chromosome notations: 2
2025/12/25 21:54:46  -A look at invalid chromosome notations: {'1.0001', '-1'}
2025/12/25 21:54:46  -Identifying non-autosomal chromosomes : X, Y, and MT ...
2025/12/25 21:54:46  -Identified  1  variants on sex chromosomes...
2025/12/25 21:54:46  -Standardizing sex chromosome notations: X to 23...
2025/12/25 21:54:46  -Valid CHR list: 1 - 25
2025/12/25 21:54:46  -Removed variants with chromosome notations not in CHR list: 5
2025/12/25 21:54:46  -A look at chromosome notations not in CHR list: {'300', '0', <NA>}
2025/12/25 21:54:46 Finished fixing chromosome notation (CHR).
2025/12/25 21:54:46 Start to fix basepair positions (POS) ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 58 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Removing thousands separator "," or underbar "_" ...
2025/12/25 21:54:46  -Trying to convert datatype for POS: string -> Int64...
2025/12/25 21:54:46  -Trying to convert datatype for POS: string -> Int64...
2025/12/25 21:54:46  -Position bound:(0 , 250,000,000)
2025/12/25 21:54:46  -Removed variants outliers: 2
2025/12/25 21:54:46  -Removed variants with bad positions: 4
2025/12/25 21:54:46 Finished fixing basepair positions (POS).
2025/12/25 21:54:46 Start to fix alleles (EA and NEA) ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 54 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Converted all bases to string datatype and UPPERCASE
2025/12/25 21:54:46  -Variants with bad EA: 1
2025/12/25 21:54:46  -Variants with bad NEA: 5
2025/12/25 21:54:46  -Variants with NA for EA or NEA: 1
2025/12/25 21:54:46  -Variants with same EA and NEA: 1
2025/12/25 21:54:46  -A look at the non-ATCG EA: {'<CN0>'} ...
2025/12/25 21:54:46  -A look at the non-ATCG NEA: {'<CN1>', 'N', '*', nan} ...
2025/12/25 21:54:46  -Removed variants with NA alleles or alleles that contain bases other than A/C/T/G: 5
2025/12/25 21:54:46  -Removed variants with same allele for EA and NEA: 1
2025/12/25 21:54:46 Finished fixing alleles (EA and NEA).
2025/12/25 21:54:46 Start to perform sanity check for statistics ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 48 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Comparison tolerance for floats: 1e-07
2025/12/25 21:54:46  -Checking if any columns are empty...
2025/12/25 21:54:46  -Checking if 0 <= N <= 2147483647 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 25,26,27 ...
2025/12/25 21:54:46   -Examples of invalid values (N): 12345700000000000,NA,-1 ...
2025/12/25 21:54:46  -Removed 3 variants with bad/na N.
2025/12/25 21:54:46  -Checking if 0 <= N_CASE <= 2147483647 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 29 ...
2025/12/25 21:54:46   -Examples of invalid values (N_CASE): -1 ...
2025/12/25 21:54:46  -Removed 1 variants with bad/na N_CASE.
2025/12/25 21:54:46  -Checking if 0 <= N_CONTROL <= 2147483647 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 28 ...
2025/12/25 21:54:46   -Examples of invalid values (N_CONTROL): -1 ...
2025/12/25 21:54:46  -Removed 1 variants with bad/na N_CONTROL.
2025/12/25 21:54:46  -Checking if -1e-07 < EAF < 1.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 31,32,33 ...
2025/12/25 21:54:46   -Examples of invalid values (EAF): 1.02,-0.01,NA ...
2025/12/25 21:54:46  -Removed 3 variants with bad/na EAF.
2025/12/25 21:54:46  -Checking if -1e-07 < CHISQ < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 38,39 ...
2025/12/25 21:54:46   -Examples of invalid values (CHISQ): -0.01,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na CHISQ.
2025/12/25 21:54:46  -Checking if -9999.0000001 < Z < 9999.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 40,41 ...
2025/12/25 21:54:46   -Examples of invalid values (Z): NA,999999.0 ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na Z.
2025/12/25 21:54:46  -Checking if -1e-07 < P < 1.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 48,49,50 ...
2025/12/25 21:54:46   -Examples of invalid values (P): 1.1,-0.01,NA ...
2025/12/25 21:54:46  -Removed 3 variants with bad/na P.
2025/12/25 21:54:46  -Checking if -1e-07 < MLOG10P < 99999.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 52,53 ...
2025/12/25 21:54:46   -Examples of invalid values (MLOG10P): -0.1,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na MLOG10P.
2025/12/25 21:54:46  -Checking if -100.0000001 < BETA < 100.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 34,35 ...
2025/12/25 21:54:46   -Examples of invalid values (BETA): 99999.0,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na BETA.
2025/12/25 21:54:46  -Checking if -1e-07 < SE < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 37 ...
2025/12/25 21:54:46   -Examples of invalid values (SE): NA ...
2025/12/25 21:54:46  -Removed 1 variants with bad/na SE.
2025/12/25 21:54:46  -Checking if -1e-07 < OR < 100.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 42,43 ...
2025/12/25 21:54:46   -Examples of invalid values (OR): 999999.0,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na OR.
2025/12/25 21:54:46  -Checking if -1e-07 < OR_95L < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 44,45 ...
2025/12/25 21:54:46   -Examples of invalid values (OR_95L): -0.01,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na OR_95L.
2025/12/25 21:54:46  -Checking if -1e-07 < OR_95U < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 46,47 ...
2025/12/25 21:54:46   -Examples of invalid values (OR_95U): -0.01,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na OR_95U.
2025/12/25 21:54:46  -Checking STATUS and converting STATUS to Int64....
2025/12/25 21:54:46  -Removed 26 variants with bad statistics in total.
2025/12/25 21:54:46  -Data types for each column:
2025/12/25 21:54:46 Finished sanity check for statistics.
2025/12/25 21:54:46 Start to check data consistency across columns ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 22 x 21 ; Memory usage: 0.00 MB
2025/12/25 21:54:46  -Tolerance: 0.001 (Relative) and 0.001 (Absolute)
2025/12/25 21:54:46  -Checking if BETA/SE-derived-MLOG10P is consistent with MLOG10P...
2025/12/25 21:54:46   -Potentially inconsistent (likely due to rounding): 1 variant(s)
2025/12/25 21:54:46   -Variant SNPID with max difference: 1:1_G_A with -0.007524575787453358
2025/12/25 21:54:46  -Checking if BETA/SE-derived-P is consistent with P...
2025/12/25 21:54:46   -Variants with inconsistent values were not detected.
2025/12/25 21:54:46  -Checking if MLOG10P-derived-P is consistent with P...
2025/12/25 21:54:46   -Variants with inconsistent values were not detected.
2025/12/25 21:54:46  -Checking if N is consistent with N_CASE + N_CONTROL ...
2025/12/25 21:54:46   -Potentially inconsistent: 1 variant(s)
2025/12/25 21:54:46   -Variant SNPID with max difference: 30 with 10000
2025/12/25 21:54:46  -Note: Minor differences are typically due to rounding.
2025/12/25 21:54:46   -If the max difference is greater than expected, please check your original sumstats.
2025/12/25 21:54:46 Finished checking data consistency across columns.
2025/12/25 21:54:46 Start to normalize indels ...(v4.0.0)
2025/12/25 21:54:46  -Number of variants to check:1
2025/12/25 21:54:46  -Chunk size:3000000
2025/12/25 21:54:46  -Processing in chunks:0 
2025/12/25 21:54:46 Finished normalizing indels.
2025/12/25 21:54:46 Start to remove duplicated/multiallelic variants ...(v4.0.0)
2025/12/25 21:54:46  -Removing mode:dm
2025/12/25 21:54:46 Start to sort the sumstats using P ...
2025/12/25 21:54:46 Start to remove duplicated variants based on snpid ...(v4.0.0)
2025/12/25 21:54:46  -Which variant to keep:  first
2025/12/25 21:54:46  -Removed variants based on SNPID: 2
2025/12/25 21:54:46 Start to remove duplicated variants based on CHR,POS,EA and NEA ...
2025/12/25 21:54:46  -Current Dataframe shape : 20 x 21 ; Memory usage: 0.00 MB
2025/12/25 21:54:46  -Which variant to keep:  first
2025/12/25 21:54:46  -Removed variants based on CHR,POS,EA and NEA: 1
2025/12/25 21:54:46 Start to remove multiallelic variants based on chr:pos ...
2025/12/25 21:54:46  -Current Dataframe shape : 19 x 21 ; Memory usage: 0.00 MB
2025/12/25 21:54:46  -Which variant to keep:  first
2025/12/25 21:54:46  -Removed variants multiallelic variants: 1
2025/12/25 21:54:46  -Removed variants in total: 4
2025/12/25 21:54:46  -Sort the coordinates based on CHR and POS...
2025/12/25 21:54:46 Finished removing duplicated/multiallelic variants.
2025/12/25 21:54:46 Start to sort the genome coordinates ...(v4.0.0)
2025/12/25 21:54:46 Finished sorting coordinates.
2025/12/25 21:54:46 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:54:46  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,CHISQ,Z,P,MLOG10P,DIRECTION,N,N_CASE,N_CONTROL,NOTE
2025/12/25 21:54:46 Finished reordering the columns.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fef521b7980>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1_G_A | 1 | 1 | A | G | 9960099 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:2 | 1 | 2 | T | TAA | 9980399 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:3_T_GA | 1 | 3 | T | GA | 9960399 | 0.996 | 0.0603 | 0.0103 |
| 3 | 3 | 1 | 6 | A | G | 9980099 | 0.996 | 0.0603 | 0.0103 |
| 4 | 4 | 1 | 7 | A | G | 9980099 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 13 | --++ | 160000 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |
| 14 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |
| 15 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |
| 16 | --++ | 160000 | 120000 | 40000 | POS with separator |  |  |  |  |
| 17 | --++ | 160000 | 120000 | 40000 | Sex chromosomes |  |  |  |  |

*[18 rows x 21 columns]*
```

## Separate functions

```python
#reload
mysumstats = gl.Sumstats("../0_sample_data/toy_data/dirty_sumstats.tsv",fmt="gwaslab",other=["NOTE"])
```

**stdout:**
```
2025/12/25 21:54:46 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:54:46 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:54:46 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/12/25 21:54:46 Start to load format from formatbook....
2025/12/25 21:54:46  -gwaslab format meta info:
2025/12/25 21:54:46   - format_name  : gwaslab
2025/12/25 21:54:46   - format_source  : https://cloufield.github.io/gwaslab/
2025/12/25 21:54:46   - format_version  : 20231220_v4
2025/12/25 21:54:46 Start to initialize gl.Sumstats from file :../0_sample_data/toy_data/dirty_sumstats.tsv
2025/12/25 21:54:46  -Reading columns          : NOTE,NEA,POS,N_CASE,N_CONTROL,DIRECTION,MLOG10P,Z,OR_95U,OR,P,SE,SNPID,EAF,EA,BETA,N,CHISQ,CHR,OR_95L
2025/12/25 21:54:46  -Renaming columns to      : NOTE,NEA,POS,N_CASE,N_CONTROL,DIRECTION,MLOG10P,Z,OR_95U,OR,P,SE,SNPID,EAF,EA,BETA,N,CHISQ,CHR,OR_95L
2025/12/25 21:54:46  -Current Dataframe shape : 63  x  20
2025/12/25 21:54:46  -Initiating a status column: STATUS ...
2025/12/25 21:54:46 #WARNING! Version of genomic coordinates is unknown...
2025/12/25 21:54:46 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:54:46  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,CHISQ,Z,P,MLOG10P,DIRECTION,N,N_CASE,N_CONTROL,NOTE
2025/12/25 21:54:46 Finished reordering the columns.
2025/12/25 21:54:46  -Trying to convert datatype for CHR: string -> Int64...Failed...
2025/12/25 21:54:46  -Trying to convert datatype for POS: object -> int64...Failed...
2025/12/25 21:54:46  -Trying to convert datatype for N: float64 -> Int64...Failed...
2025/12/25 21:54:46  -Column  : SNPID  CHR    POS    EA       NEA      STATUS EAF     BETA    SE      OR      OR_95U  OR_95L  CHISQ   Z       P       MLOG10P DIRECTION N       N_CASE N_CONTROL NOTE  
2025/12/25 21:54:46  -DType   : object string object category category int64  float64 float64 float64 float64 float64 float64 float64 float64 float64 float64 object    float64 int64  int64     object
2025/12/25 21:54:46  -Verified: T      F      F      T        T        T      T       T       T       T       T       T       T       T       T       T       T         F       T      T         NA    
2025/12/25 21:54:46 #WARNING! Columns with possibly incompatible dtypes: CHR,POS,N
2025/12/25 21:54:46 #WARNING! Consider using Sumstats.fix_chr() to fix CHR dtype
2025/12/25 21:54:46 #WARNING! Consider using Sumstats.fix_pos() to fix POS dtype
2025/12/25 21:54:46  -Current Dataframe memory usage: 0.01 MB
2025/12/25 21:54:46 Finished loading data successfully!
```

### fix id

```python
mysumstats.fix_id(fixsep=True)
```

**stdout:**
```
2025/12/25 21:54:46 Start to check SNPID/rsID ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 63 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Checking SNPID data type...
2025/12/25 21:54:46  -Converted datatype for SNPID: object -> string
2025/12/25 21:54:46  -Checking NA strings :na,NA,Na,Nan,NaN,<NA>,null,NULL,#N/A,#VALUE!,N/A,n/a,missing,
2025/12/25 21:54:46  -Checking if SNPID contains NA strings...
2025/12/25 21:54:46  -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)
2025/12/25 21:54:46  -Replacing separators in SNPID with ":" ...
2025/12/25 21:54:46 Finished checking SNPID/rsID.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID CHR POS | EA | NEA | STATUS | EAF | BETA | SE | OR | ... | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9969999 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:1:A:G | 1 | 1 | A | G | 9969999 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:1:A:G | 1 | 1 | A | G | 9969999 | 0.996 | 0.0603 | 0.0103 |
| 3 | 1:2 | 1 | 2 | T | C | 9989999 | 0.996 | 0.0603 | 0.0103 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9989999 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 58 | 160000.0 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |  |
| 59 | 160000.0 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |  |
| 60 | 160000.0 | 120000 | 40000 | MLOG10P missing |  |  |  |  |  |
| 61 | 160000.0 | 120000 | 40000 | Clean sumstats |  |  |  |  |  |
| 62 | 160000.0 | 120000 | 40000 | Clean sumstats |  |  |  |  |  |

*[63 rows x 21 columns]*
```

### fix chromosome

```python
mysumstats.fix_chr(remove=True)
```

**stdout:**
```
2025/12/25 21:54:46 Start to fix chromosome notation (CHR) ...(v4.0.0)
2025/12/25 21:54:46  -Checking CHR data type...
2025/12/25 21:54:46  -Variants with standardized chromosome notation: 56
2025/12/25 21:54:46  -Variants with fixable chromosome notations: 4
2025/12/25 21:54:46  -Variants with NA chromosome notations: 1
2025/12/25 21:54:46  -Variants with invalid chromosome notations: 2
2025/12/25 21:54:46  -A look at invalid chromosome notations: {'1.0001', '-1'}
2025/12/25 21:54:46  -Identifying non-autosomal chromosomes : X, Y, and MT ...
2025/12/25 21:54:46  -Identified  1  variants on sex chromosomes...
2025/12/25 21:54:46  -Standardizing sex chromosome notations: X to 23...
2025/12/25 21:54:46  -Valid CHR list: 1 - 25
2025/12/25 21:54:46  -Removed variants with chromosome notations not in CHR list: 5
2025/12/25 21:54:46  -A look at chromosome notations not in CHR list: {'300', '0', <NA>}
2025/12/25 21:54:46 Finished fixing chromosome notation (CHR).
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9965999 | 0.004 | 0.0603 |
| 1 | 1:1:A:G | 1 | 1 | A | G | 9965999 | 0.996 | 0.0603 |
| 2 | 1:1:A:G | 1 | 1 | A | G | 9965999 | 0.996 | 0.0603 |
| 3 | 1:2 | 1 | 2 | T | C | 9985999 | 0.996 | 0.0603 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9985999 | 0.996 | 0.0603 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 58 | MLOG10P out of range |  |  |  |  |  |  |  |
| 59 | MLOG10P out of range |  |  |  |  |  |  |  |
| 60 | MLOG10P missing |  |  |  |  |  |  |  |
| 61 | Clean sumstats |  |  |  |  |  |  |  |
| 62 | Clean sumstats |  |  |  |  |  |  |  |

*[58 rows x 21 columns]*
```

### fix position

```python
mysumstats.fix_pos(remove=True)
```

**stdout:**
```
2025/12/25 21:54:46 Start to fix basepair positions (POS) ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 58 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Removing thousands separator "," or underbar "_" ...
2025/12/25 21:54:46  -Trying to convert datatype for POS: string -> Int64...
2025/12/25 21:54:46  -Trying to convert datatype for POS: string -> Int64...
2025/12/25 21:54:46  -Position bound:(0 , 250,000,000)
2025/12/25 21:54:46  -Removed variants outliers: 2
2025/12/25 21:54:46  -Removed variants with bad positions: 4
2025/12/25 21:54:46 Finished fixing basepair positions (POS).
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9960999 | 0.004 | 0.0603 |
| 1 | 1:1:A:G | 1 | 1 | A | G | 9960999 | 0.996 | 0.0603 |
| 2 | 1:1:A:G | 1 | 1 | A | G | 9960999 | 0.996 | 0.0603 |
| 3 | 1:2 | 1 | 2 | T | C | 9980999 | 0.996 | 0.0603 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9980999 | 0.996 | 0.0603 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 58 | MLOG10P out of range |  |  |  |  |  |  |  |
| 59 | MLOG10P out of range |  |  |  |  |  |  |  |
| 60 | MLOG10P missing |  |  |  |  |  |  |  |
| 61 | Clean sumstats |  |  |  |  |  |  |  |
| 62 | Clean sumstats |  |  |  |  |  |  |  |

*[54 rows x 21 columns]*
```

### fix allele

```python
mysumstats.fix_allele(remove=True)
```

**stdout:**
```
2025/12/25 21:54:46 Start to fix alleles (EA and NEA) ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 54 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Converted all bases to string datatype and UPPERCASE
2025/12/25 21:54:46  -Variants with bad EA: 1
2025/12/25 21:54:46  -Variants with bad NEA: 5
2025/12/25 21:54:46  -Variants with NA for EA or NEA: 1
2025/12/25 21:54:46  -Variants with same EA and NEA: 1
2025/12/25 21:54:46  -A look at the non-ATCG EA: {'<CN0>'} ...
2025/12/25 21:54:46  -A look at the non-ATCG NEA: {'<CN1>', 'N', '*', nan} ...
2025/12/25 21:54:46  -Removed variants with NA alleles or alleles that contain bases other than A/C/T/G: 5
2025/12/25 21:54:46  -Removed variants with same allele for EA and NEA: 1
2025/12/25 21:54:46 Finished fixing alleles (EA and NEA).
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9960099 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:1:A:G | 1 | 1 | A | G | 9960099 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:1:A:G | 1 | 1 | A | G | 9960099 | 0.996 | 0.0603 | 0.0103 |
| 3 | 1:2 | 1 | 2 | T | C | 9980099 | 0.996 | 0.0603 | 0.0103 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9980399 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 58 | MLOG10P out of range |  |  |  |  |  |  |  |  |
| 59 | MLOG10P out of range |  |  |  |  |  |  |  |  |
| 60 | MLOG10P missing |  |  |  |  |  |  |  |  |
| 61 | Clean sumstats |  |  |  |  |  |  |  |  |
| 62 | Clean sumstats |  |  |  |  |  |  |  |  |

*[48 rows x 21 columns]*
```

### sanity check for statistics

```python
mysumstats.check_sanity()
```

**stdout:**
```
2025/12/25 21:54:46 Start to perform sanity check for statistics ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 48 x 21 ; Memory usage: 0.01 MB
2025/12/25 21:54:46  -Comparison tolerance for floats: 1e-07
2025/12/25 21:54:46  -Checking if any columns are empty...
2025/12/25 21:54:46  -Checking if 0 <= N <= 2147483647 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 25,26,27 ...
2025/12/25 21:54:46   -Examples of invalid values (N): 12345700000000000,NA,-1 ...
2025/12/25 21:54:46  -Removed 3 variants with bad/na N.
2025/12/25 21:54:46  -Checking if 0 <= N_CASE <= 2147483647 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 29 ...
2025/12/25 21:54:46   -Examples of invalid values (N_CASE): -1 ...
2025/12/25 21:54:46  -Removed 1 variants with bad/na N_CASE.
2025/12/25 21:54:46  -Checking if 0 <= N_CONTROL <= 2147483647 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 28 ...
2025/12/25 21:54:46   -Examples of invalid values (N_CONTROL): -1 ...
2025/12/25 21:54:46  -Removed 1 variants with bad/na N_CONTROL.
2025/12/25 21:54:46  -Checking if -1e-07 < EAF < 1.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 31,32,33 ...
2025/12/25 21:54:46   -Examples of invalid values (EAF): 1.02,-0.01,NA ...
2025/12/25 21:54:46  -Removed 3 variants with bad/na EAF.
2025/12/25 21:54:46  -Checking if -1e-07 < CHISQ < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 38,39 ...
2025/12/25 21:54:46   -Examples of invalid values (CHISQ): -0.01,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na CHISQ.
2025/12/25 21:54:46  -Checking if -9999.0000001 < Z < 9999.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 40,41 ...
2025/12/25 21:54:46   -Examples of invalid values (Z): NA,999999.0 ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na Z.
2025/12/25 21:54:46  -Checking if -1e-07 < P < 1.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 48,49,50 ...
2025/12/25 21:54:46   -Examples of invalid values (P): 1.1,-0.01,NA ...
2025/12/25 21:54:46  -Removed 3 variants with bad/na P.
2025/12/25 21:54:46  -Checking if -1e-07 < MLOG10P < 99999.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 52,53 ...
2025/12/25 21:54:46   -Examples of invalid values (MLOG10P): -0.1,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na MLOG10P.
2025/12/25 21:54:46  -Checking if -100.0000001 < BETA < 100.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 34,35 ...
2025/12/25 21:54:46   -Examples of invalid values (BETA): 99999.0,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na BETA.
2025/12/25 21:54:46  -Checking if -1e-07 < SE < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 37 ...
2025/12/25 21:54:46   -Examples of invalid values (SE): NA ...
2025/12/25 21:54:46  -Removed 1 variants with bad/na SE.
2025/12/25 21:54:46  -Checking if -1e-07 < OR < 100.0000001 ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 42,43 ...
2025/12/25 21:54:46   -Examples of invalid values (OR): 999999.0,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na OR.
2025/12/25 21:54:46  -Checking if -1e-07 < OR_95L < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 44,45 ...
2025/12/25 21:54:46   -Examples of invalid values (OR_95L): -0.01,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na OR_95L.
2025/12/25 21:54:46  -Checking if -1e-07 < OR_95U < inf ...
2025/12/25 21:54:46   -Examples of invalid variants(SNPID): 46,47 ...
2025/12/25 21:54:46   -Examples of invalid values (OR_95U): -0.01,NA ...
2025/12/25 21:54:46  -Removed 2 variants with bad/na OR_95U.
2025/12/25 21:54:46  -Checking STATUS and converting STATUS to Int64....
2025/12/25 21:54:46  -Removed 26 variants with bad statistics in total.
2025/12/25 21:54:46  -Data types for each column:
2025/12/25 21:54:46 Finished sanity check for statistics.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9960099 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:1:A:G | 1 | 1 | A | G | 9960099 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:1:A:G | 1 | 1 | A | G | 9960099 | 0.996 | 0.0603 | 0.0103 |
| 3 | 1:2 | 1 | 2 | T | C | 9980099 | 0.996 | 0.0603 | 0.0103 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9980399 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 37 | --++ | 160000 | 130000 | 40000 | N!=N_CONTROL +N_CASE |  |  |  |  |
| 43 | --++ | 160000 | 120000 | 40000 | SE out of range |  |  |  |  |
| 58 | --++ | 160000 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |
| 61 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |
| 62 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |

*[22 rows x 21 columns]*
```

### check data consistency

```python
mysumstats.check_data_consistency()
```

**stdout:**
```
2025/12/25 21:54:46 Start to check data consistency across columns ...(v4.0.0)
2025/12/25 21:54:46  -Current Dataframe shape : 22 x 21 ; Memory usage: 0.00 MB
2025/12/25 21:54:46  -Tolerance: 0.001 (Relative) and 0.001 (Absolute)
2025/12/25 21:54:46  -Checking if BETA/SE-derived-MLOG10P is consistent with MLOG10P...
2025/12/25 21:54:46   -Potentially inconsistent (likely due to rounding): 1 variant(s)
2025/12/25 21:54:46   -Variant SNPID with max difference: 1:1:G:A with -0.007524575787453358
2025/12/25 21:54:46  -Checking if BETA/SE-derived-P is consistent with P...
2025/12/25 21:54:46   -Variants with inconsistent values were not detected.
2025/12/25 21:54:46  -Checking if MLOG10P-derived-P is consistent with P...
2025/12/25 21:54:46   -Variants with inconsistent values were not detected.
2025/12/25 21:54:46  -Checking if N is consistent with N_CASE + N_CONTROL ...
2025/12/25 21:54:46   -Potentially inconsistent: 1 variant(s)
2025/12/25 21:54:46   -Variant SNPID with max difference: 30 with 10000
2025/12/25 21:54:46  -Note: Minor differences are typically due to rounding.
2025/12/25 21:54:46   -If the max difference is greater than expected, please check your original sumstats.
2025/12/25 21:54:46 Finished checking data consistency across columns.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

### normalize variants

```python
mysumstats.normalize_allele()
```

**stdout:**
```
2025/12/25 21:54:46 Start to normalize indels ...(v4.0.0)
2025/12/25 21:54:46  -Number of variants to check:1
2025/12/25 21:54:46  -Chunk size:3000000
2025/12/25 21:54:46  -Processing in chunks:0 
2025/12/25 21:54:46 Finished normalizing indels.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9960099 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:1:A:G | 1 | 1 | A | G | 9960099 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:1:A:G | 1 | 1 | A | G | 9960099 | 0.996 | 0.0603 | 0.0103 |
| 3 | 1:2 | 1 | 2 | T | C | 9980099 | 0.996 | 0.0603 | 0.0103 |
| 4 | 1:2 | 1 | 2 | T | TAA | 9980399 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 37 | --++ | 160000 | 130000 | 40000 | N!=N_CONTROL +N_CASE |  |  |  |  |
| 43 | --++ | 160000 | 120000 | 40000 | SE out of range |  |  |  |  |
| 58 | --++ | 160000 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |
| 61 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |
| 62 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |

*[22 rows x 21 columns]*
```

### remove duplicated / multiallelic variants

```python
mysumstats.remove_dup(mode="md")
```

**stdout:**
```
2025/12/25 21:54:46 Start to remove duplicated/multiallelic variants ...(v4.0.0)
2025/12/25 21:54:46  -Removing mode:md
2025/12/25 21:54:46 Start to sort the sumstats using P ...
2025/12/25 21:54:46 Start to remove duplicated variants based on snpid ...(v4.0.0)
2025/12/25 21:54:46  -Which variant to keep:  first
2025/12/25 21:54:46  -Removed variants based on SNPID: 2
2025/12/25 21:54:46 Start to remove duplicated variants based on CHR,POS,EA and NEA ...
2025/12/25 21:54:46  -Current Dataframe shape : 20 x 21 ; Memory usage: 0.00 MB
2025/12/25 21:54:46  -Which variant to keep:  first
2025/12/25 21:54:46  -Removed variants based on CHR,POS,EA and NEA: 1
2025/12/25 21:54:46 Start to remove multiallelic variants based on chr:pos ...
2025/12/25 21:54:46  -Current Dataframe shape : 19 x 21 ; Memory usage: 0.00 MB
2025/12/25 21:54:46  -Which variant to keep:  first
2025/12/25 21:54:46  -Removed variants multiallelic variants: 1
2025/12/25 21:54:46  -Removed variants in total: 4
2025/12/25 21:54:46  -Sort the coordinates based on CHR and POS...
2025/12/25 21:54:46 Finished removing duplicated/multiallelic variants.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:1:G:A | 1 | 1 | A | G | 9960099 | 0.004 | 0.0603 | 0.0103 |
| 1 | 1:2 | 1 | 2 | T | TAA | 9980399 | 0.996 | 0.0603 | 0.0103 |
| 2 | 1:3:T:GA | 1 | 3 | T | GA | 9960399 | 0.996 | 0.0603 | 0.0103 |
| 3 | 3 | 1 | 6 | A | G | 9980099 | 0.996 | 0.0603 | 0.0103 |
| 4 | 4 | 1 | 7 | A | G | 9980099 | 0.996 | 0.0603 | 0.0103 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 13 | --++ | 160000 | 120000 | 40000 | MLOG10P out of range |  |  |  |  |
| 14 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |
| 15 | --++ | 160000 | 120000 | 40000 | Clean sumstats |  |  |  |  |
| 16 | --++ | 160000 | 120000 | 40000 | POS with separator |  |  |  |  |
| 17 | --++ | 160000 | 120000 | 40000 | Sex chromosomes |  |  |  |  |

*[18 rows x 21 columns]*
```

### sort genome coordinate

```python
mysumstats.sort_coordinate()
```

**stdout:**
```
2025/12/25 21:54:46 Start to sort the genome coordinates ...(v4.0.0)
2025/12/25 21:54:46 Finished sorting coordinates.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```

### sort column

```python
mysumstats.sort_column()
```

**stdout:**
```
2025/12/25 21:54:46 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:54:46  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,CHISQ,Z,P,MLOG10P,DIRECTION,N,N_CASE,N_CONTROL,NOTE
2025/12/25 21:54:46 Finished reordering the columns.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7fefb5f106b0>
```
