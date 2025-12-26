# Data conversion

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2025/12/26 11:37:56 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/26 11:37:56 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/26 11:37:56 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Loading sample data

```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             neaf="Frq",
             beta="BETA",
             se="SE",nrows=5,verbose=False)
mysumstats.basic_check(verbose=False)
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 |

*[5 rows x 9 columns]*

## BETA -> OR

```python
mysumstats.fill_data(to_fill=["OR"])
```

**stdout:**
```
2025/12/26 11:37:56 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:56  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE     
2025/12/26 11:37:56  -DType   : string Int64 Int64 category category Int64  float32 float64 float64
2025/12/26 11:37:56  -Verified: T      T     T     T        T        T      T       T       T      
2025/12/26 11:37:56  -Target columns: ['OR']
2025/12/26 11:37:56     Filling OR from BETA...
2025/12/26 11:37:56     Filling OR_95L/OR_95U from BETA/SE...
2025/12/26 11:37:56   [Round 1] Filled: ['OR']
2025/12/26 11:37:56   [Round 1] All columns filled!
2025/12/26 11:37:56  -Successfully filled all requested columns.
2025/12/26 11:37:56 Finished filling data using existing columns.
2025/12/26 11:37:56 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:56  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L
2025/12/26 11:37:56 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | OR | OR_95U | OR_95L |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 | 0.928950 | 1.220815 | 0.706863 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 | 1.076484 | 1.414702 | 0.819125 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 | 1.050220 | 1.336790 | 0.825083 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 | 1.021528 | 1.062159 | 0.982452 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 | 1.017349 | 1.048935 | 0.986714 |

*[5 rows x 12 columns]*

## OR -> BETA

```python
mysumstats.data.drop(labels=["BETA","SE"],axis=1,inplace=True)
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | OR | OR_95U | OR_95L |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | 0.928950 | 1.220815 | 0.706863 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 1.076484 | 1.414702 | 0.819125 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 1.050220 | 1.336790 | 0.825083 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 1.021528 | 1.062159 | 0.982452 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 1.017349 | 1.048935 | 0.986714 |

*[5 rows x 10 columns]*

```python
mysumstats.fill_data(to_fill=["BETA","SE"])
```

**stdout:**
```
2025/12/26 11:37:56 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:56  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     OR      OR_95U  OR_95L 
2025/12/26 11:37:56  -DType   : string Int64 Int64 category category Int64  float32 float64 float64 float64
2025/12/26 11:37:56  -Verified: T      T     T     T        T        T      T       T       T       T      
2025/12/26 11:37:56  -Target columns: ['BETA', 'SE']
2025/12/26 11:37:56     Filling SE from OR/OR_95U...
2025/12/26 11:37:56     Filling BETA from OR...
2025/12/26 11:37:56   [Round 1] Filled: ['SE', 'BETA']
2025/12/26 11:37:56   [Round 1] All columns filled!
2025/12/26 11:37:56  -Successfully filled all requested columns.
2025/12/26 11:37:56 Finished filling data using existing columns.
2025/12/26 11:37:56 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:56  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L
2025/12/26 11:37:56 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | OR | OR_95U | OR_95L |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 | 0.928950 | 1.220815 | 0.706863 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 | 1.076484 | 1.414702 | 0.819125 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 | 1.050220 | 1.336790 | 0.825083 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 | 1.021528 | 1.062159 | 0.982452 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 | 1.017349 | 1.048935 | 0.986714 |

*[5 rows x 12 columns]*

## BETA/SE -> Z

```python
mysumstats.fill_data(to_fill=["Z"])
```

**stdout:**
```
2025/12/26 11:37:57 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:57  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      OR      OR_95U  OR_95L 
2025/12/26 11:37:57  -DType   : string Int64 Int64 category category Int64  float32 float64 float64 float64 float64 float64
2025/12/26 11:37:57  -Verified: T      T     T     T        T        T      T       T       T       T       T       T      
2025/12/26 11:37:57  -Target columns: ['Z']
2025/12/26 11:37:57     Filling Z from BETA/SE...
2025/12/26 11:37:57   [Round 1] Filled: ['Z']
2025/12/26 11:37:57   [Round 1] All columns filled!
2025/12/26 11:37:57  -Successfully filled all requested columns.
2025/12/26 11:37:57 Finished filling data using existing columns.
2025/12/26 11:37:57 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:57  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,Z
2025/12/26 11:37:57 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | OR | OR_95U | OR_95L | Z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 | 0.928950 | 1.220815 | 0.706863 | -0.528694 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 | 1.076484 | 1.414702 | 0.819125 | 0.528694 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 | 1.050220 | 1.336790 | 0.825083 | 0.398050 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 | 1.021528 | 1.062159 | 0.982452 | 1.070352 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 | 1.017349 | 1.048935 | 0.986714 | 1.102564 |

*[5 rows x 13 columns]*

## P -> MLOG10P

```python
mysumstats.fill_data(to_fill=["MLOG10P"])
```

**stdout:**
```
2025/12/26 11:37:57 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:57  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      OR      OR_95U  OR_95L  Z      
2025/12/26 11:37:57  -DType   : string Int64 Int64 category category Int64  float32 float64 float64 float64 float64 float64 float64
2025/12/26 11:37:57  -Verified: T      T     T     T        T        T      T       T       T       T       T       T       T      
2025/12/26 11:37:57  -Target columns: ['MLOG10P']
2025/12/26 11:37:57     Filling MLOG10P from Z (extreme)...
2025/12/26 11:37:57   [Round 1] Filled: ['MLOG10P']
2025/12/26 11:37:57   [Round 1] All columns filled!
2025/12/26 11:37:57  -Successfully filled all requested columns.
2025/12/26 11:37:57 Finished filling data using existing columns.
2025/12/26 11:37:57 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:57  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,Z,MLOG10P
2025/12/26 11:37:57 Finished reordering the columns.
```

## MLOG10P -> P 

```python
mysumstats.fill_data(to_fill=["P"])
```

**stdout:**
```
2025/12/26 11:37:57 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:57  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      OR      OR_95U  OR_95L  Z       MLOG10P
2025/12/26 11:37:57  -DType   : string Int64 Int64 category category Int64  float32 float64 float64 float64 float64 float64 float64 float64
2025/12/26 11:37:57  -Verified: T      T     T     T        T        T      T       T       T       T       T       T       T       T      
2025/12/26 11:37:57  -Target columns: ['P']
2025/12/26 11:37:57     Filling P from MLOG10P...
2025/12/26 11:37:57   [Round 1] Filled: ['P']
2025/12/26 11:37:57   [Round 1] All columns filled!
2025/12/26 11:37:57  -Successfully filled all requested columns.
2025/12/26 11:37:57 Finished filling data using existing columns.
2025/12/26 11:37:57 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:57  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,OR,OR_95U,OR_95L,Z,P,MLOG10P
2025/12/26 11:37:57 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | OR | OR_95U | OR_95L | Z | P | MLOG10P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 | 0.928950 | 1.220815 | 0.706863 | -0.528694 | 0.597017 | 0.224013 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 | 1.076484 | 1.414702 | 0.819125 | 0.528694 | 0.597017 | 0.224013 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 | 1.050220 | 1.336790 | 0.825083 | 0.398050 | 0.690593 | 0.160778 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 | 1.021528 | 1.062159 | 0.982452 | 1.070352 | 0.284461 | 0.545977 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 | 1.017349 | 1.048935 | 0.986714 | 1.102564 | 0.270217 | 0.568288 |

*[5 rows x 15 columns]*

## EAF -> MAF

```python
mysumstats.fill_data(to_fill=["MAF"])
```

**stdout:**
```
2025/12/26 11:37:57 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:57  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      OR      OR_95U  OR_95L  Z       P       MLOG10P
2025/12/26 11:37:57  -DType   : string Int64 Int64 category category Int64  float32 float64 float64 float64 float64 float64 float64 float64 float64
2025/12/26 11:37:57  -Verified: T      T     T     T        T        T      T       T       T       T       T       T       T       T       T      
2025/12/26 11:37:57  -Target columns: ['MAF']
2025/12/26 11:37:57     Filling MAF from EAF...
2025/12/26 11:37:57   [Round 1] Filled: ['MAF']
2025/12/26 11:37:57   [Round 1] All columns filled!
2025/12/26 11:37:57  -Successfully filled all requested columns.
2025/12/26 11:37:57 Finished filling data using existing columns.
2025/12/26 11:37:57 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:57  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,MAF,BETA,SE,OR,OR_95U,OR_95L,Z,P,MLOG10P
2025/12/26 11:37:57 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | MAF | BETA | SE | OR | OR_95U | OR_95L | Z | P | MLOG10P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | 0.0040 | -0.0737 | 0.1394 | 0.928950 | 1.220815 | 0.706863 | -0.528694 | 0.597017 | 0.224013 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0040 | 0.0737 | 0.1394 | 1.076484 | 1.414702 | 0.819125 | 0.528694 | 0.597017 | 0.224013 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0051 | 0.0490 | 0.1231 | 1.050220 | 1.336790 | 0.825083 | 0.398050 | 0.690593 | 0.160778 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.1626 | 0.0213 | 0.0199 | 1.021528 | 1.062159 | 0.982452 | 1.070352 | 0.284461 | 0.545977 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.1407 | 0.0172 | 0.0156 | 1.017349 | 1.048935 | 0.986714 | 1.102564 | 0.270217 | 0.568288 |

*[5 rows x 16 columns]*

## Simulation of extreme P values

```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             beta="BETA",
             se="SE",nrows=5, verbose=False)
# simulate some extreme P values by shrinking the SE
mysumstats.data["SE"] = mysumstats.data["SE"]/100
mysumstats.data
```

| SNPID | CHR | POS | STATUS | BETA | SE |
| --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | 9999999 | -0.0737 | 0.001394 |
| 1:725933_A_G | 1 | 725933 | 9999999 | 0.0737 | 0.001394 |
| 1:737801_T_C | 1 | 737801 | 9999999 | 0.0490 | 0.001231 |
| 1:749963_T_TAA | 1 | 749963 | 9999999 | 0.0213 | 0.000199 |
| 1:751343_T_A | 1 | 751343 | 9999999 | 0.0172 | 0.000156 |

*[5 rows x 6 columns]*

## Limited precision of float64

For P < 1e-308, they become 0 due to limnited precision of float64

```python
mysumstats.fill_data(to_fill=["Z","P"])
```

**stdout:**
```
2025/12/26 11:37:57 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:57  -Column  : SNPID  CHR   POS   STATUS BETA    SE     
2025/12/26 11:37:57  -DType   : object Int64 int64 int64  float64 float64
2025/12/26 11:37:57  -Verified: T      T     T     T      T       T      
2025/12/26 11:37:57  -Target columns: ['Z', 'P']
2025/12/26 11:37:57     Filling Z from BETA/SE...
2025/12/26 11:37:57   [Round 1] Filled: ['Z']
2025/12/26 11:37:57   [Round 1] Remaining: ['P']
2025/12/26 11:37:57     Filling P from Z...
2025/12/26 11:37:57   [Round 2] Filled: ['P']
2025/12/26 11:37:57   [Round 2] All columns filled!
2025/12/26 11:37:57  -Successfully filled all requested columns.
2025/12/26 11:37:57 Finished filling data using existing columns.
2025/12/26 11:37:57 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:57  -Reordering columns to    : SNPID,CHR,POS,STATUS,BETA,SE,Z,P
2025/12/26 11:37:57 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | STATUS | BETA | SE | Z | P |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | 9999999 | -0.0737 | 0.001394 | -52.869440 | 0.0 |
| 1:725933_A_G | 1 | 725933 | 9999999 | 0.0737 | 0.001394 | 52.869440 | 0.0 |
| 1:737801_T_C | 1 | 737801 | 9999999 | 0.0490 | 0.001231 | 39.805037 | 0.0 |
| 1:749963_T_TAA | 1 | 749963 | 9999999 | 0.0213 | 0.000199 | 107.035176 | 0.0 |
| 1:751343_T_A | 1 | 751343 | 9999999 | 0.0172 | 0.000156 | 110.256410 | 0.0 |

*[5 rows x 8 columns]*

## Recalculate MLOG10P with extreme P value mode

```python
mysumstats.fill_data(to_fill=["MLOG10P"],extreme=True)
```

**stdout:**
```
2025/12/26 11:37:57 Start filling data using existing columns...v4.0.0
2025/12/26 11:37:57  -Column  : SNPID  CHR   POS   STATUS BETA    SE      Z       P      
2025/12/26 11:37:57  -DType   : object Int64 int64 int64  float64 float64 float64 float64
2025/12/26 11:37:57  -Verified: T      T     T     T      T       T       T       T      
2025/12/26 11:37:57  -Target columns: ['MLOG10P']
2025/12/26 11:37:57     Filling MLOG10P from Z (extreme)...
2025/12/26 11:37:57   [Round 1] Filled: ['MLOG10P']
2025/12/26 11:37:57   [Round 1] All columns filled!
2025/12/26 11:37:57  -Successfully filled all requested columns.
2025/12/26 11:37:57 Finished filling data using existing columns.
2025/12/26 11:37:57 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:37:57  -Reordering columns to    : SNPID,CHR,POS,STATUS,BETA,SE,Z,P,MLOG10P
2025/12/26 11:37:57 Finished reordering the columns.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | STATUS | BETA | SE | Z | P | MLOG10P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | 9999999 | -0.0737 | 0.001394 | -52.869440 | 0.0 | 608.786553 |
| 1:725933_A_G | 1 | 725933 | 9999999 | 0.0737 | 0.001394 | 52.869440 | 0.0 | 608.786553 |
| 1:737801_T_C | 1 | 737801 | 9999999 | 0.0490 | 0.001231 | 39.805037 | 0.0 | 345.755249 |
| 1:749963_T_TAA | 1 | 749963 | 9999999 | 0.0213 | 0.000199 | 107.035176 | 0.0 | 2489.881261 |
| 1:751343_T_A | 1 | 751343 | 9999999 | 0.0172 | 0.000156 | 110.256410 | 0.0 | 2641.885723 |

*[5 rows x 9 columns]*

## Calculate Per-SNP r2

```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             neaf="Frq",
             beta="BETA",n=170000,
             se="SE",nrows=5,verbose=False)
mysumstats.basic_check(verbose=False)
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | N |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 | 170000 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 | 170000 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 | 170000 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 | 170000 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 | 170000 |

*[5 rows x 10 columns]*

```python
mysumstats.get_per_snp_r2()
```

**stdout:**
```
2025/12/26 11:37:57 Start to calculate per-SNP heritibility ...(v4.0.0)
2025/12/26 11:37:57  -Calculating per-SNP rsq by 2 * (BETA**2) * AF * (1-AF) / Var(y)...
2025/12/26 11:37:57  -Var(y) is provided: 1...
2025/12/26 11:37:57  -Calculating F-statistic: F = [(N-k-1)/k] * (r2/1-r2)... where k = 1
2025/12/26 11:37:57  -For r2, SNPR2 is used.
2025/12/26 11:37:57 Finished calculating per-SNP heritability!
2025/12/26 11:37:57  -Current Dataframe shape : 5 x 13 ; Memory usage: 0.00 MB
2025/12/26 11:37:57 Finished calculating per-SNP heritibility.
```

```python
mysumstats.data
```

| SNPID | CHR | POS | EA | NEA | STATUS | EAF | BETA | SE | N | _VAR(BETAX) | SNPR2 | F |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 | -0.0737 | 0.1394 | 170000 | 0.000043 | 0.000043 | 7.357797 |
| 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 | 0.0737 | 0.1394 | 170000 | 0.000043 | 0.000043 | 7.357782 |
| 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 | 0.0490 | 0.1231 | 170000 | 0.000024 | 0.000024 | 4.142153 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 | 0.0213 | 0.0199 | 170000 | 0.000124 | 0.000124 | 21.005844 |
| 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 | 0.0172 | 0.0156 | 170000 | 0.000072 | 0.000072 | 12.161878 |

*[5 rows x 13 columns]*
