# Data conversion

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```python
2024/12/22 22:33:17 GWASLab v3.5.4 https://cloufield.github.io/gwaslab/
2024/12/22 22:33:17 (C) 2022-2024, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com
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

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | STATUS |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 |
```

## BETA -> OR

```python
mysumstats.fill_data(to_fill=["OR"])
```

**stdout:**
```python
2024/12/22 22:33:32 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:32  -Column  : SNPID  CHR   POS   EA       NEA      EAF     BETA    SE      STATUS  
2024/12/22 22:33:32  -DType   : string Int64 Int64 category category float32 float64 float64 category
2024/12/22 22:33:32  -Verified: T      T     T     T        T        T       T       T       T       
2024/12/22 22:33:32  -Overwrite mode:  False
2024/12/22 22:33:32   -Skipping columns:  []
2024/12/22 22:33:32  -Filling columns:  ['OR']
2024/12/22 22:33:32   - Filling Columns iteratively...
2024/12/22 22:33:32   - Filling OR using BETA column...
2024/12/22 22:33:32   - Filling OR_95L/OR_95U using BETA/SE columns...
2024/12/22 22:33:32 Finished filling data using existing columns.
2024/12/22 22:33:32 Start to reorder the columns...v3.5.4
2024/12/22 22:33:32  -Current Dataframe shape : 5 x 12 ; Memory usage: 21.47 MB
2024/12/22 22:33:32  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,OR,OR_95L,OR_95U,STATUS
2024/12/22 22:33:32 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | OR | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 | 0.928950 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 | 0.1394 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 | 0.1231 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 | 0.0199 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 | 0.0156 |
|  | OR_95L | OR_95U | STATUS |  |  |  |  |  |
| 0 | 0.706863 | 1.220815 | 9960099 |  |  |  |  |  |
| 1 | 0.819125 | 1.414702 | 9960099 |  |  |  |  |  |
| 2 | 0.825083 | 1.336790 | 9960099 |  |  |  |  |  |
| 3 | 0.982452 | 1.062159 | 9960399 |  |  |  |  |  |
| 4 | 0.986714 | 1.048935 | 9960099 |  |  |  |  |  |
```

## OR -> BETA

```python
mysumstats.data.drop(labels=["BETA","SE"],axis=1,inplace=True)
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | OR | OR_95L | OR_95U | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 | 0.928950 | 0.706863 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 1.076484 | 0.819125 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 1.050220 | 0.825083 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 1.021528 | 0.982452 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 1.017349 | 0.986714 |
|  | STATUS |  |  |  |  |  |  |  |
| 0 | 9960099 |  |  |  |  |  |  |  |
| 1 | 9960099 |  |  |  |  |  |  |  |
| 2 | 9960099 |  |  |  |  |  |  |  |
| 3 | 9960399 |  |  |  |  |  |  |  |
| 4 | 9960099 |  |  |  |  |  |  |  |
```

```python
mysumstats.fill_data(to_fill=["BETA","SE"])
```

**stdout:**
```python
2024/12/22 22:33:32 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:32  -Column  : SNPID  CHR   POS   EA       NEA      EAF     OR      OR_95L  OR_95U  STATUS  
2024/12/22 22:33:32  -DType   : string Int64 Int64 category category float32 float64 float64 float64 category
2024/12/22 22:33:32  -Verified: T      T     T     T        T        T       T       T       T       T       
2024/12/22 22:33:32  -Overwrite mode:  False
2024/12/22 22:33:32   -Skipping columns:  []
2024/12/22 22:33:32  -Filling columns:  ['BETA', 'SE']
2024/12/22 22:33:32   - Filling Columns iteratively...
2024/12/22 22:33:32   - Filling BETA value using OR column...
2024/12/22 22:33:32   - Filling SE value using OR/OR_95U column...
2024/12/22 22:33:32 Finished filling data using existing columns.
2024/12/22 22:33:32 Start to reorder the columns...v3.5.4
2024/12/22 22:33:32  -Current Dataframe shape : 5 x 12 ; Memory usage: 21.47 MB
2024/12/22 22:33:32  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,OR,OR_95L,OR_95U,STATUS
2024/12/22 22:33:32 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | OR | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 | 0.928950 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 | 0.1394 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 | 0.1231 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 | 0.0199 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 | 0.0156 |
|  | OR_95L | OR_95U | STATUS |  |  |  |  |  |
| 0 | 0.706863 | 1.220815 | 9960099 |  |  |  |  |  |
| 1 | 0.819125 | 1.414702 | 9960099 |  |  |  |  |  |
| 2 | 0.825083 | 1.336790 | 9960099 |  |  |  |  |  |
| 3 | 0.982452 | 1.062159 | 9960399 |  |  |  |  |  |
| 4 | 0.986714 | 1.048935 | 9960099 |  |  |  |  |  |
```

## BETA/SE -> Z

```python
mysumstats.fill_data(to_fill=["Z"])
```

**stdout:**
```python
2024/12/22 22:33:32 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:32  -Column  : SNPID  CHR   POS   EA       NEA      EAF     BETA    SE      OR      OR_95L  OR_95U  STATUS  
2024/12/22 22:33:32  -DType   : string Int64 Int64 category category float32 float64 float64 float64 float64 float64 category
2024/12/22 22:33:32  -Verified: T      T     T     T        T        T       T       T       T       T       T       T       
2024/12/22 22:33:32  -Overwrite mode:  False
2024/12/22 22:33:32   -Skipping columns:  []
2024/12/22 22:33:32  -Filling columns:  ['Z']
2024/12/22 22:33:32   - Filling Columns iteratively...
2024/12/22 22:33:32   - Filling Z using BETA/SE column...
2024/12/22 22:33:32 Finished filling data using existing columns.
2024/12/22 22:33:32 Start to reorder the columns...v3.5.4
2024/12/22 22:33:32  -Current Dataframe shape : 5 x 13 ; Memory usage: 21.47 MB
2024/12/22 22:33:32  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,Z,OR,OR_95L,OR_95U,STATUS
2024/12/22 22:33:32 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | Z | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 -0.528694 |  |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 | 0.1394 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 | 0.1231 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 | 0.0199 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 | 0.0156 |
|  | OR | OR_95L | OR_95U | STATUS |  |  |  |  |
| 0 | 0.928950 | 0.706863 | 1.220815 | 9960099 |  |  |  |  |
| 1 | 1.076484 | 0.819125 | 1.414702 | 9960099 |  |  |  |  |
| 2 | 1.050220 | 0.825083 | 1.336790 | 9960099 |  |  |  |  |
| 3 | 1.021528 | 0.982452 | 1.062159 | 9960399 |  |  |  |  |
| 4 | 1.017349 | 0.986714 | 1.048935 | 9960099 |  |  |  |  |
```

## P -> MLOG10P

```python
mysumstats.fill_data(to_fill=["MLOG10P"])
```

**stdout:**
```python
2024/12/22 22:33:32 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:32  -Column  : SNPID  CHR   POS   EA       NEA      EAF     BETA    SE      Z       OR      OR_95L  OR_95U  STATUS  
2024/12/22 22:33:32  -DType   : string Int64 Int64 category category float32 float64 float64 float64 float64 float64 float64 category
2024/12/22 22:33:32  -Verified: T      T     T     T        T        T       T       T       T       T       T       T       T       
2024/12/22 22:33:32  -Overwrite mode:  False
2024/12/22 22:33:32   -Skipping columns:  []
2024/12/22 22:33:32  -Filling columns:  ['MLOG10P']
2024/12/22 22:33:32   - Filling Columns iteratively...
2024/12/22 22:33:32   - Filling P value using Z column...
2024/12/22 22:33:32   - Filling MLOG10P using P column...
2024/12/22 22:33:32 Finished filling data using existing columns.
2024/12/22 22:33:32 Start to reorder the columns...v3.5.4
2024/12/22 22:33:32  -Current Dataframe shape : 5 x 15 ; Memory usage: 21.47 MB
2024/12/22 22:33:32  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,Z,P,MLOG10P,OR,OR_95L,OR_95U,STATUS
2024/12/22 22:33:32 Finished reordering the columns.
```

## MLOG10P -> P 

```python
mysumstats.fill_data(to_fill=["P"])
```

**stdout:**
```python
2024/12/22 22:33:33 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:33  -Column  : SNPID  CHR   POS   EA       NEA      EAF     BETA    SE      Z       P       MLOG10P OR      OR_95L  OR_95U  STATUS  
2024/12/22 22:33:33  -DType   : string Int64 Int64 category category float32 float64 float64 float64 float64 float64 float64 float64 float64 category
2024/12/22 22:33:33  -Verified: T      T     T     T        T        T       T       T       T       T       T       T       T       T       T       
2024/12/22 22:33:33  -Overwrite mode:  False
2024/12/22 22:33:33   -Skipping columns:  ['P']
2024/12/22 22:33:33  -No available columns to fill. Skipping.
2024/12/22 22:33:33 Finished filling data using existing columns.
2024/12/22 22:33:33 Start to reorder the columns...v3.5.4
2024/12/22 22:33:33  -Current Dataframe shape : 5 x 15 ; Memory usage: 21.47 MB
2024/12/22 22:33:33  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,Z,P,MLOG10P,OR,OR_95L,OR_95U,STATUS
2024/12/22 22:33:33 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | Z | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 -0.528694 |  |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 | 0.1394 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 | 0.1231 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 | 0.0199 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 | 0.0156 |
|  | P | MLOG10P | OR | OR_95L | OR_95U | STATUS |  |  |
| 0 | 0.597017 | 0.224013 | 0.928950 | 0.706863 | 1.220815 | 9960099 |  |  |
| 1 | 0.597017 | 0.224013 | 1.076484 | 0.819125 | 1.414702 | 9960099 |  |  |
| 2 | 0.690593 | 0.160778 | 1.050220 | 0.825083 | 1.336790 | 9960099 |  |  |
| 3 | 0.284461 | 0.545977 | 1.021528 | 0.982452 | 1.062159 | 9960399 |  |  |
| 4 | 0.270217 | 0.568288 | 1.017349 | 0.986714 | 1.048935 | 9960099 |  |  |
```

## EAF -> MAF

```python
mysumstats.fill_data(to_fill=["MAF"])
```

**stdout:**
```python
2024/12/22 22:33:33 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:33  -Column  : SNPID  CHR   POS   EA       NEA      EAF     BETA    SE      Z       P       MLOG10P OR      OR_95L  OR_95U  STATUS  
2024/12/22 22:33:33  -DType   : string Int64 Int64 category category float32 float64 float64 float64 float64 float64 float64 float64 float64 category
2024/12/22 22:33:33  -Verified: T      T     T     T        T        T       T       T       T       T       T       T       T       T       T       
2024/12/22 22:33:33  -Overwrite mode:  False
2024/12/22 22:33:33   -Skipping columns:  []
2024/12/22 22:33:33  -Filling columns:  ['MAF']
2024/12/22 22:33:33   - Filling Columns iteratively...
2024/12/22 22:33:33   - Filling MAF using EAF column...
2024/12/22 22:33:33 Finished filling data using existing columns.
2024/12/22 22:33:33 Start to reorder the columns...v3.5.4
2024/12/22 22:33:33  -Current Dataframe shape : 5 x 16 ; Memory usage: 21.47 MB
2024/12/22 22:33:33  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,MAF,BETA,SE,Z,P,MLOG10P,OR,OR_95L,OR_95U,STATUS
2024/12/22 22:33:33 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | MAF | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 | 0.0040 -0.0737 | 0.1394 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0040 | 0.0737 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0051 | 0.0490 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.1626 | 0.0213 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.1407 | 0.0172 |
|  | Z | P | MLOG10P | OR | OR_95L | OR_95U | STATUS |  |
| 0 -0.528694 | 0.597017 | 0.224013 | 0.928950 | 0.706863 | 1.220815 | 9960099 |  |  |
| 1 | 0.528694 | 0.597017 | 0.224013 | 1.076484 | 0.819125 | 1.414702 | 9960099 |  |
| 2 | 0.398050 | 0.690593 | 0.160778 | 1.050220 | 0.825083 | 1.336790 | 9960099 |  |
| 3 | 1.070352 | 0.284461 | 0.545977 | 1.021528 | 0.982452 | 1.062159 | 9960399 |  |
| 4 | 1.102564 | 0.270217 | 0.568288 | 1.017349 | 0.986714 | 1.048935 | 9960099 |  |
```

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

```python
| SNPID CHR | POS | BETA | SE | STATUS |
| --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 -0.0737 | 0.001394 |
| 1 | 1:725933_A_G | 1 | 725933 | 0.0737 |
| 2 | 1:737801_T_C | 1 | 737801 | 0.0490 |
| 3 | 1:749963_T_TAA | 1 | 749963 | 0.0213 |
| 4 | 1:751343_T_A | 1 | 751343 | 0.0172 |
```

## Limited precision of float64

For P < 1e-308, they become 0 due to limnited precision of float64

```python
mysumstats.fill_data(to_fill=["Z","P"])
```

**stdout:**
```python
2024/12/22 22:33:34 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:34  -Column  : SNPID  CHR    POS   BETA    SE      STATUS  
2024/12/22 22:33:34  -DType   : object string int64 float64 float64 category
2024/12/22 22:33:34  -Verified: T      F      T     T       T       T       
2024/12/22 22:33:34  #WARNING! Columns with possibly incompatible dtypes: CHR
2024/12/22 22:33:34  -Overwrite mode:  False
2024/12/22 22:33:34   -Skipping columns:  []
2024/12/22 22:33:34  -Filling columns:  ['Z', 'P']
2024/12/22 22:33:34   - Filling Columns iteratively...
2024/12/22 22:33:34   - Filling Z using BETA/SE column...
2024/12/22 22:33:34   - Filling P value using Z column...
2024/12/22 22:33:34 Finished filling data using existing columns.
2024/12/22 22:33:34 Start to reorder the columns...v3.5.4
2024/12/22 22:33:34  -Current Dataframe shape : 5 x 8 ; Memory usage: 21.47 MB
2024/12/22 22:33:34  -Reordering columns to    : SNPID,CHR,POS,BETA,SE,Z,P,STATUS
2024/12/22 22:33:34 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID CHR | POS | BETA | SE | Z | P | STATUS |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 -0.0737 | 0.001394 | -52.869440 | 0.0 |
| 1 | 1:725933_A_G | 1 | 725933 | 0.0737 | 0.001394 | 52.869440 |
| 2 | 1:737801_T_C | 1 | 737801 | 0.0490 | 0.001231 | 39.805037 |
| 3 | 1:749963_T_TAA | 1 | 749963 | 0.0213 | 0.000199 | 107.035176 |
| 4 | 1:751343_T_A | 1 | 751343 | 0.0172 | 0.000156 | 110.256410 |
```

## Recalculate MLOG10P with extreme P value mode

```python
mysumstats.fill_data(to_fill=["MLOG10P"],extreme=True)
```

**stdout:**
```python
2024/12/22 22:33:34 Start filling data using existing columns...v3.5.4
2024/12/22 22:33:34  -Column  : SNPID  CHR    POS   BETA    SE      Z       P       STATUS  
2024/12/22 22:33:34  -DType   : object string int64 float64 float64 float64 float64 category
2024/12/22 22:33:34  -Verified: T      F      T     T       T       T       T       T       
2024/12/22 22:33:34  #WARNING! Columns with possibly incompatible dtypes: CHR
2024/12/22 22:33:34  -Overwrite mode:  False
2024/12/22 22:33:34   -Skipping columns:  []
2024/12/22 22:33:34  -Filling columns:  ['MLOG10P']
2024/12/22 22:33:34   - Filling Columns iteratively...
2024/12/22 22:33:34   - Filling MLOG10P using Z column...
2024/12/22 22:33:34 Finished filling data using existing columns.
2024/12/22 22:33:34 Start to reorder the columns...v3.5.4
2024/12/22 22:33:34  -Current Dataframe shape : 5 x 11 ; Memory usage: 21.47 MB
2024/12/22 22:33:34  -Reordering columns to    : SNPID,CHR,POS,BETA,SE,Z,P,MLOG10P,STATUS,P_MANTISSA,P_EXPONENT
2024/12/22 22:33:34 Finished reordering the columns.
```

```python
mysumstats.data
```

```python
| SNPID CHR | POS | BETA | SE | Z | P | MLOG10P | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 -0.0737 | 0.001394 | -52.869440 | 0.0 | 608.786553 |
| 1 | 1:725933_A_G | 1 | 725933 | 0.0737 | 0.001394 | 52.869440 | 0.0 |
| 2 | 1:737801_T_C | 1 | 737801 | 0.0490 | 0.001231 | 39.805037 | 0.0 |
| 3 | 1:749963_T_TAA | 1 | 749963 | 0.0213 | 0.000199 | 107.035176 | 0.0 |
| 4 | 1:751343_T_A | 1 | 751343 | 0.0172 | 0.000156 | 110.256410 | 0.0 |
|  | STATUS | P_MANTISSA | P_EXPONENT |  |  |  |  |
| 0 | 9999999 | 1.634734 | -609.0 |  |  |  |  |
| 1 | 9999999 | 1.634734 | -609.0 |  |  |  |  |
| 2 | 9999999 | 1.756915 | -346.0 |  |  |  |  |
| 3 | 9999999 | 1.314436 | -2490.0 |  |  |  |  |
| 4 | 9999999 | 1.300999 | -2642.0 |  |  |  |  |
```

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

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | N | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 | 170000 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 | 0.1394 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 | 0.1231 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 | 0.0199 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 | 0.0156 |
|  | STATUS |  |  |  |  |  |  |  |
| 0 | 9960099 |  |  |  |  |  |  |  |
| 1 | 9960099 |  |  |  |  |  |  |  |
| 2 | 9960099 |  |  |  |  |  |  |  |
| 3 | 9960399 |  |  |  |  |  |  |  |
| 4 | 9960099 |  |  |  |  |  |  |  |
```

```python
mysumstats.get_per_snp_r2()
```

**stdout:**
```python
2024/12/22 22:33:49 Start to calculate per-SNP heritibility...
2024/12/22 22:33:49  -Calculating per-SNP rsq by 2 * (BETA**2) * AF * (1-AF) / Var(y)...
2024/12/22 22:33:49  -Var(y) is provided: 1...
2024/12/22 22:33:49  -Calculating F-statistic: F = [(N-k-1)/k] * (r2/1-r2)... where k = 1
2024/12/22 22:33:49  -For r2, SNPR2 is used.
2024/12/22 22:33:49 Finished calculating per-SNP heritability!
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | EAF | BETA | SE | N | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 0.9960 -0.0737 | 0.1394 | 170000 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 0.0040 | 0.0737 | 0.1394 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 0.0051 | 0.0490 | 0.1231 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 0.8374 | 0.0213 | 0.0199 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 0.8593 | 0.0172 | 0.0156 |
|  | STATUS | _VAR(BETAX) | SNPR2 | F |  |  |  |  |
| 0 | 9960099 | 0.000043 | 0.000043 | 7.357797 |  |  |  |  |
| 1 | 9960099 | 0.000043 | 0.000043 | 7.357782 |  |  |  |  |
| 2 | 9960099 | 0.000024 | 0.000024 | 4.142153 |  |  |  |  |
| 3 | 9960399 | 0.000124 | 0.000124 | 21.005844 |  |  |  |  |
| 4 | 9960099 | 0.000072 | 0.000072 | 12.161878 |  |  |  |  |
```
