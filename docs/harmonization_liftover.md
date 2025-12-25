# Liftover

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```python
2025/12/25 22:07:57 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 22:07:57 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 22:07:57 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Load sample data

```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             neaf="Frq",
             beta="BETA",
             se="SE",
             p="P",
             nrows=5000,
             verbose=False)
mysumstats.basic_check(verbose=False)
```

```python
<gwaslab.g_Sumstats.Sumstats at 0x7fb73278db20>
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 725932 | G | A | 9960099 | 0.9960 -0.0737 |
| 1 | 1:725933_A_G | 1 | 725933 | G | A | 9960099 | 0.0040 |
| 2 | 1:737801_T_C | 1 | 737801 | C | T | 9960099 | 0.0051 |
| 3 | 1:749963_T_TAA | 1 | 749963 | TAA | T | 9960399 | 0.8374 |
| 4 | 1:751343_T_A | 1 | 751343 | T | A | 9960099 | 0.8593 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 4995 | 0.0528 | 0.4172 |  |  |  |  |  |
| 4996 | 0.0754 | 0.1015 |  |  |  |  |  |
| 4997 | 0.1271 | 0.2371 |  |  |  |  |  |
| 4998 | 0.0152 | 0.2977 |  |  |  |  |  |
| 4999 | 0.0152 | 0.3028 |  |  |  |  |  |

*[5000 rows x 10 columns]*
```

## Liftover

```python
mysumstats.liftover(from_build="19", to_build="38")
```

**stdout:**
```python
2025/12/25 22:07:59 Start to perform liftover ...(v4.0.0)
2025/12/25 22:07:59  -Using built-in chain file: /home/yunye/anaconda3/envs/py312/lib/python3.12/site-packages/gwaslab/data/chains/hg19ToHg38.over.chain.gz
2025/12/25 22:07:59  -Converting variants with status code xxx0xxx: 5,000
2025/12/25 22:07:59  -Target build: 38
2025/12/25 22:07:59  -Input positions are 1-based
2025/12/25 22:07:59  -Output positions will be 1-based
2025/12/25 22:07:59  -Mapped: 4992 variants
2025/12/25 22:07:59  -Unmapped: 8 variants
2025/12/25 22:07:59  -Examples of unmapped variants:
2025/12/25 22:07:59    SNPID=1:1581586_T_C | CHR=1 | POS=1581586 | STATUS=9960099
2025/12/25 22:07:59    SNPID=1:1584117_G_C | CHR=1 | POS=1584117 | STATUS=9960099
2025/12/25 22:07:59    SNPID=1:1584218_T_C | CHR=1 | POS=1584218 | STATUS=9960099
2025/12/25 22:07:59    SNPID=1:1585257_A_G | CHR=1 | POS=1585257 | STATUS=9960099
2025/12/25 22:07:59    SNPID=1:1585313_CACAA_C | CHR=1 | POS=1585313 | STATUS=9960399
2025/12/25 22:07:59  -Removed 8 unmapped variants
2025/12/25 22:07:59 Start to fix chromosome notation (CHR) ...(v4.0.0)
2025/12/25 22:07:59  -Current Dataframe shape : 4992 x 10 ; Memory usage: 0.36 MB
2025/12/25 22:07:59  -Checking CHR data type...
2025/12/25 22:07:59  -Variants with standardized chromosome notation: 4992
2025/12/25 22:07:59  -All CHR are already fixed...
2025/12/25 22:07:59 Finished fixing chromosome notation (CHR).
2025/12/25 22:07:59 Start to fix basepair positions (POS) ...(v4.0.0)
2025/12/25 22:07:59  -Trying to convert datatype for POS: Int64 -> Int64...
2025/12/25 22:07:59  -Position bound:(0 , 250,000,000)
2025/12/25 22:07:59  -Removed variants outliers: 0
2025/12/25 22:07:59  -Removed variants with bad positions: 0
2025/12/25 22:07:59 Finished fixing basepair positions (POS).
2025/12/25 22:07:59 Finished liftover.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 790552 | G | A | 3860099 | 0.9960 -0.0737 |
| 1 | 1:725933_A_G | 1 | 790553 | G | A | 3860099 | 0.0040 |
| 2 | 1:737801_T_C | 1 | 802421 | C | T | 3860099 | 0.0051 |
| 3 | 1:749963_T_TAA | 1 | 814583 | TAA | T | 3860399 | 0.8374 |
| 4 | 1:751343_T_A | 1 | 815963 | T | A | 3860099 | 0.8593 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 4995 | 0.0528 | 0.4172 |  |  |  |  |  |
| 4996 | 0.0754 | 0.1015 |  |  |  |  |  |
| 4997 | 0.1271 | 0.2371 |  |  |  |  |  |
| 4998 | 0.0152 | 0.2977 |  |  |  |  |  |
| 4999 | 0.0152 | 0.3028 |  |  |  |  |  |

*[4992 rows x 10 columns]*
```

# Liftover using user-provided chain

```python
#https://github.com/marbl/CHM13
```

```python
! wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
```

**stdout:**
```python
--2025-12-25 22:08:00--  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
Resolving s3-us-west-2.amazonaws.com (s3-us-west-2.amazonaws.com)... 52.92.233.24, 52.92.131.48, 52.92.233.136, ...
Connecting to s3-us-west-2.amazonaws.com (s3-us-west-2.amazonaws.com)|52.92.233.24|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 6288201 (6.0M) [binary/octet-stream]
Saving to: ‘grch38-chm13v2.chain.4’

grch38-chm13v2.chai 100%[===================>]   6.00M  3.44MB/s    in 1.7s    

2025-12-25 22:08:02 (3.44 MB/s) - ‘grch38-chm13v2.chain.4’ saved [6288201/6288201]
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 790552 | G | A | 3860099 | 0.9960 -0.0737 |
| 1 | 1:725933_A_G | 1 | 790553 | G | A | 3860099 | 0.0040 |
| 2 | 1:737801_T_C | 1 | 802421 | C | T | 3860099 | 0.0051 |
| 3 | 1:749963_T_TAA | 1 | 814583 | TAA | T | 3860399 | 0.8374 |
| 4 | 1:751343_T_A | 1 | 815963 | T | A | 3860099 | 0.8593 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 4995 | 0.0528 | 0.4172 |  |  |  |  |  |
| 4996 | 0.0754 | 0.1015 |  |  |  |  |  |
| 4997 | 0.1271 | 0.2371 |  |  |  |  |  |
| 4998 | 0.0152 | 0.2977 |  |  |  |  |  |
| 4999 | 0.0152 | 0.3028 |  |  |  |  |  |

*[4992 rows x 10 columns]*
```

```python
mysumstats.liftover(from_build="38", to_build="13",chain_path="./grch38-chm13v2.chain")
```

**stdout:**
```python
2025/12/25 22:08:02 Start to perform liftover ...(v4.0.0)
2025/12/25 22:08:02  -Using provided chain file: ./grch38-chm13v2.chain
2025/12/25 22:08:02  -Converting variants with status code xxx0xxx: 4,992
2025/12/25 22:08:02  -Target build: 13
2025/12/25 22:08:02  -Input positions are 1-based
2025/12/25 22:08:02  -Output positions will be 1-based
2025/12/25 22:08:05  -Mapped: 4947 variants
2025/12/25 22:08:05  -Unmapped: 45 variants
2025/12/25 22:08:05  -Examples of unmapped variants:
2025/12/25 22:08:05    SNPID=1:844476_T_G | CHR=1 | POS=909096 | STATUS=3860099
2025/12/25 22:08:05    SNPID=1:860795_CTGCG_C | CHR=1 | POS=925415 | STATUS=3860399
2025/12/25 22:08:05    SNPID=1:871683_G_A | CHR=1 | POS=936303 | STATUS=3860099
2025/12/25 22:08:05    SNPID=1:928373_A_G | CHR=1 | POS=992993 | STATUS=3860099
2025/12/25 22:08:05    SNPID=1:930892_T_C | CHR=1 | POS=995512 | STATUS=3860099
2025/12/25 22:08:05  -Removed 45 unmapped variants
2025/12/25 22:08:05 Start to fix chromosome notation (CHR) ...(v4.0.0)
2025/12/25 22:08:05  -Current Dataframe shape : 4947 x 10 ; Memory usage: 0.36 MB
2025/12/25 22:08:05  -Checking CHR data type...
2025/12/25 22:08:05  -Variants with standardized chromosome notation: 4947
2025/12/25 22:08:05  -All CHR are already fixed...
2025/12/25 22:08:05 Finished fixing chromosome notation (CHR).
2025/12/25 22:08:05 Start to fix basepair positions (POS) ...(v4.0.0)
2025/12/25 22:08:05  -Trying to convert datatype for POS: Int64 -> Int64...
2025/12/25 22:08:05  -Position bound:(0 , 250,000,000)
2025/12/25 22:08:05  -Removed variants outliers: 0
2025/12/25 22:08:05  -Removed variants with bad positions: 0
2025/12/25 22:08:05 Finished fixing basepair positions (POS).
2025/12/25 22:08:05 Finished liftover.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:725932_G_A | 1 | 219468 | G | A | 1360099 | 0.9960 -0.0737 |
| 1 | 1:725933_A_G | 1 | 219469 | G | A | 1360099 | 0.0040 |
| 2 | 1:737801_T_C | 1 | 231327 | C | T | 1360099 | 0.0051 |
| 3 | 1:749963_T_TAA | 1 | 243588 | TAA | T | 1360399 | 0.8374 |
| 4 | 1:751343_T_A | 1 | 244969 | T | A | 1360099 | 0.8593 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 4995 | 0.0528 | 0.4172 |  |  |  |  |  |
| 4996 | 0.0754 | 0.1015 |  |  |  |  |  |
| 4997 | 0.1271 | 0.2371 |  |  |  |  |  |
| 4998 | 0.0152 | 0.2977 |  |  |  |  |  |
| 4999 | 0.0152 | 0.3028 |  |  |  |  |  |

*[4947 rows x 10 columns]*
```
