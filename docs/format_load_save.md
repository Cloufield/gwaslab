# Input and output sumstats

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2025/12/25 21:52:10 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:52:10 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:52:10 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

# Input

## Loading data

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
             direction="Dir",
             build="19",
             n="N", verbose=False)

# select just 1000 variants for example
mysumstats.random_variants(n=1000, inplace=True, random_state=123,verbose=False)

# basic_check
mysumstats.basic_check(verbose=False)
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7f10939b73b0>
```

```python
mysumstats.data
```

```
| SNPID | CHR | POS | EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:2005486_C_T | 1 | 2005486 | C | T | 1960099 | 0.9863 -0.0969 |
| 1 | 1:2247939_AAGG_A | 1 | 2247939 | AAGG | A | 1960399 | 0.9966 |
| 2 | 1:3741853_G_A | 1 | 3741853 | G | A | 1960099 | 0.8849 -0.0375 |
| 3 | 1:5017526_G_A | 1 | 5017526 | G | A | 1960099 | 0.9822 |
| 4 | 1:5843475_C_T | 1 | 5843475 | C | T | 1960099 | 0.9857 -0.0011 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 995 | 0.0071 | 0.061670 | +--- | 191764 |  |  |  |
| 996 | 0.0363 | 0.295800 | ++++ | 191764 |  |  |  |
| 997 | 0.1671 | 0.691400 | ---+ | 191764 |  |  |  |
| 998 | 0.0078 | 0.635100 | ++-+ | 191764 |  |  |  |
| 999 | 0.0082 | 0.622200 | -++- | 191764 |  |  |  |

*[1000 rows x 12 columns]*
```

## Load and filtering by pattern

### Load and filtering by snp pattern

```python
mysumstats_snp2 = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",fmt="auto",
             snpid_pat="^2:123",
             build="19",
             n="N", verbose=True)
mysumstats_snp2.data
```

**stdout:**
```
2025/12/25 21:52:37 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:52:37 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:52:37 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/12/25 21:52:37 Start to load format from formatbook....
2025/12/25 21:52:37  -auto format meta info:
2025/12/25 21:52:37   - format_name  : auto
2025/12/25 21:52:37   - format_separator  : \t
2025/12/25 21:52:37   - format_na  : #NA
2025/12/25 21:52:37   - format_assumption  : Note: auto-detection assumes A1=EA; Alt=EA and Frq=EAF
2025/12/25 21:52:37   - format_version  :  20230328
2025/12/25 21:52:37   - Auto-detection mode. Note: auto-detection assumes A1=EA; Alt=EA and Frq=EAF...
2025/12/25 21:52:37   - Header conversion source: https://github.com/Cloufield/formatbook/blob/main/formats/auto.json
2025/12/25 21:52:37 Start to initialize gl.Sumstats from file :../0_sample_data/t2d_bbj.txt.gz
2025/12/25 21:52:37  -Columns used to filter variants: SNP
2025/12/25 21:52:37  -Loading only variants with pattern :  ^2:123 ...
2025/12/25 21:52:53  -Loaded 5722 variants with pattern : ^2:123 ...
2025/12/25 21:52:53  -Reading columns          : POS,ALT,N,CHR,BETA,SNP,REF,Frq,P,SE
2025/12/25 21:52:53  -Renaming columns to      : POS,EA,N,CHR,BETA,SNPID,NEA,EAF,P,SE
2025/12/25 21:52:53  -Current Dataframe shape : 5722  x  10
2025/12/25 21:52:53  -Initiating a status column: STATUS ...
2025/12/25 21:52:53  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 21:52:53 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:52:53  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,P,N
2025/12/25 21:52:53 Finished reordering the columns.
2025/12/25 21:52:53  -Trying to convert datatype for CHR: string -> Int64...Success
2025/12/25 21:52:53  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      P       N    
2025/12/25 21:52:53  -DType   : object Int64 int64 category category int64  float64 float64 float64 float64 int64
2025/12/25 21:52:53  -Verified: T      T     T     T        T        T      T       T       T       T       T    
2025/12/25 21:52:53  -Current Dataframe memory usage: 0.46 MB
2025/12/25 21:52:53 Finished loading data successfully!
```

```
| SNPID | CHR | POS | EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 967172 | 2:12320_AAT_A | 2 | 12320 | AAT | A | 1999999 | 0.0572 -0.0158 |
| 967173 | 2:12371_G_C | 2 | 12371 | G | C | 1999999 | 0.0350 -0.0190 |
| 967529 | 2:123233_G_A | 2 | 123233 | G | A | 1999999 | 0.8736 |
| 967530 | 2:123332_C_G | 2 | 123332 | G | C | 1999999 | 0.1264 -0.0191 |
| 967531 | 2:123554_C_A | 2 | 123554 | C | A | 1999999 | 0.4748 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 1490418 | 0.0144 | 0.4834 | 191764 |  |  |  |  |
| 1490419 | 0.0119 | 0.4062 | 191764 |  |  |  |  |
| 1490420 | 0.0222 | 0.7132 | 191764 |  |  |  |  |
| 1490421 | 0.0144 | 0.4760 | 191764 |  |  |  |  |
| 1490422 | 0.1654 | 0.4861 | 191764 |  |  |  |  |

*[5722 rows x 11 columns]*
```

### Load and filtering by chr pattern

```python
mysumstats_chr22 = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",fmt="auto",
             chrom_pat="^22",
             build="19",
             n="N", verbose=True)
mysumstats_chr22.data
```

**stdout:**
```
2025/12/25 21:52:53 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:52:53 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:52:53 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/12/25 21:52:53 Start to load format from formatbook....
2025/12/25 21:52:53  -auto format meta info:
2025/12/25 21:52:53   - format_name  : auto
2025/12/25 21:52:53   - format_separator  : \t
2025/12/25 21:52:53   - format_na  : #NA
2025/12/25 21:52:53   - format_assumption  : Note: auto-detection assumes A1=EA; Alt=EA and Frq=EAF
2025/12/25 21:52:53   - format_version  :  20230328
2025/12/25 21:52:53   - Auto-detection mode. Note: auto-detection assumes A1=EA; Alt=EA and Frq=EAF...
2025/12/25 21:52:53   - Header conversion source: https://github.com/Cloufield/formatbook/blob/main/formats/auto.json
2025/12/25 21:52:53 Start to initialize gl.Sumstats from file :../0_sample_data/t2d_bbj.txt.gz
2025/12/25 21:52:53  -Columns used to filter variants: CHR
2025/12/25 21:52:53  -Loading only variants on chromosome with pattern : ^22 ...
2025/12/25 21:53:08  -Loaded 157050 variants on chromosome with pattern :^22 ...
2025/12/25 21:53:08  -Reading columns          : POS,ALT,N,CHR,BETA,SNP,REF,Frq,P,SE
2025/12/25 21:53:08  -Renaming columns to      : POS,EA,N,CHR,BETA,SNPID,NEA,EAF,P,SE
2025/12/25 21:53:08  -Current Dataframe shape : 157050  x  10
2025/12/25 21:53:08  -Initiating a status column: STATUS ...
2025/12/25 21:53:08  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 21:53:08 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:53:08  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,P,N
2025/12/25 21:53:08 Finished reordering the columns.
2025/12/25 21:53:08  -Trying to convert datatype for CHR: string -> Int64...Success
2025/12/25 21:53:09  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      P       N    
2025/12/25 21:53:09  -DType   : object Int64 int64 category category int64  float64 float64 float64 float64 int64
2025/12/25 21:53:09  -Verified: T      T     T     T        T        T      T       T       T       T       T    
2025/12/25 21:53:09  -Current Dataframe memory usage: 12.93 MB
2025/12/25 21:53:09 Finished loading data successfully!
```

**stderr:**
```
/home/yunye/anaconda3/envs/py312/lib/python3.12/site-packages/gwaslab/io/io_preformat_input.py:860: FutureWarning: The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. To retain the old behavior, exclude the relevant entries before the concat operation.
  sumstats_filtered = pd.concat([chunk[chunk[chunk_chrom].str.match(chrom_pat, case=False,na=False) ] for chunk in sumstats_iter])
/home/yunye/anaconda3/envs/py312/lib/python3.12/site-packages/gwaslab/io/io_preformat_input.py:860: FutureWarning: The behavior of DataFrame concatenation with empty or all-NA entries is deprecated. In a future version, this will no longer exclude empty or all-NA columns when determining the result dtypes. To retain the old behavior, exclude the relevant entries before the concat operation.
  sumstats_filtered = pd.concat([chunk[chunk[chunk_chrom].str.match(chrom_pat, case=False,na=False) ] for chunk in sumstats_iter])
```

```
| SNPID | CHR | POS EA NEA | STATUS | EAF | BETA | \ |
| --- | --- | --- | --- | --- | --- | --- |
| 12071920 | 22:16847963_A_G | 22 | 16847963 | G | A | 1999999 |
| 12071921 | 22:16848015_C_G | 22 | 16848015 | G | C | 1999999 |
| 12071922 | 22:16848470_A_G | 22 | 16848470 | G | A | 1999999 |
| 12071923 | 22:16848520_A_T | 22 | 16848520 | T | A | 1999999 |
| 12071924 | 22:16849105_A_G | 22 | 16849105 | G | A | 1999999 |
| ... | ... | ... | ... | ... | ... | ... |
| 12228965 | 0.2269 | 0.7762 | 166718 |  |  |  |
| 12228966 | 0.0766 | 0.4780 | 166718 |  |  |  |
| 12228967 | 0.0720 | 0.6212 | 191764 |  |  |  |
| 12228968 | 0.0331 | 0.9077 | 166718 |  |  |  |
| 12228969 | 0.2318 | 0.2051 | 166718 |  |  |  |

*[157050 rows x 11 columns]*
```

# Output

## general output

```python
mysumstats.to_format("./mysumstats",fmt="gwaslab",xymt_number=True)
```

**stdout:**
```
2025/12/25 21:53:09 Start to convert the output sumstats in:  gwaslab  format
2025/12/25 21:53:09  -Formatting statistics ...
2025/12/25 21:53:09  -Float statistics formats:
2025/12/25 21:53:09   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:09   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:09  -Start outputting sumstats in gwaslab format...
2025/12/25 21:53:09  -gwaslab format will be loaded...
2025/12/25 21:53:09  -gwaslab format meta info:
2025/12/25 21:53:09   - format_name  : gwaslab
2025/12/25 21:53:09   - format_source  : https://cloufield.github.io/gwaslab/
2025/12/25 21:53:09   - format_version  : 20231220_v4
2025/12/25 21:53:09  -Output path: ./mysumstats.gwaslab.tsv.gz
2025/12/25 21:53:09  -Output columns: SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,P,DIRECTION,N
2025/12/25 21:53:09  -Writing sumstats to: ./mysumstats.gwaslab.tsv.gz...
2025/12/25 21:53:09  -Saving log file to: ./mysumstats.gwaslab.log
2025/12/25 21:53:09 Finished outputting successfully!
```

```python
!zcat mysumstats.gwaslab.tsv.gz | head
```

**stdout:**
```
SNPID	CHR	POS	EA	NEA	STATUS	EAF	BETA	SE	P	DIRECTION	N
1:2005486_C_T	1	2005486	C	T	1960099	0.9863	-0.0969	0.0471	3.9820e-02	+---	191764
1:2247939_AAGG_A	1	2247939	AAGG	A	1960399	0.9966	0.0330	0.1249	7.9190e-01	++--	191764
1:3741853_G_A	1	3741853	G	A	1960099	0.8849	-0.0375	0.0142	8.2820e-03	----	191764
1:5017526_G_A	1	5017526	G	A	1960099	0.9822	0.0126	0.0373	7.3620e-01	+-++	191764
1:5843475_C_T	1	5843475	C	T	1960099	0.9857	-0.0011	0.0433	9.8010e-01	--++	191764
1:9405103_T_C	1	9405103	C	T	1960099	0.0021	-0.0729	0.1516	6.3050e-01	+---	191764
1:9443411_G_A	1	9443411	G	A	1960099	0.9916	0.0362	0.0532	4.9690e-01	+-++	191764
1:12866348_G_C	1	12866348	G	C	1960099	0.9728	-0.0352	0.0431	4.1450e-01	---+	191764
1:14466316_A_G	1	14466316	G	A	1960099	0.6942	-0.0042	0.0096	6.6360e-01	--+-	191764

gzip: stdout: Broken pipe
```

```python
!zcat mysumstats.gwaslab.tsv.gz | tail
```

**stdout:**
```
X:121023171_A_ACTT	23	121023171	ACTT	A	1960399	0.5117	0.0096	0.0071	1.7390e-01	++++	191764
X:134838698_C_CTA	23	134838698	C	CTA	1960399	0.9355	0.0060	0.0145	6.8080e-01	++-+	191764
X:135939006_G_T	23	135939006	G	T	1960099	0.2068	-0.0037	0.0085	6.6300e-01	-+-+	191764
X:136020644_C_T	23	136020644	C	T	1960099	0.8756	0.0103	0.0106	3.3370e-01	+++-	191764
X:138148816_C_T	23	138148816	C	T	1960099	0.7842	0.0089	0.0088	3.1330e-01	++++	191764
X:139318378_A_G	23	139318378	G	A	1960099	0.5252	-0.0132	0.0071	6.1670e-02	+---	191764
X:144540038_G_T	23	144540038	G	T	1960099	0.9866	0.0379	0.0363	2.9580e-01	++++	191764
X:145299627_C_T	23	145299627	C	T	1960099	0.9984	-0.0663	0.1671	6.9140e-01	---+	191764
X:146441317_G_A	23	146441317	G	A	1960099	0.7345	0.0037	0.0078	6.3510e-01	++-+	191764
X:152025052_A_G	23	152025052	G	A	1960099	0.2417	0.0041	0.0082	6.2220e-01	-++-	191764
```

## output each chromosome to a single file

```python
mysumstats.to_format("./mysumstats.@",fmt="gwaslab",xymt_number=True)
```

**stdout:**
```
2025/12/25 21:53:09 Start to convert the output sumstats in:  gwaslab  format
2025/12/25 21:53:09  -Formatting statistics ...
2025/12/25 21:53:09  -Float statistics formats:
2025/12/25 21:53:09   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:09   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:09  -Start outputting sumstats in gwaslab format...
2025/12/25 21:53:09  -gwaslab format will be loaded...
2025/12/25 21:53:09  -gwaslab format meta info:
2025/12/25 21:53:09   - format_name  : gwaslab
2025/12/25 21:53:09   - format_source  : https://cloufield.github.io/gwaslab/
2025/12/25 21:53:09   - format_version  : 20231220_v4
2025/12/25 21:53:09  -Output path: ./mysumstats.@.gwaslab.tsv.gz
2025/12/25 21:53:09  -Output columns: SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,P,DIRECTION,N
2025/12/25 21:53:09  -Writing sumstats to: ./mysumstats.@.gwaslab.tsv.gz...
2025/12/25 21:53:09   -@ detected: writing each chromosome to a single file...
2025/12/25 21:53:09   -Chromosomes:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]...
2025/12/25 21:53:09  -Saving log file to: ./mysumstats.@.gwaslab.log
2025/12/25 21:53:09 Finished outputting successfully!
```

## load sumstats for each chromosome

```python
mysumstats = gl.Sumstats("./mysumstats.@.gwaslab.tsv.gz",
             fmt="gwaslab", verbose=True)
```

**stdout:**
```
2025/12/25 21:53:09 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/25 21:53:09 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/25 21:53:09 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/12/25 21:53:09 Start to load format from formatbook....
2025/12/25 21:53:09  -gwaslab format meta info:
2025/12/25 21:53:09   - format_name  : gwaslab
2025/12/25 21:53:09   - format_source  : https://cloufield.github.io/gwaslab/
2025/12/25 21:53:09   - format_version  : 20231220_v4
2025/12/25 21:53:09  -Detected @ in path: load sumstats by each chromosome...
2025/12/25 21:53:09  -Chromosomes detected: 1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,3,4,5,6,7,8,9
2025/12/25 21:53:09 Start to initialize gl.Sumstats from files with pattern :./mysumstats.@.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.1.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.10.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.11.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.12.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.13.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.14.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.15.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.16.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.17.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.18.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.19.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.2.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.20.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.21.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.22.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.23.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.3.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.4.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.5.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.6.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.7.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.8.gwaslab.tsv.gz
2025/12/25 21:53:09  -Loading:./mysumstats.9.gwaslab.tsv.gz
2025/12/25 21:53:09  -Merging sumstats for chromosomes: 1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,23,3,4,5,6,7,8,9
2025/12/25 21:53:09  -Reading columns          : POS,N,CHR,DIRECTION,BETA,EA,P,EAF,SNPID,STATUS,NEA,SE
2025/12/25 21:53:09  -Renaming columns to      : POS,N,CHR,DIRECTION,BETA,EA,P,EAF,SNPID,STATUS,NEA,SE
2025/12/25 21:53:09  -Current Dataframe shape : 1000  x  12
2025/12/25 21:53:09  -Initiating a status column: STATUS ...
2025/12/25 21:53:09 #WARNING! Version of genomic coordinates is unknown...
2025/12/25 21:53:09 Start to reorder the columns ...(v4.0.0)
2025/12/25 21:53:09  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,P,DIRECTION,N
2025/12/25 21:53:09 Finished reordering the columns.
2025/12/25 21:53:09  -Trying to convert datatype for CHR: string -> Int64...Success
2025/12/25 21:53:09  -Column  : SNPID  CHR   POS   EA       NEA      STATUS EAF     BETA    SE      P       DIRECTION N    
2025/12/25 21:53:09  -DType   : object Int64 int64 category category int64  float64 float64 float64 float64 object    int64
2025/12/25 21:53:09  -Verified: T      T     T     T        T        T      T       T       T       T       T         T    
2025/12/25 21:53:09  -Current Dataframe memory usage: 0.08 MB
2025/12/25 21:53:09 Finished loading data successfully!
```

## Check available formats

List the formats that GWASLab supports

```python
formats = gl.list_formats()
```

**stdout:**
```
2025/12/25 21:53:09 Available formats: auto,auto_0,auto_1,auto_2,auto_neaf,auto_ref,bolt_lmm,ccgwas,cojo,fastgwa,genomicsem,gwascatalog,gwascatalog_hm,gwaslab,ldsc,mesusie,metal,mrmega,mtag,pgscatalog,pgscatalog_hm,pheweb,plink,plink2,plink2_firth,plink2_linear,plink2_logistic,plink_assoc,plink_bim,plink_dosage,plink_fam,plink_fisher,plink_linear,plink_logistic,plink_psam,plink_pvar,popcorn,regenie,regenie_gene,saige,ssf,template,vcf
```

Check the contents of the specified format

```python
ssf_format_dict = gl.check_format("ssf")
```

**stdout:**
```
2025/12/25 21:53:09 Available formats:2025/12/25 21:53:09 meta_data2025/12/25 21:53:09 format_dict2025/12/25 21:53:09 
2025/12/25 21:53:09 {'format_name': 'ssf', 'format_source': 'https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1.full', 'format_cite_name': 'GWAS-SSF v0.1', 'format_separator': '\t', 'format_na': '#NA', 'format_comment': None, 'format_col_order': ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'odds_ratio', 'hazard_ratio', 'standard_error', 'effect_allele_frequency', 'p_value', 'neg_log_10_p_value', 'ci_upper', 'ci_lower', 'rsid', 'variant_id', 'info', 'ref_allele', 'n'], 'format_version': 20230328}2025/12/25 21:53:09 {'variant_id': 'SNPID', 'rsid': 'rsID', 'chromosome': 'CHR', 'base_pair_location': 'POS', 'other_allele': 'NEA', 'effect_allele': 'EA', 'effect_allele_frequency': 'EAF', 'n': 'N', 'beta': 'BETA', 'standard_error': 'SE', 'p_value': 'P', 'neg_log_10_p_value': 'MLOG10P', 'info': 'INFO', 'odds_ratio': 'OR', 'hazard_ratio': 'HR', 'ci_lower': 'OR_95L', 'ci_upper': 'OR_95U'}
```

## Formatting and saving

### get ready for submission to gwas catalog (GWAS-ssf format)

- `fmt`: specify the output format
- `ssfmeta`: if True, output the meta file
- `md5sum`: if True, create a file with the md5sum of the output sumstats

```python
mysumstats.to_format("./mysumstats", fmt="ssf", ssfmeta=True, md5sum=True)
```

**stdout:**
```
2025/12/25 21:53:09 Start to convert the output sumstats in:  ssf  format
2025/12/25 21:53:09  -Formatting statistics ...
2025/12/25 21:53:09  -Float statistics formats:
2025/12/25 21:53:09   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:09   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:09  -Replacing SNPID separator from ":" to "_"...
2025/12/25 21:53:09  -Start outputting sumstats in ssf format...
2025/12/25 21:53:09  -ssf format will be loaded...
2025/12/25 21:53:09  -ssf format meta info:
2025/12/25 21:53:09   - format_name  : ssf
2025/12/25 21:53:09   - format_source  : https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1.full
2025/12/25 21:53:09   - format_cite_name  : GWAS-SSF v0.1
2025/12/25 21:53:09   - format_separator  : \t
2025/12/25 21:53:09   - format_na  : #NA
2025/12/25 21:53:09   - format_col_order  : chromosome,base_pair_location,effect_allele,other_allele,beta,odds_ratio,hazard_ratio,standard_error,effect_allele_frequency,p_value,neg_log_10_p_value,ci_upper,ci_lower,rsid,variant_id,info,ref_allele,n
2025/12/25 21:53:09   - format_version  :  20230328
2025/12/25 21:53:09  -gwaslab to ssf format dictionary:
2025/12/25 21:53:09   - gwaslab keys: SNPID,rsID,CHR,POS,NEA,EA,EAF,N,BETA,SE,P,MLOG10P,INFO,OR,HR,OR_95L,OR_95U
2025/12/25 21:53:09   - ssf values: variant_id,rsid,chromosome,base_pair_location,other_allele,effect_allele,effect_allele_frequency,n,beta,standard_error,p_value,neg_log_10_p_value,info,odds_ratio,hazard_ratio,ci_lower,ci_upper
2025/12/25 21:53:09  -Output path: ./mysumstats.ssf.tsv.gz
2025/12/25 21:53:09  -Output columns: chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,effect_allele_frequency,p_value,variant_id,n
2025/12/25 21:53:09  -Writing sumstats to: ./mysumstats.ssf.tsv.gz...
2025/12/25 21:53:09  -md5sum hashing for the file: ./mysumstats.ssf.tsv.gz
2025/12/25 21:53:09  -md5sum path: ./mysumstats.ssf.tsv.gz.md5sum
2025/12/25 21:53:09  -md5sum: e8e9b4cf01d7b0166d0cd9a208a2808f
2025/12/25 21:53:09  -Exporting SSF-style meta data to ./mysumstats.ssf.tsv.gz.ssf.tsv-meta.yaml
2025/12/25 21:53:09  -Saving log file to: ./mysumstats.ssf.log
2025/12/25 21:53:09 Finished outputting successfully!
2025/12/25 21:53:09 Validating SSF format output...
2025/12/25 21:53:09 Using built-in SSF validator (gwas-ssf CLI not available)...
2025/12/25 21:53:09 Validating file extension...
2025/12/25 21:53:09 ✓ File extension OK
2025/12/25 21:53:09 Validating column order...
2025/12/25 21:53:09 ✓ Column order OK
2025/12/25 21:53:09 Validating chromosomes...
2025/12/25 21:53:09 ✓ Chromosomes OK (All autosomes exist. Optional chromosomes ['24', '25'] do not exist.)
2025/12/25 21:53:09 Validating minimum row count...
2025/12/25 21:53:09 #WARNING! ✗ SSF validation failed: The file has fewer than the minimum rows required: 1000 < 100000
2025/12/25 21:53:09 #WARNING!   - The file has fewer than the minimum rows required: 1000 < 100000
```

```python
!zcat mysumstats.ssf.tsv.gz | head
```

**stdout:**
```
chromosome	base_pair_location	effect_allele	other_allele	beta	standard_error	effect_allele_frequency	p_value	variant_id	n
1	2005486	C	T	-0.0969	0.0471	0.9863	3.9820e-02	1_2005486_C_T	191764
1	2247939	AAGG	A	0.0330	0.1249	0.9966	7.9190e-01	1_2247939_AAGG_A	191764
1	3741853	G	A	-0.0375	0.0142	0.8849	8.2820e-03	1_3741853_G_A	191764
1	5017526	G	A	0.0126	0.0373	0.9822	7.3620e-01	1_5017526_G_A	191764
1	5843475	C	T	-0.0011	0.0433	0.9857	9.8010e-01	1_5843475_C_T	191764
1	9405103	C	T	-0.0729	0.1516	0.0021	6.3050e-01	1_9405103_T_C	191764
1	9443411	G	A	0.0362	0.0532	0.9916	4.9690e-01	1_9443411_G_A	191764
1	12866348	G	C	-0.0352	0.0431	0.9728	4.1450e-01	1_12866348_G_C	191764
1	14466316	G	A	-0.0042	0.0096	0.6942	6.6360e-01	1_14466316_A_G	191764
```

```python
!head mysumstats.ssf.tsv.gz.md5sum
```

**stdout:**
```
e8e9b4cf01d7b0166d0cd9a208a2808f
```

```python
!head ./mysumstats.ssf.tsv.gz.ssf.tsv-meta.yaml
```

**stdout:**
```
coordinate_system: 1-based
data_file_md5sum: e8e9b4cf01d7b0166d0cd9a208a2808f
data_file_name: ./mysumstats.ssf.tsv.gz
date_last_modified: 2025-12-25-21:53:09
file_type: GWAS-SSF v0.1
genome_assembly: Unknown
genotyping_technology: Unknown
gwas_id: Unknown
gwaslab:
  basic_check:
```

### ldsc default format

- `hapmap3`: if True, only output hapmap3 SNPs
- `exclude_hla`: if True, exclude variants in HLA region from output

```python
mysumstats.to_format("./mysumstats",fmt="ldsc",hapmap3=True,exclude_hla=True,build="19")
```

**stdout:**
```
2025/12/25 21:53:10 Start to convert the output sumstats in:  ldsc  format
2025/12/25 21:53:10 Start to exclude variants in HLA regions ...(v4.0.0)
2025/12/25 21:53:10  -Current Dataframe shape : 1000 x 12 ; Memory usage: 0.08 MB
2025/12/25 21:53:10  -Excluded 3 variants in HLA region (chr6: 25000000-34000000 )...
2025/12/25 21:53:10  -Filtered out variants: 3
2025/12/25 21:53:10  -Current Dataframe shape : 997 x 12 ; Memory usage: 0.09 MB
2025/12/25 21:53:10 Finished filtering variants.
2025/12/25 21:53:10  -Genomic coordinates are based on GRCh37/hg19...
2025/12/25 21:53:10 Start to extract HapMap3 SNPs ...(v4.0.0)
2025/12/25 21:53:10  -Current Dataframe shape : 997 x 12 ; Memory usage: 0.09 MB
2025/12/25 21:53:10  -Loading Hapmap3 variants from built-in datasets...
2025/12/25 21:53:11  -Since rsID not in sumstats, CHR:POS( build 19) will be used for matching...
2025/12/25 21:53:12  -Checking if alleles are same...
2025/12/25 21:53:12  -Variants with macthed alleles: 81
2025/12/25 21:53:12  -Raw input contains 81 Hapmap3 variants based on CHR:POS...
2025/12/25 21:53:12  -Current Dataframe shape : 997 x 13 ; Memory usage: 0.10 MB
2025/12/25 21:53:12 Finished extracting HapMap3 SNPs.
2025/12/25 21:53:12  -Extract 81 variants in Hapmap3 datasets for build 19.
2025/12/25 21:53:12  -Formatting statistics ...
2025/12/25 21:53:12  -Float statistics formats:
2025/12/25 21:53:12   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:12   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:12  -Start outputting sumstats in ldsc format...
2025/12/25 21:53:12  -ldsc format will be loaded...
2025/12/25 21:53:12  -ldsc format meta info:
2025/12/25 21:53:12   - format_name  : ldsc
2025/12/25 21:53:12   - format_source  : https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
2025/12/25 21:53:12   - format_source2  : https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
2025/12/25 21:53:12   - format_version  :  20150306
2025/12/25 21:53:12  -gwaslab to ldsc format dictionary:
2025/12/25 21:53:12   - gwaslab keys: rsID,NEA,EA,EAF,N,BETA,P,Z,INFO,OR,CHR,POS
2025/12/25 21:53:12   - ldsc values: SNP,A2,A1,Frq,N,Beta,P,Z,INFO,OR,CHR,POS
2025/12/25 21:53:12  -Output path: ./mysumstats.hapmap3.noMHC.ldsc.tsv.gz
2025/12/25 21:53:12  -Output columns: CHR,POS,A1,A2,Frq,Beta,P,N,SNP
2025/12/25 21:53:12  -Writing sumstats to: ./mysumstats.hapmap3.noMHC.ldsc.tsv.gz...
2025/12/25 21:53:12  -Saving log file to: ./mysumstats.hapmap3.noMHC.ldsc.log
2025/12/25 21:53:12 Finished outputting successfully!
```

```python
!zcat ./mysumstats.hapmap3.noMHC.ldsc.tsv.gz | head
```

**stdout:**
```
CHR	POS	A1	A2	Frq	Beta	P	N	SNP
1	14900419	G	A	0.3952	0.0144	1.3750e-01	191764	rs6703840
1	19593199	C	T	0.1323	-0.0127	3.2570e-01	191764	rs7527253
1	35282297	G	A	0.5434	0.0041	6.4190e-01	191764	rs1407135
1	66001402	C	T	0.2103	-0.0148	1.7720e-01	191764	rs1171261
1	83510491	G	A	0.0025	0.0378	6.9800e-01	191764	rs2022427
1	166110693	C	T	0.8627	0.0286	2.5250e-02	191764	rs4656480
1	175886511	G	A	0.1828	-0.0141	2.2480e-01	191764	rs6656281
1	181612041	C	T	0.9603	0.0135	5.5050e-01	191764	rs199955
1	196329362	C	T	0.0301	0.0300	2.5060e-01	191764	rs11801881
```

### vcf

- `bgzip` : if True, bgzip the output vcf/bed
- `tabix` : if True, index the bgzipped file with tabix

```python
mysumstats.to_format("./mysumstats",fmt="vcf",bgzip=True,tabix=True,build="19")
```

**stdout:**
```
2025/12/25 21:53:12 Start to convert the output sumstats in:  vcf  format
2025/12/25 21:53:12  -Formatting statistics ...
2025/12/25 21:53:12  -Float statistics formats:
2025/12/25 21:53:12   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:12   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:12  -Start outputting sumstats in vcf format...
2025/12/25 21:53:12  -vcf format will be loaded...
2025/12/25 21:53:12  -vcf format meta info:
2025/12/25 21:53:12   - format_name  : vcf
2025/12/25 21:53:12   - format_source  : https://github.com/MRCIEU/gwas-vcf-specification/tree/1.0.0
2025/12/25 21:53:12   - format_version  :  20220923
2025/12/25 21:53:12   - format_citation  : Lyon, M.S., Andrews, S.J., Elsworth, B. et al. The variant call format provides efficient and robust storage of GWAS summary statistics. Genome Biol 22, 32 (2021). https://doi.org/10.1186/s13059-020-02248-0
2025/12/25 21:53:12   - format_fixed  : #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
2025/12/25 21:53:12   - format_format  : ID,SS,ES,SE,LP,SI,EZ
2025/12/25 21:53:12  -gwaslab to vcf format dictionary:
2025/12/25 21:53:12   - gwaslab keys: rsID,CHR,POS,NEA,EA,N,EAF,BETA,SE,MLOG10P,INFO,Z
2025/12/25 21:53:12   - vcf values: ID,#CHROM,POS,REF,ALT,SS,AF,ES,SE,LP,SI,EZ
2025/12/25 21:53:12  -Creating VCF file header...
2025/12/25 21:53:12   -VCF header contig build:19
2025/12/25 21:53:12   -ID:Study_1
2025/12/25 21:53:12   -StudyType:Unknown
2025/12/25 21:53:12   -TotalVariants:1000
2025/12/25 21:53:12   -HarmonisedVariants:0
2025/12/25 21:53:12   -VariantsNotHarmonised:1000
2025/12/25 21:53:12   -SwitchedAlleles:0
2025/12/25 21:53:12  -Writing sumstats to: ./mysumstats.vcf...
2025/12/25 21:53:12  -bgzip compressing : ./mysumstats.vcf.gz...
2025/12/25 21:53:12  -tabix indexing : : ./mysumstats.vcf.gz.tbi...
2025/12/25 21:53:12  -Saving log file to: ./mysumstats.vcf.log
2025/12/25 21:53:12 Finished outputting successfully!
```

### parquet

```python
!mkdir -p mysumstats
```

```python
mysumstats.to_format("./mysumstats",
                     fmt="gwaslab",
                     tab_fmt="parquet",
                     to_tabular_kwargs={"partition_cols":["CHR"]})
```

**stdout:**
```
2025/12/25 21:53:12 Start to convert the output sumstats in:  gwaslab  format
2025/12/25 21:53:12  -Start outputting sumstats in gwaslab format...
2025/12/25 21:53:12  -gwaslab format will be loaded...
2025/12/25 21:53:12  -gwaslab format meta info:
2025/12/25 21:53:12   - format_name  : gwaslab
2025/12/25 21:53:12   - format_source  : https://cloufield.github.io/gwaslab/
2025/12/25 21:53:12   - format_version  : 20231220_v4
2025/12/25 21:53:12  -Output path: ./mysumstats.gwaslab.parquet
2025/12/25 21:53:12  -Output columns: SNPID,CHR,POS,EA,NEA,STATUS,EAF,BETA,SE,P,DIRECTION,N
2025/12/25 21:53:12  -Writing sumstats to: ./mysumstats.gwaslab.parquet...
2025/12/25 21:53:12   -Partitioned parquet will be written to directory: ./mysumstats.gwaslab
2025/12/25 21:53:13  -Saving log file to: ./mysumstats.gwaslab.log
2025/12/25 21:53:13 Finished outputting successfully!
```

## For annotation

### convert to bed format

```python
mysumstats.to_format("./mysumstats",fmt="bed")
```

**stdout:**
```
2025/12/25 21:53:13 Start to convert the output sumstats in:  bed  format
2025/12/25 21:53:13  -Formatting statistics ...
2025/12/25 21:53:13  -Float statistics formats:
2025/12/25 21:53:13   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:13   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:13  -Start outputting sumstats in bed format...
2025/12/25 21:53:13  -Number of SNPs : 920
2025/12/25 21:53:13  -Number of Insertions : 52
2025/12/25 21:53:13  -Number of Deletions : 28
2025/12/25 21:53:13  -formatting to 0-based bed-like file...
2025/12/25 21:53:13  -format description: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
2025/12/25 21:53:13  -Adjusting positions in format-specific manner..
2025/12/25 21:53:13  -Output columns: CHR,START,END,NEA/EA,STRAND,SNPID
2025/12/25 21:53:13  -Writing sumstats to: ./mysumstats.bed...
2025/12/25 21:53:13  -Saving log file to: ./mysumstats.bed.log
2025/12/25 21:53:13 Finished outputting successfully!
```

```python
!cat mysumstats.bed | head
```

**stdout:**
```
1	2005485	2005486	T/C	+	1:2005486_C_T
1	2247939	2247939	-/AGG	+	1:2247939_AAGG_A
1	3741852	3741853	A/G	+	1:3741853_G_A
1	5017525	5017526	A/G	+	1:5017526_G_A
1	5843474	5843475	T/C	+	1:5843475_C_T
1	9405102	9405103	T/C	+	1:9405103_T_C
1	9443410	9443411	A/G	+	1:9443411_G_A
1	12866347	12866348	C/G	+	1:12866348_G_C
1	14466315	14466316	A/G	+	1:14466316_A_G
1	14900418	14900419	A/G	+	1:14900419_A_G
```

### convert to vep default format

```python
mysumstats.to_format("./mysumstats",fmt="vep")
```

**stdout:**
```
2025/12/25 21:53:13 Start to convert the output sumstats in:  vep  format
2025/12/25 21:53:13  -Formatting statistics ...
2025/12/25 21:53:13  -Float statistics formats:
2025/12/25 21:53:13   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:13   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:13  -Start outputting sumstats in vep format...
2025/12/25 21:53:13  -Number of SNPs : 920
2025/12/25 21:53:13  -Number of Insertions : 52
2025/12/25 21:53:13  -Number of Deletions : 28
2025/12/25 21:53:13  -formatting to 1-based bed-like file (for vep)...
2025/12/25 21:53:13  -format description: http://asia.ensembl.org/info/docs/tools/vep/vep_formats.html
2025/12/25 21:53:13  -Adjusting positions in format-specific manner..
2025/12/25 21:53:13  -Output columns: CHR,START,END,NEA/EA,STRAND,SNPID
2025/12/25 21:53:13  -Writing sumstats to: ./mysumstats.vep...
2025/12/25 21:53:13  -Saving log file to: ./mysumstats.vep.log
2025/12/25 21:53:13 Finished outputting successfully!
```

```python
!cat mysumstats.vep | head
```

**stdout:**
```
1	2005486	2005486	T/C	+	1:2005486_C_T
1	2247940	2247939	-/AGG	+	1:2247939_AAGG_A
1	3741853	3741853	A/G	+	1:3741853_G_A
1	5017526	5017526	A/G	+	1:5017526_G_A
1	5843475	5843475	T/C	+	1:5843475_C_T
1	9405103	9405103	T/C	+	1:9405103_T_C
1	9443411	9443411	A/G	+	1:9443411_G_A
1	12866348	12866348	C/G	+	1:12866348_G_C
1	14466316	14466316	A/G	+	1:14466316_A_G
1	14900419	14900419	A/G	+	1:14900419_A_G
```

### convert to annovar default input format

```python
mysumstats.to_format("./mysumstats",fmt="annovar")
```

**stdout:**
```
2025/12/25 21:53:13 Start to convert the output sumstats in:  annovar  format
2025/12/25 21:53:13  -Formatting statistics ...
2025/12/25 21:53:13  -Float statistics formats:
2025/12/25 21:53:13   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:13   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:13  -Start outputting sumstats in annovar format...
2025/12/25 21:53:13  -Number of SNPs : 920
2025/12/25 21:53:13  -Number of Insertions : 52
2025/12/25 21:53:13  -Number of Deletions : 28
2025/12/25 21:53:13  -formatting to 1-based bed-like file...
2025/12/25 21:53:13  -format description: https://annovar.openbioinformatics.org/en/latest/user-guide/input/
2025/12/25 21:53:13  -Adjusting positions in format-specific manner..
2025/12/25 21:53:13  -Output columns: CHR,START,END,NEA_out,EA_out,SNPID
2025/12/25 21:53:13  -Writing sumstats to: ./mysumstats.annovar...
2025/12/25 21:53:13  -Saving log file to: ./mysumstats.annovar.log
2025/12/25 21:53:13 Finished outputting successfully!
```

```python
!cat mysumstats.annovar | head
```

**stdout:**
```
1	2005486	2005486	T	C	1:2005486_C_T
1	2247939	2247939	-	AGG	1:2247939_AAGG_A
1	3741853	3741853	A	G	1:3741853_G_A
1	5017526	5017526	A	G	1:5017526_G_A
1	5843475	5843475	T	C	1:5843475_C_T
1	9405103	9405103	T	C	1:9405103_T_C
1	9443411	9443411	A	G	1:9443411_G_A
1	12866348	12866348	C	G	1:12866348_G_C
1	14466316	14466316	A	G	1:14466316_A_G
1	14900419	14900419	A	G	1:14900419_A_G
```

## Filter and then output

```python
mysumstats.filter_value("EAF >0.05 and EAF < 0.95").to_format("./mysumstats_maf005", fmt="ssf", ssfmeta=True, md5sum=True)
```

**stdout:**
```
2025/12/25 21:53:13 Start to filter variants by condition... ...(v4.0.0)
2025/12/25 21:53:13  -Current Dataframe shape : 1000 x 12 ; Memory usage: 0.08 MB
2025/12/25 21:53:13  -Expression: EAF >0.05 and EAF < 0.95
2025/12/25 21:53:13  -Filtered out variants: 483
2025/12/25 21:53:13  -Current Dataframe shape : 517 x 12 ; Memory usage: 0.05 MB
2025/12/25 21:53:13 Finished filtering variants.
2025/12/25 21:53:13 Start to convert the output sumstats in:  ssf  format
2025/12/25 21:53:13  -Formatting statistics ...
2025/12/25 21:53:13  -Float statistics formats:
2025/12/25 21:53:13   - Columns       : ['EAF', 'BETA', 'SE', 'P']
2025/12/25 21:53:13   - Output formats: ['{:.4g}', '{:.4f}', '{:.4f}', '{:.4e}']
2025/12/25 21:53:13  -Replacing SNPID separator from ":" to "_"...
2025/12/25 21:53:13  -Start outputting sumstats in ssf format...
2025/12/25 21:53:13  -ssf format will be loaded...
2025/12/25 21:53:13  -ssf format meta info:
2025/12/25 21:53:13   - format_name  : ssf
2025/12/25 21:53:13   - format_source  : https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1.full
2025/12/25 21:53:13   - format_cite_name  : GWAS-SSF v0.1
2025/12/25 21:53:13   - format_separator  : \t
2025/12/25 21:53:13   - format_na  : #NA
2025/12/25 21:53:13   - format_col_order  : chromosome,base_pair_location,effect_allele,other_allele,beta,odds_ratio,hazard_ratio,standard_error,effect_allele_frequency,p_value,neg_log_10_p_value,ci_upper,ci_lower,rsid,variant_id,info,ref_allele,n
2025/12/25 21:53:13   - format_version  :  20230328
2025/12/25 21:53:13  -gwaslab to ssf format dictionary:
2025/12/25 21:53:13   - gwaslab keys: SNPID,rsID,CHR,POS,NEA,EA,EAF,N,BETA,SE,P,MLOG10P,INFO,OR,HR,OR_95L,OR_95U
2025/12/25 21:53:13   - ssf values: variant_id,rsid,chromosome,base_pair_location,other_allele,effect_allele,effect_allele_frequency,n,beta,standard_error,p_value,neg_log_10_p_value,info,odds_ratio,hazard_ratio,ci_lower,ci_upper
2025/12/25 21:53:13  -Output path: ./mysumstats_maf005.ssf.tsv.gz
2025/12/25 21:53:13  -Output columns: chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,effect_allele_frequency,p_value,variant_id,n
2025/12/25 21:53:13  -Writing sumstats to: ./mysumstats_maf005.ssf.tsv.gz...
2025/12/25 21:53:13  -md5sum hashing for the file: ./mysumstats_maf005.ssf.tsv.gz
2025/12/25 21:53:13  -md5sum path: ./mysumstats_maf005.ssf.tsv.gz.md5sum
2025/12/25 21:53:13  -md5sum: 0370f9f0dcb5360b8d1f13460fba43c0
2025/12/25 21:53:13  -Exporting SSF-style meta data to ./mysumstats_maf005.ssf.tsv.gz.ssf.tsv-meta.yaml
2025/12/25 21:53:13  -Saving log file to: ./mysumstats_maf005.ssf.log
2025/12/25 21:53:13 Finished outputting successfully!
2025/12/25 21:53:13 Validating SSF format output...
2025/12/25 21:53:13 Using built-in SSF validator (gwas-ssf CLI not available)...
2025/12/25 21:53:13 Validating file extension...
2025/12/25 21:53:13 ✓ File extension OK
2025/12/25 21:53:13 Validating column order...
2025/12/25 21:53:13 ✓ Column order OK
2025/12/25 21:53:13 Validating chromosomes...
2025/12/25 21:53:13 ✓ Chromosomes OK (All autosomes exist. Optional chromosomes ['24', '25'] do not exist.)
2025/12/25 21:53:13 Validating minimum row count...
2025/12/25 21:53:13 #WARNING! ✗ SSF validation failed: The file has fewer than the minimum rows required: 517 < 100000
2025/12/25 21:53:13 #WARNING!   - The file has fewer than the minimum rows required: 517 < 100000
```
