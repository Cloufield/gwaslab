# Clumping by calling PLINK

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```python
2025/05/11 19:15:52 GWASLab v3.6.3 https://cloufield.github.io/gwaslab/
2025/05/11 19:15:52 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/05/11 19:15:52 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Load sample data and perform QC

```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             se="SE",
             p="P",
             nrows=1000000  # Load only the first 1 million lines for demonstration purposes
             )
mysumstats.basic_check(verbose=False)
mysumstats.fix_id(fixsep=True)
```

**stdout:**
```python
2025/05/11 19:13:45 GWASLab v3.6.2 https://cloufield.github.io/gwaslab/
2025/05/11 19:13:45 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/05/11 19:13:45 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/05/11 19:13:45 Start to initialize gl.Sumstats from file :../0_sample_data/t2d_bbj.txt.gz
2025/05/11 19:13:46  -Reading columns          : SE,SNP,P,CHR,POS
2025/05/11 19:13:46  -Renaming columns to      : SE,SNPID,P,CHR,POS
2025/05/11 19:13:46  -Current Dataframe shape : 1000000  x  5
2025/05/11 19:13:46  -Initiating a status column: STATUS ...
2025/05/11 19:13:46  #WARNING! Version of genomic coordinates is unknown...
2025/05/11 19:13:47 Start to reorder the columns...v3.6.2
2025/05/11 19:13:47  -Current Dataframe shape : 1000000 x 6 ; Memory usage: 63.43 MB
2025/05/11 19:13:47  -Reordering columns to    : SNPID,CHR,POS,SE,P,STATUS
2025/05/11 19:13:47 Finished reordering the columns.
2025/05/11 19:13:47  -Trying to convert datatype for CHR: string -> Int64...Int64
2025/05/11 19:13:47  -Column  : SNPID  CHR   POS   SE      P       STATUS  
2025/05/11 19:13:47  -DType   : object Int64 int64 float64 float64 category
2025/05/11 19:13:47  -Verified: T      T     T     T       T       T       
2025/05/11 19:13:47  -Current Dataframe memory usage: 64.38 MB
2025/05/11 19:13:47 Finished loading data successfully!
2025/05/11 19:13:58  #WARNING! Necessary columns for .fix_allele() were not detected:EA,NEA
2025/05/11 19:14:00  #WARNING! Necessary columns for .normalize() were not detected:EA,NEA
2025/05/11 19:14:01 Start to check SNPID/rsID...v3.6.2
2025/05/11 19:14:01  -Current Dataframe shape : 1000000 x 6 ; Memory usage: 65.33 MB
2025/05/11 19:14:01  -Checking SNPID data type...
2025/05/11 19:14:01  -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)
2025/05/11 19:14:04  -Replacing [_-] in SNPID with ":" ...
2025/05/11 19:14:05 Finished checking SNPID/rsID.
```

```python
mysumstats.data
```

```python
| SNPID | CHR | POS | SE | P | STATUS |
| --- | --- | --- | --- | --- | --- |
| 0 | 1:725932:G:A | 1 | 725932 | 0.1394 | 0.59700 |
| 1 | 1:725933:A:G | 1 | 725933 | 0.1394 | 0.59730 |
| 2 | 1:737801:T:C | 1 | 737801 | 0.1231 | 0.69080 |
| 3 | 1:749963:T:TAA | 1 | 749963 | 0.0199 | 0.28460 |
| 4 | 1:751343:T:A | 1 | 751343 | 0.0156 | 0.27050 |
| ... | ... | ... | ... | ... | ... |
| 999995 | 2:6347639:C:A | 2 | 6347639 | 0.0111 | 0.15100 |
| 999996 | 2:6347694:G:C | 2 | 6347694 | 0.0111 | 0.17250 |
| 999997 | 2:6348478:G:A | 2 | 6348478 | 0.0111 | 0.17160 |
| 999998 | 2:6348490:G:C | 2 | 6348490 | 0.2032 | 0.11750 |
| 999999 | 2:6348754:A:C | 2 | 6348754 | 0.0126 | 0.01846 |

*[1000000 rows x 6 columns]*
```

## Run clumping by calling plink2

GWASLab will extract the chromosomes with lead variants and prepare input files for running clumping by calling PLINK2

- VCF: reference LD VCF file (it will be converted to bfile)
- bfile: PLINK bfile prefix
- pfile: PLINK2 pfile prefix

Default parameters:
    
- clump_p1=5e-8
- clump_p2=5e-8
- clump_r2=0.01
- clump_kb=250

Note:
- plink2 need to be avaialble in your environment, or you can specify the path for PLINK2.

```python
mysumstats.clump(   plink2="plink2", # default
                    clump_r2=0.1,
                    clump_p2=1e-6,
                    vcf=gl.get_path("1kg_eas_hg19"))
```

**stdout:**
```python
2025/05/11 19:14:05 Start to perfrom clumping...v3.6.2
2025/05/11 19:14:05  -Current Dataframe shape : 1000000 x 6 ; Memory usage: 65.33 MB
2025/05/11 19:14:05 Start to perform clumping...
2025/05/11 19:14:05  -Clumping parameters for PLINK2:
2025/05/11 19:14:05   -clump_p1 : 5e-08...
2025/05/11 19:14:05   -clump_p2 : 1e-06...
2025/05/11 19:14:05   -clump_kb : 250...
2025/05/11 19:14:05   -clump_r2 : 0.1...
2025/05/11 19:14:05  -Clumping will be performed using P
2025/05/11 19:14:05  -Significant variants on CHR:  [1, 2]
2025/05/11 19:14:05  -Processing VCF : /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz...
2025/05/11 19:14:05  -PLINK version: PLINK v2.00a5.9LM AVX2 AMD (12 Dec 2023)
2025/05/11 19:14:05   -Processing VCF for CHR 1...
2025/05/11 19:14:05   -Plink bfile for CHR 1 exists: /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.1. Skipping...
2025/05/11 19:14:05   -Processing VCF for CHR 2...
2025/05/11 19:14:05   -Plink bfile for CHR 2 exists: /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2. Skipping...
2025/05/11 19:14:05  -Processing sumstats for CHR 1...
2025/05/11 19:14:12  -Variants in reference file: 6471908...
2025/05/11 19:14:12  -Variants in sumstats: 624...
2025/05/11 19:14:12  -Variants available in both reference and sumstats: 624...
2025/05/11 19:14:12  -Processing sumstats for CHR 2...
2025/05/11 19:14:20  -Variants in reference file: 7085910...
2025/05/11 19:14:20  -Variants in sumstats: 233...
2025/05/11 19:14:20  -Variants available in both reference and sumstats: 233...
2025/05/11 19:14:20  -Performing clumping for CHR 1...
2025/05/11 19:14:20  -PLINK version: PLINK v2.00a5.9LM AVX2 AMD (12 Dec 2023)
2025/05/11 19:14:21  -Saved results for CHR 1 to : clumping_plink2.1.clumps
2025/05/11 19:14:21  -Performing clumping for CHR 2...
2025/05/11 19:14:21  -PLINK version: PLINK v2.00a5.9LM AVX2 AMD (12 Dec 2023)
2025/05/11 19:14:23  -Saved results for CHR 2 to : clumping_plink2.2.clumps
2025/05/11 19:14:23 Finished clumping.
2025/05/11 19:14:23 Finished clumping.
```

if using multiple bfile or pfile, @ can be used to indicate each chromosome.

```python
#mysumstats.clump(   plink2="plink2", # default
#                    clump_r2=0.1,
#                    clump_p2=1e-6,
#                    bfile="/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.@")
```

### clump results are stored in mysumstats.clumps

```python
mysumstats.clumps["clumps"]
```

```python
| SNPID | CHR | POS | SE | P | STATUS |
| --- | --- | --- | --- | --- | --- |
| 96739 | 1:22068326:A:G | 1 | 22068326 | 0.0103 | 1.629000e-09 |
| 213860 | 1:51103268:T:C | 1 | 51103268 | 0.0120 | 2.519000e-11 |
| 214500 | 1:51401146:CT:C | 1 | 51401146 | 0.0145 | 6.090000e-10 |
| 534095 | 1:154309595:TA:T | 1 | 154309595 | 0.0166 | 3.289000e-08 |
| 969974 | 2:640986:CACAT:C | 2 | 640986 | 0.0150 | 2.665000e-10 |
```

```python
mysumstats.clumps["clumps_raw"]
```

```python
| CHR | POS | SNPID | P | TOTAL | NONSIG | S0.05 | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 1 | 22068326 | 1:22068326:A:G | 1.629000e-09 | 89 | 0 | 0 |
| 0 | 1 | 51103268 | 1:51103268:T:C | 2.519000e-11 | 347 | 0 | 0 |
| 1 | 1 | 51401146 | 1:51401146:CT:C | 6.090000e-10 | 134 | 0 | 0 |
| 3 | 1 | 154309595 | 1:154309595:TA:T | 3.289000e-08 | 5 | 0 | 0 |
| 4 | 2 | 640986 | 2:640986:CACAT:C | 2.665000e-10 | 232 | 0 | 0 |
|  | S0.01 | S0.001 | S0.0001 | SP2 |  |  |  |
| 2 | 0 | 0 | 89 | 1:21982702:A:G,1:21983254:C:T,1:21983611:C:T,1... |  |  |  |
| 0 | 0 | 0 | 347 | 1:50854197:G:C,1:50854654:A:G,1:50855083:A:G,1... |  |  |  |
| 1 | 0 | 0 | 134 | 1:51353720:C:G,1:51355823:A:G,1:51356091:G:T,1... |  |  |  |
| 3 | 0 | 0 | 5 | 1:154264194:G:A,1:154321382:C:T,1:154345344:G:... |  |  |  |
| 4 | 0 | 0 | 232 | 2:601905:T:G,2:605176:T:A,2:609177:A:G,2:61060... |  |  |  |
```

```python
print(mysumstats.clumps["plink_log"])
```

**stdout:**
```python
PLINK v2.00a5.9LM AVX2 AMD (12 Dec 2023)       www.cog-genomics.org/plink/2.0/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to clumping_plink2.1.log.
Options in effect:
  --bfile /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.1
  --chr 1
  --clump _gwaslab_tmp.1.SNPIDP
  --clump-id-field SNPID
  --clump-kb 250
  --clump-p-field P
  --clump-p1 5e-08
  --clump-p2 1e-06
  --clump-r2 0.1
  --out clumping_plink2.1
  --threads 1

Start time: Sun May 11 19:14:20 2025
31093 MiB RAM detected, ~10647 available; reserving 10583 MiB for main
workspace.
Using 1 compute thread.
504 samples (0 females, 0 males, 504 ambiguous; 504 founders) loaded from
/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.1.fam.
6471908 variants loaded from
/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.1.bim.
Note: No phenotype data present.

--clump: 0/312 index candidates processed.
--clump: 4 clumps formed from 312 index candidates.  
Results written to clumping_plink2.1.clumps .
End time: Sun May 11 19:14:21 2025

PLINK v2.00a5.9LM AVX2 AMD (12 Dec 2023)       www.cog-genomics.org/plink/2.0/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to clumping_plink2.2.log.
Options in effect:
  --bfile /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2
  --chr 2
  --clump _gwaslab_tmp.2.SNPIDP
  --clump-id-field SNPID
  --clump-kb 250
  --clump-p-field P
  --clump-p1 5e-08
  --clump-p2 1e-06
  --clump-r2 0.1
  --out clumping_plink2.2
  --threads 1

Start time: Sun May 11 19:14:22 2025
31093 MiB RAM detected, ~10619 available; reserving 10555 MiB for main
workspace.
Using 1 compute thread.
504 samples (0 females, 0 males, 504 ambiguous; 504 founders) loaded from
/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2.fam.
7085910 variants loaded from
/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2.bim.
Note: No phenotype data present.

--clump: 0/231 index candidates processed.
--clump: 1 clump formed from 231 index candidates.  
Results written to clumping_plink2.2.clumps .
End time: Sun May 11 19:14:23 2025
```
