# Clumping by calling PLINK

!!! example
    ```python
    import gwaslab as gl
    ```

!!! example
    ```python
    gl.show_version()
    ```

**stdout:**
```
2025/12/26 11:38:17 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/26 11:38:17 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/26 11:38:17 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Load sample data and perform QC

!!! example
    ```python
    mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
                 snpid="SNP",
                 chrom="CHR",
                 pos="POS",
                 ea="ALT",nea="REF",
                 se="SE",
                 p="P",
                 nrows=1000000  # Load only the first 1 million lines for demonstration purposes
                 )
    mysumstats.basic_check(verbose=False)
    mysumstats.fix_id(fixsep=True)
    ```

**stdout:**
```
2025/12/26 11:38:17 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/26 11:38:17 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/26 11:38:17 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
2025/12/26 11:38:17 Start to initialize gl.Sumstats from file :../0_sample_data/t2d_bbj.txt.gz
2025/12/26 11:38:18  -Reading columns          : P,POS,SE,SNP,ALT,REF,CHR
2025/12/26 11:38:18  -Renaming columns to      : P,POS,SE,SNPID,EA,NEA,CHR
2025/12/26 11:38:18  -Current Dataframe shape : 1000000  x  7
2025/12/26 11:38:18  -Initiating a status column: STATUS ...
2025/12/26 11:38:18 #WARNING! Version of genomic coordinates is unknown...
2025/12/26 11:38:18 Start to reorder the columns ...(v4.0.0)
2025/12/26 11:38:18  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,STATUS,SE,P
2025/12/26 11:38:18 Finished reordering the columns.
2025/12/26 11:38:18  -Trying to convert datatype for CHR: string -> Int64...Success
2025/12/26 11:38:18  -Column  : SNPID  CHR   POS   EA       NEA      STATUS SE      P      
2025/12/26 11:38:18  -DType   : object Int64 int64 category category int64  float64 float64
2025/12/26 11:38:18  -Verified: T      T     T     T        T        T      T       T      
2025/12/26 11:38:18  -Current Dataframe memory usage: 50.93 MB
2025/12/26 11:38:18 Finished loading data successfully!
2025/12/26 11:38:19 Start to check SNPID/rsID ...(v4.0.0)
2025/12/26 11:38:19  -Checking SNPID data type...
2025/12/26 11:38:19  -Checking NA strings :na,NA,Na,Nan,NaN,<NA>,null,NULL,#N/A,#VALUE!,N/A,n/a,missing,
2025/12/26 11:38:19  -Checking if SNPID contains NA strings...
2025/12/26 11:38:19  -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)
2025/12/26 11:38:20  -Replacing separators in SNPID with ":" ...
2025/12/26 11:38:20 Finished checking SNPID/rsID.
```

| SNPID | CHR | POS | EA | NEA | STATUS | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932:G:A | 1 | 725932 | G | A | 9960099 | 0.1394 | 0.59700 |
| 1:725933:A:G | 1 | 725933 | G | A | 9960099 | 0.1394 | 0.59730 |
| 1:737801:T:C | 1 | 737801 | C | T | 9960099 | 0.1231 | 0.69080 |
| 1:749963:T:TAA | 1 | 749963 | TAA | T | 9960399 | 0.0199 | 0.28460 |
| 1:751343:T:A | 1 | 751343 | T | A | 9960099 | 0.0156 | 0.27050 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 2:6347639:C:A | 2 | 6347639 | C | A | 9960099 | 0.0111 | 0.15100 |
| 2:6347694:G:C | 2 | 6347694 | G | C | 9960099 | 0.0111 | 0.17250 |
| 2:6348478:G:A | 2 | 6348478 | G | A | 9960099 | 0.0111 | 0.17160 |
| 2:6348490:G:C | 2 | 6348490 | G | C | 9960099 | 0.2032 | 0.11750 |
| 2:6348754:A:C | 2 | 6348754 | C | A | 9960099 | 0.0126 | 0.01846 |

!!! example
    ```python
    mysumstats.data
    ```

| SNPID | CHR | POS | EA | NEA | STATUS | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932:G:A | 1 | 725932 | G | A | 9960099 | 0.1394 | 0.59700 |
| 1:725933:A:G | 1 | 725933 | G | A | 9960099 | 0.1394 | 0.59730 |
| 1:737801:T:C | 1 | 737801 | C | T | 9960099 | 0.1231 | 0.69080 |
| 1:749963:T:TAA | 1 | 749963 | TAA | T | 9960399 | 0.0199 | 0.28460 |
| 1:751343:T:A | 1 | 751343 | T | A | 9960099 | 0.0156 | 0.27050 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 2:6347639:C:A | 2 | 6347639 | C | A | 9960099 | 0.0111 | 0.15100 |
| 2:6347694:G:C | 2 | 6347694 | G | C | 9960099 | 0.0111 | 0.17250 |
| 2:6348478:G:A | 2 | 6348478 | G | A | 9960099 | 0.0111 | 0.17160 |
| 2:6348490:G:C | 2 | 6348490 | G | C | 9960099 | 0.2032 | 0.11750 |
| 2:6348754:A:C | 2 | 6348754 | C | A | 9960099 | 0.0126 | 0.01846 |

*[1000000 rows x 8 columns]*

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

!!! example
    ```python
    mysumstats.clump(   plink2="plink2", # default
                        clump_r2=0.1,
                        clump_p2=1e-6,
                        vcf=gl.get_path("1kg_eas_hg19"))
    ```

**stdout:**
```
2025/12/26 11:38:20  -Number of threads/cores to use: 1
2025/12/26 11:38:20 Start to perfrom clumping ...(v4.0.0)
2025/12/26 11:38:20 Dumpping dataframe to :  ./125110533172752
2025/12/26 11:38:20 Start to perform clumping...
2025/12/26 11:38:20  -Clumping parameters for PLINK2:
2025/12/26 11:38:20   -clump_p1 : 5e-08...
2025/12/26 11:38:20   -clump_p2 : 1e-06...
2025/12/26 11:38:20   -clump_kb : 250...
2025/12/26 11:38:20   -clump_r2 : 0.1...
2025/12/26 11:38:20  -Clumping will be performed using P
2025/12/26 11:38:20  -Significant variants on CHR:  [1, 2]
2025/12/26 11:38:20  -Processing VCF : /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz...
2025/12/26 11:38:20  -PLINK version: PLINK v2.0.0-a.7LM AVX2 Intel (7 Jul 2025)
2025/12/26 11:38:20   -Processing VCF for CHR 1...
2025/12/26 11:38:20   -Plink bfile for CHR 1 exists: /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.1. Skipping...
2025/12/26 11:38:24    -Variants in ref file: 6471908
2025/12/26 11:38:24   -Processing VCF for CHR 2...
2025/12/26 11:38:24   -Plink bfile for CHR 2 exists: /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2. Skipping...
2025/12/26 11:38:28    -Variants in ref file: 7085912
2025/12/26 11:38:30  -Total variants in reference BIM: 13557820...
2025/12/26 11:38:33  -Filtered reference BIM from 13557820 to 859 variants matching sumstats CHR/POS...
2025/12/26 11:38:33  -Matching sumstats with reference BIM using CHR, POS, EA, NEA...
2025/12/26 11:38:33  -Matched 857 variants using CHR, POS, EA, NEA...
2025/12/26 11:38:33  -Created temporary directory: /tmp/gwaslab_clump_7t7vpq9m_125110533172752_
2025/12/26 11:38:33  -Processing sumstats for CHR 1...
2025/12/26 11:38:38  -Variants in reference file: 6471908...
2025/12/26 11:38:38  -Variants in sumstats: 624...
2025/12/26 11:38:38  -Variants available in both reference and sumstats: 624...
2025/12/26 11:38:39  -Processing sumstats for CHR 2...
2025/12/26 11:38:44  -Variants in reference file: 7085912...
2025/12/26 11:38:44  -Variants in sumstats: 233...
2025/12/26 11:38:44  -Variants available in both reference and sumstats: 233...
2025/12/26 11:38:46  -Performing clumping for CHR 1...
2025/12/26 11:38:46  -PLINK version: PLINK v2.0.0-a.7LM AVX2 Intel (7 Jul 2025)
2025/12/26 11:38:47  -Saved results for CHR 1 to : clumping_plink2.1.clumps
2025/12/26 11:38:47  -Performing clumping for CHR 2...
2025/12/26 11:38:47  -PLINK version: PLINK v2.0.0-a.7LM AVX2 Intel (7 Jul 2025)
2025/12/26 11:38:48  -Saved results for CHR 2 to : clumping_plink2.2.clumps
2025/12/26 11:38:48  -Cleaned up temporary directory: /tmp/gwaslab_clump_7t7vpq9m_125110533172752_
2025/12/26 11:38:48  -Mapped BIM SNPIDs back to original SNPIDs in clumping results...
2025/12/26 11:38:48 Loaded dataframe back from :  ./125110533172752
2025/12/26 11:38:48  -Cleaned up 2 additional file(s) after successful reload...
2025/12/26 11:38:49 Finished clumping.
```

if using multiple bfile or pfile, @ can be used to indicate each chromosome.

!!! example
    ```python
    #mysumstats.clump(   plink2="plink2", # default
    #                    clump_r2=0.1,
    #                    clump_p2=1e-6,
    #                    bfile="/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.@")
    ```

### clump results are stored in mysumstats.clumps

!!! example
    ```python
    mysumstats.clumps["clumps"]
    ```

| SNPID | CHR | POS | EA | NEA | STATUS | SE | P | SNPID_bim |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326:A:G | 1 | 22068326 | G | A | 9960099 | 0.0103 | 1.629000e-09 | 1:22068326:A:G |
| 1:51103268:T:C | 1 | 51103268 | C | T | 9960099 | 0.0120 | 2.519000e-11 | 1:51103268:T:C |
| 1:51401146:CT:C | 1 | 51401146 | C | CT | 9960399 | 0.0145 | 6.090000e-10 | 1:51401146:CT:C |
| 1:154309595:TA:T | 1 | 154309595 | TA | T | 9960399 | 0.0166 | 3.289000e-08 | 1:154309595:TA:T |
| 2:640986:CACAT:C | 2 | 640986 | C | CACAT | 9960399 | 0.0150 | 2.665000e-10 | 2:640986:CACAT:C |

!!! example
    ```python
    mysumstats.clumps["clumps_raw"]
    ```

| CHR | POS | SNPID | P | TOTAL | NONSIG | S0.05 | S0.01 | S0.001 | S0.0001 | SP2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 22068326 | 1:22068326:A:G | 1.629000e-09 | 89 | 0 | 0 | 0 | 0 | 89 | 1:21982702:A:G,1:21983254:C:T,1:21983611:C:T,1... |
| 1 | 51103268 | 1:51103268:T:C | 2.519000e-11 | 347 | 0 | 0 | 0 | 0 | 347 | 1:50854197:G:C,1:50854654:A:G,1:50855083:A:G,1... |
| 1 | 51401146 | 1:51401146:CT:C | 6.090000e-10 | 134 | 0 | 0 | 0 | 0 | 134 | 1:51353720:C:G,1:51355823:A:G,1:51356091:G:T,1... |
| 1 | 154309595 | 1:154309595:TA:T | 3.289000e-08 | 5 | 0 | 0 | 0 | 0 | 5 | 1:154264194:G:A,1:154321382:C:T,1:154345344:G:... |
| 2 | 640986 | 2:640986:CACAT:C | 2.665000e-10 | 232 | 0 | 0 | 0 | 0 | 232 | 2:601905:T:G,2:605176:T:A,2:609177:A:G,2:61060... |

!!! example
    ```python
    print(mysumstats.clumps["plink_log"])
    ```

**stdout:**
```
PLINK v2.0.0-a.7LM AVX2 Intel (7 Jul 2025)         cog-genomics.org/plink/2.0/
(C) 2005-2025 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to clumping_plink2.1.log.
Options in effect:
  --bfile /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.1
  --chr 1
  --clump /tmp/gwaslab_clump_7t7vpq9m_125110533172752_/chr1.SNPIDP
  --clump-id-field SNPID
  --clump-kb 250
  --clump-p-field P
  --clump-p1 5e-08
  --clump-p2 1e-06
  --clump-r2 0.1
  --out clumping_plink2.1
  --threads 1

Start time: Fri Dec 26 11:38:46 2025
48175 MiB RAM detected, ~38678 available; reserving 24087 MiB for main
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
End time: Fri Dec 26 11:38:47 2025

PLINK v2.0.0-a.7LM AVX2 Intel (7 Jul 2025)         cog-genomics.org/plink/2.0/
(C) 2005-2025 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to clumping_plink2.2.log.
Options in effect:
  --bfile /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2
  --chr 2
  --clump /tmp/gwaslab_clump_7t7vpq9m_125110533172752_/chr2.SNPIDP
  --clump-id-field SNPID
  --clump-kb 250
  --clump-p-field P
  --clump-p1 5e-08
  --clump-p2 1e-06
  --clump-r2 0.1
  --out clumping_plink2.2
  --threads 1

Start time: Fri Dec 26 11:38:47 2025
48175 MiB RAM detected, ~38701 available; reserving 24087 MiB for main
workspace.
Using 1 compute thread.
504 samples (0 females, 0 males, 504 ambiguous; 504 founders) loaded from
/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2.fam.
7085912 variants loaded from
/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.2.bim.
Note: No phenotype data present.

--clump: 0/231 index candidates processed.
--clump: 1 clump formed from 231 index candidates.  
Results written to clumping_plink2.2.clumps .
End time: Fri Dec 26 11:38:48 2025
```
