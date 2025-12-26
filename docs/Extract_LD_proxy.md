# Extract LD Proxy Variants

This notebook demonstrates how to find LD (Linkage Disequilibrium) proxy variants for SNPs using a VCF reference panel.


## Setup

Import GWASLab and check version.


```python
import gwaslab as gl
gl.show_version()
```

**stdout:**
```
2025/12/26 12:09:10 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/26 12:09:10 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/26 12:09:10 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Load Summary Statistics

Load and prepare the summary statistics data.


```python
mysumstats = gl.Sumstats("../0_sample_data/t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",            
             beta="BETA",
             se="SE",
             p="P", 
             build="19",
             verbose=False, 
             nrows=1000000)
mysumstats.basic_check(verbose=False)
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:725932_G_A | 1 | 725932 | G | A | 1960099 | -0.0737 | 0.1394 | 0.59700 |
| 1:725933_A_G | 1 | 725933 | G | A | 1960099 | 0.0737 | 0.1394 | 0.59730 |
| 1:737801_T_C | 1 | 737801 | C | T | 1960099 | 0.0490 | 0.1231 | 0.69080 |
| 1:749963_T_TAA | 1 | 749963 | TAA | T | 1960399 | 0.0213 | 0.0199 | 0.28460 |
| 1:751343_T_A | 1 | 751343 | T | A | 1960099 | 0.0172 | 0.0156 | 0.27050 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 2:6347639_C_A | 2 | 6347639 | C | A | 1960099 | -0.0159 | 0.0111 | 0.15100 |
| 2:6347694_G_C | 2 | 6347694 | G | C | 1960099 | -0.0151 | 0.0111 | 0.17250 |
| 2:6348478_G_A | 2 | 6348478 | G | A | 1960099 | -0.0152 | 0.0111 | 0.17160 |
| 2:6348490_G_C | 2 | 6348490 | G | C | 1960099 | -0.3180 | 0.2032 | 0.11750 |
| 2:6348754_A_C | 2 | 6348754 | C | A | 1960099 | 0.0297 | 0.0126 | 0.01846 |

## Extract Lead Variants

Identify lead variants from the summary statistics.


```python
mysumstats.get_lead()
```

**stdout:**
```
2025/12/26 12:09:14  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 12:09:14 Start to extract lead variants ...(v4.0.0)
2025/12/26 12:09:14  -Processing 1000000 variants...
2025/12/26 12:09:14  -Significance threshold : 5e-08
2025/12/26 12:09:14  -Sliding window size: 500  kb
2025/12/26 12:09:14  -Using P for extracting lead variants...
2025/12/26 12:09:14  -Found 543 significant variants in total...
2025/12/26 12:09:14  -Identified 4 lead variants!
2025/12/26 12:09:14 Finished extracting lead variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 |
| 1:51103268_T_C | 1 | 51103268 | C | T | 1960099 | -0.0802 | 0.0120 | 2.519000e-11 |
| 1:154309595_TA_T | 1 | 154309595 | TA | T | 1960399 | -0.0915 | 0.0166 | 3.289000e-08 |
| 2:640986_CACAT_C | 2 | 640986 | C | CACAT | 1960399 | -0.0946 | 0.0150 | 2.665000e-10 |

## Find LD Proxies (Basic)

Find LD proxy variants that are present in both the summary statistics and the VCF reference panel.


```python
mysumstats.get_proxy(snplist=["1:22068326_A_G"], 
                        vcf_path="/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz")
```

**stdout:**
```
2025/12/26 12:09:32 Start to find LD proxies for variants ...(v4.0.0)
2025/12/26 12:09:32 Start to load reference genotype...
2025/12/26 12:09:32  -reference vcf path : /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
2025/12/26 12:09:32  -Checking chromosome notations in VCF/BCF files...
2025/12/26 12:09:32  -Checking prefix for chromosomes in VCF/BCF files...
2025/12/26 12:09:32  -No prefix for chromosomes in the VCF/BCF files.
2025/12/26 12:09:32  -1 available variants for LD proxy checking... 
2025/12/26 12:09:32  -Normalized region: (CHR=1, START=21568326, END=22568326)
2025/12/26 12:09:32  -Extract 4950 variants in flanking region of 1:22068326_A_G for checking: 1:21568326-22568326
2025/12/26 12:09:34   -Retrieving index...
2025/12/26 12:09:34   -Ref variants in the region: 31267
2025/12/26 12:09:34   -Matching variants using POS, NEA, EA ...
2025/12/26 12:09:34   -Matched variants in sumstats and vcf:4947 
2025/12/26 12:09:34   -Calculating Rsq...
2025/12/26 12:09:34   -Variants in LD with 1:22068326_A_G (RSQ > 0.8): 21
2025/12/26 12:09:34   -Top Proxy for 1:22068326_A_G is found: 1:22046558_A_C (LD RSQ= 0.8692399263381958)
2025/12/26 12:09:34 Finished loading reference genotype successfully!
2025/12/26 12:09:34 Finished finding LD proxies for variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P | RSQ | LD_REF_VARIANT |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22061876_T_TA | 1 | 22061876 | TA | T | 1960399 | 0.0614 | 0.0103 | 2.522000e-09 | 1.000000 | 1:22068326_A_G |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 | 1.000000 | 1:22068326_A_G |
| 1:22069036_CA_C | 1 | 22069036 | C | CA | 1960399 | 0.0615 | 0.0103 | 2.371000e-09 | 1.000000 | 1:22068326_A_G |
| 1:22060267_T_C | 1 | 22060267 | C | T | 1960099 | 0.0612 | 0.0103 | 2.789000e-09 | 0.994905 | 1:22068326_A_G |
| 1:22046611_G_C | 1 | 22046611 | G | C | 1960099 | -0.0612 | 0.0103 | 2.783000e-09 | 0.989802 | 1:22068326_A_G |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 1:22105568_C_T | 1 | 22105568 | C | T | 1960099 | -0.0585 | 0.0104 | 2.018000e-08 | 0.949362 | 1:22068326_A_G |
| 1:22106875_A_G | 1 | 22106875 | G | A | 1960099 | 0.0592 | 0.0104 | 1.422000e-08 | 0.949362 | 1:22068326_A_G |
| 1:22048891_C_G | 1 | 22048891 | G | C | 1960099 | 0.0603 | 0.0103 | 4.706000e-09 | 0.886967 | 1:22068326_A_G |
| 1:22046558_A_C | 1 | 22046558 | C | A | 1960099 | 0.0607 | 0.0103 | 3.941000e-09 | 0.869240 | 1:22068326_A_G |
| 1:22120540_C_T | 1 | 22120540 | C | T | 1960099 | -0.0531 | 0.0112 | 2.124000e-06 | 0.804480 | 1:22068326_A_G |

## Find LD Proxies (Include All Variants)

Find LD proxy variants including those from the VCF reference panel that are not in the summary statistics. Set `include_all=True` to include all proxy variants from the VCF.

```python
mysumstats.get_proxy(snplist=["1:22068326_A_G"], include_all =True,
                        vcf_path="/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz")
```

**stdout:**
```
2025/12/26 12:09:14 Start to find LD proxies for variants ...(v4.0.0)
2025/12/26 12:09:14 Start to load reference genotype...
2025/12/26 12:09:14  -reference vcf path : /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
2025/12/26 12:09:14  -Checking chromosome notations in VCF/BCF files...
2025/12/26 12:09:14  -Checking prefix for chromosomes in VCF/BCF files...
2025/12/26 12:09:14  -No prefix for chromosomes in the VCF/BCF files.
2025/12/26 12:09:14  -1 available variants for LD proxy checking... 
2025/12/26 12:09:14  -Normalized region: (CHR=1, START=21568326, END=22568326)
2025/12/26 12:09:14  -Extract 4950 variants in flanking region of 1:22068326_A_G for checking: 1:21568326-22568326
2025/12/26 12:09:16   -Retrieving index...
2025/12/26 12:09:16   -Ref variants in the region: 31267
2025/12/26 12:09:16   -Matching variants using POS, NEA, EA ...
2025/12/26 12:09:16   -Matched variants in sumstats and vcf:4947 
2025/12/26 12:09:16   -Calculating Rsq...
2025/12/26 12:09:16   -Variants in LD with 1:22068326_A_G (RSQ > 0.8): 21
2025/12/26 12:09:16   -Top Proxy for 1:22068326_A_G is found: 1:22046558_A_C (LD RSQ= 0.8692399263381958)
2025/12/26 12:09:18   -Retrieving index...
2025/12/26 12:09:18   -Ref variants in the region: 31267
2025/12/26 12:09:18   -Calculating Rsq...
2025/12/26 12:09:20   -Found 20 proxy variants from VCF not in sumstats (RSQ > 0.8)
2025/12/26 12:09:20 Finished loading reference genotype successfully!
2025/12/26 12:09:20 Finished finding LD proxies for variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P | RSQ | LD_REF_VARIANT |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22061876:T:TA | 1 | 22061876 | TA | T | <NA> | NaN | NaN | NaN | 1.000000 | 1:22068326_A_G |
| 1:22061876_T_TA | 1 | 22061876 | TA | T | 1960399 | 0.0614 | 0.0103 | 2.522000e-09 | 1.000000 | 1:22068326_A_G |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 | 1.000000 | 1:22068326_A_G |
| 1:22069036_CA_C | 1 | 22069036 | C | CA | 1960399 | 0.0615 | 0.0103 | 2.371000e-09 | 1.000000 | 1:22068326_A_G |
| 1:22069036:CA:C | 1 | 22069036 | C | CA | <NA> | NaN | NaN | NaN | 1.000000 | 1:22068326_A_G |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 1:22048891_C_G | 1 | 22048891 | G | C | 1960099 | 0.0603 | 0.0103 | 4.706000e-09 | 0.886967 | 1:22068326_A_G |
| 1:22046558:A:C | 1 | 22046558 | C | A | <NA> | NaN | NaN | NaN | 0.869240 | 1:22068326_A_G |
| 1:22046558_A_C | 1 | 22046558 | C | A | 1960099 | 0.0607 | 0.0103 | 3.941000e-09 | 0.869240 | 1:22068326_A_G |
| 1:22120540_C_T | 1 | 22120540 | C | T | 1960099 | -0.0531 | 0.0112 | 2.124000e-06 | 0.804480 | 1:22068326_A_G |
| 1:22120540:C:T | 1 | 22120540 | T | C | <NA> | NaN | NaN | NaN | 0.804480 | 1:22068326_A_G |
