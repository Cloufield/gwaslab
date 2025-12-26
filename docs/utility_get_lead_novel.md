# Lead and novel variants

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2025/12/26 11:36:20 GWASLab v4.0.0 https://cloufield.github.io/gwaslab/
2025/12/26 11:36:20 (C) 2022-2025, Yunye He, Kamatani Lab, GPL-3.0 license, gwaslab@gmail.com
2025/12/26 11:36:20 Python version: 3.12.0 | packaged by conda-forge | (main, Oct  3 2023, 08:43:22) [GCC 12.3.0]
```

## Load sample data

Use only first 1000000 variants as example

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

## Get lead variants

GWASLab will use MLOG10P first by default. If MLOG10P is not avaiable, it will look for P column.

```python
mysumstats.get_lead()
```

**stdout:**
```
2025/12/26 11:36:24  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:24 Start to extract lead variants ...(v4.0.0)
2025/12/26 11:36:24  -Processing 1000000 variants...
2025/12/26 11:36:24  -Significance threshold : 5e-08
2025/12/26 11:36:24  -Sliding window size: 500  kb
2025/12/26 11:36:24  -Using P for extracting lead variants...
2025/12/26 11:36:24  -Found 543 significant variants in total...
2025/12/26 11:36:24  -Identified 4 lead variants!
2025/12/26 11:36:24 Finished extracting lead variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 |
| 1:51103268_T_C | 1 | 51103268 | C | T | 1960099 | -0.0802 | 0.0120 | 2.519000e-11 |
| 1:154309595_TA_T | 1 | 154309595 | TA | T | 1960399 | -0.0915 | 0.0166 | 3.289000e-08 |
| 2:640986_CACAT_C | 2 | 640986 | C | CACAT | 1960399 | -0.0946 | 0.0150 | 2.665000e-10 |

## Get lead variants with gene name annotation

```python
mysumstats.get_lead(anno=True)
```

**stdout:**
```
2025/12/26 11:36:24  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:24 Start to extract lead variants ...(v4.0.0)
2025/12/26 11:36:24  -Processing 1000000 variants...
2025/12/26 11:36:24  -Significance threshold : 5e-08
2025/12/26 11:36:24  -Sliding window size: 500  kb
2025/12/26 11:36:24  -Using P for extracting lead variants...
2025/12/26 11:36:24  -Found 543 significant variants in total...
2025/12/26 11:36:24  -Identified 4 lead variants!
2025/12/26 11:36:24  -Annotating variants using references:ensembl
2025/12/26 11:36:24  -Annotating variants using references based on genome build:19
2025/12/26 11:36:24  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:24 Start to annotate variants with nearest gene name(s) ...(v4.0.0)
2025/12/26 11:36:24  -Current Dataframe shape : 4 x 10 ; Memory usage: 0.64 MB
2025/12/26 11:36:24  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:24  -Assigning Gene name using ensembl_hg19_gtf for protein coding genes
2025/12/26 11:36:24  -Loading and indexing GTF file...
2025/12/26 11:36:26  -Processing variants by chromosome...
2025/12/26 11:36:26    -Processed 3 variants on chromosome 1
2025/12/26 11:36:26    -Processed 1 variants on chromosome 2
2025/12/26 11:36:26 Finished annotating variants with nearest gene name(s) successfully!.
2025/12/26 11:36:26 Finished extracting lead variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P | LOCATION | GENE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 | 0 | USP48 |
| 1:51103268_T_C | 1 | 51103268 | C | T | 1960099 | -0.0802 | 0.0120 | 2.519000e-11 | 0 | FAF1 |
| 1:154309595_TA_T | 1 | 154309595 | TA | T | 1960399 | -0.0915 | 0.0166 | 3.289000e-08 | 0 | ATP8B2 |
| 2:640986_CACAT_C | 2 | 640986 | C | CACAT | 1960399 | -0.0946 | 0.0150 | 2.665000e-10 | -26349 | TMEM18 |

## Different window sizes

```python
mysumstats.get_lead(windowsizekb=1000)
```

**stdout:**
```
2025/12/26 11:36:28  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:28 Start to extract lead variants ...(v4.0.0)
2025/12/26 11:36:28  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 60.72 MB
2025/12/26 11:36:28  -Processing 1000000 variants...
2025/12/26 11:36:28  -Significance threshold : 5e-08
2025/12/26 11:36:28  -Sliding window size: 1000  kb
2025/12/26 11:36:29  -Using P for extracting lead variants...
2025/12/26 11:36:29  -Found 543 significant variants in total...
2025/12/26 11:36:29  -Identified 4 lead variants!
2025/12/26 11:36:29 Finished extracting lead variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 |
| 1:51103268_T_C | 1 | 51103268 | C | T | 1960099 | -0.0802 | 0.0120 | 2.519000e-11 |
| 1:154309595_TA_T | 1 | 154309595 | TA | T | 1960399 | -0.0915 | 0.0166 | 3.289000e-08 |
| 2:640986_CACAT_C | 2 | 640986 | C | CACAT | 1960399 | -0.0946 | 0.0150 | 2.665000e-10 |

## Different thresholds

```python
mysumstats.get_lead(sig_level=1e-10)
```

**stdout:**
```
2025/12/26 11:36:30  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:30 Start to extract lead variants ...(v4.0.0)
2025/12/26 11:36:30  -Processing 1000000 variants...
2025/12/26 11:36:30  -Significance threshold : 1e-10
2025/12/26 11:36:30  -Sliding window size: 500  kb
2025/12/26 11:36:30  -Using P for extracting lead variants...
2025/12/26 11:36:30  -Found 1 significant variants in total...
2025/12/26 11:36:30  -Identified 1 lead variants!
2025/12/26 11:36:30 Finished extracting lead variants.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:51103268_T_C | 1 | 51103268 | C | T | 1960099 | -0.0802 | 0.012 | 2.519000e-11 |

## Check if novel against a tabular file

```python
novel = mysumstats.get_novel(known="../0_sample_data//toy_data/known_loci.txt")
```

**stdout:**
```
2025/12/26 11:36:36  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:36 Start to check if lead variants are known ...(v4.0.0)
2025/12/26 11:36:36  -Genomic coordinates are based on GRCh37/hg19...
2025/12/26 11:36:36 Start to extract lead variants ...(v4.0.0)
2025/12/26 11:36:36  -Processing 1000000 variants...
2025/12/26 11:36:36  -Significance threshold : 5e-08
2025/12/26 11:36:36  -Sliding window size: 500  kb
2025/12/26 11:36:36  -Using P for extracting lead variants...
2025/12/26 11:36:36  -Found 543 significant variants in total...
2025/12/26 11:36:36  -Identified 4 lead variants!
2025/12/26 11:36:36 Finished extracting lead variants.
2025/12/26 11:36:36  -Lead variants in known loci: 2
2025/12/26 11:36:36  -Checking the minimum distance between identified lead variants and provided known variants...
2025/12/26 11:36:36  -Identified  2  known vairants in current sumstats...
2025/12/26 11:36:36  -Identified  2  novel vairants in current sumstats...
2025/12/26 11:36:36 Finished checking if lead variants are known.
```

```python
novel
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P | DISTANCE_TO_KNOWN | KNOWN_ID | NOVEL | LOCATION_OF_KNOWN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326_A_G | 1 | 22068326 | G | A | 1960099 | 0.0621 | 0.0103 | 1.629000e-09 | 29034942.0 | 1:51103268 | True | Upstream |
| 1:51103268_T_C | 1 | 51103268 | C | T | 1960099 | -0.0802 | 0.0120 | 2.519000e-11 | 0.0 | 1:51103268 | False | Same |
| 1:154309595_TA_T | 1 | 154309595 | TA | T | 1960399 | -0.0915 | 0.0166 | 3.289000e-08 | 0.0 | 1:154309595 | False | Same |
| 2:640986_CACAT_C | 2 | 640986 | C | CACAT | 1960399 | -0.0946 | 0.0150 | 2.665000e-10 | NaN | <NA> | True | NoneOnThisChr |

## Ccheck against GWAS Catalog using EFO ID

```python
mysumstats.liftover(to_build="38")
```

**stdout:**
```
2025/12/26 11:37:18 Start to perform liftover ...(v4.0.0)
2025/12/26 11:37:18  -Using built-in chain file: /home/yunye/anaconda3/envs/py312/lib/python3.12/site-packages/gwaslab/data/chains/hg19ToHg38.over.chain.gz
2025/12/26 11:37:18  -Converting variants with status code xxx0xxx: 1,000,000
2025/12/26 11:37:18  -Target build: 38
2025/12/26 11:37:18  -Input positions are 1-based
2025/12/26 11:37:18  -Output positions will be 1-based
2025/12/26 11:37:19  -Chromosome mismatches detected: 52 variants (treated as unmapped)
2025/12/26 11:37:19  -Examples of chromosome mismatches:
2025/12/26 11:37:19    SNPID=1:13127868_C_A | CHR=1 | POS=13127868 | CHR_LIFT=1_KI270766v1_alt | POS_LIFT=5939 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:146705000_G_A | CHR=1 | POS=146705000 | CHR_LIFT=15 | POS_LIFT=39894340 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:147758522_T_C | CHR=1 | POS=147758522 | CHR_LIFT=8 | POS_LIFT=38661844 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:223725692_T_C | CHR=1 | POS=223725692 | CHR_LIFT=9 | POS_LIFT=136033409 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:223725695_T_G | CHR=1 | POS=223725695 | CHR_LIFT=9 | POS_LIFT=136033412 | STATUS=1960099
2025/12/26 11:37:19  -Mapped: 999381 variants
2025/12/26 11:37:19  -Unmapped: 619 variants
2025/12/26 11:37:19  -Examples of unmapped variants:
2025/12/26 11:37:19    SNPID=1:1581586_T_C | CHR=1 | POS=1581586 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:1584117_G_C | CHR=1 | POS=1584117 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:1584218_T_C | CHR=1 | POS=1584218 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:1585257_A_G | CHR=1 | POS=1585257 | STATUS=1960099
2025/12/26 11:37:19    SNPID=1:1585313_CACAA_C | CHR=1 | POS=1585313 | STATUS=1960399
2025/12/26 11:37:20  -Removed 619 unmapped variants
2025/12/26 11:37:20 Start to fix chromosome notation (CHR) ...(v4.0.0)
2025/12/26 11:37:20  -Current Dataframe shape : 999381 x 9 ; Memory usage: 67.36 MB
2025/12/26 11:37:20  -Checking CHR data type...
2025/12/26 11:37:20  -Variants with standardized chromosome notation: 999381
2025/12/26 11:37:20  -All CHR are already fixed...
2025/12/26 11:37:20 Finished fixing chromosome notation (CHR).
2025/12/26 11:37:20 Start to fix basepair positions (POS) ...(v4.0.0)
2025/12/26 11:37:20  -Trying to convert datatype for POS: Int64 -> Int64...
2025/12/26 11:37:20  -Position bound:(0 , 250,000,000)
2025/12/26 11:37:20  -Removed variants outliers: 0
2025/12/26 11:37:20  -Removed variants with bad positions: 0
2025/12/26 11:37:20 Finished fixing basepair positions (POS).
2025/12/26 11:37:20 Finished liftover.
```

```python
# EFO ID can be found on gwas catalog
mysumstats.get_novel(efo="MONDO_0005148")
```

**stdout:**
```
2025/12/26 11:37:25  -Genomic coordinates are based on GRCh38/hg38...
2025/12/26 11:37:25 Start to check if lead variants are known ...(v4.0.0)
2025/12/26 11:37:25  -Genomic coordinates are based on GRCh38/hg38...
2025/12/26 11:37:25 Start to extract lead variants ...(v4.0.0)
2025/12/26 11:37:25  -Processing 999381 variants...
2025/12/26 11:37:25  -Significance threshold : 5e-08
2025/12/26 11:37:25  -Sliding window size: 500  kb
2025/12/26 11:37:25  -Using P for extracting lead variants...
2025/12/26 11:37:25  -Found 543 significant variants in total...
2025/12/26 11:37:25  -Identified 4 lead variants!
2025/12/26 11:37:25 Finished extracting lead variants.
2025/12/26 11:37:25  -Genomic coordinates are based on GRCh38/hg38...
2025/12/26 11:37:25  -Genomic coordinates are based on GRCh38/hg38...
2025/12/26 11:37:25  -Sumstats build matches target build
2025/12/26 11:37:25 Start to retrieve data using EFO: MONDO_0005148...
2025/12/26 11:37:25  -Querying GWAS Catalog API v2 for trait: MONDO_0005148...
2025/12/26 11:37:25  -MONDO ID detected, looking up EFO equivalent...
2025/12/26 11:37:25 Requesting: https://www.ebi.ac.uk/gwas/rest/api/v2/efoTraits/MONDO_0005148
2025/12/26 11:37:27 #WARNING! Error: Received status code 404
2025/12/26 11:37:27 #WARNING! URL: https://www.ebi.ac.uk/gwas/rest/api/v2/efoTraits/MONDO_0005148
2025/12/26 11:37:27 #WARNING! Response: {"timestamp":"2025-12-26T02:37:28.946+00:00","status":404,"error":"Not Found","path":"/gwas/rest/api/v2/efoTraits/MONDO_0005148"}
2025/12/26 11:37:27  -Warning: Could not retrieve trait info for MONDO ID. Trying directly...
2025/12/26 11:37:27 Requesting: https://www.ebi.ac.uk/gwas/rest/api/v2/associations
2025/12/26 11:37:27 Parameters: {'efo_trait': 'MONDO_0005148', 'sort': 'p_value', 'direction': 'asc', 'size': 200, 'page': 0}
2025/12/26 11:37:27 Status code: 200
2025/12/26 11:37:27 #WARNING! No more data or an error occurred. Stopping.
2025/12/26 11:37:27 Retrieved 0 total items
2025/12/26 11:37:27  -No results from API v2, falling back to old API v1 for MONDO ID...
2025/12/26 11:37:27 Start to retrieve data from GWASCatalog ...(v4.0.0)
2025/12/26 11:37:27 searching cache in : ./
2025/12/26 11:37:27  -Cache not found for False... Downloading from GWASCatalog...
2025/12/26 11:37:27  -Requesting (GET) trait information through the GWASCatalog API...
2025/12/26 11:37:27  -EFO trait api: https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005148
2025/12/26 11:37:28  -Status code: 200
2025/12/26 11:37:28  -Trait Name: type 2 diabetes mellitus
2025/12/26 11:37:28  -Trait URL: http://purl.obolibrary.org/obo/MONDO_0005148
2025/12/26 11:37:28  -Requesting (GET) GWAS associations through the GWASCatalog API...
2025/12/26 11:37:28  -associationsByTraitSummary API: https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005148/associations?projection=associationByEfoTrait
2025/12/26 11:37:28  -Note: this step might take a while...
2025/12/26 11:42:13  -Status code 200 OK: Retrieved data from GWASCatalog successffully ...
2025/12/26 11:42:13  -Loading json ...
2025/12/26 11:42:14  -Saving json to: ./GWASCatalog_MONDO_0005148_associationsByTraitSummary_text_20251226.json ...
2025/12/26 11:42:20  -Parsing json ...
2025/12/26 11:42:20  -Number of reported associations for MONDO_0005148 in GWASCatalog: 8794
2025/12/26 11:42:21  -Loading retrieved data into gwaslab Sumstats object ...
2025/12/26 11:42:21 Finished retrieving data from GWASCatalog...
2025/12/26 11:42:21 Finished retrieving data from GWASCatalog.
2025/12/26 11:42:21  -Retrieved 6927 variants using old API v1
2025/12/26 11:42:21  -Removed 2512 duplicate variants (kept 4415 unique)
2025/12/26 11:42:21  -Processed 4415 unique variants from GWAS Catalog
2025/12/26 11:42:21  -Retrieved 4415 associations from GWAS catalog.
2025/12/26 11:42:21  -Lead variants in known loci: 4415
2025/12/26 11:42:21  -Checking the minimum distance between identified lead variants and provided known variants...
2025/12/26 11:42:21  -Identified  4  known vairants in current sumstats...
2025/12/26 11:42:21  -Identified  0  novel vairants in current sumstats...
2025/12/26 11:42:21 Finished checking if lead variants are known.
```

| SNPID | CHR | POS | EA | NEA | STATUS | BETA | SE | P | DISTANCE_TO_KNOWN | KNOWN_ID | KNOWN_PUBMED_ID | KNOWN_AUTHOR | NOVEL | LOCATION_OF_KNOWN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:22068326_A_G | 1 | 21741833 | G | A | 3860099 | 0.0621 | 0.0103 | 1.629000e-09 | 0 | rs1825307 | 30718926 | Suzuki K | False | Same |
| 1:51103268_T_C | 1 | 50637596 | C | T | 3860099 | -0.0802 | 0.0120 | 2.519000e-11 | 0 | rs12031188 | 30718926 | Suzuki K | False | Same |
| 1:154309595_TA_T | 1 | 154337119 | TA | T | 3860399 | -0.0915 | 0.0166 | 3.289000e-08 | 1 | rs68062313 | 30718926 | Suzuki K | False | Upstream |
| 2:640986_CACAT_C | 2 | 640986 | C | CACAT | 3860399 | -0.0946 | 0.0150 | 2.665000e-10 | -1931 | rs7564708 | 34594039 | Sakaue S | False | Downstream |
