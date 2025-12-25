# Lead and novel variants

```python
import gwaslab as gl
```

```python
gl.show_version()
```

**stdout:**
```
2024/12/22 22:40:34 GWASLab v3.5.4 https://cloufield.github.io/gwaslab/
2024/12/22 22:40:34 (C) 2022-2024, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com
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

## Get lead variants

GWASLab will use MLOG10P first by default. If MLOG10P is not avaiable, it will look for P column.

```python
mysumstats.get_lead()
```

**stdout:**
```
2024/12/22 22:40:58 Start to extract lead variants...v3.5.4
2024/12/22 22:40:58  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:40:58  -Processing 1000000 variants...
2024/12/22 22:40:58  -Significance threshold : 5e-08
2024/12/22 22:40:58  -Sliding window size: 500  kb
2024/12/22 22:40:58  -Using P for extracting lead variants...
2024/12/22 22:40:58  -Found 543 significant variants in total...
2024/12/22 22:40:58  -Identified 4 lead variants!
2024/12/22 22:40:58 Finished extracting lead variants.
```

```
| SNPID | CHR | POS | EA | NEA | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 96739 | 1:22068326_A_G | 1 | 22068326 | G | A | 0.0621 | 0.0103 |
| 213860 | 1:51103268_T_C | 1 | 51103268 | C | T -0.0802 | 0.0120 |  |
| 534095 | 1:154309595_TA_T | 1 | 154309595 | TA | T -0.0915 | 0.0166 |  |
| 969974 | 2:640986_CACAT_C | 2 | 640986 | C | CACAT -0.0946 | 0.0150 |  |
|  | P | STATUS |  |  |  |  |  |
| 96739 | 1.629000e-09 | 1960099 |  |  |  |  |  |
| 213860 | 2.519000e-11 | 1960099 |  |  |  |  |  |
| 534095 | 3.289000e-08 | 1960399 |  |  |  |  |  |
| 969974 | 2.665000e-10 | 1960399 |  |  |  |  |  |
```

## Get lead variants with gene name annotation

```python
mysumstats.get_lead(anno=True)
```

**stdout:**
```
2024/12/22 22:40:58 Start to extract lead variants...v3.5.4
2024/12/22 22:40:58  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:40:58  -Processing 1000000 variants...
2024/12/22 22:40:58  -Significance threshold : 5e-08
2024/12/22 22:40:58  -Sliding window size: 500  kb
2024/12/22 22:40:58  -Using P for extracting lead variants...
2024/12/22 22:40:58  -Found 543 significant variants in total...
2024/12/22 22:40:58  -Identified 4 lead variants!
2024/12/22 22:40:58  -Annotating variants using references:ensembl
2024/12/22 22:40:58  -Annotating variants using references based on genome build:19
2024/12/22 22:40:58 Start to annotate variants with nearest gene name(s)...
2024/12/22 22:40:58  -Assigning Gene name using ensembl_hg19_gtf for protein coding genes
2024/12/22 22:40:58 Finished annotating variants with nearest gene name(s) successfully!
2024/12/22 22:40:58 Finished extracting lead variants.
```

```
| SNPID | CHR | POS | EA | NEA | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 96739 | 1:22068326_A_G | 1 | 22068326 | G | A | 0.0621 | 0.0103 |
| 213860 | 1:51103268_T_C | 1 | 51103268 | C | T -0.0802 | 0.0120 |  |
| 534095 | 1:154309595_TA_T | 1 | 154309595 | TA | T -0.0915 | 0.0166 |  |
| 969974 | 2:640986_CACAT_C | 2 | 640986 | C | CACAT -0.0946 | 0.0150 |  |
|  | P | STATUS | LOCATION | GENE |  |  |  |
| 96739 | 1.629000e-09 | 1960099 | 0 | USP48 |  |  |  |
| 213860 | 2.519000e-11 | 1960099 | 0 | FAF1 |  |  |  |
| 534095 | 3.289000e-08 | 1960399 | 0 | ATP8B2 |  |  |  |
| 969974 | 2.665000e-10 | 1960399 | 26349 | TMEM18 |  |  |  |
```

## Different window sizes

```python
mysumstats.get_lead(windowsizekb=1000)
```

**stdout:**
```
2024/12/22 22:40:59 Start to extract lead variants...v3.5.4
2024/12/22 22:40:59  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:40:59  -Processing 1000000 variants...
2024/12/22 22:40:59  -Significance threshold : 5e-08
2024/12/22 22:40:59  -Sliding window size: 1000  kb
2024/12/22 22:40:59  -Using P for extracting lead variants...
2024/12/22 22:40:59  -Found 543 significant variants in total...
2024/12/22 22:40:59  -Identified 4 lead variants!
2024/12/22 22:40:59 Finished extracting lead variants.
```

```
| SNPID | CHR | POS | EA | NEA | BETA | SE | \ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 96739 | 1:22068326_A_G | 1 | 22068326 | G | A | 0.0621 | 0.0103 |
| 213860 | 1:51103268_T_C | 1 | 51103268 | C | T -0.0802 | 0.0120 |  |
| 534095 | 1:154309595_TA_T | 1 | 154309595 | TA | T -0.0915 | 0.0166 |  |
| 969974 | 2:640986_CACAT_C | 2 | 640986 | C | CACAT -0.0946 | 0.0150 |  |
|  | P | STATUS |  |  |  |  |  |
| 96739 | 1.629000e-09 | 1960099 |  |  |  |  |  |
| 213860 | 2.519000e-11 | 1960099 |  |  |  |  |  |
| 534095 | 3.289000e-08 | 1960399 |  |  |  |  |  |
| 969974 | 2.665000e-10 | 1960399 |  |  |  |  |  |
```

## Different thresholds

```python
mysumstats.get_lead(sig_level=1e-10)
```

**stdout:**
```
2024/12/22 22:40:59 Start to extract lead variants...v3.5.4
2024/12/22 22:40:59  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:40:59  -Processing 1000000 variants...
2024/12/22 22:40:59  -Significance threshold : 1e-10
2024/12/22 22:40:59  -Sliding window size: 500  kb
2024/12/22 22:40:59  -Using P for extracting lead variants...
2024/12/22 22:40:59  -Found 1 significant variants in total...
2024/12/22 22:40:59  -Identified 1 lead variants!
2024/12/22 22:40:59 Finished extracting lead variants.
```

```
| SNPID | CHR | POS EA NEA | BETA | SE | P | \ |
| --- | --- | --- | --- | --- | --- | --- |
| 213860 | 1:51103268_T_C | 1 | 51103268 | C | T -0.0802 | 0.012 |
|  | STATUS |  |  |  |  |  |
| 213860 | 1960099 |  |  |  |  |  |
```

## Check if novel against a tabular file

```python
novel = mysumstats.get_novel(known="../0_sample_data//toy_data/known_loci.txt")
```

**stdout:**
```
2024/12/22 22:40:59 Start to check if lead variants are known...v3.5.4
2024/12/22 22:40:59  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:40:59 Start to extract lead variants...v3.5.4
2024/12/22 22:40:59  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:40:59  -Processing 1000000 variants...
2024/12/22 22:40:59  -Significance threshold : 5e-08
2024/12/22 22:40:59  -Sliding window size: 500  kb
2024/12/22 22:40:59  -Using P for extracting lead variants...
2024/12/22 22:40:59  -Found 543 significant variants in total...
2024/12/22 22:40:59  -Identified 4 lead variants!
2024/12/22 22:40:59 Finished extracting lead variants.
2024/12/22 22:40:59  -Lead variants in known loci: 2
2024/12/22 22:40:59  -Checking the minimum distance between identified lead variants and provided known variants...
2024/12/22 22:40:59  -Identified  2  known vairants in current sumstats...
2024/12/22 22:40:59  -Identified  2  novel vairants in current sumstats...
2024/12/22 22:40:59 Finished checking if lead variants are known.
```

```python
novel
```

```
| SNPID | CHR | POS | EA | NEA | BETA | SE | P | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:22068326_A_G | 1 | 22068326 | G | A | 0.0621 | 0.0103 | 1.629000e-09 |
| 1 | 1:51103268_T_C | 1 | 51103268 | C | T -0.0802 | 0.0120 | 2.519000e-11 |  |
| 2 | 1:154309595_TA_T | 1 | 154309595 | TA | T -0.0915 | 0.0166 | 3.289000e-08 |  |
| 3 | 2:640986_CACAT_C | 2 | 640986 | C | CACAT -0.0946 | 0.0150 | 2.665000e-10 |  |
|  | STATUS | DISTANCE_TO_KNOWN | KNOWN_ID | NOVEL LOCATION_OF_KNOWN |  |  |  |  |
| 0 | 1960099 | 29034942.0 | 1:51103268 | True | Upstream |  |  |  |
| 1 | 1960099 | 0.0 | 1:51103268 | False | Same |  |  |  |
| 2 | 1960399 | 0.0 | 1:154309595 | False | Same |  |  |  |
| 3 | 1960399 | NaN | <NA> | True | NoneOnThisChr |  |  |  |
```

## Ccheck against GWAS Catalog using EFO ID

```python
# EFO ID can be found on gwas catalog
mysumstats.get_novel(efo="MONDO_0005148")
```

**stdout:**
```
2024/12/22 22:41:00 Start to check if lead variants are known...v3.5.4
2024/12/22 22:41:00  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:41:00 Start to extract lead variants...v3.5.4
2024/12/22 22:41:00  -Current Dataframe shape : 1000000 x 9 ; Memory usage: 77.42 MB
2024/12/22 22:41:00  -Processing 1000000 variants...
2024/12/22 22:41:00  -Significance threshold : 5e-08
2024/12/22 22:41:00  -Sliding window size: 500  kb
2024/12/22 22:41:00  -Using P for extracting lead variants...
2024/12/22 22:41:00  -Found 543 significant variants in total...
2024/12/22 22:41:00  -Identified 4 lead variants!
2024/12/22 22:41:00 Finished extracting lead variants.
2024/12/22 22:41:00 Start to retrieve data using EFO: MONDO_0005148...
2024/12/22 22:41:00 Start to retrieve data from GWASCatalog...
2024/12/22 22:41:00  -Please make sure your sumstats is based on GRCh38...
2024/12/22 22:41:00  -Requesting (GET) trait information through the GWASCatalog API...
2024/12/22 22:41:00  -EFO trait api: https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005148
2024/12/22 22:41:01  -Status code: 200
2024/12/22 22:41:01  -Trait Name: type 2 diabetes mellitus
2024/12/22 22:41:01  -Trait URL: http://purl.obolibrary.org/obo/MONDO_0005148
2024/12/22 22:41:01  -Requesting (GET) GWAS associations through the GWASCatalog API...
2024/12/22 22:41:01  -associationsByTraitSummary API: https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0005148/associations?projection=associationByEfoTrait
2024/12/22 22:41:01  -Note: this step might take a while...
2024/12/22 22:50:50  -Status code 200 OK: Retrieved data from GWASCatalog successffully ...
2024/12/22 22:50:50  -Loading json ...
2024/12/22 22:50:52  -Parsing json ...
2024/12/22 22:50:52  -Number of reported associations for MONDO_0005148 in GWASCatalog: 5911
2024/12/22 22:50:52  -Loading retrieved data into gwaslab Sumstats object ...
2024/12/22 22:51:03 Finished retrieving data from GWASCatalog...
2024/12/22 22:51:03  -Retrieved 4593 associations from GWAS catalog.
2024/12/22 22:51:03  -Lead variants in known loci: 4593
2024/12/22 22:51:03  -Checking the minimum distance between identified lead variants and provided known variants...
2024/12/22 22:51:03  -Identified  4  known vairants in current sumstats...
2024/12/22 22:51:03  -Identified  0  novel vairants in current sumstats...
2024/12/22 22:51:03 Finished checking if lead variants are known.
```

```
| SNPID | CHR | POS | EA | NEA | BETA | SE | P | \ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1:22068326_A_G | 1 | 22068326 | G | A | 0.0621 | 0.0103 | 1.629000e-09 |
| 1 | 1:51103268_T_C | 1 | 51103268 | C | T -0.0802 | 0.0120 | 2.519000e-11 |  |
| 2 | 1:154309595_TA_T | 1 | 154309595 | TA | T -0.0915 | 0.0166 | 3.289000e-08 |  |
| 3 | 2:640986_CACAT_C | 2 | 640986 | C | CACAT -0.0946 | 0.0150 | 2.665000e-10 |  |
|  | STATUS | DISTANCE_TO_KNOWN | KNOWN_ID KNOWN_PUBMED_ID KNOWN_AUTHOR | NOVEL | \ |  |  |  |
| 0 | 1960099 | -279528 | rs34759301 | 34594039 | Sakaue S | False |  |  |
| 1 | 1960099 | -62054 | rs12088739 | 30054458 | Xue A | False |  |  |
| 2 | 1960399 | 12189 | rs1194606 | 32541925 | Vujkovic M | False |  |  |
| 3 | 1960399 | -1931 | rs7564708 | 34594039 | Sakaue S | False |  |  |
|  | LOCATION_OF_KNOWN |  |  |  |  |  |  |  |
| 0 | Downstream |  |  |  |  |  |  |  |
| 1 | Downstream |  |  |  |  |  |  |  |
| 2 | Upstream |  |  |  |  |  |  |  |
| 3 | Downstream |  |  |  |  |  |  |  |
```
