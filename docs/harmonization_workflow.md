# Harmonization

!!! example
    ```python
    import gwaslab as gl
    ```

## Load test data

!!! example
    ```python
    mysumstats = gl.Sumstats("/home/yunye/work/gwaslab/examples/toy_data/to_harmonize.tsv",fmt="gwaslab",other=["NOTE"])
    ```

**stdout:**
```
Fri Feb  2 19:36:24 2024 GWASLab v3.4.38 https://cloufield.github.io/gwaslab/
Fri Feb  2 19:36:24 2024 (C) 2022-2024, Yunye He, Kamatani Lab, MIT License, gwaslab@gmail.com
Fri Feb  2 19:36:24 2024 Start to load format from formatbook....
Fri Feb  2 19:36:24 2024  -gwaslab format meta info:
Fri Feb  2 19:36:24 2024   - format_name  : gwaslab
Fri Feb  2 19:36:24 2024   - format_source  : https://cloufield.github.io/gwaslab/
Fri Feb  2 19:36:24 2024   - format_version  : 20231220_v4
Fri Feb  2 19:36:24 2024 Start to initialize gl.Sumstats from file :/home/yunye/work/gwaslab/examples/toy_data/to_harmonize.tsv
Fri Feb  2 19:36:24 2024  -Reading columns          : CHR,NEA,SE,P,DIRECTION,BETA,SNPID,NOTE,N,EAF,EA,POS
Fri Feb  2 19:36:24 2024  -Renaming columns to      : CHR,NEA,SE,P,DIRECTION,BETA,SNPID,NOTE,N,EAF,EA,POS
Fri Feb  2 19:36:24 2024  -Current Dataframe shape : 15  x  12
Fri Feb  2 19:36:24 2024  -Initiating a status column: STATUS ...
Fri Feb  2 19:36:24 2024  #WARNING! Version of genomic coordinates is unknown...
Fri Feb  2 19:36:24 2024 Start to reorder the columns...v3.4.38
Fri Feb  2 19:36:24 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:24 2024  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,P,N,DIRECTION,STATUS,NOTE
Fri Feb  2 19:36:24 2024 Finished reordering the columns.
Fri Feb  2 19:36:24 2024  -Column  : SNPID  CHR    POS   EA       NEA      EAF     BETA    SE      P       N     DIRECTION STATUS   NOTE  
Fri Feb  2 19:36:24 2024  -DType   : object string int64 category category float64 float64 float64 float64 int64 object    category object
Fri Feb  2 19:36:24 2024  -Verified: T      F      T     T        T        T       T       T       T       T     T         T        NA    
Fri Feb  2 19:36:24 2024  #WARNING! Columns with possibly incompatable dtypes: CHR
Fri Feb  2 19:36:24 2024  -Current Dataframe memory usage: 19.94 MB
Fri Feb  2 19:36:24 2024 Finished loading data successfully!
```

A list of variants that needs to be harmonized. Issues are described in NOTE column.

!!! example
    ```python
    mysumstats.data
    ```

| SNPID | CHR | POS | EA | NEA | EAF | BETA | SE | P | N | DIRECTION | STATUS | NOTE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:1066403:T:C | 1 | 1066403 | C | T | 0.5276 | 0.0043 | 0.0109 | 0.68910 | 191764 | ++-+ | 9999999 | Clean |
| 1:1163266:G:A | 1 | 1163266 | G | A | 0.8355 | 0.0122 | 0.0120 | 0.30810 | 191764 | -+++ | 9999999 | Flip |
| 1:1213224:T:A | 1 | 1213224 | T | A | 0.3345 | 0.0071 | 0.0095 | 0.45280 | 191764 | +++- | 9999999 | Flip + Palindromic |
| 1:3066761:A:T | 1 | 3066761 | T | A | 0.1245 | -0.0066 | 0.0157 | 0.67420 | 191764 | -+++ | 9999999 | Palindromic + No flip |
| 1:3997271:T:TTTTA | 1 | 3997271 | TTTTA | T | 0.3840 | 0.0188 | 0.0135 | 0.16300 | 191764 | ++++ | 9999999 | Indel + Both on genome + No flip |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 3:183629306:T:TTCTC | 3 | 183629306 | TTCTC | T | 0.0082 | -0.0083 | 0.0555 | 0.88060 | 191764 | +--+ | 9999999 | Indel  + Both on genome + No information |
| 4:99731866:G:C | 4 | 99731866 | G | C | 0.3883 | 0.0052 | 0.0093 | 0.57360 | 191764 | +-++ | 9999999 | Palindromic+ Flip + Different strand (REF G, T... |
| 8:89935201:C:G | 8 | 89935201 | G | C | 0.3533 | 0.0103 | 0.0094 | 0.27560 | 191764 | ++++ | 9999999 | Palindromic + Different strand (REF A; TOMMO G... |
| X:4206466:G:C | X | 4206466 | G | C | 0.8727 | -0.0055 | 0.0104 | 0.59840 | 191764 | -+-- | 9999999 | Palindromic + Flip + No information |
| X:5053229:A:T | X | 5053229 | T | A | 0.9747 | 0.0218 | 0.0225 | 0.33180 | 191764 | ++++ | 9999999 | Palindromic + No information |

## Perform harmonization

- ref_seq : reference genome fasta file for allele alignment
- ref_infer : vcf file with allele frequency information for inferring strand and comparing allele frequency 
- ref_alt_freq : field in INFO of vcf file for alternative allele frequency

!!! example
    ```python
    mysumstats.harmonize(   basic_check=True,
                            n_cores=1,
                            ref_seq=gl.get_path("ucsc_genome_hg19"),
                            ref_infer=gl.get_path("1kg_eas_hg19"),
                            ref_alt_freq="AF")
    ```

**stdout:**
```
Fri Feb  2 19:36:24 2024 Start to check SNPID/rsID...v3.4.38
Fri Feb  2 19:36:24 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:24 2024  -Checking SNPID data type...
Fri Feb  2 19:36:24 2024  -Converting SNPID to pd.string data type...
Fri Feb  2 19:36:24 2024  -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)
Fri Feb  2 19:36:25 2024 Finished checking SNPID/rsID.
Fri Feb  2 19:36:25 2024 Start to fix chromosome notation (CHR)...v3.4.38
Fri Feb  2 19:36:25 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:25 2024  -Checking CHR data type...
Fri Feb  2 19:36:25 2024  -Variants with standardized chromosome notation: 13
Fri Feb  2 19:36:25 2024  -Variants with fixable chromosome notations: 2
Fri Feb  2 19:36:25 2024  -No unrecognized chromosome notations...
Fri Feb  2 19:36:25 2024  -Identifying non-autosomal chromosomes : X, Y, and MT ...
Fri Feb  2 19:36:25 2024  -Identified  2  variants on sex chromosomes...
Fri Feb  2 19:36:25 2024  -Standardizing sex chromosome notations: X to 23...
Fri Feb  2 19:36:27 2024 Finished fixing chromosome notation (CHR).
Fri Feb  2 19:36:27 2024 Start to fix basepair positions (POS)...v3.4.38
Fri Feb  2 19:36:27 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:27 2024  -Converting to Int64 data type ...
Fri Feb  2 19:36:28 2024  -Position bound:(0 , 250,000,000)
Fri Feb  2 19:36:28 2024  -Removed outliers: 0
Fri Feb  2 19:36:28 2024 Finished fixing basepair positions (POS).
Fri Feb  2 19:36:28 2024 Start to fix alleles (EA and NEA)...v3.4.38
Fri Feb  2 19:36:28 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:28 2024  -Converted all bases to string datatype and UPPERCASE.
Fri Feb  2 19:36:28 2024  -Variants with bad EA  : 0
Fri Feb  2 19:36:28 2024  -Variants with bad NEA : 0
Fri Feb  2 19:36:28 2024  -Variants with NA for EA or NEA: 0
Fri Feb  2 19:36:28 2024  -Variants with same EA and NEA: 0
Fri Feb  2 19:36:28 2024  -Detected 0 variants with alleles that contain bases other than A/C/T/G .
Fri Feb  2 19:36:30 2024 Finished fixing alleles (EA and NEA).
Fri Feb  2 19:36:30 2024 Start to perform sanity check for statistics...v3.4.38
Fri Feb  2 19:36:30 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:30 2024  -Comparison tolerance for floats: 1e-07
Fri Feb  2 19:36:30 2024  -Checking if 0 <= N <= 2147483647 ...
Fri Feb  2 19:36:30 2024  -Removed 0 variants with bad/na N.
Fri Feb  2 19:36:30 2024  -Checking if -1e-07 < EAF < 1.0000001 ...
Fri Feb  2 19:36:30 2024  -Removed 0 variants with bad/na EAF.
Fri Feb  2 19:36:30 2024  -Checking if -1e-07 < P < 1.0000001 ...
Fri Feb  2 19:36:30 2024  -Removed 0 variants with bad/na P.
Fri Feb  2 19:36:30 2024  -Checking if -100.0000001 < BETA < 100.0000001 ...
Fri Feb  2 19:36:30 2024  -Removed 0 variants with bad/na BETA.
Fri Feb  2 19:36:30 2024  -Checking if -1e-07 < SE < inf ...
Fri Feb  2 19:36:30 2024  -Removed 0 variants with bad/na SE.
Fri Feb  2 19:36:30 2024  -Checking STATUS and converting STATUS to categories....
Fri Feb  2 19:36:30 2024  -Removed 0 variants with bad statistics in total.
Fri Feb  2 19:36:30 2024  -Data types for each column:
Fri Feb  2 19:36:30 2024  -Column  : SNPID  CHR   POS   EA       NEA      EAF     BETA    SE      P       N     DIRECTION STATUS   NOTE  
Fri Feb  2 19:36:30 2024  -DType   : string Int64 Int64 category category float32 float64 float64 float64 Int64 object    category object
Fri Feb  2 19:36:30 2024  -Verified: T      T     T     T        T        T       T       T       T       T     T         T        NA    
Fri Feb  2 19:36:30 2024 Finished sanity check for statistics.
Fri Feb  2 19:36:30 2024 Start to normalize indels...v3.4.38
Fri Feb  2 19:36:30 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:30 2024  -No available variants to normalize..
Fri Feb  2 19:36:30 2024 Finished normalizing variants successfully!
Fri Feb  2 19:36:30 2024 Start to reorder the columns...v3.4.38
Fri Feb  2 19:36:30 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:30 2024  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,P,N,DIRECTION,STATUS,NOTE
Fri Feb  2 19:36:30 2024 Finished reordering the columns.
Fri Feb  2 19:36:31 2024 Start to check if NEA is aligned with reference sequence...v3.4.38
Fri Feb  2 19:36:31 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 19.94 MB
Fri Feb  2 19:36:31 2024  -Reference genome FASTA file: /home/yunye/.gwaslab/hg19.fa
Fri Feb  2 19:36:31 2024  -Checking records: 1  2  3  4  5  6  7  X  8  9  10  11  12  13  14  15  16  17  18  20  Y  19  22  21  M  
Fri Feb  2 19:36:47 2024  -Variants allele on given reference sequence :  6
Fri Feb  2 19:36:47 2024  -Variants flipped :  6
Fri Feb  2 19:36:47 2024   -Raw Matching rate :  80.00%
Fri Feb  2 19:36:47 2024  -Variants inferred reverse_complement :  0
Fri Feb  2 19:36:47 2024  -Variants inferred reverse_complement_flipped :  0
Fri Feb  2 19:36:47 2024  -Both allele on genome + unable to distinguish :  3
Fri Feb  2 19:36:47 2024  -Variants not on given reference sequence :  0
Fri Feb  2 19:36:47 2024 Finished checking if NEA is aligned with reference sequence.
Fri Feb  2 19:36:47 2024 Start to adjust statistics based on STATUS code...v3.4.38
Fri Feb  2 19:36:47 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 0.00 MB
Fri Feb  2 19:36:47 2024 Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...v3.4.38
Fri Feb  2 19:36:47 2024  -Flipping 6 variants...
Fri Feb  2 19:36:47 2024  -Swapping column: NEA <=> EA...
Fri Feb  2 19:36:47 2024  -Flipping column: BETA = - BETA...
Fri Feb  2 19:36:47 2024  -Flipping column: DIRECTION +-?0 <=> -+?0 ...
Fri Feb  2 19:36:47 2024  -Flipping column: EAF = 1 - EAF...
Fri Feb  2 19:36:47 2024  -Changed the status for flipped variants : xxxxx[35]x -> xxxxx[12]x
Fri Feb  2 19:36:47 2024 Finished adjusting.
Fri Feb  2 19:36:47 2024 Start to infer strand for palindromic SNPs/align indistinguishable indels...v3.4.38
Fri Feb  2 19:36:47 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 0.00 MB
Fri Feb  2 19:36:47 2024  -Number of threads/cores to use: 1
Fri Feb  2 19:36:47 2024  -Reference VCF: /home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
Fri Feb  2 19:36:47 2024  -Checking prefix for chromosomes in vcf files...
Fri Feb  2 19:36:47 2024  -No prefix for chromosomes in the VCF files.
Fri Feb  2 19:36:47 2024  -Field for alternative allele frequency in VCF INFO: AF
Fri Feb  2 19:36:47 2024  -Identified  8  palindromic SNPs...
Fri Feb  2 19:36:47 2024  -After filtering by MAF< 0.4 , 6 palindromic SNPs with unknown strand will be inferred...
Fri Feb  2 19:36:47 2024   -Non-palindromic :  2
Fri Feb  2 19:36:47 2024   -Palindromic SNPs on + strand:  2
Fri Feb  2 19:36:47 2024   -Palindromic SNPs on - strand and needed to be flipped: 2
Fri Feb  2 19:36:47 2024   -Palindromic SNPs with MAF not available to infer :  2
Fri Feb  2 19:36:47 2024   -Palindromic SNPs with no macthes or no information :  1
Fri Feb  2 19:36:47 2024  -Identified  3  indistinguishable Indels...
Fri Feb  2 19:36:47 2024  -Indistinguishable indels will be inferred from reference vcf ref and alt...
Fri Feb  2 19:36:47 2024  -DAF tolerance: 0.2
Fri Feb  2 19:36:47 2024   -Indels ea/nea match reference :  1
Fri Feb  2 19:36:47 2024   -Indels ea/nea need to be flipped :  1
Fri Feb  2 19:36:47 2024   -Indels with no macthes or no information :  1
Fri Feb  2 19:36:47 2024 Finished inferring strand for palindromic SNPs/align indistinguishable indels.
Fri Feb  2 19:36:47 2024 Start to adjust statistics based on STATUS code...v3.4.38
Fri Feb  2 19:36:47 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 0.00 MB
Fri Feb  2 19:36:47 2024 Start to flip allele-specific stats for standardized indels with status xxxx[123][67][6]: ALT->EA , REF->NEA...v3.4.38
Fri Feb  2 19:36:47 2024  -Flipping 1 variants...
Fri Feb  2 19:36:47 2024  -Swapping column: NEA <=> EA...
Fri Feb  2 19:36:47 2024  -Flipping column: BETA = - BETA...
Fri Feb  2 19:36:47 2024  -Flipping column: DIRECTION +-?0 <=> -+?0 ...
Fri Feb  2 19:36:47 2024  -Flipping column: EAF = 1 - EAF...
Fri Feb  2 19:36:47 2024  -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4
Fri Feb  2 19:36:47 2024 Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...v3.4.38
Fri Feb  2 19:36:47 2024  -Flipping 2 variants...
Fri Feb  2 19:36:47 2024  -Flipping column: BETA = - BETA...
Fri Feb  2 19:36:47 2024  -Flipping column: DIRECTION +-?0 <=> -+?0 ...
Fri Feb  2 19:36:47 2024  -Flipping column: EAF = 1 - EAF...
Fri Feb  2 19:36:47 2024  -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2
Fri Feb  2 19:36:47 2024 Finished adjusting.
```

**stdout:**
```
Fri Feb  2 19:36:48 2024 Start to sort the genome coordinates...v3.4.38
Fri Feb  2 19:36:48 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 0.00 MB
Fri Feb  2 19:36:48 2024 Finished sorting coordinates.
Fri Feb  2 19:36:48 2024 Start to reorder the columns...v3.4.38
Fri Feb  2 19:36:48 2024  -Current Dataframe shape : 15 x 13 ; Memory usage: 0.00 MB
Fri Feb  2 19:36:48 2024  -Reordering columns to    : SNPID,CHR,POS,EA,NEA,EAF,BETA,SE,P,N,DIRECTION,STATUS,NOTE
Fri Feb  2 19:36:48 2024 Finished reordering the columns.
```

```
<gwaslab.g_Sumstats.Sumstats at 0x7f6c69f86550>
```

All variants were checked and harmonized based on the reference datasets. The manipulations are reflected in STATUS column.

!!! example
    ```python
    mysumstats.data
    ```

| SNPID | CHR | POS | EA | NEA | EAF | BETA | SE | P | N | DIRECTION | STATUS | NOTE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1:1066403:T:C | 1 | 1066403 | C | T | 0.5276 | 0.0043 | 0.0109 | 0.68910 | 191764 | ++-+ | 9960000 | Clean |
| 1:1163266:G:A | 1 | 1163266 | A | G | 0.1645 | -0.0122 | 0.0120 | 0.30810 | 191764 | +--- | 9960010 | Flip |
| 1:1213224:T:A | 1 | 1213224 | A | T | 0.6655 | -0.0071 | 0.0095 | 0.45280 | 191764 | ---+ | 9960011 | Flip + Palindromic |
| 1:3066761:A:T | 1 | 3066761 | T | A | 0.1245 | -0.0066 | 0.0157 | 0.67420 | 191764 | -+++ | 9960001 | Palindromic + No flip |
| 1:3997271:T:TTTTA | 1 | 3997271 | TTTTA | T | 0.3840 | 0.0188 | 0.0135 | 0.16300 | 191764 | ++++ | 9960363 | Indel + Both on genome + No flip |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 3:183629306:T:TTCTC | 3 | 183629306 | TTCTC | T | 0.0082 | -0.0083 | 0.0555 | 0.88060 | 191764 | +--+ | 9960368 | Indel  + Both on genome + No information |
| 4:99731866:G:C | 4 | 99731866 | C | G | 0.3883 | 0.0052 | 0.0093 | 0.57360 | 191764 | +-++ | 9960012 | Palindromic+ Flip + Different strand (REF G, T... |
| 8:89935201:C:G | 8 | 89935201 | G | C | 0.6467 | -0.0103 | 0.0094 | 0.27560 | 191764 | ---- | 9960002 | Palindromic + Different strand (REF A; TOMMO G... |
| X:4206466:G:C | 23 | 4206466 | C | G | 0.1273 | 0.0055 | 0.0104 | 0.59840 | 191764 | +-++ | 9960018 | Palindromic + Flip + No information |
| X:5053229:A:T | 23 | 5053229 | T | A | 0.9747 | 0.0218 | 0.0225 | 0.33180 | 191764 | ++++ | 9960008 | Palindromic + No information |

## Check summary

!!! example
    ```python
    mysumstats.summary()
    ```

| Values | Percentage |
| --- | --- |
| 15 | <NA> |
| <NA> |  |
| <NA> |  |
| <NA> |  |
| 0 | 0.0 |
| ... | ... |
| 6.67 |  |
| 6.67 |  |
| 6.67 |  |
| 6.67 |  |
| 6.67 |  |

!!! example
    ```python
    mysumstats.lookup_status()
    ```

| Genome_Build | rsID&SNPID | CHR&POS | Stadardize&Normalize | Align | Panlidromic_SNP&Indel | Count | Percentage(%) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized SNP | Match: NEA=REF | Not_palindromic_SNPs | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized SNP | Match: NEA=REF | Palindromic+strand | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized SNP | Match: NEA=REF | Palindromic-strand_fixed | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized SNP | Match: NEA=REF | Indistinguishable | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized SNP | Match: NEA=REF | No_matching_or_no_info | 1 | 6.67 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized & normalized indel | Match: NEA=REF | Unchecked | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized & normalized indel | Flipped_fixed | Unchecked | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized & normalized indel | Both_alleles_on_ref+indistinguishable | Indel_match | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized & normalized indel | Both_alleles_on_ref+indistinguishable | Indel_flipped_fixed | 1 | 6.67 |
| Unchecked | rsid unknown & SNPID valid | CHR valid & POS valid | standardized & normalized indel | Both_alleles_on_ref+indistinguishable | No_matching_or_no_info | 1 | 6.67 |

## VCF : NCBI Sequence Identifiers

For some reference file the chromosome notation is not in the form of 1, or chr1, instead, they are using NCBI Sequence Identifiers (for example: NC_000001.11).

See https://www.ncbi.nlm.nih.gov/genbank/sequenceids/

For example: GCF_000001405.25.vcf.gz

NC_000001.10    10059   rs1570391745    C       G       .       .       RS=1570391745;dbSNPBuildID=154;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9997,0.0003425|dbGaP_PopFreq:1,0

NC_000001.10    10060   rs1639544146    C       CT      .       .       RS=1639544146;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0

NC_000001.10    10060   rs1639544159    CT      C       .       .       RS=1639544159;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=DEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0

 In such case, gwaslab provide built-in conversion table. 

!!! example
    ```python
    gl.get_number_to_NC(build="19")
    ```

```
{1: 'NC_000001.10',
 2: 'NC_000002.11',
 3: 'NC_000003.11',
 4: 'NC_000004.11',
 5: 'NC_000005.9',
 6: 'NC_000006.11',
 7: 'NC_000007.13',
 8: 'NC_000008.10',
 9: 'NC_000009.11',
 10: 'NC_000010.10',
 11: 'NC_000011.9',
 12: 'NC_000012.11',
 13: 'NC_000013.10',
 14: 'NC_000014.8',
 15: 'NC_000015.9',
 16: 'NC_000016.9',
 17: 'NC_000017.10',
 18: 'NC_000018.9',
 19: 'NC_000019.9',
 20: 'NC_000020.10',
 21: 'NC_000021.8',
 22: 'NC_000022.10',
 23: 'NC_000023.10',
 24: 'NC_000024.9',
 25: 'NC_012920.1'}
```

Specify it in assignrsid_args and inferstrand_args for the all-in-one function:

!!! example
    ```python
    mysumstats.harmonize(   basic_check=False,
                            n_cores=3,
                            ref_infer="/home/yunye/mydata/d_disk/dbsnp/GCF_000001405.25.vcf.gz",
                            inferstrand_args= {"chr_dict" : gl.get_number_to_NC(build="19")})
    ```
