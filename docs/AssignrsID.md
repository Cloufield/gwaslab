# Assigning rsID

GWASLab uses a two-step strategy (both steps are optional).

- For quick annotation, GWASLab iterates over a SNPID-rsID table and assign rsID by joining on SNPID (CHR:POS:REF:ALT) with sumstats. GWASLab provides a curated table  (1KG autosome variants). 
- For full annotation, GWASLab will query a large reference VCF file (dbSNP for example, >20GB ) by CHR, POS, NEA, EA. It will assign the ID in VCF file to sumstats if the CHR, POS and EN/NEA matches.

## Reference data

### SNPID-rsID table

GWASLab provides a download function `gl.download_ref()` and two curated tables which contains ~80M 1KG variants:

- `hg19` : `gl.download_ref("1kg_dbsnp151_hg19_auto")`
- `hg38` : `gl.download_ref("1kg_dbsnp151_hg38_auto")`

### VCF file

You can download this from dbSNP:

`hg19`

- `vcf`:https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
- `tbi`:https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

`hg38`

- `vcf`:https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz
- `tbi`:https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz.tbi


!!! note "1kg_dbsnp151_hg19_auto"
    ```
    ~/.gwaslab$ zcat 1kg_dbsnp151_hg19_auto.txt.gz |head
    SNPID   rsID    CHR     POS     NEA     EA
    1:10177:A:AC    rs367896724     1       10177   A       AC
    1:10235:T:TA    rs540431307     1       10235   T       TA
    1:10352:T:TA    rs555500075     1       10352   T       TA
    1:10505:A:T     rs548419688     1       10505   A       T
    1:10511:G:A     rs534229142     1       10511   G       A
    1:10539:C:A     rs537182016     1       10539   C       A
    1:10542:C:T     rs572818783     1       10542   C       T
    1:10579:C:A     rs538322974     1       10579   C       A
    1:10616:CCGCCGTTGCAAAGGCGCGCCG:C        rs376342519     1       10616   CCGCCGTTGCAAAGGCGCGCCG  C
    ```

!!! note "VCF file form dbSNP"
    ```
     zcat GCF_000001405.25.vcf.gz | head -100 | tail -10
    NC_000001.10    10059   rs1570391745    C       G       .       .       RS=1570391745;dbSNPBuildID=154;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9997,0.0003425|dbGaP_PopFreq:1,0
    NC_000001.10    10060   rs1639544146    C       CT      .       .       RS=1639544146;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10060   rs1639544159    CT      C       .       .       RS=1639544159;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=DEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10063   rs1010989343    A       C,G     .       .       RS=1010989343;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9928,0.004112,0.003084|Siberian:0.5,0.5,.|dbGaP_PopFreq:1,.,0
    NC_000001.10    10067   rs1489251879    T       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC      .       .       RS=1489251879;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=GnomAD:1,1.789e-05
    NC_000001.10    10067   rs1639545042    T       C       .       .       RS=1639545042;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10067   rs1639545104    TA      T       .       .       RS=1639545104;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10068   rs1639545079    A       T       .       .       RS=1639545079;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    NC_000001.10    10069   rs1570391755    A       C,G     .       .       RS=1570391755;dbSNPBuildID=154;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV;R5;GNO;FREQ=KOREAN:0.9966,.,0.003425|dbGaP_PopFreq:1,0,0
    NC_000001.10    10069   rs1639545200    A       AC      .       .       RS=1639545200;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL;R5;GNO;FREQ=dbGaP_PopFreq:1,0
    ```

## Usage

```
mysumstats.basic_check()
mysumstats.assign_rsid( 
                        ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg19_auto"),
                        ref_rsid_vcf = "/home/yunye/mydata/d_disk/dbsnp/GCF_000001405.25.vcf.gz",
                        chr_dict = gl.get_number_to_NC(build="19"),
                        n_cores = 2)
```
!!! note
    Please always run `.basic_check()` first. This will convert the data to right data type, standardize and normalize the sumstats.

## Options

- `ref_rsid_tsv` : tsv file for annotation of common used variants. Using SNPID (like 1:725932:G:A)
- `ref_rsid_vcf` : vcf file for annotation of other variants.
- `chr_dict`: a dictionary for to convert 1-25 to CHR in vcf files. For example, the notation in dbSNP vcf file is based on RefSeq (like NC_000001.10). `gwaslab` provides built-in conversion dictionaries.   `gl.get_number_to_NC(build="19")` and `gl.get_number_to_NC(build="19")`
- `n_cores`: number of cores to use.


!!! note "Conversion for RefSeq sequence"

    ```
    gl.get_number_to_NC(build="19")
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

    gl.get_number_to_NC(build="19")
    {1: 'NC_000001.11',
     2: 'NC_000002.12',
     3: 'NC_000003.12',
     4: 'NC_000004.12',
     5: 'NC_000005.10',
     6: 'NC_000006.12',
     7: 'NC_000007.14',
     8: 'NC_000008.11',
     9: 'NC_000009.12',
     10: 'NC_000010.11',
     11: 'NC_000011.10',
     12: 'NC_000012.12',
     13: 'NC_000013.11',
     14: 'NC_000014.9',
     15: 'NC_000015.10',
     16: 'NC_000016.10',
     17: 'NC_000017.11',
     18: 'NC_000018.10',
     19: 'NC_000019.10',
     20: 'NC_000020.11',
     21: 'NC_000021.9',
     22: 'NC_000022.11',
     23: 'NC_000023.11',
     24: 'NC_000024.1',
     25: 'NC_012920.1'}
    ```

## Example
!!! example
    ```
    # download ref SNPID-rsID table first
    gl.download_ref("1kg_dbsnp151_hg19_auto") 
    
    # if you want to annotate as much as possible. Please download the very large dbSNP vcf file.
    
    mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
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
                 n="N",
                 nrows=10000)
    
    # run basic_check first
    mysumstats.basic_check() 
    ```
    
    ![image](https://user-images.githubusercontent.com/40289485/211799713-032b3571-ae01-4307-909b-973cdf043e17.png)
    
    ```
    # if you SNPID is like 1:725932_G_A	, you can use fix_id to fix the separator.
    mysumstats.fix_id(fixsep=True)
    ```
    
    ![image](https://user-images.githubusercontent.com/40289485/211799561-dddf7649-09fc-4eb5-b8e3-d7301a8944d1.png)
    
    ```
    # rsID annotation
    mysumstats.assign_rsid( n_cores = 2,
                            ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg19_auto"),
                            ref_rsid_vcf ="/home/yunye/mydata/d_disk/dbsnp/GCF_000001405.25.vcf.gz",
                            chr_dict = gl.get_number_to_NC(build="19"))
    ```
    
    ```
    Wed Jan 11 20:05:38 2023 Start to annotate rsID based on chromosome and position information...
    Wed Jan 11 20:05:38 2023  -Current Dataframe shape : 10000  x  12
    Wed Jan 11 20:05:38 2023  -SNPID-rsID text file: /home/yunye/.gwaslab/1kg_dbsnp151_hg19_auto.txt.gz
    Wed Jan 11 20:05:38 2023  -10000 rsID could be possibly fixed...
    Wed Jan 11 20:05:38 2023  -Setting block size:  5000000
    Wed Jan 11 20:05:38 2023  -Loading block: 0   1   2   3   4   5   6   7   8   9   10   11   12   13   14   15   
    Wed Jan 11 20:10:16 2023  -rsID Annotation for 58 need to be fixed!
    Wed Jan 11 20:10:16 2023  -Annotated 9942 rsID successfully!
    Wed Jan 11 20:10:16 2023 Start to assign rsID using vcf...
    Wed Jan 11 20:10:16 2023  -Current Dataframe shape : 10000  x  13
    Wed Jan 11 20:10:16 2023  -CPU Cores to use : 2
    Wed Jan 11 20:10:16 2023  -Reference VCF file: /home/yunye/mydata/d_disk/dbsnp/GCF_000001405.25.vcf.gz
    Wed Jan 11 20:10:16 2023  -Assigning rsID based on chr:pos and ref:alt/alt:ref...
    Wed Jan 11 20:10:17 2023  -rsID Annotation for 1 need to be fixed!
    Wed Jan 11 20:10:17 2023  -Annotated 57 rsID successfully!
    ```
    
    As you can see, SNPID-rsID (`1kg_dbsnp151_hg19_auto`) annotated 9942 rsID and the large reference VCF file (from dbSNP) annotated additonal 57 rare rsID.
    
    ![image](https://user-images.githubusercontent.com/40289485/211800319-68f33eaa-4c48-4ba4-aa52-afb9ac145dee.png)

    
