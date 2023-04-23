# Infer Genome Build

GWASLab use the chromosome and basepair position information for Hapmap3 SNPs to infer the reference genome build for sumstats.

Reference genome build will be simply assigned based on the matching count. (Note: the results are more reliable if you have more than 10,000 variants)

Status codes (first two digits) will be changed based on the matching results.

!!! example
    ```
    mysumstats.infer_build()
    
    Wed Oct 19 11:01:01 2022  -Start to infer genome build version using hapmap3 SNPs...
    Wed Oct 19 11:01:01 2022  -Loading Hapmap3 variants data...
    Wed Oct 19 11:01:04 2022  -chr:pos will be used for matching...
    Wed Oct 19 11:01:33 2022  -Matching variants for hg19: num_hg19= 1092441
    Wed Oct 19 11:01:33 2022  -Matching variants for hg38: num_hg38= 15997
    Wed Oct 19 11:01:33 2022  -Since num_hg19>num_hg38, assigning genome build hg19...
    ```
