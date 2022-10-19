# Infer Genome Build Version
gwaslab use hapmaps chromosome and basepair position information to infer the reference genome build for sumstats.

Status codes (first two digits) will be changed based on the results

Example:
```
mysumstats.infer_build()

Wed Oct 19 11:01:01 2022  -Start to infer genome build version using hapmap3 SNPs...
Wed Oct 19 11:01:01 2022  -Loading Hapmap3 variants data...
Wed Oct 19 11:01:04 2022  -chr:pos will be used for matching...
Wed Oct 19 11:01:33 2022  -Matching variants for hg19: num_hg19= 1092441
Wed Oct 19 11:01:33 2022  -Matching variants for hg38: num_hg38= 15997
Wed Oct 19 11:01:33 2022  -Since num_hg19>num_hg38, assigning genome build hg19...
```
