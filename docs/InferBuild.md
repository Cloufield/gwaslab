# Infer Genome Build

GWASLab uses chromosome and base pair position information for Hapmap3 SNPs to infer the reference genome build for sumstats.

```
.infer_build()
```

Reference genome build will be simply assigned based on the matching count.  Status codes (first two digits) will be changed based on the matching results.

!!! note 
    The results are more reliable if you have a reasonable amount of variants.


!!! example
    ```
    mysumstats.infer_build()
    ```
