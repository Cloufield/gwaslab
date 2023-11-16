# Clumping in GWASLab by calling PLINK2

GWASLab provides a clumping function using PLINK2. 
You can get clumping results without worrying about the file format.  

# Clumping 

!!! warning "You need to install PLINK2 first and add to your environment."

```
clump(  vcf=None, 
        bfile=None,
        scaled=False, 
        out="clumping_plink2", 
        overwrite=False, 
        n_cores=2, 
        chrom=None, 
        clump_p1=5e-8, 
        clump_p2=5e-8, 
        clump_r2=0.2, 
        clump_kb=250,
        log=Log())
```

|Option|DataType|Description|Default|
|-|-|-|-|
|`vcf`|`string`|path to reference VCF file (it will be converted to plink binary format). You need to specify either `vcf` or `bfile`|-|
|`bfile`|`string`|path to PLINK bfile. You need to specify either `vcf` or `bfile`|-|
|`scaled`|`boolean`|If tru, use MLOG10P instead of P|`False`|
|`out`|`string`|output file prefix|`clumping_plink2`|
|`overwrite`|`boolean`|if True, overwrite the existing bfile when vcf path is provided|`False`|
|`n_cores`|`int`|number of cores to use|2|
|`clump_p1`|`float`|clump_p1|5e-8|
|`clump_p2`|`float`|clump_p2|5e-8|
|`clump_r2`|`float`|clump_r2|0.2|
|`clump_kb`|`float`|clump_kb|250|


!!! info "Plink2 script" 
    
    ```
    # using P
    plink2 \
        --bfile {}\
        --chr {} \
        --clump {} \
        --clump-field P \
        --clump-snp-field SNPID \
        --clump-p1 {} \
        --clump-p2 {} \
        --clump-r2 {} \
        --clump-kb {} \
        --threads {} \
        --out {}
    ```
    
    ```
    # using MLOG10P
    plink2 \
        --bfile {}\
        --chr {} \
        --clump {} \
        --clump-field P \
        --clump-snp-field SNPID \
        --clump-p1 {} \
        --clump-p2 {} \
        --clump-r2 {} \
        --clump-kb {} \
        --threads {} \
        --out {}
    ```

# Example 

!!! example

```
clumps = mysumstats.clump(vcf= "/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz")



```