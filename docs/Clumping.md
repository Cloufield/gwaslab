# Clumping in GWASLab by calling PLINK2

GWASLab provides a wrapper for clumping using PLINK2. 
You can get clumping results without worrying about the file format.  

# Clumping 

!!! warning "You need to install PLINK2 first and add to your environment."

```python
clump(  plink2="plink2",
        vcf=None, 
        bfile=None,
        pfile=None，
        scaled=False, 
        out="clumping_plink2", 
        overwrite=False, 
        threads=2, 
        chrom=None, 
        clump_p1=5e-8, 
        clump_p2=5e-8, 
        clump_r2=0.2, 
        clump_kb=250,
        log=Log())
```

|Option|DataType|Description|Default|
|-|-|-|-|
|`plink2`|`string`|path to plink2|'plink2'|
|`vcf`|`string`|path to reference VCF file (it will be converted to plink binary format). You need to specify either `vcf` ,`bfile` or `bfile`|-|
|`bfile`|`string`|path to PLINK bfile. You need to specify either `vcf` ,`bfile` or `bfile`|-|
|`pfile`|`string`|path to PLINK bfile. You need to specify either `vcf` `bfile`， or `pfile`|-|
|`scaled`|`boolean`|If True, use MLOG10P instead of P|`False`|
|`out`|`string`|output file prefix|`clumping_plink2`|
|`overwrite`|`boolean`|if True, overwrite the existing bfile when vcf path is provided|`False`|
|`threads`|`int`|number of threads to use|2|
|`memory`|`int`|number of memory to use for plink (in MB) |None|
|`clump_p1`|`float`|clump_p1|5e-8|
|`clump_p2`|`float`|clump_p2|5e-8|
|`clump_r2`|`float`|clump_r2|0.2|
|`clump_kb`|`float`|clump_kb|250|


# Example 

!!! example

```python
clumps = mysumstats.clump(vcf= "/home/yunye/.gwaslab/EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz")

```

!!! example "Brisbane plot"
    See [Clumping in gwaslab](https://cloufield.github.io/gwaslab/util_ex_clumping/)

!!! info "Plink2 script for clumping" 

    ```python
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

    or

    ```python
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
