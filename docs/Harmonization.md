# Harmonization

GWASLab provides reference-dependent harmonization functions.

## Methods summary

| Sumstats Methods       | Options                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | Description                                                                                                                       |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `.check_ref()`         | `ref_path`,<br /> `chr_dict=get_chr_to_number()`                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | Check alignment with a reference genome FASTA file                                                                                |
| `.assign_rsid()`       | `ref_rsid_tsv`,<br /> `ref_rsid_vcf`,<br /> `n_cores=1`, <br />`chunksize=5000000`, <br />`chr_dict=get_number_to_chr()`, <br />`overwrite="empty"`                                                                                                                                                                                                                                                                                                                                                                       | Annotate rsID using a reference tabular file or VCF/BCF file                                                                      |
| `.infer_strand()`      | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`remove_snp=""`,<br />`daf_tolerance=0.2`, ,<br />`mode="pi"`,<br />`n_cores=1`,<br />`remove_indel=""`                                                                                                                                                                                                                                                                                                                                            | Infer the strand of palindromic variants/indistinguishable indels using reference VCF/BCF files based on allele frequency in INFO |
| `.check_daf()`         | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`                                                                                                                                                                                                                                                                                                                                                                                                                                        | Calculate difference in allele frequencies  with a reference VCF/BCF file                                                         |
| `.flip_allele_stats()` |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | After alignment and inferring, flip the alleles and allele-specific statistics to harmonise the variants.                         |
| `.harmonize()`         | `basic_check=True`, <br /> `ref_seq=None`,<br />`ref_rsid_tsv=None`,<br />`ref_rsid_vcf=None`,<br />`ref_infer=None`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`,<br />`remove=False`,<br />`checkref_args={}`,<br />`removedup_args={}`,<br />`assignrsid_args={}`,<br />`inferstrand_args={}`,<br />`flipallelestats_args={}`,<br />`fixid_args={}`,<br />`fixchr_agrs={}`,<br />`fixpos_args={}`,<br />`fixallele_args={}`,<br />`sanitycheckstats_args={}`,<br />`normalizeallele_args={}` | All-in-one function for harmonization                                                                                             |

## Align NEA with REF in the reference genome

`.check_ref()`:  Check if NEA is aligned with the reference sequence. After checking, the tracking status code will be changed accordingly. 

| `.check_ref()` options                  | DataType   | Description                                                                               | Default                            |
|-----------------------------------------|------------|-------------------------------------------------------------------------------------------|------------------------------------|
| `ref_path`                              | `string`   | path to the reference genome FASTA file.                                                  |                                    |
| `ref_seq_mode` (avaiable since v3.4.42) | `v` or `s` | `v` for vectorized implementation (faster); `s` for single row iteration mode                     | `v` since v3.4.42 (except v3.4.43) |
| `chr_dict`                              | `dict`     | a conversion dictionary for chromosome notations in reference FASTA and those in sumstats | `gl.get_chr_to_number()`           |

!!! example
    ```python
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    ```
!!! note
    `check_ref()` only change the status code. Use [flip function](#flipping-based-on-status-code) `.flip_allele_stats()` to flip the allele-specific stats.

## Assign rsID according to CHR, POS, REF/ALT

`.assign_rsid()` : Annotated variants with rsID using a reference tsv file (1KG variants) and reference vcf file (tabix indexed, entire dbSNP).

See [https://cloufield.github.io/gwaslab/AssignrsID/](https://cloufield.github.io/gwaslab/AssignrsID/)

!!! example
    ```python
    mysumstats.assign_rsid(ref_rsid_tsv = gl.get_path("1kg_dbsnp151_hg19_auto"), 
                           ref_rsid_vcf = "/home/yunye/mydata/d_disk/dbsnp/GCF_000001405.25.vcf.gz",
                           chr_dict = gl.get_number_to_NC(build="19") 
                           n_cores=1)
    ```

- For TSV file, variants will be matched using SNPID (CHR:POS:NEA:EA) for quick assigning.
- For VCF file, GWASLab will first extract all variants in the reference file with matching CHR and POS. And then compare EA/NEA in sumstats with REF/ALT in reference vcf. When matching, it will annotate the variant in sumstats with the matching rsID in reference vcf.  

## Check palindromic SNPs or indistinguishable Indels

`.infer_strand()`:

- Infer the strand for palindromic SNPs (AT, or CG) with MAF <`maf_threshold`. 
- Checking the alignment status of indels with the REF allele in a reference VCF/BCF file and check if the allele frequencies are consistent (DAF < `daf_tolerance`).

!!! note "DAF : Difference between Effect allele frequency in sumstats and Reference ALT frequency in reference VCF/BCF file"

| `.infer_strand()` options | DataType          | Description                                                                                                                 | Default                  |
|---------------------------|-------------------|-----------------------------------------------------------------------------------------------------------------------------|--------------------------|
| `ref_infer`               | `string`          | path to the reference VCF/BCF file (index file is required).                                                                |                          |
| `ref_alt_freq`            | `string`          | the field for alternative allele frequency in INFO                                                                          |                          |
| `chr_dict`                | `dict`            | a conversion dictionary for chromosome notations in sumstats and those in reference VCF/BCF                                 | `gl.get_number_to_chr()` |
| `maf_threshold`           | `string`          | only palindromic SNPs with MAF < `maf_threshold` will be inferrred                                                          | `0.4`                    |
| `daf_tolerance`           | `string`          | only indistinguishable indels with difference in allele frequency < `daf_tolerance` will be inferrred                       | `0.2`                    |
| `remove_snp`              | ` `, `7` or `8`   | `7` remove palindromic SNPs with MAF unable to infer. `8`: remove palindromic SNPs with No infromation in reference VCF/BCF | ``                       |
| `remove_indel`            | ` ` or `8`        | `8`: indistinguishable indels with No infromation in reference VCF/BCF                                                      | ``                       |
| `mode`                    | `p`, `i`, or `pi` | `p`: infer panlindromic SNPs. `i`: infer indels.                                                                            | `pi`                     |
| `n_cores`                 | `int`             | number of CPU threads to use                                                                                                | `1`                      |
| `n_cores`                 | `int`             | number of CPU threads to use                                                                                                | `1`                      |
| `cache_options` (available since v3.4.42)           | `dict`            | options for using cache to speed up this step                                                                                                     | `None`                   |


| `cache_options` | DataType                                  | Description                                                                                                                                                                                                                                             | Default |
|-----------------|-------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `cache_manager` | CacheManager object or `None`               | If any between cache_loader and cache_process is not None, or use_cache is True, a CacheManager object will be created automatically.                                                                                                                   |         |
| `trust_cache`   | `bool`                                      | Whether to completely trust the cache or not. Trusting the cache means that any key not found inside the cache will be considered as a missing value even in the VCF file.                                                                              | True    |
| `cache_loader`  | Object with a get_cache() method or `None`  | Object with a get_cache() method or None.                                                                                                                                                                                                               |         |
| `cache_process` | Object with an apply_fn() method or `None` | Object with an apply_fn() method or None.                                                                                                                                                                                                               |         |
| `use_cache`     | `bool`                                      | If any of the cache_manager, cache_loader or cache_process is not None, this will be set to True automatically. If set to True and all between cache_manager, cache_loader and cache_process are None, the cache will be loaded (or built) on the spot. | False   |

!!! example
    ```python
    mysumstats.infer_strand()
    mysumstats.flip_allele_stats()

    # using cache to speed up this method
    mysumstats.infer_strand(cache_options={"use_cache":True})
    mysumstats.flip_allele_stats()
    ```
    
!!! note
    `infer_strand()` only change the status code. Use [filp function](#flipping-based-on-status-code) `.flip_allele_stats()` to filp the allele-specific stats.


## Check the difference in allele frequency

`.check_daf()` : check the allele frequency discrepancy with a reference vcf. Please make sure your sumstats are already harmonized, and the variants in reference VCF are also aligned. gwaslab will retrieve information only for matched variants (CHR, POS, EA-ALT, and NEA-REF).

| `.check_daf()` options | DataType | Description                                                                                 | Default                  |
|------------------------|----------|---------------------------------------------------------------------------------------------|--------------------------|
| `ref_infer`            | `string` | path to the reference VCF/BCF file (index file is required).                                |                          |
| `ref_alt_freq`         | `string` | the field for alternative allele frequency in INFO                                          |                          |
| `chr_dict`             | `dict`   | a conversion dictionary for chromosome notations in sumstats and those in reference VCF/BCF | `gl.get_number_to_chr()` |
| `n_cores`              | `int`    | number of CPU threads to use                                                                | `1`                      |

!!! example
    ```
    mysumstats.check_af(ref_infer=gl.get_path("1kg_eas_hg19"), 
                        ref_alt_freq="AF",
                        n_cores=2)
    ```

DAF : Difference between Effect allele frequency and Reference ALT frequency
EAF: Effect allele frequency
RAF: Reference ALT allele frequency

You may want to check the allele frequency discrepancy with a reference VCF. Just specify the path and the right allele frequency for your target ancestry in INFO field.

## Allele frequency correlation plot

GWASlab will simply calculate DAF = AF-EAF - AF-ALT , and store the results in DAF column. DAF can then be used for plotting (`.plot_daf()`) or filter variants.

!!! example
    ```python
    mysumstats.plot_daf(threshold=0.12)
    ```
![image](https://github.com/Cloufield/gwaslab/assets/40289485/0c607470-bbb6-4f11-93fe-038a53f6eebb)

## Flip allele-specific statistics

`.flip_allele_stats()` :  Flip allele-specific statistics to harmonize the variants based on the tracking status code in STATUS. 

!!! example
    ```python 
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    
    mysumstats.infer_strand()
    mysumstats.flip_allele_stats()
    ```


## Assign CHR and POS according to rsID and reference data

```python
mysumstats.rsid_to_chrpos()  

mysumstats.rsid_to_chrpos2()  
```

