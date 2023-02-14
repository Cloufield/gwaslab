# Harmonization

GWASLab provides reference-dependent harmonization functions.

See examples [here.](https://cloufield.github.io/gwaslab/harmonization_workflow/)

## Methods summary

| Sumstats Methods| Options| Description |
|-|-|-|
| `.check_ref()` | `ref_path`,<br /> `chr_dict=get_chr_to_number()`     | Check alignment with a reference sequence                                  |
| `.assign_rsid()` | `ref_rsid_tsv`,<br /> `ref_rsid_vcf`,<br /> `n_cores=1`, <br />`chunksize=5000000`, <br />`chr_dict=get_number_to_chr()`, <br />`overwrite="empty"` | Annotate rsid using a reference vcf file                                   |
| `.infer_strand()`      | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`remove_snp=""`,<br />`mode="pi"`,<br />`n_cores=1`,<br />`remove_indel=""` | Infer the strand of a variant using reference vcf file with EAF in INFO    |
| `.check_daf()`         | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`                   | Calculate difference in allele frequencies                                 |
| `.flip_allele_stats()` |                                                       | After alignment and inferring, flip the alleles to harmonise the variants. |
|`.harmonize()`|`basic_check=True`, <br /> `ref_seq=None`,<br />`ref_rsid_tsv=None`,<br />`ref_rsid_vcf=None`,<br />`ref_infer=None`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`,<br />`remove=False`,<br />`checkref_args={}`,<br />`removedup_args={}`,<br />`assignrsid_args={}`,<br />`inferstrand_args={}`,<br />`flipallelestats_args={}`,<br />`fixid_args={}`,<br />`fixchr_agrs={}`,<br />`fixpos_args={}`,<br />`fixallele_args={}`,<br />`sanitycheckstats_args={}`,<br />`normalizeallele_args={}` |all-in-one function for harmonization|

## Align NEA with REF in reference genome

`.check_ref()`:  Check if NEA is aligned with the reference sequence. After checking, the tracking status code will be changed accordingly. 

!!! example
    ```python
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    ```
!!! note
    `check_ref()` only change the status code. Use [filp function](#flipping-based-on-status-code) `.flip_allele_stats()` to filp the allele-specific stats.

## Assign rsID according to CHR, POS, REF/ALT

`.assign_rsid()` : Annotated variants with rsID using a reference tsv file (1KG variants) and reference vcf file (tabix indexd, entire dbSNP).

!!! example
    ```python
    mysumstats.assign_rsid(ref_rsid_tsv, ref_rsid_vcf, n_cores=1)
    ```

- For tsv file, variants will be matched using SNPID (CHR:POS:NEA:EA) for quick assigning.
- For VCF file, Gwaslab will first extract all variants in reference file with matching CHR and POS. And then comapre EA/NEA in sumstats with REF/ALT in reference vcf. When matching, it will annotate the vairant in sumstats with the matching rsID in reference vcf.  

## Check panlidromic SNPs or undistingushable Indels

`.infer_strand()`:

- Infer the strand for palindromic SNPs (AT, or CG), the default threshlod is 0.40. 
- Checking the alignment status of indels with the REF allele in a reference vcf file.

!!! example
    ```python
    mysumstats.infer_strand()
    mysumstats.flip_allele_stats()
    ```
    
!!! note
    `infer_strand()` only change the status code. Use [filp function](#flipping-based-on-status-code) `.flip_allele_stats()` to filp the allele-specific stats.


## Check difference in allele frequency

`.check_daf()` : check the allele frequency discrepancy with a reference vcf.

!!! example
    ```
    mysumstats.check_daf()
    ```

You may want to check the allele frequency discrepancy with a reference vcf. Just specify the path and the right allele frequency for you target ancestry in INFO field.

GWASlab will simply calculate DAF = AF-EAF - AF-ALT , and store the results in DAF column. DAF can then be used for plotting (`.plot_daf()`) or filter variants.

!!! example
    ```python
    mysumstats.plot_daf()
    ```

## Flipping based on status code

`.flip_allele_stats()` :  Flip allele-specific statistics to harmonise the variants based on the tracking status code. 

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
```

