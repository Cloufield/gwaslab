# Harmonization

GWASLab provides reference-dependent harmonization functions.

## Methods summary

| Sumstats Methods       | Options                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | Description                                                                |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| `.check_ref()`         | `ref_path`,<br /> `chr_dict=get_chr_to_number()`                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | Check alignment with a reference sequence                                  |
| `.assign_rsid()`       | `ref_rsid_tsv`,<br /> `ref_rsid_vcf`,<br /> `n_cores=1`, <br />`chunksize=5000000`, <br />`chr_dict=get_number_to_chr()`, <br />`overwrite="empty"`                                                                                                                                                                                                                                                                                                                                                                       | Annotate rsid using a reference vcf file                                   |
| `.infer_strand()`      | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`remove_snp=""`,<br />`mode="pi"`,<br />`n_cores=1`,<br />`remove_indel=""`                                                                                                                                                                                                                                                                                                                                                                        | Infer the strand of a variant using reference vcf file with EAF in INFO    |
| `.check_daf()`         | `ref_infer`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`                                                                                                                                                                                                                                                                                                                                                                                                                                        | Calculate difference in allele frequencies                                 |
| `.flip_allele_stats()` |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | After alignment and inferring, flip the alleles to harmonise the variants. |
| `.harmonize()`         | `basic_check=True`, <br /> `ref_seq=None`,<br />`ref_rsid_tsv=None`,<br />`ref_rsid_vcf=None`,<br />`ref_infer=None`,<br />`ref_alt_freq=None`,<br />`maf_threshold=0.40`,<br />`n_cores=1`,<br />`remove=False`,<br />`checkref_args={}`,<br />`removedup_args={}`,<br />`assignrsid_args={}`,<br />`inferstrand_args={}`,<br />`flipallelestats_args={}`,<br />`fixid_args={}`,<br />`fixchr_agrs={}`,<br />`fixpos_args={}`,<br />`fixallele_args={}`,<br />`sanitycheckstats_args={}`,<br />`normalizeallele_args={}` | all-in-one function for harmonization                                      |

## Align NEA with REF in the reference genome

`.check_ref()`:  Check if NEA is aligned with the reference sequence. After checking, the tracking status code will be changed accordingly. 

!!! example
    ```python
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    ```
!!! note
    `check_ref()` only change the status code. Use [flip function](#flipping-based-on-status-code) `.flip_allele_stats()` to flip the allele-specific stats.

## Assign rsID according to CHR, POS, REF/ALT

`.assign_rsid()` : Annotated variants with rsID using a reference tsv file (1KG variants) and reference vcf file (tabix indexed, entire dbSNP).

!!! example
    ```python
    mysumstats.assign_rsid(ref_rsid_tsv, ref_rsid_vcf, n_cores=1)
    ```

- For tsv file, variants will be matched using SNPID (CHR:POS:NEA:EA) for quick assigning.
- For VCF file, GWASLab will first extract all variants in the reference file with matching CHR and POS. And then compare EA/NEA in sumstats with REF/ALT in reference vcf. When matching, it will annotate the variant in sumstats with the matching rsID in reference vcf.  

## Check palindromic SNPs or indistinguishable Indels

`.infer_strand()`:

- Infer the strand for palindromic SNPs (AT, or CG), the default threshold is 0.40. 
- Checking the alignment status of indels with the REF allele in a reference vcf file.

!!! example
    ```python
    mysumstats.infer_strand()
    mysumstats.flip_allele_stats()
    ```
    
!!! note
    `infer_strand()` only change the status code. Use [filp function](#flipping-based-on-status-code) `.flip_allele_stats()` to filp the allele-specific stats.


## Check the difference in allele frequency

`.check_daf()` : check the allele frequency discrepancy with a reference vcf. Please make sure your sumstats are already harmonized, and the variants in reference VCF are also aligned. gwaslab will retrieve information only for matched variants (CHR, POS, EA-ALT, and NEA-REF).

`ref_infer`: reference VCF file path.
`ref_alt_freq`:  allele frequency for ALT in the INFO field of reference VCF file.
`n_cores`: number of cores to use.

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
```

