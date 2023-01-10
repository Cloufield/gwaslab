# Harmonization

gwaslab provides reference-dependent harmonization functions.

See examples [here.](https://cloufield.github.io/gwaslab/harmonization_workflow/)

## Methods summary

| Sumstats Methods| Options| Description |
|-|-|-|
| `.check_ref()` | `ref_path`, `chr_dict=get_chr_to_number()`     | Check alignment with a reference sequence                                  |
| `.assign_rsid()` | `ref_rsid_tsv`, `ref_rsid_vcf`, `n_cores=1`, `chunksize=5000000`, `chr_dict=get_number_to_chr()`, `overwrite="empty"` | Annotate rsid using a reference vcf file                                   |
| `.infer_strand()`      | `ref_infer`,`ref_alt_freq=None`,`maf_threshold=0.40`,`remove_snp=""`,`mode="pi"`,`n_cores=1`,`remove_indel=""` | Infer the strand of a variant using reference vcf file with EAF in INFO    |
| `.check_daf()`         | `ref_infer`,`ref_alt_freq=None`,`maf_threshold=0.40`,`n_cores=1`                   | Calculate difference in allele frequencies                                 |
| `.flip_allele_stats()` |                                                       | After alignment and inferring, flip the alleles to harmonise the variants. |
| `.liftover()`          | `n_cores=1`,`from_build="19"`, `to_build="38"`              | Perform liftover for POS                                                          |

## Align NEA with REF in reference genome

`.check_ref()`
Check if NEA is aligned with the reference sequence. After checking, the tracking status code will be changed accordingly. 

!!! example
    ```python
    mysumstats.check_ref(ref_path="ref_genome.fa")
    mysumstats.flip_allele_stats()
    ```
!!! note
    `check_ref()` only change the status code. Use [filp function](#flipping-based-on-status-code) `.flip_allele_stats()` to filp the allele-specific stats.

## Assign rsID according to CHR, POS, REF/ALT

`.assign_rsid()`

!!! example
    ```python
    mysumstats.assign_rsid(ref_rsid_tsv, ref_rsid_vcf, n_cores=1)
    ```

Annotated variants with rsID using a reference tsv file (1KG variants) and reference vcf file (tabix indexd, entire dbSNP).

- For tsv file, variants will be matched using SNPID (CHR:POS:NEA:EA) for quick assigning.
- For VCF file, Gwaslab will first extract all variants in reference file with matching CHR and POS. And then comapre EA/NEA in sumstats with REF/ALT in reference vcf. When matching, it will annotate the vairant in sumstats with the matching rsID in reference vcf.  

## Check panlidromic SNPs or undistingushable Indels

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

`.check_daf()`

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

`.flip_allele_stats()`

Flip allele-specific statistics to harmonise the variants based on the tracking status code

!!! example
    ```python
    mysumstats.flip_allele_stats()
    ```

## Liftover

Perform liftover for POS (based on liftover [GitHub - jeremymcrae/liftover: liftover for python, made fast with cython](https://github.com/jeremymcrae/liftover))

!!!example
    ```python
    mysumstats.liftover(n_cores=1,from_build="19", to_build="38")
    ```
    
!!! note
    gwaslab will only liftover basepair positions.  

## Assign CHR and POS according to rsID and reference data

```python
mysumstats.rsid_to_chrpos()     #single thread
mysumstats.rsid_to_chrpos2()    #multithread
```
(to be tested)
