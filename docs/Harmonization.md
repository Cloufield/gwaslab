# Harmonization

gwaslab provides reference-dependent harmonization functions.

See examples [here.](https://cloufield.github.io/gwaslab/harmonization_workflow/)

## Methods summary
Outdated!!!
| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.check_ref()`         | ref_path                                              | check alignment with a reference sequence                                  |
| `.rsid_to_chrpos()`    | path, n_cores                                         | use rsid to fill CHR and POS                                               |
| `.rsid_to_chrpos2()`   | path                                                  | use rsid to fill CHR and POS (muilti-thread, need hd5 file)                |
| `.assign_rsid()`       | path                                                  | annotate rsid using a reference vcf file                                   |
| `.infer_strand()`      | ref_infer="" , ref_alt_freq=None,  maf_threshold=0.43 | infer the strand of a variant using reference vcf file with EAF in INFO    |
| `.check_daf()`         | ref_infer="" , ref_alt_freq=None,                     | calculate difference in allele frequencies                                 |
| `.flip_allele_stats()` |                                                       | After alignment and inferring, flip the alleles to harmonise the variants. |
| `.liftover()`          | n_cores=1,from_build="19", to_build="38"              | perform liftover                                                           |
|                        |                                                       |                                                                            |

## Align NEA with REF in reference genome

```python
mysumstats.check_ref(ref_path="ref_genome.fa")
```

(! Only changing the status code) 

Check if NEA is aligned with the reference sequence. After checking, the tracking status code will be changed accordingly.   

## Assign CHR and POS according to rsID and reference data

```python
mysumstats.rsid_to_chrpos()     #single thread
mysumstats.rsid_to_chrpos2()    #multithread
```

(to be tested)

## Assign rsID according to CHR, POS, REF/ALT

```python
mysumstats.assign_rsid(path="reference.vcf.gz")
```

Annotated variants with rsID using a reference vcf file (tabix indexd).

Gwaslab will first extract all variants in reference file with matching CHR and POS. And then comapre EA/NEA in sumstats with REF/ALT in reference vcf. When matching, it will annotate the vairant in sumstats with the matching rsID in reference vcf.  

## Check panlidromic SNPs or undistingushable Indels

```python
mysumstats.infer_strand()
```

(! Only changing the status code)

1. Infer the strand for palindromic SNPs (AT, or CG), the default threshlod is 0.43. 
   
   1. make sure specify the right allele frequency for you target ancestry in INFO field.

2. Checking the alignment status of  indels with the REF allele in reference vcf file.

## Check difference in allele frequency

```
mysumstats.check_daf()
```

You may want to check the allele frequency discrepancy with a reference vcf. Just specify the path and the right allele frequency for you target ancestry in INFO field.

GWASlab will simply calculate DAF = AF-EAF - AF-ALT , and store the results in DAF column. DAF can then be used to plot or filter variants.

```python
mysumstats.plot_daf()
```

## Flipping based on status code

```python
mysumstats.flip_allele_stats()
```

Flip allele-specific statistics to harmonise the variants based on the tracking status code

## Liftover

```python
mysumstats.liftover(n_cores=1,from_build="19", to_build="38")
```

Perform liftover for POS (based on liftover [GitHub - jeremymcrae/liftover: liftover for python, made fast with cython](https://github.com/jeremymcrae/liftover))
