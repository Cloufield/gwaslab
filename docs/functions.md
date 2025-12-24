# gwaslab.Sumstats methods


## gl.Sumstats


##  Standardization

| Sumstats Methods      | Options                                                      | Description                                                         |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| `.fix_id()`           | `fixchrpos=False`, <br/>`fixid=False`, <br/>`fixsep=False`,<br/>`overwrite=False`,<br/>`forcefixid=False` | check and  fix rsID or SNPID(chr:pos:ref:alt), or use snpid to fix CHR and POS |
| `.fix_chr()`          | `remove=False`                                               | standardize chromosome notation                                                 |
| `.fix_pos()`          | `remove=False`                                               | standardize basepair position notation and filter out bad values              |
| `.fix_allele()`       | `remove=False`                                               | standardize base notation to ATCG                                              |
| `.normalize_allele()` | `threads=1`                                                  | normalize indels (only support ATA:AA -> AT:A but not -:T)                     |
| `.sort_coordinate()`  |                                                              | sort the variant coordinates                                                   |
| `.sort_column()`  |                                                              | sort the column order to GWASLab default                                                   |
| `.basic_check()`  |                                                              | all in one function                                              |

##  QC and filtering

| Sumstats Methods  | Options                  | Description                                                             |
| ----------------- | ------------------------ | ----------------------------------------------------------------------- |
| `.check_sanity()` |  `n=(0,float("Inf"))`, <br/>`eaf=(0,1)`, <br/>`mac=(5,float("Inf"))`, <br/>`chisq=(0,float("Inf"))`, <br/>`p=(5e-300,1)`, <br/>`mlog10p=(0,float("Inf"))`, <br/>`beta=(-10,10)`, <br/>`z=(-37.5,37.5)`, <br/>`se=(0,float("Inf"))`, <br/>`OR=(-10,10)` , <br/>`OR_95L=(0,float("Inf"))`, <br/>`OR_95U=(0,float("Inf"))`, <br/>`info=(0,float("Inf"))`   | sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N... |
| `.remove_dup()`   |  `mode="md"`, <br/>` keep='first'`, <br/>`keep_col="P"`, <br/>`remove_na=False` | remove duplicated, multiallelic or NA variants |
| `.filter_value()`    |  expr     |    filter in variants based on expr                                                                    |
| `.filter_in()`    |  lt, gt, eq, inplace     |    filter in variants based on given threshold                                                                      |
| `.filter_out()`   |  lt, gt, eq, inplace     |       filter out variants based on given threshold                                                                      |
| `.filter_region_in()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                         |      filter in variants in the specified region defined by a bed file                                                                   |
| `.filter_region_out()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                        |      filter out variants in the specified region defined by a bed file                                                                  |

##  Harmonization

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.check_ref()`         | ref_seq                                              | check alignment with a reference sequence                                  |
| `.rsid_to_chrpos()`    | path, threads                                         | use rsid to fill CHR and POS                                               |
| `.rsid_to_chrpos2()`   | path                                                  | use rsid to fill CHR and POS (muilti-thread, need hd5 file)                |
| `.assign_rsid()`       | path                                                  | annotate rsid using a reference vcf file                                   |
| `.infer_strand()`      | ref_infer="" , ref_alt_freq=None,  maf_threshold=0.40, daf_tolerance=0.20 | infer the strand of palindromic SNPs and indels using reference vcf file with allele frequency in INFO    |
| `.infer_strand2()`     | path/vcf_path/tsv_path, maf_threshold=0.40, ref_maf_threshold=0.40, daf_tolerance=0.20 | optimized strand inference using pre-annotated RAF from lookup table (faster for large datasets) |
| `.check_af()`         | ref_infer, ref_alt_freq=None, maf_threshold=0.40 | calculate difference in allele frequencies (DAF) between sumstats EAF and reference VCF ALT frequency |
| `.infer_af()`          | ref_infer, ref_alt_freq=None | infer effect allele frequency (EAF) in sumstats using reference VCF ALT frequency |
| `.flip_allele_stats()` |                                                       | After alignment and inferring, flip the alleles to harmonise the variants. |
| `.liftover()`          | threads=1,from_build="19", to_build="38"              | perform liftover                                                           |
| `.infer_build()`       | verbose=True                                           | infer genome build (hg19/hg38) from CHR:POS coordinates                     |
| `.set_build()`         | build, verbose=True                                    | manually set the genome build                                              |

##  Data Conversion

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.fill_data()`         | to_fill, extreme=False                                 | convert equivalent statistics (P, MLOG10P, Z, CHISQ, BETA, SE, OR, etc.)  |

##  Visualization

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.plot_mqq()`          | mode="m", ...                                          | create Manhattan and Q-Q plots                                             |
| `.plot_mqq(mode="r")`  | ...                                                    | create regional plot                                                       |
| `.plot_mqq(mode="b")`  | ...                                                    | create Brisbane plot (signal density plot)                                 |
| `.plot_trumpet()`      | ...                                                    | create Trumpet plot                                                         |
| `.plot_effect()`       | ...                                                    | create effect size comparison plot                                         |
| `.plot_daf()`          | ref_path, ...                                          | plot difference in allele frequencies                                      |

##  Format and Output

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.to_format()`         | path, fmt, tab_fmt, cols, extract, exclude, ...       | output sumstats in specified format (plink, ldsc, saige, etc.)             |
| `.to_pickle()`         | path="~/mysumstats.pickle", overwrite=False           | save Sumstats object to pickle file                                        |

##  Utilities

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.get_lead()`          | scaled=False, use_p=False, windowsizekb=500, ...      | extract lead variants from significant loci                                |
| `.get_novel()`         | ref_path, ...                                         | extract novel variants not present in reference                            |
| `.summary()`           |                                                        | print summary statistics of the sumstats                                   |
| `.lookup_status()`     | status="STATUS"                                        | lookup and display STATUS code information                                 |
| `.copy()`              |                                                        | create a copy of the Sumstats object                                       |

##  Clumping

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.clump()`             | plink2, vcf, bfile, pfile, scaled, threads, ...       | perform clumping using PLINK2                                              |

##  LD Score Regression

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.estimate_h2_by_ldsc()` | ldscore_path, maf_path, n, ...                        | estimate heritability using LD score regression                            |
| `.estimate_rg_by_ldsc()` | sumstats2, ldscore_path, maf_path, ...                | estimate genetic correlation using LD score regression                     |
| `.load_ldsc_log()`     | log_path                                               | load LDSC log file and parse results                                       |




# Standalone functions

| Function               | Category      | Description                                                                |
| ---------------------- | ------------- | -------------------------------------------------------------------------- |
| `gl.plot_miami2()`     | Visualization | create Miami plot from two sumstats files                                 |
| `gl.h2_obs_to_liab()`  | Conversion    | convert heritability from observed scale to liability scale               |
