# gwaslab.Sumstats methods


## gl.Sumstats


##  Standardization

| Sumstats Methods      | Options                                                      | Description                                                         |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------ |
| `.fix_id()`           | `fixchrpos=False`, <br/>`fixid=False`, <br/>`fixsep=False`,<br/>`overwrite=False`,<br/>`forcefixid=False` | check and  fix rsID or SNPID(chr:pos:ref:alt), or use snpid to fix CHR and POS |
| `.fix_CHR()`          | `remove=False`                                               | standardize chromosome notation                                                 |
| `.fix_POS()`          | `remove=False`                                               | standardize basepair position notation and filter out bad values              |
| `.fix_allele()`       | `remove=False`                                               | standardize base notation to ATCG                                              |
| `.normalize_allele()` | `threads=1`                                                  | normalize indels (only support ATA:AA -> AT:A but not -:T)                     |
| `.sort_coordinate()`  |                                                              | sort the variant coordinates                                                   |
| `.sort_column()`  |                                                              | sort the column order to GWASLab default                                                   |
| `.basic_check()`  |                                                              | all in one function                                              |

##  QC and filtering

| Sumstats Methods  | Options                  | Description                                                             |
| ----------------- | ------------------------ | ----------------------------------------------------------------------- |
| `.check_sanity()` |  `n=(0,float("Inf"))`, <br/>`eaf=(0,1)`, <br/>`mac=(5,float("Inf"))`, <br/>`chisq=(0,float("Inf"))`, <br/>`p=(5e-300,1)`, <br/>`mlog10p=(0,float("Inf"))`, <br/>`beta=(-10,10)`, <br/>`z=(-37.5,37.5)`, <br/>`se=(0,float("Inf"))`, <br/>`OR=(-10,10)` , <br/>`OR_95L=(0,float("Inf"))`, <br/>`OR_95U=(0,float("Inf"))`, <br/>`info=(0,float("Inf"))`   | sanity check for statistics including BETA, SE, Z, CHISQ, EAF, OR, N... |
| `.remove_dup()`   |  `mode="md"`, <br/>` keep='first'`, <br/>`keep_col="P"`, <br/>`remove=False` | remove duplicated, multiallelic or NA variants |
| `.filter_value()`    |  expr     |    filter in variants based on expr                                                                    |
| `.filter_in()`    |  lt, gt, eq, inplace     |    filter in variants based on given threshold                                                                      |
| `.filter_out()`   |  lt, gt, eq, inplace     |       filter out variants based on given threshold                                                                      |
| `.filter_region_in()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                         |      filter in variants in the specified region defined by a bed file                                                                   |
| `.filter_region_out()`   | `path` , <br/> `inplace=True` , <br/>`high_ld=False`, <br/> `build="19"`                        |      filter out variants in the specified region defined by a bed file                                                                  |

##  Harmonization

| Sumstats Methods       | Options                                               | Description                                                                |
| ---------------------- | ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `.check_ref()`         | ref_path                                              | check alignment with a reference sequence                                  |
| `.rsid_to_chrpos()`    | path, threads                                         | use rsid to fill CHR and POS                                               |
| `.rsid_to_chrpos2()`   | path                                                  | use rsid to fill CHR and POS (muilti-thread, need hd5 file)                |
| `.assign_rsid()`       | path                                                  | annotate rsid using a reference vcf file                                   |
| `.infer_strand()`      | ref_infer="" , ref_alt_freq=None,  maf_threshold=0.43 | infer the strand of a variant using reference vcf file with EAF in INFO    |
| `.check_daf()`         | ref_infer="" , ref_alt_freq=None,                     | calculate difference in allele frequencies                                 |
| `.flip_allele_stats()` |                                                       | After alignment and inferring, flip the alleles to harmonise the variants. |
| `.liftover()`          | threads=1,from_build="19", to_build="38"              | perform liftover                                                           |




# Standalone functions
| function | catagory | description |
|-|-|-|
||||
||||
