# Update Logs

## 3.4.9
- update variant matching criteria for regional plot

## 3.4.8 - 2023/04/21
- fixed regional plot lead variant line error
- added tutorial for v3.4
- fixed anno_alias priority

## 3.4.7 - 2023/04/20
- fixed font_family for annotation
- get_lead can now use mlo10p to extract lead variants 
- added sig_level_lead in plot_mqq

## 3.4.6 - 2023/04/18
- updated yticklabel fontsize
- updated drop_chr_start

## 3.4.5 - 2023/04/06
- added rr_lim : input a tuple like (0,100) or "max"

## 3.4.4 - 2023/04/04
- fixed error in mqqplot (chr > 26)

## 3.4.3 - 2023/03/29
- added jagged y axis
- added cut_log

## 3.4.2 - 2023/03/28
- fixed y tick labels
- added font_family
- added ylabels

## 3.4.1 - 2023/03/28
- added sc_linewidth
- fixed use_rank
- fixed ystep

## 3.4.0
- restructured plot functions
- added dtype conversion for input pd.DataFrame
- reimplemented gtfparse and revised requirements
- update reference datasets
- fixed bugs when no variants were selected for mqqplot
- fixed bugs in `gl.check_downloaded_ref()`
- added suggestive significance line


## 3.3.24 - 2023/02/01
- update bugs in __init__.py

## 3.3.23 - 2023/01/31
- update annotation arrow style `anno_style` for `.plot_mqq`
- fixed effect size comparison bugs (added sorting)
- added allele check for effect size comparison
- update large number selection algorithm
- fixed dtype errors in `fix_chr` and `fix_pos`
- added perSNPh2 `.get_per_snp_r2()` and F statistics
- implemented save in Miami plots
- fixed annotation error in Miami plots

## 3.3.22 - 2023/01/27
- added checks for duplicates and NAs in compare_effect()

## 3.3.21 - 2023/01/27

- implement package info gl.show_version()
- fixed rsid_to_chrpos()

## 3.3.20 

- updated get_density() to calculate the signal density for sumstats
- implemented winner's curse correction for effect size comparison


## 3.3.19

- updated fontsize options for plot_mqq(). Added anno_fontsize, title_fontsize.
- updated default values and optimized methods for remove_dup(), fix_chr(), basic_check(),check_sanity().
- added dump_pickle() load_pickle() to save half-finished sumstats object.


## 3.3.18

- updated config and downloading system


## 3.3.17

- added xtcik_chr_dict
- fixed bugs in miami plot


## 3.3.16

- added hg38 recombination rate file
- fixed bugs in get_novel

## 3.3.15

- fixed bugs for reading gtf files

## 3.3.14

- updated requirements for dependencies 
    - pandas>=1.3,<1.5
    - pyensembl==2.2.3
- support customized gtf/vcf/recombination_rate files
- fixed bugs for regional plot
- calculate_gc()
- fill MAF : `.fill_data(to_fill=["MAF"])`

## 3.3.13

- specified python engine for query 

## 3.3.12

- fixed bugs for matplotlib v3.6.x
- added method chain for filter_xxx functions
- updated requirements for dependencies 
    - pySAM>=0.18.1,<0.20
    - matplotlib>=3.5
    - pyensembl>=2.2.3

## 3.3.11

- updated requirements for dependencies

## 3.3.10 

- updated bugs for mqqplot

## v3.3.9 

- updated download system

## v3.3.8

- included recombination data

## v3.3.7 

- updated packaging methods. Now when installing gwaslab, pip will install all dependencies as well.

## v3.3.6 - 2022/11/05

- added download function: 
    - now you can download reference files from predefined list via gwaslab
    - `gl.check_available_ref()` : list available reference files
    - `gl.check_downloaded_ref()` : list downloaded reference files
    - `gl.download_ref(name)` : download reference files
    - `gl.remove_file(name)` : remove the local reference files
    - `gl.get_path(name)` : get the local path for the reference data name
- implemented parsing gwas-vcf (`fmt="vcf"`)
- implemented `Sumstats.filter_value(expr)`
- fixed bugs for check_allele
- optimized functions for sorting columns
- removed outdated codes in Sumstats

## v3.3.5 - 2022/11/02 

- added `filter_value`
- integrate `gwascatalog` to `get_novel`
- optimized `remove_dup`
- fixed bugs

## v3.3.4 - 2022/10/31

- added gwascatalog_trait()
- optimized check_sanity()
- optimized the logic for removing duplicated and multiallelic variants 
- added update_formatbook()
- added functions to read vcf.gz `gl.Sumstats("myvcf.vcf.gz",fmt="vcf")`
- gwaslab is now able to read chromosome-separated files
- fixed bugs

## v3.3.3

- added Miami plot
- added Brisbane plot
- updated tutorials 
- mqq plot annotation: new customization options

## v3.3.2 - 2022/10/21

- added forcefixid for fix_id()
- fixed bugs for plotting gene tracks

## v3.3.1 

- extract novel loci given a list of known lead variants
- fixed bugs in fill_data()
- fixed path for hapmap3 snps for infer_build()

## v3.3.0 - 2022/10/18

- added forest plot
- fixed options for mqqplot
- supported vcf

## v3.2.0

- incorporated pyensembl and scikit-allel.
- get_lead() : support automatic gene name annotation (using pyensembl)
- to_format():
	- support common sumstats formats 
	- support 1-based bed-like formats for VEP
	- support 0-based bed-like formats
- manhattan plot:
	- optimized plotting logic
	- annotate gene names
	- added regional plot feature using a user-provided reference panel
- comparison effect plot: 
	- fix using OR


## v3.1.0 

- implemented formatbook: easily import sumstats and output sumstats in certain formats (support for commonly used formats including ldsc, plink, plink2, gwas-ssf, saige, regenie, fastgwa, metal, mrmega, pgscatalog, pgscatalog_hm, gwascatalog, gwascatalog_hm and gwaslab)
- added `.filter_region_in/out` using bed files (or in-built regions like high-ld or hla)
- implemented `.summay()` methods.
- optimized rsID annotation pipeline. Support annotation using curated chr:pos:ref:alt - rsID tsv for quick annotation.
- changed some datatypes and optimized memory usage.
- replaced pyVCF with pySAM
