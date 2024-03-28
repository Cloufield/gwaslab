# Update Logs

## 3.4.42 20240328
- fast implementation of check_ref and to_format. (credit to @sup3rgiu Mr. Andrea)
- added highligh and pinpoint for plot_trumpet()
- fixed typo (credit to @sup3rgiu Mr. Andrea)
- changed `**args` to `**kwargs` (credit to @sup3rgiu Mr. Andrea)
- implemented munge-like filters `munge=True`  for ldsc in gwaslab

## 3.4.41 20240219
- fixed error in region_ref_second

## 3.4.40 20240215
- fixed color issue for regional plots. The color assigned to each variants is actually the color for the lower LD r2 category. For example, variants with LD>0.8 will be colored with the color for 0.8>LD>0.6.
- integrated LDSC (partitioned h2/h2-cts)
- added wc_correction for get_lead()
- updated log

## 3.4.39 20240207
- integrated LDSC h2/rg
- updated LICENSE from MIT to GPL-3.0 license
- added check_novel_set / check_cis 
- added filter_snp/palindromic/indel and filter_hapmap3
- updated log system
- fixed bug when gene name is empty in GTF in regional plot
- fixed error in harmonization when there is no palindromic SNPS or indels to check.
- fixed bug when saving plot as pdf with matplotlib>3.6
- fixed typos

## 3.4.38 20240203
- fixed seaborn version
- replaced array_split
- fixed warnings due to pandas upgrade 
- renamed get_flanking to filter_flanking
- fixed error in cbar for regional plot
- restructured log system
- restructured sanity checking
- restructured stats flipping
- updated tutorial and examples
- fixed error in log

## 3.4.37 20240129
- added data consistency check
- added "s","r","n" mode for remove_dup()
- updated fix_id
- added sample data
- added datatype check for rsID and SNPID
- added memory usage check 
- added extreme P value check

## 3.4.36 20240126
- removed statsmodels
- updated functions in basic_check()
- added a test example for basic_check()
- fixed fontsize and font_family errors in plot_mqq() for ax4
- added new options for plot_mqq(mode="r"); cbar_fontsize, cbar_font_family, cbar_title
- added new options for plot_mqq(mode="r"); track_n, track_n_offset, track_fontsize_ratio, track_exon_ratio, track_text_offset, track_font_family 

## 3.4.35 20240123
- fixed version number

## 3.4.34 20240123
- get_lead will use MLOG10P first instead of P (`scaled` is deprecated in `get_lead`)
- added datatype check for fill_data
- added datatype verification
- updated sanity check default values (BETA, OR, HR)

## 3.4.33 20240122
- fix_allele() now prints out allele information 
- updated sanity check default values (added tolerance for floats)
- fixed error in loading pickle created by older versions
- fixed bug in flipping OR and HR
- added version information to the start line of each checking function
- updated reference VCF
- updated clump() / to_finemapping() / run_susie_rss() (beta)
- updated Manhattan-like plotting system 
- updated datatype for certain statistics
- updated file naming system (beta)
- added stacked mqq plot (beta)
- plot_miami2() can iteratively call plot_mqq to create miami plot, which supports more functions (beta)

## 3.4.32 20231118
- fixed bug in sanity_check : N = N_CASE + N_CONTROL
- fixed the logic chain for normalization
- added extra log in compare_effect
- fixed a few typos (credit to @gmauro Mr. Gianmauro Cuccuru)
- added clump() for beta testing
- added to_finemapping() for beta testing
- added run_susie_rss() for beta testing

## 3.4.31 20231115
- updated version requirements

## 3.4.30 20231103
- updated random seed range (0,2^32-1).
- fixed error for suggestive_sig_line in miami plot.
- loosened version requirements for python, pandas and matplotlib.  

## 3.4.29 20231003
- set `scatter_kwargs={"rasterized":True}` as default for miami plot

## 3.4.28 20231003
- fixed a bug in `compare_effect`

## 3.4.27 20231002
- update reference book
- added extra parameters for tabix_index()
- update random variants

## 3.4.26 20230911
- added `infer_af()`
- updated reference datasets

## 3.4.24 20230821
- supported `highlight` and `pinpoint` for multiple sets of loci and variants in `.plot_mqq()`
- added `overwrite` option for `gl.download_ref()`

## 3.4.24 20230821
- update references

## 3.4.23 20230817
- fixed bug in `get_lead()`. In some rare cases, it was not counted when the last variant is a new lead variant.

## 3.4.22 20230803
- added plot_power() and plot_power_x()
- added support for multiple EFO IDs
- updated built-in formatbook
- fixed regex error in read_ldsc for numbers in the format of 1e-01
- fixed typos

## 3.4.21 2023/07/18
- fix error in gc calculation for mqq plot when `expected_min_mlog10p` is not 0 and `stratified=True`.

## 3.4.20 2023/07/15

- reimplement and unified the module for saving figures
- added trumpet plot document page   
- added CHR range check in `.fix_chr()` and `.plot_mqq()` (remove variants with CHR<=0) (#42)
- fixed error in column headers for plot_mqq "b" mode (#40)
- fixed many typos.. (#41)
- fixed error in plot_mqq logging (warning for genome build) 
- fixed error in annotation for miami plot 

## 3.4.19 2023/06/29
- added annotation using custom column for miami plot
- updated alogorithm for extracting lead variants (using scaled P)

## 3.4.18 2023/06/28
- fixed bugs for miami plot

## 3.4.17 2023/06/27
- added `expected_min_mlog10p` for qq mode in `plot_mqq()`
- added trumpet plot
- added support for gwaslab Sumstats object for gl.compare_effect()
- added saving options for gl.compare_effect()
- fixed bug for `is_q=False` in gl.compare_effect().

## 3.4.16 2023/06/22
- fixed bug for qq mode in `plot_mqq()`

## 3.4.15 2023/06/20
- [LDSC-rg genetic correlation heatmap](https://cloufield.github.io/gwaslab/GeneticCorrelation/)
- [Allele frequency correlation plot](https://cloufield.github.io/gwaslab/AlleleFrequency/)
- Miami plot using gl.Sumstats Object pickle files
- Auto-check for VCF chromosome prefix (chr1 or 1)
- fill_data() is now implemented iteratively
- Downloaded files auto detection.

## 3.4.14 2023/06/09
- fixed bug

## 3.4.13 2023/06/09
- added two-reference-variant mode for regional plot
- auto-check for genome build version
- added `additional_line` and `additional_line_color` for plot_mqq
- fixed bugs for sig_level_lead

## 3.4.12 2023/05/30
- added highlight_anno_args in plot_mqq()
- added GWAS-SSF style metadata yaml output support
- added region_ref for regional plots
- updated default formatbook

## 3.4.11 - 2023/04/28
- fixed anno=None in miami plot

## 3.4.10 - 2023/04/27
- fixed variant_id sep for ssf
- update tutorials

## 3.4.9 - 2023/04/25
- update variant matching criteria for regional plot
- fixed default values for mqqplot
- added qq_scatter_args
- added connection timeout, status code check and md5sum verification for download_ref

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
- added jagged y-axis
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
