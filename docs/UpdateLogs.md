# Update Logs

# v4.0.8 20260204

- Fixed `plot_miami2()` build version not being passed to underlying `_mqqplot()` function: now correctly displays genome build (e.g., GRCh38/hg38) instead of defaulting to hg19
- Added `build` parameter to `plot_miami2()`: auto-detects from Sumstats objects if not provided
- Added build consistency check in `plot_miami2()`: raises `ValueError` when sumstats1 and sumstats2 have different genome builds
- Fixed `AttributeError` when `anno_alias=None` in annotation functions: added defensive null check
- Added annotation parameters to `plot_miami2()`: `anno_alias`, `anno_d`, `anno_kwargs`, `anno_style`, `anno_fontsize`, `anno_adjust`
- Fixed handling of chrX in `read_gtf()`

# v4.0.7 20260128

- Fixed `get_novel()` TypeError (wrong argument name for `_get_sig()`)
- `get_novel()`: when `build="19"`, liftover to hg38 first; queries and results use hg38
- `get_novel()`: added `show_child_traits` (default True); `use_cache` and `cache_dir` now work when using `efo`
- `get_novel()`: `efo` accepts EFO IDs, MONDO IDs, or trait names; lists may mix them (e.g. `['coffee consumption','MONDO_0004247','EFO_0004330']`)
- GWAS Catalog: MONDO used as `efo_id` (no lookup); v2 JSON cache; retry with minimal params on 500; legacy fallback honors `show_child_traits` and cache options
- Docs: ExtractNovel (Caching), utility_get_lead_novel, GWASCatalogAPI
- Tests: `get_novel()` with `known=` (DataFrame/path), no live API

# v4.0.6 20260124

- Fixed TypeError in regional plot gene track: filtered out genes with None names from `uniq_gene_region` in `process_gtf()` function to prevent concatenation errors when creating gene annotations
- Added logging for gene filtering: reports how many genes were filtered out due to missing names
- Fixed regional plot lead variant hiding: properly excludes lead variant from main Manhattan scatter plot while displaying it as pinpoint marker

# v4.0.5 20260123

- Fixed regex flag incompatibility: resolved `ValueError: ASCII and UNICODE flags are incompatible` error in `qc_pattern.py` by removing `re.ASCII` from `FLAGS` constant
- Updated pandas version constraint: set pandas version to `<=2.3.3` in `pyproject.toml`
- Added `io_bedpe` module for BEDPE file format support
- Enhanced meta-analysis functionality with improved parameter handling

# v4.0.4 20260104

- Added new visualization function: `plot_ld_block()` for plotting LD (linkage disequilibrium) matrices as 45°-rotated inverted triangles, supporting both standalone mode and integration with regional plots
- Added `ld_block=True` option to `plot_mqq(mode="r")` for regional plots: enables LD block visualization below the regional plot with automatic alignment and position bar integration
- Added new visualization feature: `ld_link=True` option to `plot_mqq(mode="r")` for regional plots: draws lines connecting variant pairs with high LD (r²), with colors matching LD categories and optional significance filtering
- LD link visualization features: line colors match variant LD categories using `region_ld_threshold` and `region_ld_colors`, transparency scales with LD value, and supports filtering by significance threshold via `ld_link_sig_level`
- Updated log messages in Manhattan plot annotation: changed "significant variants" to "lead variants" and added scaled threshold information to improve clarity

# v4.0.3 20260104

- Added new visualization: `plot_phenogram()` for chromosome-wide association visualization with cytoband annotations
- Enhanced BED file support: added `io_ucsc_bed.py` and `io_plink.py` for comprehensive BED file reading and processing
- Added ChromosomeMapper integration: unified chromosome format conversion system in `bd_chromosome_mapper.py` and `bd_common_data.py`
- Enhanced Polars support: major refactoring of Polars-based QC functions (`qc_check_datatype_polars.py`, `qc_fix_sumstats_polars.py`)
- Improved `fill_data()`: enhanced allele frequency inference and data filling capabilities with expanded test coverage
- Enhanced `filter_value()`: improved filtering logic with better expression parsing and performance optimizations
- Added cytoband data files: included hg19 and hg38 cytoband annotation files for visualization
- Enhanced regional plot: improved stacked regional plot functionality with better visualization options
- Improved GTF parsing: enhanced gene annotation and region extraction capabilities
- Enhanced QC decorators: improved error handling and logging in QC functions
- Comprehensive test suite: added extensive tests for `fill_data()`, BED filtering, and `fix_id()` functionality
- Updated CLI: improved command-line interface functionality
- Documentation updates: enhanced visualization and index documentation

# v4.0.2 

- Fixed log message accuracy issues: < or <=
- Fixed log message accuracy issues: corrected typos ("macthes" → "matches", "will will" → "will") and improved grammar in log messages across harmonization and QC modules
- Fixed bugs in `_parallelize_check_af()` and `_parallelize_infer_af()`: added missing `else` clause for `force=True` case to prevent NameError
- Fixed bug in `_fix_ID()`: corrected redundant column check `(nea in sumstats.columns) and (nea in sumstats.columns)` to properly check both `nea` and `ea` columns
- Fixed log message placement in `_infer_strand_with_annotation()`: moved renaming log message inside conditional block to only log when renaming actually occurs
- Added type hints to improve code documentation and IDE support
- Added missing dependencies to pyproject.toml: requests, pyyaml, bitarray, platformdirs

# v4.0.1 

- Added ChromosomeMapper with numeric middle layer architecture: unified chromosome format conversion system supporting numeric, string, chr-prefixed, and NCBI RefSeq formats across 12+ species
- ChromosomeMapper supports automatic format detection from sumstats and reference files (VCF, FASTA, GTF, Chain files)
- ChromosomeMapper can auto-detect genome build from NCBI RefSeq notation in reference files (hg19/hg38)
- Removed legacy chr_dict in favor of ChromosomeMapper
- Fixed infer_strand bug: added FREQ_COMPARISON_EPSILON constant (1e-6) for floating-point precision in allele frequency comparisons to handle edge cases at threshold boundaries
- Added comprehensive test cases for floating-point precision edge cases in strand inference (e.g., EAF=0.4000003453, VCF AF=0.4000000059604645)
- Enhanced P-value consistency checking in `check_data_consistency()`: now uses fold change as primary metric for better interpretability with small P-values (reports "10.23x" instead of absolute difference)
- Improved handling of P-values spanning many orders of magnitude and NaN values
- Added comprehensive test suite for P-value consistency checks
- Created ChromosomeHandling.md documentation covering architecture, formats, species support, and usage examples

# v4.0.0 20251207

## Major Changes

### Performance Improvements
- Reimplemented core functions with optimized algorithms and data structures
- Refactored GWAS Catalog API v2 integration with improved CHR/POS extraction and error handling
- Optimized status code processing using integer arithmetic
- Improved memory usage and processing speed for large datasets

### API Consistency
- Standardized parameter names across functions (e.g., unified `inplace` parameter usage)
- Consistent naming conventions throughout the codebase

### Code Quality & Testing
- Improved code organization and separation of concerns
- Comprehensive test suite with expanded coverage across core functionality

### Visualization Parameter Management System (VizParamsManager)
- New centralized parameter management system for all plotting functions
- Registry-based configuration with inheritance support (`viz_aux_params.txt`, `viz_aux_params_registry.txt`)
- Advanced parameter filtering with whitelist, numeric suffix support, and nested dict filtering
- Context-specific defaults and wildcard expansion support

### Documentation Updates
- Enhanced documentation for GWAS Catalog API, pickle operations, status codes, QC/filtering, standardization, and statistics sanity check
- Added practical examples and improved organization throughout


# v3.6.12 - v3.6.16
- temp bug fix

# 3.6.12 20251118
- added infer_strand2 and assign_rsid2 for testing
- revised some docstrings (ongoing)

# 3.6.11 20251113
- added ref_maf_threshold for inferring strand
- revised some docstrings (ongoing)

# 3.6.10 20251101
- fixed taf update_args
- fixed typos in log
- fixed legend in trumpet plot
- fixed quotation mark error for py39 in viz_plot_associations.py

# 3.6.9 20250902
- fixed wrong version numbers

# 3.6.8 20250902
- added `qq_xlabels` and `qq_xlim` for `.plot_mqq(m="qq")`
- added `highlight_lim` and `highlight_lim_mode` to accept different upper and lower bounds of the loci for highlighting in `.plot_mqq`
- fixed errors in `.get_associations()`

# 3.6.7 20250827
- fixed error for `compare_effect()` when creating legend title
- fixed error for changing STATUS in some cases
- fixed error in `.remove_dup` when columns are missing
- added checking for NA strings in `fix_id()`
- added path manager
- added support for NEA from formatbook
- rearranged package structure
- Make int columns inclusive of range ends and float columns exclusive of range ends. (credit to @jhchung)

# 3.6.6 20250704
- fixed error in to_format() when CHR is not required.

# 3.6.5 20250625
- fixed scaled in plot_miami2

# 3.6.4 20250624
- fixed x tick labels in plot_miami2
- fixed xpad in plot_mqq

# 3.6.3 202505
- fixed pfile in `.clump()` 
- fixed typos 

# 3.6.2 20250509
- updated range for filtering OR and HR `(exp(-100), exp(100))`
- fixed typos

# 3.6.1 20250506
- updated dependencies

# 3.6.0 20250504
- fixed title_pad and add anno_args_single for specific labels for '.plot_mqq()'
- added infer ancesry
- updated main tutorial

# 3.5.8 20250424
- fixed bug for align when sumstats were not sorted
- added build check for `get_novel` when using efo ids.
- added auto-saving for API responses from GWASCatalog.

# 3.5.7 20250307
- updated Python version requirements (3.9 - 3.12)
- updated font family for mqq-related functions (default: Arial)
- adjusted default region_hspace for plot_stacked_mqq

# 3.5.6 20250306
- fixed error when region_recombination=False in regional plots
- updated cbar legends in regional plots
- fixed bugs for fig_args in `plot_miami2()`
- optimized gene track plotting
- added search and get_proxy
- increased mlog10p upper bound to 99999

# 3.5.5 20250102
- update 2025 footer
- added support for parquet
- added quick datatype conversion when loading sumstats
- fixed N datatype to np.int64 instead of pd.Int64 for ldsc
- updated requirement for matplotlib version
- organized examples
- added plot pipcs (under development)

# 3.5.4 20241218
-`gl.plot_miami2`: fixed error for `titles_pad`

# 3.5.3 20241217
-`to_format()`: added parquet format. 
-`to_format()`: added separate output for each chromosome  
-`to_format()`: updated log
-`plot_mqq(mode="r")`: added `region_legend_marker=True`
-`plot_mqq(mode="r")`: fixed legend maker size 
-`plot_mqq(mode="r")`: fixed coloring for single reference variant
-`plot_mqq()`: added `anno_xshift=0`
-`plot_stacked_mqq()`: fixed error when saving as pdf

# 3.5.2 20241203
- fixed size error for variants with very low MLOG10P when plotting region plots.
- fixed error when lead variant was not available in stacked regional plot.
- fixed legend title error in trumpet plot.
- updated functions for automatically extracting kwargs.

# 3.5.2 20241203
- added marker to indicate reference variant to LD colorbar in regional plots
- fixed errors when plotting y ticks in `plot_mqq()`
- added additional fix_id() in harmonization workflow (credit to @joshchiou)

# 3.5.1 20241120
- fixed errors due to liftover version change 
- restructured `compare_effect()`
- fixed error when getting hapmap3 variants

# 3.5.0 20241029
- added `plot_gwheatmap()`
- added genename annotation in `compare_effect()`
- fixed annotation error in `compare_effect()`
- fixed is_q error in `compare_effect()`
- fixed error for calculating arm length for annotation.
- removed some outdated code

# 3.4.49 2024
- fixed the suffix for yaml in `to_format()`
- implemented `.flip_snpid` and `.strip_snpid`
- implemented support for chm13 `13` in `.liftover`
- updated reference url list (added chain files)
- implemented `chrom_pat` and `snpid_pat` filters for `gl.Sumstats()` 
- fixed bug in `.plot_miami2`
- updated reference url

## 3.4.48 20240822
- fixed bug for `is_q_mc` in `compare_effect()`
- added `anno_height` for `plot_mqq()`
- added `xtight` for `plot_mqq()`
- fixed bug for `xpad` in `plot_mqq()`
- supported baseline model for ldsc outputs 

## 3.4.47 20240703
- updated stacked regional plot
- supported user-provided plink path
- fixed a bug in column order for `to_format()`

## 3.4.46 20240624
- added version requirements for numpy (<2) and matplotlib (<3.9)

## 3.4.45 20240509
- changed python version requirements to >=3.9, <3.11
- updated the version of pysam to v0.22.1
- fixed bug for status code when flipping statistics

## 3.4.44 20240424
- fixed bug when POS > sequence length for check_ref()
- vectorized normalize_allele() 

## 3.4.43 20240403
- [Added cache to speed up strand inference](https://github.com/Cloufield/gwaslab/pull/86) (credit to @sup3rgiu Mr. Andrea)
- fixed error for .plot_mqq(m="qq")
- Added h5py==3.10.0 to dependencies.

## 3.4.42 20240328
- fast implementation of check_ref and to_format. (credit to @sup3rgiu Mr. Andrea)
- added highlight and pinpoint for plot_trumpet()
- fixed typo (credit to @sup3rgiu Mr. Andrea)
- changed `**args` to `**kwargs` (credit to @sup3rgiu Mr. Andrea)
- implemented munge-like filters `munge=True` for ldsc in GWASLab

## 3.4.41 20240219
- fixed error in region_ref_second

## 3.4.40 20240215
- fixed color issue for regional plots. The color assigned to each variant is actually the color for the lower LD r2 category. For example, variants with LD>0.8 will be colored with the color for 0.8>LD>0.6.
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
- added support for GWASLab Sumstats object for gl.compare_effect()
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

- updated packaging methods. Now when installing GWASLab, pip will install all dependencies as well.

## v3.3.6 - 2022/11/05

- added download function: 
    - now you can download reference files from predefined list via GWASLab
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
- GWASLab is now able to read chromosome-separated files
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
