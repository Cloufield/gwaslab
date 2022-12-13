# Update Logs
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
## v3.3.6 -20221105
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

## v3.3.5 -20221102 
- added `filter_value`
- integrate `gwascatalog` to `get_novel`
- optimized `remove_dup`
- fixed bugs

## v3.3.4 -20221031
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

## v3.3.2 - 20221021
- added forcefixid for fix_id()
- fixed bugs for plotting gene tracks

## v3.3.1 
- extract novel loci given a list of known lead variants
- fixed bugs in fill_data()
- fixed path for hapmap3 snps for infer_build()

## v3.3.0 2022/10/18
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
