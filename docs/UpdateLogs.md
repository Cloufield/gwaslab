# Update Logs

## v3.3.5 coming soon

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
