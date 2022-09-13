# Update Logs

## v3.2.0

Incorporated pyensembl and scikit-allel.

- Get_lead() : support automatic gene name annotation (using pyensembl)

- to_format():
	- support common sumstats formats 
	- support 1-based bed-like formats for VEP
	- support 0-based bed-like formats

- Manhattan plot:
	- optimized plotting logic
	- annotate gene names
	- added regional plot feature using a user-provided reference panel

- Comparison effect plot: 
	- fix using OR


## v3.1.0 

1. implemented formatbook: easily import sumstats and output sumstats in certain formats (support for commonly used formats including ldsc, plink, plink2, gwas-ssf, saige, regenie, fastgwa, metal, mrmega, pgscatalog, pgscatalog_hm, gwascatalog, gwascatalog_hm and gwaslab)

2. added `.filter_region_in/out` using bed files (or in-built regions like high-ld or hla)

3. implemented `.summay()` methods.

4. optimized rsID annotation pipeline. Support annotation using curated chr:pos:ref:alt - rsID tsv for quick annotation.

5. changed some datatypes and optimized memory usage.

6. replaced pyVCF with pySAM