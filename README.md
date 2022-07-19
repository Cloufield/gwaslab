![gwaslab_logo](https://cloufield.github.io/gwaslab/images/index_logo.jpg)
A collection of handy python scripts for GWAS. 
Just want to make life eaiser and save myself from repetitive work.

For usage, please check GWASLab document at https://cloufield.github.io/gwaslab/ .

## What you can do with gwaslab:
### Standardization, Normalization & Harmonization

- CHR and POS notation standardization

- Variant POS and allele normalization

- Genome build : Liftover

- Reference allele alignment using a reference genome sequence

- rsID assignment based on CHR, POS, REF and ALT

- CHR POS assignment based on rsID using a reference VCF

- Palindromic SNPs and indels strand inference using a reference VCF

- Check allele frequency discrepancy using a reference VCF

### Quality control, Value conversion & Filtering

- Statistics sanity check

- Equivalent statistics conversion
  
  - BETA/SE , OR/OR_95L/OR_95U
  
  - P, Z, CHISQ, MLOG10

- Extract/exclude hapmap3 variants 

- Extract/exclude MHC variants

- Filtering values.

### Visualization

- Mqq plot : Manhattan plot and QQ plot side by side

- Heatmap : ldsc-rg genetic correlation matrix

- Scatter Plot : variant effect size comparison with sumstats

- Scatter Plot : allele frequency comparison 

### Other Utilities

- Read ldsc h2 or rg outputs directly as DataFrames

- Extract lead SNPs given a window size

- Logging : keep a complete record of manipulations from raw data to munged data

- Formating GWAS sumstats in certain formats
  
  - LDSC / MAGMA / METAL / MR-MEGA / FUMA
...




## Requirements
- Python >= 3
- pyVCF >= 0.6.8
- Biopython >= 1.79
- liftover >= 1.1.13
- pandas >= 1.2.4
- numpy >= 1.21.2
- matplotlib>3.5
- seaborn >= 0.11.1
- scipy >= 1.6.2
- adjustText

## Install
```
pip install gwaslab
```
Current version: 3.0.0

# Usage

For usage, please check GWASLab document at https://cloufield.github.io/gwaslab/ .

# Update Log
- 3.0.0 first complete version
- 1.0.0 implemented Sumstats object

- 0.0.5 - 0.0.6
- added  compare_effect, read_ldsc 

- 0.0.4  
  -  added mqqplot feature
  -  fixed gtesig algorithm
  -  recreated mplot and qqplot


For more information: 
https://gwaslab.com/
