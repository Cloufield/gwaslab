# GWASLab

<img width="600" alt="image" src="https://github.com/user-attachments/assets/109262c6-c870-4078-94b5-66cf8c6b13c4" />

![badge](https://img.shields.io/pypi/v/gwaslab)
[![Downloads](https://static.pepy.tech/badge/gwaslab)](https://pepy.tech/project/gwaslab)
![badge_pip](https://img.shields.io/pypi/dm/gwaslab)
![badge_commit_m](https://img.shields.io/github/commit-activity/m/Cloufield/gwaslab)

* A handy Python-based toolkit for handling GWAS summary statistics (sumstats).
* Each process is modularized and can be customized to your needs.
* Sumstats-specific manipulations are designed as methods of a Python object, `gwaslab.Sumstats`.

## Installation

### install via pip

The latest version of GWASLab now supports Python 3.9, 3.10, 3.11, and 3.12.

```bash
pip install gwaslab
```

### install in conda environment

Create a Python 3.9, 3.10, 3.11 or 3.12 environment and install gwaslab using pip:

```bash
conda env create -n gwaslab -c conda-forge python=3.12

conda activate gwaslab

pip install gwaslab
```

or create a new environment using yml file [environment.yml](https://github.com/Cloufield/gwaslab/blob/main/environment.yml)

```bash
conda env create -n gwaslab -f environment.yml
```

### install using docker (deprecated)

A docker file is available [here](https://github.com/Cloufield/gwaslab/blob/main/docker/Dockerfile) for building local images.

## Quick start

```python

import gwaslab as gl

# load plink2 output
mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")

# or load sumstats with auto mode (auto-detecting commonly used headers) 
# assuming ALT/A1 is EA, and frq is EAF
mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="auto")

# or you can specify the columns:
mysumstats = gl.Sumstats("sumstats.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             eaf="Frq",
             beta="BETA",
             se="SE",
             p="P",
             direction="Dir",
             n="N",
             build="19")

# manhattan and qq plot
mysumstats.plot_mqq()
...
```

## Documentation and tutorials

Documentation and tutorials for GWASLab are avaiable at [here](https://cloufield.github.io/gwaslab/).

## Functions

### Loading and Formatting

- Loading sumstats by simply specifying the software name or format name, or specifying each column name.
- Converting GWAS sumstats to specific formats:
  - LDSC / MAGMA / METAL / PLINK / SAIGE / REGENIE / MR-MEGA / GWAS-SSF / FUMA / GWAS-VCF / BED... 
  - [check available formats](https://github.com/Cloufield/formatbook)
- Optional filtering of variants in commonly used genomic regions: Hapmap3 SNPs / High-LD regions / MHC region 
 
### Standardization & Normalization

- Variant ID standardization
- CHR and POS notation standardization
- Variant POS and allele normalization
- Genome build : Inference and Liftover 

### Quality control, Value conversion & Filtering

- Statistics sanity check
- Extreme value removal
- Equivalent statistics conversion
    - BETA/SE , OR/OR_95L/OR_95U
    - P, Z, CHISQ, MLOG10P
- Customizable value filtering

###  Harmonization

- rsID assignment based on CHR, POS, and REF/ALT
- CHR POS assignment based on rsID using a reference text file
- Palindromic SNPs and indels strand inference using a reference VCF
- Check allele frequency discrepancy using a reference VCF
- Reference allele alignment using a reference genome sequence FASTA file

### Visualization

- Mqq plot: Manhattan plot, QQ plot or MQQ plot (with a bunch of customizable features including auto-annotate nearest gene names)
- Miami plot: mirrored Manhattan plot
- Brisbane plot:  GWAS hits density plot
- Regional plot: GWAS regional plot
- Genetic correlation heatmap: ldsc-rg genetic correlation matrix
- Scatter plot: variant effect size comparison
- Scatter plot: allele frequency comparison 
- Scatter plot: trumpet plot (plot of MAF and effect size with power lines)

#### Visualization Examples

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/233836639-34b03c47-5a59-4fd4-9677-5e13b02aab15.png">
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197393168-e3e7076f-2801-4d66-9526-80778d44f3da.png">
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197463243-89352749-f882-418d-907d-27530fd4e922.png">
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

### Other Utilities

- Read ldsc h2 or rg outputs directly as DataFrames (auto-parsing).
- Extract lead variants given a sliding window size.
- Extract novel loci given a list of known lead variants / or known loci obtained from GWAS Catalog.
- Logging: keep a complete record of manipulations applied to the sumstats.
- Sumstats summary: give you a quick overview of the sumstats. 
- ...

## Issues

- GWASLab is currently under active development, with frequent updates.
- Note: Known issues are documented at https://cloufield.github.io/gwaslab/KnownIssues/.

## How to cite
- GWASLab preprint: He, Y., Koido, M., Shimmori, Y., Kamatani, Y. (2023). GWASLab: a Python package for processing and visualizing GWAS summary statistics. Preprint at Jxiv, 2023-5. https://doi.org/10.51094/jxiv.370

## Sample data used for tutorial
- Sample GWAS data used in GWASLab is obtained from: http://jenger.riken.jp/ (Suzuki, Ken, et al. "Identification of 28 new susceptibility loci for type 2 diabetes in the Japanese population." Nature genetics 51.3 (2019): 379-386.).

## Acknowledgement

Thanks to @sup3rgiu, @soumickmj and @gmauro for their contributions to the source codes.

## Contacts
* Github: [https://github.com/Cloufield/gwaslab](https://github.com/Cloufield/gwaslab)
* Blog (in Chinese): [https://gwaslab.com/](https://gwaslab.com/)
* Email: gwaslab@gmail.com
* Stats: [https://pypistats.org/packages/gwaslab](https://pypistats.org/packages/gwaslab)
