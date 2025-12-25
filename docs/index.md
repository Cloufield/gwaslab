# GWASLab

<img width="600" alt="image" src="https://github.com/user-attachments/assets/109262c6-c870-4078-94b5-66cf8c6b13c4" />

![badge](https://img.shields.io/pypi/v/gwaslab)
[![Downloads](https://static.pepy.tech/badge/gwaslab)](https://pepy.tech/project/gwaslab)
![badge_pip](https://img.shields.io/pypi/dm/gwaslab)
![badge_commit_m](https://img.shields.io/github/commit-activity/m/Cloufield/gwaslab)

* A handy Python-based toolkit for handling GWAS summary statistics (sumstats).
* Each process is modularized and can be customized to your needs.
* Sumstats-specific manipulations are designed as methods of a Python object, `gwaslab.Sumstats`.

!!! info "GWASLab v4.0.0"
    GWASLab v4.0.0 introduces major improvements including:
    
    - **Command Line Interface (CLI)**: Process sumstats directly from the command line
    - **Performance Improvements**: Optimized algorithms and data structures for faster processing
    - **Visualization Parameter Management**: Centralized parameter management system for all plotting functions
    - **API Consistency**: Standardized parameter names across functions
    - **Enhanced Documentation**: Comprehensive documentation and tutorials

## Installation

### Install via pip

The latest version of GWASLab now supports Python 3.9, 3.10, 3.11, and 3.12.

!!! tip "Recommended: Python 3.12"
    We recommend using **Python 3.12** for the best performance and compatibility with GWASLab v4.0.0.

```
pip install gwaslab
```

### Install via uv

[uv](https://github.com/astral-sh/uv) is a fast Python package installer and resolver written in Rust.

```
uv pip install gwaslab
```

### Install in conda environment

Create a Python 3.9, 3.10, 3.11 or 3.12 environment and install gwaslab using pip. We recommend Python 3.12:

```
conda create -n gwaslab -c conda-forge python=3.12

conda activate gwaslab

pip install gwaslab
```

or create a new environment using yml file [environment.yml](https://github.com/Cloufield/gwaslab/blob/main/environment.yml)

```
conda env create -n gwaslab -f environment.yml -c conda-forge
```

### Install using docker (deprecated)

A docker file is available [here](https://github.com/Cloufield/gwaslab/blob/main/docker/Dockerfile) for building local images.

## Quick Start

### Python API

```
import gwaslab as gl

# load plink2 output
mysumstats = gl.Sumstats("sumstats.txt.gz", fmt="plink2")

# or load sumstats with auto mode (auto-detecting commonly used headers) 
# assuming ALT/A1 is **EA**, and frq is **EAF**
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

# basic quality control
mysumstats.basic_check()

# manhattan and qq plot
mysumstats.plot_mqq()

# extract lead variants
lead_variants = mysumstats.get_lead()
```

### Command Line Interface (CLI) - New in v4.0.0

```
# Show version
gwaslab version

# Basic QC and output
gwaslab --input sumstats.tsv --fmt auto --qc --out cleaned.tsv --to-fmt gwaslab

# Harmonization with reference
gwaslab --input sumstats.tsv --fmt auto --ref-seq ref.fasta --harmonize --out harmonized.tsv --to-fmt gwaslab

# Format conversion
gwaslab --input sumstats.tsv --fmt gwaslab --out sumstats.ldsc --to-fmt ldsc
```

See [CLI documentation](CLI.md) for more details.

## Documentation and Tutorials

Documentation and tutorials for GWASLab are available at [https://cloufield.github.io/gwaslab/](https://cloufield.github.io/gwaslab/).

## Features

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

### Quality Control, Value Conversion & Filtering

- Statistics sanity check
- Extreme value removal
- Equivalent statistics conversion
    - BETA/SE , OR/OR_95L/OR_95U
    - P, Z, CHISQ, MLOG10P
- Customizable value filtering

### Harmonization

- rsID assignment based on CHR, POS, and REF/ALT (with sweep mode for large datasets)
- CHR POS assignment based on rsID using a reference text file
- Palindromic SNPs and indels strand inference using a reference VCF
- Check allele frequency discrepancy using a reference VCF
- Reference allele alignment using a reference genome sequence FASTA file

### Visualization

- **MQQ plot**: Manhattan plot, QQ plot or MQQ plot (with a bunch of customizable features including auto-annotate nearest gene names)
- **Miami plot**: Mirrored Manhattan plot
- **Brisbane plot**: GWAS hits density plot
- **Regional plot**: GWAS regional plot with LD information
- **Genetic correlation heatmap**: LDSC-rg genetic correlation matrix
- **Scatter plots**: 
  - Variant effect size comparison
  - Allele frequency comparison 
  - Trumpet plot (plot of MAF and effect size with power lines)

#### Visualization Examples

<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/233836639-34b03c47-5a59-4fd4-9677-5e13b02aab15.png">
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197393168-e3e7076f-2801-4d66-9526-80778d44f3da.png">
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197463243-89352749-f882-418d-907d-27530fd4e922.png">
<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

### Other Utilities

- Read LDSC h2 or rg outputs directly as DataFrames (auto-parsing).
- Extract lead variants given a sliding window size.
- Extract novel loci given a list of known lead variants or known loci obtained from GWAS Catalog.
- Clumping and fine-mapping utilities.
- LD score regression (LDSC) integration for heritability estimation.
- Logging: keep a complete record of manipulations applied to the sumstats.
- Sumstats summary: give you a quick overview of the sumstats.
- ...

## Command Line Interface (CLI) - New in v4.0.0

GWASLab now provides a unified command-line interface for processing GWAS summary statistics. The CLI supports:

- Quality control (QC)
- Harmonization
- Format conversion
- Batch processing
- Various output formatting options

See [CLI documentation](CLI.md) for complete details and examples.

## Performance Improvements in v4.0.0

- **Optimized algorithms**: Core functions have been reimplemented with optimized algorithms and data structures
- **Sweep mode**: New sweep mode for rsID assignment and strand inference, significantly faster for large datasets
- **Efficient data processing**: Improved memory usage and processing speed for large datasets
- **Better parallelization**: Enhanced multi-threading support

## Issues

- GWASLab is currently under active development, with frequent updates.
- Note: Known issues are documented at https://cloufield.github.io/gwaslab/KnownIssues/.

## How to Cite

- GWASLab preprint: He, Y., Koido, M., Shimmori, Y., Kamatani, Y. (2023). GWASLab: a Python package for processing and visualizing GWAS summary statistics. Preprint at Jxiv, 2023-5. https://doi.org/10.51094/jxiv.370

## Sample Data Used for Tutorial

- Sample GWAS data used in GWASLab is obtained from: http://jenger.riken.jp/ (Suzuki, Ken, et al. "Identification of 28 new susceptibility loci for type 2 diabetes in the Japanese population." Nature genetics 51.3 (2019): 379-386.).

## Acknowledgement

Thanks to @sup3rgiu, @soumickmj and @gmauro for their contributions to the source codes.

## Contacts

* Github: [https://github.com/Cloufield/gwaslab](https://github.com/Cloufield/gwaslab)
* Blog (in Chinese): [https://gwaslab.com/](https://gwaslab.com/)
* Email: gwaslab@gmail.com
* Stats: [https://pypistats.org/packages/gwaslab](https://pypistats.org/packages/gwaslab)
