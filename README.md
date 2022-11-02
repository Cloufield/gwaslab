<img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197167760-5f761f5e-5856-4b27-a540-8b9cd90bdadb.png">

# gwaslab

![badge](https://img.shields.io/badge/release-v3.3.4-brightgreen)
![badge_pip](https://img.shields.io/pypi/dm/gwaslab)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fcloufield.github.io%2Fgwaslab%2F&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

* A simple python toolkit for handling GWAS sumstats.
* Each process is modularized and can be customized to your needs.
* Sumstats-specific manipulations are designed as methods of a python object, `gwaslab.Sumstats`.

Please check GWASLab document at [https://cloufield.github.io/gwaslab/](https://cloufield.github.io/gwaslab/)


## Install

```
pip install gwaslab==3.3.4
```


```python
import gwaslab as gl
# load plink2 output
mysumstats = gl.Sumstats("t2d_bbj.txt.gz", fmt="plink2")

# or you can specify the columns:
mysumstats = gl.Sumstats("t2d_bbj.txt.gz",
             snpid="SNP",
             chrom="CHR",
             pos="POS",
             ea="ALT",
             nea="REF",
             neaf="Frq",
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

## Functions

### Loading and Formatting

- Loading sumstats by simply specifying the software name or format name.
- Converting GWAS sumstats to specific formats:
  - LDSC / MAGMA / METAL / MR-MEGA / FUMA / VCF / BED... 
  - [check available formats](https://github.com/Cloufield/formatbook)
- Optional Filtering of variants in Hapmap3 SNPs / High-LD regions / HLA region 
 
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
    - P, Z, CHISQ, MLOG10
- Customized value filtering

###  Harmonization

- rsID assignment based on CHR, POS, and REF/ALT
- CHR POS assignment based on rsID using a reference text file
- Palindromic SNPs and indels strand inference using a reference VCF
- Check allele frequency discrepancy using a reference VCF
- Reference allele alignment using a reference genome sequence FASTA file

### Visualization

- Mqq plot : Manhattan plot , QQ plot or MQQ plot (with a bunch of customizable features including auto-annotate nearest gene names)
- Miami plot : Manhattan plot
- Brisbane plot:  GWAS hits density plot
- Regional plot : GWAS regional plot
- Heatmap : ldsc-rg genetic correlation matrix
- Scatter Plot : variant effect size comparison with sumstats
- Scatter Plot : allele frequency comparison 
- Forest Plot : forest plots for meta-analysis of SNPs

### Visualization Examples

<img width="700" alt="image" src="https://user-images.githubusercontent.com/40289485/195526882-dff70593-752d-4672-901e-0b3cea3d8cda.png"><img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197393168-e3e7076f-2801-4d66-9526-80778d44f3da.png"><img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197463243-89352749-f882-418d-907d-27530fd4e922.png"><img width="600" alt="image" src="https://user-images.githubusercontent.com/40289485/197126045-b1c55adf-3391-4c3d-b2f6-eaeac7c26024.png">

<img width="400" alt="image" src="https://user-images.githubusercontent.com/40289485/195526481-df060ad5-dc61-4e35-ab37-3ea45ed00618.png">

### Other Utilities

- Read ldsc h2 or rg outputs directly as DataFrames (auto-parsing).
- Extract lead variants given a sliding window size.
- Extract novel loci given a list of known lead variants.
- Logging : keep a complete record of manipulations from raw data to munged data.
- Sumstats summary function: know your data better. 



## Requirements:

- Python >= 3.6
- pySAM 
- pyensembl
- scikit-allel
- Biopython >= 1.79
- liftover >= 1.1.13
- pandas >= 1.2.4
- numpy >= 1.21.2
- matplotlib >=3.5
- seaborn >=0.11.1
- scipy >=1.6.2
- statsmodels > =0.13
- adjustText

# Contacts
* Github: [https://github.com/Cloufield/gwaslab](https://github.com/Cloufield/gwaslab)
* Blog (in Chinese): [https://gwaslab.com/](https://gwaslab.com/)
* Email: gwaslab@gmail.com
* Stats: [https://pypistats.org/packages/gwaslab](https://pypistats.org/packages/gwaslab)

