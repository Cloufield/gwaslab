# pyensembl Reference file

https://github.com/openvax/pyensembl

Install pyensembl and download reference:

```
# install  pyensembl if not
pip install pyensembl

# syntax for download reference for pyensembl
pyensembl install --release <list of Ensembl release numbers> --species <species-name>
```

For gwaslab, please run the following commands:

ensembl release 75 : hg19

ensembl release 76 : hg38 

```
pyensembl install --release 75 76 --species human
```

gwaslab could use ensembl reference data to annotate lead SNPs with the nearest gene name.



# Process Reference file

## 1000 Genome

Download:

[Index of /vol1/ftp/release/20130502/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

After downloading the raw vcf, we need to normalize the variants, split multiallelic variants, rename the variant and remove duplicates. 

Also, we need to extract out target  ancestry and recalculate the allele frequency. 

Create a sample list for EAS samples

```bash
 awk '$3=="EAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel >EAS.sample
```

1000 genome:

```bash
#!/bin/bash

for chr in {1..22}
do
    bcftools view -S EAS.sample ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
        bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta | \
        bcftools annotate -x ID,INFO -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm --rm-dup both | \
        bcftools +fill-tags -Oz  -- -t AF  \
          > EAS.chr"${chr}".split_norm_af.vcf.gz
    tabix -p vcf EAS.chr"${chr}".split_norm_af.vcf.gz
done
```

# dbsnp

rsID database.

dbsnp v151 (hg19): https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz

latest release: [Index of /snp/latest_release/VCF](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/)

# gnomad

Allele frequency for major ancestries and rsID

gnomAD v2 & gnomAD v2 liftover & gnomAD v3:   [gnomAD](https://gnomad.broadinstitute.org/downloads)    


