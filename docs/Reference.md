# Reference data for handling Sumstats
## Reference Library for variants annotation : pyensembl

Install pyensembl and download reference:

```
# install  pyensembl if not
pip install pyensembl

# syntax for download reference for pyensembl
pyensembl install --release <list of Ensembl release numbers> --species <species-name>
```
GWASlab uses:
- ensembl release 75 : hg19
- ensembl release 107 : hg38 
- NCBI refseq GRCh37
- NCBI refseq GRCh38

Currently, gwaslab could use ensembl or refseq gtf reference data to annotate lead SNPs with the nearest gene name.
For details about pyensembl, please check [https://github.com/openvax/pyensembl](https://github.com/openvax/pyensembl)

## Reference genome sequence for variant allele alignment

- GRCh37 / hg19 : [ucsc_hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)
- GRCh38 / hg38 : [ucsc_hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

For details about reference genome, please check [https://cloufield.github.io/CTGCatalog/Reference_data_Genome_README/](https://cloufield.github.io/CTGCatalog/Reference_data_Genome_README/)

## Processed Reference files for harmonization
### 1000 Genome Project(hg19)
Download:
[Index of /vol1/ftp/release/20130502/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

- After downloading the raw vcf, we need to normalize the variants, split multiallelic variants, rename the variant and remove duplicates. 
- Also, we need to extract out target  ancestry and recalculate the allele frequency. 
- Create a sample list for EAS samples

```bash
 awk '$3=="EAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel >EAS.sample
```

Process the 1000 genome vcf file for EAS sample:

- Extract EAS sample
- Split multiallelic variants
- Normalize variants
- Rename variants as CHR:POS:REF:ALT
- Remove duplicated variants
- Calculate Alternative allele frequency
- Tabix index

### Sample code:
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
    echo "EAS.chr"${chr}".split_norm_af.vcf.gz" >>concat_list.txt 
done

bcftools concat -a -d both -f concat_list.txt -Ob | bcftools sort -Oz  > EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
tabix -p vcf EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
```

## dbsnp: database for annoattion of rsID

- dbsnp v151 (hg19): https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
- latest release: [Index of /snp/latest_release/VCF](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/)

## gnomad

- Allele frequency for major ancestries and rsID
- gnomAD v2 & gnomAD v2 liftover & gnomAD v3:   [gnomAD](https://gnomad.broadinstitute.org/downloads)    
