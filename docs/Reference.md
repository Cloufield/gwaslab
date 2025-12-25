# Reference data for handling Sumstats
## Reference genome sequence for variant allele alignment

- GRCh37 / hg19 : [ucsc_hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)
- GRCh38 / hg38 : [ucsc_hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)

For details about reference genome, please check [https://cloufield.github.io/CTGCatalog/Reference_data_Genome_README/](https://cloufield.github.io/CTGCatalog/Reference_data_Genome_README/)

## Processed Reference files for harmonization
### 1000 Genome Project(hg19)
Download:
[Index of /vol1/ftp/release/20130502/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

Process the 1000 genome vcf file for EAS sample:

- Extract EAS sample
- Split multiallelic variants
- Normalize variants
- Rename variants as CHR:POS:REF:ALT
- Remove duplicated variants
- Calculate Alternative allele frequency
- Tabix index

### Sample code:
```python
#!/bin/bash
# extract EAS sample ID
awk '$3=="EAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel >EAS.sample

# process
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

# merge
bcftools concat -a -d both -f concat_list.txt -Ob | bcftools sort -Oz  > EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
tabix -p vcf EAS.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz
```
The processed dataset can then be used for:

- infer palindromic SNPs / indels
- check allele frequency difference 
- regional plot

## rsID conversion table for quick annotation
Since GWASLab will check chr:pos and ea/nea to match rsID, it would take a little bit longer if we only use vcf. 

But we can use a pre-annotated conversion table for common SNPs, and then annotate the rest of SNPs using large VCF file from dbSNP. 
```python
for i in range(1,23):
    mysumstats = gl.Sumstats("./EAS.chr"+str(i)+".split_norm_af.vcf.gz",snpid="ID",fmt="vcf")
    mysumstats.harmonize( basic_check=True,
                        ref_seq="./reference_genome/hg19/human_g1k_v37_decoy.fasta",
                        ref_rsid_vcf="./All_20180423.vcf.gz", 
                        threads=4)
    mysumstats.data.loc[:,["SNPID","rsID","CHR","POS","NEA","EA"]].to_csv("./1kg_af_dbsnp151."+str(i)+".txt.gz","\t",index=None)
```

In terminal, combine the files:
```python
# get header
zcat 1kg_af_dbsnp151.1.txt.gz | head -1 > 1kg_af_dbsnp151_auto.txt

for i in $(seq 1 22)
do
# get complete SNPID-rsID pairs
zcat 1kg_af_dbsnp151.${i}.txt.gz | awk -v FS="\t" -v OFS="\t" 'NR>1 && $2!="" {print $0}' >>1kg_af_dbsnp151_auto.txt
done
```

## dbsnp: database for annotation of rsID

- dbsnp v151 (GRCh37): [https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz)
- dbsnp v151 (GRCh37): [https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.tbi)
- dbsnp v151 (GRCh38): [https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF//00-All.vcf.gz](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF//00-All.vcf.gz)
- dbsnp v151 (GRCh38): [https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF//00-All.vcf.gz.tbi](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF//00-All.vcf.gz.tbi)
- latest release: [Index of /snp/latest_release/VCF](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/)

## gnomad

- Allele frequency for major ancestries and rsID
- gnomAD v2 & gnomAD v2 liftover & gnomAD v3:   [gnomAD](https://gnomad.broadinstitute.org/downloads)    


## Reference Library for variants annotation 

GWASLab uses:

- ensembl release 87 (hg19): https://ftp.ensembl.org/pub/grch37/release-109/gtf/homo_sapiens/
- ensembl release 109 (hg38):  https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/
- NCBI refseq GRCh37 : https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz
- NCBI refseq GRCh38 : https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz

Currently, GWASLab could use either ensembl or refseq gtf reference data to annotate lead SNPs with the nearest gene name.
