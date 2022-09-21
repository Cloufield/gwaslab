#!/bin/zsh
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz
gzcat GRCh38_latest_genomic.gtf.gz |\
grep 'gene_biotype "protein_coding"' |\
sed 's/; gene "/; gene_name "/g' |\
gzip >GRCh38_latest_genomic.protein_coding.gtf.gz
