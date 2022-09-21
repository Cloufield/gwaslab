#!/bin/zsh
gzcat Homo_sapiens.GRCh37.75.gtf.gz|\
grep -E 'processed_transcript|protein_coding|_gene' |\
gzip >Homo_sapiens.GRCh37.75.protein_coding.gtf.gz 

