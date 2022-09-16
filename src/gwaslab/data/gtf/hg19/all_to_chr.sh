#!/bin/zsh
for i in $(seq 1 22) X Y MT
do
gzcat hg19.refGene.gtf.gz | awk -v FS='\t' -v i=$i '$1=="chr"i{print $0}' > hg19.refGene.chr${i}.gtf
gzip hg19.refGene.chr${i}.gtf
done
