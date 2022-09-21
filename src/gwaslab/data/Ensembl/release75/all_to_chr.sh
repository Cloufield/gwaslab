#!/bin/bash
for i in $(seq 1 22) X Y MT
do
zcat Homo_sapiens.GRCh37.75.gtf.gz | awk -v FS='\t' -v i=$i '$1==i{print $0}' | gzip > Homo_sapiens.GRCh37.75.chr${i}.gtf.gz
done
