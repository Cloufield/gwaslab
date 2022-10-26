#!/bin/zsh
for i in "X"
#$(seq 1 22)
do
gzcat ../t2d_bbj.txt.gz |awk -v i=${i} 'NR==1 || $2==i {print $0}' >t2d_bbj_chr${i}.txt
done
