import pandas as pd
import os

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

for i in range(0,75):
    path = "/home/heyunye/mydata/reference/snp/organisms/human_9606_b151_GRCh37p13/VCF/db151_rsn_"+str(i)+".vcf"
    print(i)
    if is_non_zero_file(path) == True:
        df = pd.read_csv(path,"\s+",header= None,dtype={0:"int",1:"str",2:"int"})
        df.rename(columns={0:"rsn",1:"CHR",2:"POS"},inplace=True)
        df.set_index("rsn",inplace=True)
        df.to_hdf("rsid_chrpos_db151_hg19_uniq.h5", key="part"+str(i))