import pandas as pd
import os
import numpy as np
from gwaslab.g_Log import Log

def process_ref_vcf(vcf, directory=None, chr_dict=None, group_size=20000000,complevel=9,chunksize=20000000,log=Log()):
    #load vcf
    log.write("Start processing VCF files:")
    log.write(" -Reference VCF path:{}".format(vcf))
    log.write(" -Output group size:{}".format(group_size))
    log.write(" -Compression level:{}".format(complevel))
    log.write(" -Loading chunksize:{}".format(chunksize))

    if directory is None:
        directory="./"

    elif directory[-1] == "/":
        directory = directory.rstrip('/')
    
    h5_path = "{}/rsID_CHR_POS_groups_{}.h5".format(directory,int(group_size))
    log_path = "{}/rsID_CHR_POS_groups_{}.log".format(directory,int(group_size))
    log.write(" -HDF5 Output path: {}".format(h5_path))
    log.write(" -Log output path: {}".format(log_path))
    df = pd.read_table(vcf,comment="#",usecols=[0,1,2],header=None,chunksize=chunksize)

    
    log.write(" -Processing chunk: ",end="")
    
    for index,chunk in enumerate(df):
        log.write(index,end=" ",show_time=False)
        chunk = chunk.rename(columns={0:"CHR",1:"POS",2:"rsn"})
        if chr_dict is not None:
            chunk["CHR"] = chunk["CHR"].map(chr_dict)
        
        chunk["rsn"] = chunk["rsn"].str.strip("rs")
        chunk = chunk.dropna()
        chunk = chunk.drop_duplicates(subset="rsn")
        chunk = chunk.astype("int64")
        
        if len(chunk)>0:
            chunk["group"]=chunk["rsn"]//group_size
            for i in chunk["group"].unique():
                chunk.loc[chunk["group"]==i,["CHR","POS","rsn"]].to_hdf(h5_path,
                                                                        key="group_"+str(i),
                                                                        append=True,
                                                                        index=None,
                                                                        dropna=True,
                                                                        format="table",
                                                                        complevel=complevel)
    log.write("Processing finished!")
    log.save(log_path, verbose=False)