import pandas as pd
import os
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished

def process_vcf_to_hfd5(vcf, 
                    directory=None, 
                    chr_dict=None, 
                    group_size=20000000,
                    complevel=9,
                    chunksize=20000000,
                    log=Log(),
                    verbose=True):
    
    #load vcf
    log.write("Start to process VCF file to HDF5:", verbose=verbose)
    log.write(" -Reference VCF path:{}".format(vcf), verbose=verbose)
    log.write(" -Output group size:{}".format(group_size), verbose=verbose)
    log.write(" -Compression level:{}".format(complevel), verbose=verbose)
    log.write(" -Loading chunksize:{}".format(chunksize), verbose=verbose)

    vcf_file_name = os.path.basename(vcf)
    vcf_dir_path = os.path.dirname(vcf)
    
    if directory is None:
        directory = vcf_dir_path
    elif directory[-1] == "/":
        directory = directory.rstrip('/')
    
    h5_path = "{}/{}.rsID_CHR_POS_groups_{}.h5".format(directory,vcf_file_name,int(group_size))
    log_path = "{}/{}.rsID_CHR_POS_groups_{}.log".format(directory,vcf_file_name, int(group_size))
    log.write(" -HDF5 Output path: {}".format(h5_path), verbose=verbose)
    log.write(" -Log output path: {}".format(log_path), verbose=verbose)
    df = pd.read_table(vcf,comment="#",usecols=[0,1,2],header=None,chunksize=chunksize)

    log.write(" -Processing chunk: ",end="", verbose=verbose)
    
    for index,chunk in enumerate(df):
        log.write(index,end=" ",show_time=False, verbose=verbose)
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
    log.write("Processing finished!", verbose=verbose)
    log.save(log_path, verbose=verbose)