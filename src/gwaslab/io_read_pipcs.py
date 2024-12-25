import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.qc_check_datatype import check_datatype
from gwaslab.qc_check_datatype import check_dataframe_memory_usage

def _read_pipcs(data, output_prefix, log=Log(),verbose=True):
    log.write("Start to load PIP and CREDIBLE_SET_INDEX from file...",verbose=verbose)
    log.write(" -File:{}.pipcs".format(output_prefix),verbose=verbose)

    pipcs = pd.read_csv("{}.pipcs".format(output_prefix))
    
    log.write(" -Merging CHR and POS from main dataframe...",verbose=verbose)
    pipcs = _merge_chrpos(data,pipcs)

    log.write(" -Current pipcs Dataframe shape :",len(pipcs)," x ", len(pipcs.columns),verbose=verbose) 
    check_datatype(pipcs,log=log,verbose=verbose)
    check_dataframe_memory_usage(pipcs,log=log,verbose=verbose)
    log.write("Finished loading PIP and CREDIBLE_SET_INDEX from file!",verbose=verbose)
    return pipcs

def _merge_chrpos(data,pipcs):
    df = pd.merge(pipcs, data,on="SNPID",how="left")
    return df