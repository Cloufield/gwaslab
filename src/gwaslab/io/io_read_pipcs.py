import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_memory_usage
import re
import os

def _read_pipcs(data_or_dataframe, 
                output_prefix, 
                study=None, 
                group=None,
                studie_names=None,
                log=Log(),
                verbose=True,
                **readcsv_kwargs):
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(data_or_dataframe, pd.DataFrame):
        data = data_or_dataframe
    else:
        data = data_or_dataframe.data
    
    log.write("Start to load PIP and CREDIBLE_SET_INDEX from file...",verbose=verbose)
    log.write(" -File:{}".format(output_prefix),verbose=verbose)
    
    if "@" in output_prefix:
        log.write(" -Detected @ in path: load all matching pipcs files ...",verbose=verbose)
        pipcs_path_list = []
        pipcs_loci_list = []

        dirname = os.path.dirname(output_prefix)
        files = os.listdir(dirname)
        target_file_name = os.path.basename(output_prefix).replace('@','([\w:_]+)')
        for file in files:
            if re.search(target_file_name, file) is not None:
                pipcs_path_list.append(dirname+"/"+file)
                pipcs_loci_list.append(re.search(target_file_name, file)[1])

        pipcs_single_list=[]
        for index,pipcs_path in enumerate(pipcs_path_list):
            log.write(" -Loading {}:".format(pipcs_loci_list[index]) + pipcs_path)
            pipcs_single = pd.read_csv(pipcs_path,**readcsv_kwargs)
            if "LOCUS" not in pipcs_single.columns:
                pipcs_single["LOCUS"]=pipcs_loci_list[index]
            pipcs_single_list.append(pipcs_single)

        pipcs = pd.concat(pipcs_single_list, axis=0, ignore_index=True) 
    else:
        pipcs = pd.read_csv("{}".format(output_prefix),**readcsv_kwargs)
    
    if "CHR" not in pipcs.columns:
        log.write(" -Merging CHR and POS from main dataframe...",verbose=verbose)
        pipcs = _merge_chrpos(data,pipcs)

    pipcs = pipcs.rename(columns={
        "cs":"CREDIBLE_SET_INDEX",
        "variable_prob":"PIP",
        "variable":"N_SNP"
    })

    log.write(" -Current pipcs Dataframe shape :",len(pipcs)," x ", len(pipcs.columns),verbose=verbose) 
    
    if group is not None:
        pipcs["GROUP"] = group
    if study is not None:
        pipcs["STUDY"] = study

    pipcs = _process_pip(pipcs, group, studie_names)

    check_datatype(pipcs,log=log,verbose=verbose)
    check_dataframe_memory_usage(pipcs,log=log,verbose=verbose)
    log.write("Finished loading PIP and CREDIBLE_SET_INDEX from file!",verbose=verbose)
    return pipcs

def _merge_chrpos(data,pipcs):
    df = pd.merge(pipcs, data,on="SNPID",how="left")
    return df

def _process_pip(pipcs, group, studie_names):
    if group is not None and "PIP" not in pipcs.columns:
        pipcs["PIP"] = pipcs[studie_names]

        for i in pipcs["CS_CATEGORY"].dropna().unique():
            print(i)
            pipcs.loc[pipcs["CS_CATEGORY"]==i,"PIP"] = pipcs.loc[pipcs["CS_CATEGORY"]==i,i]
    return pipcs