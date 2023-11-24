import pandas as pd
from gwaslab.bd_common_data import get_formats_list
from gwaslab.g_Log import Log
from gwaslab.bd_common_data import get_format_dict

def _read_tabular(path, fmt, **args):
    
    # default
    load_args_dict = {"sep":"\t",
                      "header":None}
    
    # if specified by user
    if len(args)>0:
        load_args_dict = args
    
    # load format
    meta_data, rename_dictionary = get_format_dict(fmt)
    
    if "format_separator" in meta_data and "sep" not in args:
        load_args_dict["sep"] = meta_data["format_separator"]
    
    if "format_comment" in meta_data and "comment" not in args:
        if  meta_data["format_comment"] is not None:    
            load_args_dict["comment"] = meta_data["format_comment"]

    if "format_header" in meta_data and "header" not in args:
        load_args_dict["header"] = meta_data["format_header"]

    if "format_na" in meta_data and "na_values" not in args:
        if  meta_data["format_na"] is not None:    
            load_args_dict["na_values"] = meta_data["format_na"]

    #######################################################################################
    df = pd.read_csv(path, **load_args_dict)
    #######################################################################################

    # configure header
    if "format_header" in meta_data:
        if meta_data["format_header"] is None:
            num_to_name =  {}
            for key,value in rename_dictionary.items():
                num_to_name[int(key)] = value

    # rename columns
    if "format_header" in meta_data:
        if meta_data["format_header"] is None:
            df = df.rename(columns=num_to_name)
        else:
            df = df.rename(columns=rename_dictionary)
    else:
        df = df.rename(columns=rename_dictionary)

    # convert datatype
    if "format_datatype" in meta_data:
        df = df.astype(meta_data["format_datatype"])

    return df

def read_bim(path):
    df = _read_tabular(path,fmt="plink_bim")
    return df

def read_fam(path):
    df = _read_tabular(path,fmt="plink_fam")
    return df

def read_psam(path):
    df = _read_tabular(path,fmt="plink_psam")
    return df

def read_pvar(path):
    df = _read_tabular(path,fmt="plink_pvar")
    return df

def read_bgen_sample(path):
    df = _read_tabular(path,fmt="bgen_sample")
    return df