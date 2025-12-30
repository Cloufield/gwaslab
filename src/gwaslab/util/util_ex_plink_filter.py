from typing import Tuple, Any, List
import pandas as pd

def _run_plink_filter(filter_flag: str, out_prefix: str) -> None:
    '''
    run plink filter functions to generate bim and fam files
    
    Returns:
        bim as pd.DataFrame: SNPID,CHR,POS,NEA,EA
        fam as pd.DataFrame: FID,IID
    '''
    plink_script='''
    plink2 \
    {} \
    --make-just-bim \
    --make-just-fam \ 
    --out {}
    '''.format(filter_flag, out_prefix)

def _plink2_filter_to_flag(tmpdir: str = "./", **kwargs: Any) -> Tuple[str, List[str]]:
    combined_flag=""
    temp_file_list=[]
    for flag_with_underbar,value in kwargs.items():
        if isinstance(value, pd.DataFrame) or isinstance(value, pd.Series):
            formated_flag, temp_file = _process_df_to_file(flag_with_underbar=flag_with_underbar,
                                                           df=value,
                                                           tmpdir=tmpdir)
            temp_file_list.append(temp_file)
        else:
            formated_flag = "{} {} ".format( flag_with_underbar.replace("_","-"),
                                             value)
        combined_flag += formated_flag
    return combined_flag,temp_file_list

def _process_df_to_file(flag_with_underbar: str, df: pd.DataFrame, tmpdir: str) -> Tuple[str, str]:
    temp_path ="{}/memory_address_{}.{}".format(tmpdir.rstrip("/"),
                                                id(df),
                                                flag_with_underbar.replace("_",""))
    formated_flag = "{} {} ".format( flag_with_underbar.replace("_","-"),
                                     temp_path)
    df.to_csv(temp_path,sep="\t",index=None,header=None)
    return formated_flag, temp_path
    #value.to_csv(temp_path, sep="\t", index=None, header=None)

