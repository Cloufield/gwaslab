import pandas as pd

def _plink2_filter_to_flag(tmpdir="./",**args):
    combined_flag=""
    temp_file_list=[]
    for flag_with_underbar,value in args.items():
        if isinstance(value, pd.DataFrame) or isinstance(value, pd.Series):
            formated_flag, temp_file = _process_df_to_file(flag_with_underbar=flag_with_underbar,df=value,tmpdir=tmpdir)
            temp_file_list.append(temp_file)
        else:
            formated_flag = "{} {} ".format( flag_with_underbar.replace("_","-"), value)
        combined_flag += formated_flag
    return combined_flag,temp_file_list

def _process_df_to_file(flag_with_underbar, df, tmpdir):
    temp_path ="{}/memory_address_{}.{}".format(tmpdir.rstrip("/"),id(df),flag_with_underbar.replace("_",""))
    formated_flag = "{} {} ".format( flag_with_underbar.replace("_","-"), temp_path)
    df.to_csv(temp_path,sep="\t",index=None,header=None)
    return formated_flag, temp_path
    #value.to_csv(temp_path, sep="\t", index=None, header=None)
