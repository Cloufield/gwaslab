import pandas as pd
from os import path
import json

def get_chr_NC_dict(build,inverse=False):
    #https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
    if build =="19":
        dic={
        "1":"NC_000001.10",
        "2":"NC_000002.11",
        "3":"NC_000003.11",
        "4":"NC_000004.11",
        "5":"NC_000005.9",
        "6":"NC_000006.11",
        "7":"NC_000007.13",
        "8":"NC_000008.10",
        "9":"NC_000009.11",
        "10":"NC_000010.10",
        "11":"NC_000011.9",
        "12":"NC_000012.11",
        "13":"NC_000013.10",
        "14":"NC_000014.8",
        "15":"NC_000015.9",
        "16":"NC_000016.9",
        "17":"NC_000017.10",
        "18":"NC_000018.9",
        "19":"NC_000019.9",
        "20":"NC_000020.10",
        "21":"NC_000021.8",
        "22":"NC_000022.10",
        "X":"NC_000023.10",
        "Y":"NC_000024.9"}
    elif build=="39":
        dic={
        "1":"NC_000001.11",
        "2":"NC_000002.12",
        "3":"NC_000003.12",
        "4":"NC_000004.12",
        "5":"NC_000005.10",
        "6":"NC_000006.12",
        "7":"NC_000007.14",
        "8":"NC_000008.11",
        "9":"NC_000009.12",
        "10":"NC_000010.11",
        "11":"NC_000011.10",
        "12":"NC_000012.12",
        "13":"NC_000013.11",
        "14":"NC_000014.9",
        "15":"NC_000015.10",
        "16":"NC_000016.10",
        "17":"NC_000017.11",
        "18":"NC_000018.10",
        "19":"NC_000019.10",
        "20":"NC_000020.11",
        "21":"NC_000021.9",
        "22":"NC_000022.11",
        "X":"NC_000023.11",
        "Y":"NC_000024.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic

def get_chr_list():
    chrom_list=["1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18",
            "19","20","21","22","23","24","25","X","Y","M","MT"]
    return chrom_list

def get_chr_to_number(out_chr=False,xymt=["X","Y","MT"]):
    if out_chr is True:
        dic= {str(i):str(i) for i in range(1,23)}
        dic[xymt[0]]="23"
        dic[xymt[1]]="24"
        dic[xymt[2]]="25"
    
    else:
        dic= {str(i):i for i in range(1,23)}
        dic[xymt[0]]=23
        dic[xymt[1]]=24
        dic[xymt[2]]=25
    return dic

def get_number_to_chr(in_chr=False,xymt=["X","Y","MT"]):
    if in_chr is True:
        dic= {str(i):str(i) for i in range(1,23)}
        dic["23"]=xymt[0]
        dic["24"]=xymt[1]
        dic["25"]=xymt[2]
    else:
        dic= {i:str(i) for i in range(1,23)}
        dic[23]=xymt[0]
        dic[24]=xymt[1]
        dic[25]=xymt[2]
    return dic


    
###################################################################################################################    
def get_high_ld(build="19"):
    if build=="19":
        data_path =  path.dirname(__file__) + '/data/high_ld/high_ld_hla_hg19.bed'
    elif build=="38":
        data_path =  path.dirname(__file__) + '/data/high_ld/high_ld_hla_hg38.bed'
    return data_path

def get_format_dict(fmt,inverse=False):
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    dicts = json.load(open(data_path))
    dic_meta = dicts[fmt]["meta_data"]
    dic_dict = dicts[fmt]["format_dict"]
    if inverse is True:
        inv_dic = {v: k for k, v in dic_dict.items()}
        return dic_meta,inv_dic
    return dic_meta, dic_dict

####################################################################################################################
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    