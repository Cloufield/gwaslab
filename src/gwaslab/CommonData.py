import pandas as pd
from os import path

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
            "19","20","21","22","23","24","25"]
    return chrom_list

def get_chr_to_number(out_chr=False):
    if out_chr:
        dic= {str(i):str(i) for i in range(1,23)}
        dic["X"]="23"
        dic["Y"]="24"
        dic["MT"]="25"
    
    else:
        dic= {str(i):i for i in range(1,23)}
        dic["X"]=23
        dic["Y"]=24
        dic["MT"]=25
    return dic

def get_number_to_chr(in_chr=False):
    if in_chr:
        dic= {str(i):str(i) for i in range(1,23)}
        dic["23"]="X"
        dic["24"]="Y"
        dic["25"]="MT"
    else:
        dic= {i:str(i) for i in range(1,23)}
        dic[23]="X"
        dic[24]="Y"
        dic[25]="MT"
    return dic

def lookup_status():
    status_dic={
    "19xxx":"",
    "19xxx":"",
    "19xxx":"",
    "19xxx":"",
    "19xxx":"",
    "19xxx":"",
    "19xxx":"",
    "19xxx":""
    }
    
    
def get_high_ld(build="19"):
    if build=="19":
        data_path =  path.dirname(__file__) + '/data/high-ld/high-ld-hla-hg19.bed'
    elif build=="38":
        data_path =  path.dirname(__file__) + '/data/high-ld/high-ld-hla-hg38.bed'
    return data_path

def get_format_dict(fmt,inverse=False):
    if fmt.lower()=="plink":
         #https://www.cog-genomics.org/plink/1.9/formats#assoc_linear
         #INFO
        dic={
        "SNP":"SNPID",
        "CHR":"CHR",
        "BP":"POS",
        "A1":"EA",
        "A2":"NEA",
        "NMISS":"N",
        "FRQ":"EAF",
        "BETA":"BETA",
        "SE":"SE",
        "P":"P",
        "STAT":"Z",
        "INFO":"INFO",
        "OR":"OR",
        "L95":"OR_95L",
        "U95":"OR_95U"
        }
    if fmt.lower()=="plink2":
        #https://www.cog-genomics.org/plink/2.0/formats#glm_logistic
        #https://www.cog-genomics.org/plink/2.0/formats#glm_linear
         dic={
        "ID":"SNPID",
        "CHROM":"CHR",
        "POS":"POS",
        "REF":"NEA",
        "ALT":"EA",
        "OBS_CT":"N",
        "BETA":"BETA",
        "SE":"SE",
        "P":"P",
        "[LOG(OR)_]SE":"OR_SE",
        "Z_[OR_F_]STAT":"Z",
        "MACH_R2":"INFO",
        "OR":"OR",
        "L95":"OR_95L",
        "U95":"OR_95U"
        }
    if fmt.lower()=="ssf":
        #https://www.biorxiv.org/content/10.1101/2022.07.15.500230v1.full
        dic={
        "variant_id":"SNPID",
        "rsid":"rsID",
        "chromosome":"CHR",
        "bas_pair_location":"POS",
        "other_allele":"NEA",
        "effect_allele":"EA",
        "n":"N",
        "beta":"BETA",
        "standard_error":"SE",
        "p_value":"P",
        "info":"INFO",
        "odds_ratio":"OR",
        "ci_lower":"OR_95L",
        "ci_upper":"OR_95U"
        }
    if fmt.lower()=="saige":
        dic={}
    if fmt.lower()=="regenie":
        dic={}
    if fmt.lower()=="fastgwa":
        #https://yanglab.westlake.edu.cn/software/gcta/#fastGWA
        # CHR     SNP     POS     A1      A2      N      AF1     BETA    SE      P
        #https://yanglab.westlake.edu.cn/software/gcta/#fastGWA-GLMM
        # CHR    SNP    POS    A1    A2    N    AF1    T    SE_T    P_noSPA    BETA    SE    P    CONVERGE
        dic={
        "SNP":"SNPID",
        "CHR":"CHR",
        "POS":"POS",
        "A2":"NEA",
        "A1":"EA",
        "AF1":"EAF",
        "BETA":"BETA",
        "SE":"SE",
        "P":"P"
        }
    if fmt.lower()=="gcta":
        dic={}
    if fmt.lower()=="metal":
        dic={
        "MarkerName":"SNPID",
        "Allele1":"NEA",
        "Allele2":"EA",
        "Freq1":"EAF",
        "Effect":"BETA",
        "StdErr":"SE",
        "P-value":"P",
        "Direction":"DIRECTION"
        }
    if fmt.lower()=="mrmega":
        dic={}
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    