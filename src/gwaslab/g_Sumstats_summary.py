import pandas as pd
import time
import numpy as np
def summarize(insumstats,
              snpid="SNPID",
              rsid="rsID",
              eaf="EAF",
              p="P",
              n="N",
              status="STATUS",
             ):
    #print("Start summarizing the current sumstats...")
    ###############################################################################
    cols=[]
    for i in [snpid,rsid,eaf,p,n,status]:
        if i in insumstats.columns:
            cols.append(i)
    sumstats= insumstats.loc[:,cols].copy()
    ###############################################################################
    numeric_cols=[]
    output = {}
    meta_dic={}
    meta_dic["Row_num"] = str(insumstats.shape[0])
    meta_dic["Column_num"] = str(insumstats.shape[1])
    meta_dic["Column_names"] = ",".join(list(insumstats.columns.values))
    meta_dic["Last_checked_time"]=str(time.ctime(time.time()))
    output["META"]=meta_dic
    ###############################################################################
    ##check missing
    missing_dict={}
    missing = sumstats.isna()
    if missing.any(axis=None):
        missing_dict["Missing_total"]=sum(missing.any(axis=1))
        for i in missing.columns:
            if sum(missing[i])>0:
                missing_dict["Missing_"+i]=sum(missing[i])
    else:
        missing_dict["Missing_total"]=0
    output["MISSING"]=missing_dict
    numeric_cols.append("MISSING")
    ###############################################################################
    ##check maf 
    maf_dic={}
    if eaf in sumstats.columns:
        maf = pd.to_numeric(sumstats[eaf], errors='coerce')
        maf[maf>0.5] = 1 - maf[maf>0.5]
        maf_dic["Common"] =  sum(maf>0.05)
        maf_dic["Low_frequency"] = sum((maf>0.01 )& (maf<0.05))
        maf_dic["Rare"] = sum(maf<0.01)
        output["MAF"]=maf_dic
        numeric_cols.append("MAF")
    ##############################################################################
    ##check p values
    if p in sumstats.columns:
        p_dic={}
        p_dic["Minimum"] = sumstats[p].min()
        p_dic["Significant"] = sum(sumstats[p]<5e-8)
        p_dic["Suggestive"] = sum(sumstats[p]<5e-6)
        output["P"]=p_dic
        numeric_cols.append("P")
    ##############################################################################
    ##check status

    if status in sumstats.columns:
        id_to_use = "uniq_index"
        sumstats[id_to_use] = range(len(sumstats))
        status_summary = sum_status(id_to_use,sumstats[[id_to_use,status]])
        sumstats.drop(columns='uniq_index',inplace=True)
        status_dic = {}
        for index,row in status_summary.iterrows():
            status_dic[str(index)]=row[0]
        output["STATUS"]=status_dic
        numeric_cols.append("STATUS")
    df = pd.DataFrame.from_dict({(i,j): output[i][j] 
                           for i in output.keys() 
                           for j in output[i].keys()},
                           orient='index',columns=["Values"],dtype="string")
    index = pd.MultiIndex.from_tuples(df.index.values, names=["Category", "Items"])
    df = pd.DataFrame(df,index=index)
    if len(numeric_cols)>0:
        df.loc[numeric_cols,"Percentage"] = np.around(pd.to_numeric(df.loc[numeric_cols,"Values"],errors='coerce')*100 / int(sumstats.shape[0]),2)
    #df = pd.DataFrame.from_dict(output,orient="index",columns=["Values"],dtype="string").sort_index()
    ###############################################################################
    return df 

def sum_status(id_to_use, sumstats):
        results = sumstats.groupby("STATUS").count()
        results = results.loc[results[id_to_use]>0,:].sort_values(id_to_use,ascending=False)
        return results
    
def lookupstatus(status):
    uniq_status = status.unique()
    status_dic_12={
    "13":"CHM13",
    "19":"hg19",
    "38":"hg38",
    "97":"Unmapped",
    "98":"Unknown",
    "99":"Unchecked"
        }  
    status_dic_3={
    "0":"rsid valid & SNPID valid",
    "1":"rsid valid & SNPID invalid",
    "2":"rsid invalid & SNPID valid",
    "3":"rsid invalid & SNPID invalid",
    "4":"rsid valid & SNPID valid",
    "5":"rsid valid & SNPID unknown",
    "6":"rsid unknown & SNPID valid",
    "7":"rsid invalid & SNPID unknown",
    "8":"rsid unknown & SNPID invalid",
    "9":"Unchecked"
        }  
    status_dic_4={
    "0":"CHR valid & POS valid",
    "2":"CHR invalid & POS invalid",
    "3":"CHR invalid & POS valid",
    "4":"CHR valid & POS invalid",
    "5":"CHR valid & POS unknown",
    "6":"CHR unknown & POS valid",
    "7":"CHR invalid & POS unknown",
    "8":"CHR unknown & POS invalid",
    "9":"Unchecked"
        }  
    status_dic_5={
    "0":"standardized SNP",
    "1":"standardized & normalized insertion",
    "2":"standardized & normalized deletion",
    "3":"standardized & normalized indel",
    "4":"standardized indel",
    "5":"indistinguishable or not normalized allele",
    "6":"invalid allele notation",
    "7":"Unknown",
    "9":"Unchecked"
        }  
    status_dic_6={
    "0":"Match: NEA=REF",
    "1":"Flipped_fixed",
    "2":"Reverse_complementary_fixed",
    "3":"Flipped",
    "4":"Reverse_complementary",
    "5":"Reverse_complementary+Flipped",
    "6":"Both_alleles_on_ref+indistinguishable",
    "8":"Not_on_reference_genome",
    "9":"Unchecked"
        }  
    status_dic_7={
    "0":"Not_palindromic_SNPs",
    "1":"Palindromic+strand",
    "2":"Palindromic-strand_fixed",
    "3":"Indel_match",
    "4":"Indel_flipped_fixed",
    "5":"Palindromic-strand",
    "6":"Indel_flipped",
    "7":"Indistinguishable",
    "8":"No_matching_or_no_info",
    "9":"Unchecked"
        }   
    status_dic={}
    for i in uniq_status:
        count = sum(status==i)
        percentage =  np.around(count*100 / len(status),2)
        description = [status_dic_12[i[:2]] , status_dic_3[i[2]] , status_dic_4[i[3]], status_dic_5[i[4]], status_dic_6[i[5]], status_dic_7[i[6]],count,percentage]
        status_dic[i] = description
    df = pd.DataFrame.from_dict({i: status_dic[i] 
                           for i in status_dic.keys()},
                           orient='index',columns=["Genome_Build","rsID&SNPID","CHR&POS","Stadardize&Normalize","Align","Panlidromic_SNP&Indel","Count","Percentage(%)"], dtype="string")
    df = df.sort_index()
    return df
        
