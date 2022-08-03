import pandas as pd

def summarize(sumstats,
              snpid="SNPID",
              rsid="rsID",
              eaf="EAF",
              p="P",
              n="N",
              status="STATUS",
             ):
    #print("Start summarizing the current sumstats...")
    output = {}
    meta_dic={}
    meta_dic["Row_num"] = str(sumstats.shape[0])
    meta_dic["Column_num"] = str(sumstats.shape[1])
    meta_dic["Column_names"] = ",".join(list(sumstats.columns.values))
    output["META"]=meta_dic
    ###############################################################################
    ##check missing
    missing_dict={}
    missing =  sumstats.isna()
    if missing.any(axis=None) is True:
        missing_dict["Missing_total"]=sum(missing.any(axis=0).any())
        for i in missing.columns:
            if sum(missing[i])>0:
                missing_dict["Missing_"+i]=sum(missing[i])
    else:
        missing_dict["Missing_total"]=0
    output["MISSING"]=missing_dict
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
    ##############################################################################
    ##check p values
    if p in sumstats.columns:
        p_dic={}
        p_dic["Minimum"] = sumstats[p].min()
        p_dic["Significant"] = sum(sumstats[p]<5e-8)
        p_dic["Suggestive"] = sum(sumstats[p]<5e-6)
        output["P"]=p_dic
    ##############################################################################
    ##check status
    id_to_use = None
    if rsid in sumstats.columns:
        id_to_use=rsid
    elif snpid in sumstats.columns:
        id_to_use=snpid
    if status in sumstats.columns and id_to_use in sumstats.columns:
        status_summary = sum_status(id_to_use,sumstats[[id_to_use,status]])
    if status in sumstats:
        status_dic = {}
        for index,row in status_summary.iterrows():
            status_dic[str(index)]=row[0]
        output["STATUS"]=status_dic
    
    df = pd.DataFrame.from_dict({(i,j): output[i][j] 
                           for i in output.keys() 
                           for j in output[i].keys()},
                           orient='index',columns=["Values"],dtype="string")
    index = pd.MultiIndex.from_tuples(df.index.values, names=["Category", "Items"])
    df = pd.DataFrame(df,index=index)
    #df = pd.DataFrame.from_dict(output,orient="index",columns=["Values"],dtype="string").sort_index()
    ###############################################################################
    return df 

def sum_status(id_to_use, sumstats):
        results = sumstats.groupby("STATUS").count()
        results = results.loc[results[id_to_use]>0,:].sort_values(id_to_use,ascending=False)
        return results