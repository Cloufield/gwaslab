import numpy as np
from gwaslab.util_in_filter_value import inferbuild
from gwaslab.g_Log import Log
import time

def _update_meta(meta, sumstats, object="Sumstats",log=Log(), verbose=True):  
    
    meta["gwaslab"]["variants"]["variant_number"] = len(sumstats)
    
    if "CHR" in sumstats.columns:
        meta["gwaslab"]["variants"]["number_of_chromosomes"] = len(sumstats["CHR"].unique())
    
    if  meta["gwaslab"]["gwaslab_object"]=="gwaslab.Sumstats":      
        if "P" in sumstats.columns:
            meta["gwaslab"]["variants"]["min_P"]=np.nanmin(sumstats["P"])
        if "EAF" in sumstats.columns:
            meta["gwaslab"]["variants"]["min_minor_allele_freq"]=min (np.min(sumstats["EAF"]) , 1- np.max(sumstats["EAF"]))
        if "N" in sumstats.columns:
            meta["gwaslab"]["samples"]["sample_size"] = int(sumstats["N"].max())
            meta["gwaslab"]["samples"]["sample_size_median"] = sumstats["N"].median()
            meta["gwaslab"]["samples"]["sample_size_min"] = int(sumstats["N"].min())
    
    if  meta["gwaslab"]["gwaslab_object"]=="gwaslab.SumstatsMulti" or meta["gwaslab"]["gwaslab_object"]=="gwaslab.SumstatsPair":   
        nstudy = meta["gwaslab"]['number_of_studies']
        for i in range(nstudy):
            i_form_1 = i + 1
            meta["gwaslab"]["variants"][i_form_1]=dict()
            meta["gwaslab"]["samples"][i_form_1] =dict()

            if "P_{}".format(i_form_1) in sumstats.columns:
                p = "P_{}".format(i_form_1)
                
                meta["gwaslab"]["variants"][i_form_1]["min_P"]= np.nanmin(sumstats[p])
            if "N_{}".format(i_form_1) in sumstats.columns:
                n = "N_{}".format(i_form_1)
                meta["gwaslab"]["samples"][i_form_1]["sample_size"] = int(sumstats[n].max())
                meta["gwaslab"]["samples"][i_form_1]["sample_size_median"] = sumstats[n].median()
                meta["gwaslab"]["samples"][i_form_1]["sample_size_min"] = int(sumstats[n].min())
            if "EAF_{}".format(i_form_1) in sumstats.columns:
                eaf="EAF_{}".format(i_form_1)
                meta["gwaslab"]["variants"][i_form_1]["min_minor_allele_freq"]=min (np.min(sumstats[eaf]) , 1- np.max(sumstats[eaf]))

    if meta["gwaslab"]["genome_build"] == "99":
        _, meta["gwaslab"]["genome_build"] = inferbuild(sumstats, change_status=False, log=log, verbose=verbose)
    
    meta["date_last_modified"] = str(time.strftime('%Y/%m/%d'))
    
    return meta