from typing import TYPE_CHECKING, Union, List, Dict, Set
import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
import re

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

@with_logging(
        start_to_msg="initialize gl.SumstatsSet",
        finished_msg="initializing gl.SumstatsSet"
)
def _extract_variant(variant_set: Union[List[str], List[List[Union[int, float]]], Set[str]], sumstats_dic: Dict[str, 'Sumstats'], log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    
    combined = pd.DataFrame()
    for key, sumstats_gls in sumstats_dic.items():
        log.write(" -{} : {}".format(key, sumstats_gls), verbose=verbose)

    for key, sumstats_gls in sumstats_dic.items():
        
        sumstats_single = sumstats_gls.data
        
        # create a boolean col with FALSE 
        is_extract = sumstats_single["SNPID"]!=sumstats_single["SNPID"]
        
        for variant in variant_set:
            
            if pd.api.types.is_list_like(variant):

                chrom=variant[0]
                pos=variant[1]

                is_extract = is_extract | ((sumstats_single["POS"] == pos ) &(sumstats_single["CHR"] == chrom))
            elif pd.api.types.is_string_dtype(type(variant)):
                
                is_extract = is_extract | (sumstats_single["SNPID"] == variant)

                a= re.search(r'^(chr|Chr|CHR)?(\d+)[:_-](\d+)[:_-][ATCG]+[:_-][ATCG]+$', variant, flags=0)
                if a is not None:
                    chrom=int(a[2])
                    pos=int(a[3])
                    is_extract = is_extract | ((sumstats_single["POS"] == pos ) &(sumstats_single["CHR"] == chrom))

        to_extract =  sumstats_single.loc[is_extract,:].copy()
        log.write(" -Extracted {} variants from {}".format(len(to_extract), key),verbose=verbose)
        to_extract["STUDY"] = key
        
        to_extract_cols=["STUDY"]

        default_cols=["SNPID","EA","NEA","CHR","POS","BETA","SE","P","MLOG10P","EAF","MAF","STATUS"]
        
        for i in default_cols:
            if i in sumstats_single.columns:
                to_extract_cols.append(i)

        combined = pd.concat([combined, to_extract[to_extract_cols]], ignore_index=True)
    return combined
