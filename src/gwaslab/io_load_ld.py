
import scipy.sparse as sparse
import numpy as np
import pandas as pd
from gwaslab.hm_casting import _merge_mold_with_sumstats_by_chrpos

def _load_ld_matrix(path, 
                    fmt="npz", 
                    if_square=False, 
                    if_add_T=False):
    
    if fmt == "npz":
        r_matrix = sparse.load_npz(path).toarray()
        if if_add_T==True:
            r_matrix += r_matrix.T
        if if_square==True:
            r_matrix = np.power(r_matrix,2)
    return r_matrix
    
def _load_ld_map(path, 
                 snpid="SNPID", 
                 chrom="CHR", 
                 pos="POS", 
                 ref="EA", 
                 alt="NEA",
                 **args):
    
    ld_map = pd.read_csv(path,**args)
    ld_map = ld_map.rename(columns={alt:'EA', ref:'NEA', pos:'POS', chrom:'CHR', snpid:'SNPID'}, errors='ignore')

    return ld_map

def _merge_ld_map_with_sumstats(sumstats, 
                                ld_map,
                                r_matrix):
    merged_sumstats = _merge_mold_with_sumstats_by_chrpos(ld_map, sumstats, add_raw_index=True, suffixes=("_1","_2"))
    
    return merged_sumstats

def _extract_variants(merged_sumstats, r_matrix):
    avaiable_index= merged_sumstats["_INDEX_1"].values 
    reduced_r_matrix = r_matrix[np.ix_(avaiable_index, avaiable_index)]
    return reduced_r_matrix
