import pandas as pd
import numpy as np
from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm_casting import _merge_mold_with_sumstats
from gwaslab.hm_casting import _align_with_mold

class SumstatsT( ):
    def __init__(self, sumstatsObject):
        
        print(isinstance(sumstatsObject, Sumstats))

        self.snp_info_cols = []
        self.stats_cols =[]
        self.other_cols=[]

        for i in sumstatsObject.data.columns:
            if i in ["SNPID","rsID","CHR","POS","EA","NEA","STATUS"]:
                self.snp_info_cols.append(i)
            elif i in ["BETA","SE","P","N","EAF","Z","T","OR","OR95L","OR95U"]:
                self.stats_cols.append(i)
            else:
                self.other_cols.append(i)

        self.snp_info = sumstatsObject.data.loc[:,self.snp_info_cols]

        self.snp_info = self.snp_info.rename(columns={"EA":"EA_MOLD","NEA":"NEA_MOLD"})
    
    def cast(self, sumstatsObject):

        molded_sumstats = _merge_mold_with_sumstats(self.snp_info, sumstatsObject.data)

        return molded_sumstats

            
        
