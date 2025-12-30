import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm.hm_casting import _align_with_mold
from gwaslab.hm.hm_casting import _fill_missing_columns
from gwaslab.hm.hm_casting import _check_daf
from gwaslab.hm.hm_casting import _assign_warning_code

from gwaslab.qc.qc_fix_sumstats import _flip_allele_stats


class SumstatsT( ):
    def __init__(self, sumstatsObject: 'Sumstats') -> None:
        
        if not isinstance(sumstatsObject, Sumstats):
            raise ValueError("Please provide a GWASLab Sumstats Object.")

        self.snp_info_cols = []
        self.stats_cols =[]
        self.other_cols=[]

        for i in sumstatsObject.data.columns:
            if i in ["SNPID","rsID","CHR","POS","EA","NEA","EAF","STATUS"]:
                self.snp_info_cols.append(i)
            elif i in ["BETA","SE","P","MLOG10P","N","Z","OR","OR95L","OR95U","MAF"]:
                self.stats_cols.append(i)
            else:
                self.other_cols.append(i)

        self.snp_info = sumstatsObject.data.loc[:,self.snp_info_cols]

        self.snp_info = self.snp_info.rename(columns={"EA":"EA_MOLD","NEA":"NEA_MOLD","EAF":"EAF_MOLD"})
    
    def cast(
        self, 
        sumstatsObject: 'Sumstats', 
        threshold: float = 0.2, 
        verbose: bool = True,
        windowsizeb: int = 10, 
        ref_path: Optional[str] = None
    ) -> None:

        molded_sumstats = _merge_mold_with_sumstats_by_chrpos(self.snp_info, sumstats_or_dataframe=sumstatsObject.data, log=sumstatsObject.log, verbose=verbose, windowsizeb=windowsizeb,ref_path=ref_path)

        molded_sumstats = _align_with_mold(molded_sumstats, log=sumstatsObject.log, verbose=verbose)
        
        molded_sumstats = _flip_allele_stats(molded_sumstats, log=sumstatsObject.log, verbose=verbose)

        molded_sumstats = _fill_missing_columns(molded_sumstats, self.stats_cols, log=sumstatsObject.log, verbose=verbose)

        
        molded_sumstats = _check_daf(molded_sumstats, log=sumstatsObject.log, verbose=verbose)

        molded_sumstats = _assign_warning_code(molded_sumstats, threshold=threshold, log=sumstatsObject.log, verbose=verbose)

        sumstatsObject.data = molded_sumstats

            
        
