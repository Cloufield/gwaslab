import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm_casting import _merge_mold_with_sumstats
from gwaslab.hm_casting import _align_with_mold
from gwaslab.hm_casting import _fill_missing_columns
from gwaslab.hm_casting import _check_daf
from gwaslab.hm_casting import _assign_warning_code
from gwaslab.qc_fix_sumstats import flipallelestats
from gwaslab.hm_casting import _renaming_cols

class SumstatsPair( ):
    def __init__(self, sumstatsObject1, sumstatsObject2, suffixes = ("_1","_2") ):
        
        if not isinstance(sumstatsObject1, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats Object #1.")
        if not isinstance(sumstatsObject2, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats Object #2.")

        self.snp_info_cols = []
        self.stats_cols =[]
        self.other_cols=[]
        self.log = Log()

        for i in sumstatsObject1.data.columns:
            if i in ["SNPID","rsID","CHR","POS","EA","NEA","EAF","STATUS"]:
                self.snp_info_cols.append(i)
            elif i in ["BETA","SE","P","MLOG10P","N","Z","OR","OR95L","OR95U","MAF"]:
                self.stats_cols.append(i)
            else:
                self.other_cols.append(i)

        self.data = sumstatsObject1.data.loc[:,self.snp_info_cols + self.stats_cols]

        self.data = self.data.rename(columns={"EA":"EA_1","NEA":"NEA_1","EAF":"EAF_1"})

        self.data = self.data.rename(columns={i:i + suffixes[0] for i in self.stats_cols})

        self._merge_two_sumstats(sumstatsObject2, suffixes=suffixes)
    
    def _merge_two_sumstats(self, sumstatsObject2, threshold=0.2, verbose=True,windowsizeb=10, ref_path=None,suffixes=("_1","_2")):

        molded_sumstats = _merge_mold_with_sumstats(self.data, 
                                                    sumstatsObject2.data, 
                                                    log=self.log,
                                                    verbose=verbose,
                                                    mold_suffix=suffixes[0])

        molded_sumstats = _align_with_mold(molded_sumstats, log=self.log, verbose=verbose,mold_suffix=suffixes[0])
        
        molded_sumstats = flipallelestats(molded_sumstats, log=self.log, verbose=verbose)

        molded_sumstats = _fill_missing_columns(molded_sumstats, self.stats_cols, log=self.log, verbose=verbose)

        molded_sumstats = _renaming_cols(molded_sumstats, self.stats_cols, log=self.log, verbose=verbose, suffix=suffixes[1])
        
        self.data = molded_sumstats