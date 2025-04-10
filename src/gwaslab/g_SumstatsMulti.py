import pandas as pd
import numpy as np
import copy
import gc
from gwaslab.util_in_filter_value import filtervalues
from gwaslab.g_Log import Log
from math import floor
from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm_casting import _align_with_mold
from gwaslab.hm_casting import _fill_missing_columns
from gwaslab.hm_casting import _check_daf
from gwaslab.hm_casting import _assign_warning_code
from gwaslab.qc_fix_sumstats import flipallelestats
from gwaslab.qc_check_datatype import check_datatype
from gwaslab.qc_check_datatype import check_dataframe_shape
from gwaslab.hm_casting import _renaming_cols
from gwaslab.hm_casting import _sort_pair_cols
from gwaslab.util_ex_calculate_ldmatrix import tofinemapping
from gwaslab.util_ex_run_coloc import _run_coloc_susie
from gwaslab.viz_plot_miamiplot2 import plot_miami2
from gwaslab.viz_plot_compare_af import  plotdaf
from gwaslab.util_ex_run_2samplemr import _run_two_sample_mr
from gwaslab.util_ex_run_clumping import _clump
from gwaslab.util_ex_ldproxyfinder import _extract_with_ld_proxy
from gwaslab.g_headers import _get_headers
from gwaslab.util_ex_match_ldmatrix import tofinemapping_m
from gwaslab.util_ex_run_mesusie import _run_mesusie
from gwaslab.util_in_meta import meta_analyze_multi

class SumstatsMulti( ):
    def __init__(self, 
                 sumstatsObjects, 
                 study=None, 
                 verbose=True ):
        
        for i,sumstatsObject in enumerate(sumstatsObjects):
            if not isinstance(sumstatsObject, Sumstats):
                raise ValueError("Please provide GWASLab Sumstats Object #{}.".format(i+1))
        
        self.nstudy = len(sumstatsObjects)
        self.study_name = "Group1" 
        
        self.snp_info_cols = dict()
        self.stats_cols =  dict()
        self.other_cols= dict()

        self.log = Log()

        self.log.write( "Start to create SumstatsMulti object..." )

        for i,sumstatsObject in enumerate(sumstatsObjects):
            self.log.write( " -Checking sumstats Object #{}...".format(i+1), verbose=verbose)
            check_datatype(sumstatsObject.data, log=self.log, verbose=verbose)
            check_dataframe_shape(sumstats=sumstatsObject.data, 
                            log=self.log, 
                            verbose=verbose)
            
            self.snp_info_cols[i] = list()
            self.stats_cols[i] = list()
            self.other_cols[i] = list()

            for col in sumstatsObject.data.columns:

                if col in _get_headers(mode="info"):
                    # extract SNP info columns from sumstats1
                    self.snp_info_cols[i].append(col)
                elif col in _get_headers(mode="stats"):
                    self.stats_cols[i].append(col)
                else:
                    self.other_cols[i].append(col)


        self.log.write( " -Variant Info columns: {}".format(self.snp_info_cols[0]) , verbose=verbose)
        for i in range(len(sumstatsObjects)):
            self.log.write( " -Sumstats #{} variant statistics columns: {}".format(i+1, self.stats_cols[i]) , verbose=verbose)
            self.log.write( " -Sumstats #{} other columns: {}".format(i+1, self.other_cols[i]) , verbose=verbose)
        
        #for i,sumstatsObject in enumerate(sumstatsObjects):
        #    sumstatsObject.data["_RAW_INDEX_{}".format(i+1)] = range(len(sumstatsObject.data))

        # extract only info and stats cols
        self.data = sumstatsObjects[0].data
        
        #rename with _1
        self.data = self.data.rename(columns={"EA":"EA_1","NEA":"NEA_1"})
        self.data = self.data.rename(columns={i:i + "_1" for i in self.stats_cols[0]})
        self.data = self.data.rename(columns={i:i + "_1" for i in self.other_cols[0]})

        for i, sumstatsObject in enumerate(sumstatsObjects):
            if i >0:
                self.data = self._merge_two_sumstats(sumstatsObject,i=i)

    def _merge_two_sumstats(self, 
                            sumstatsObject2, 
                            threshold=0.2, 
                            verbose=True,
                            windowsizeb=10, 
                            ref_path=None,
                            i=0,
                            suffixes=("_1","_2")):

        # _1 _2 
        # add suffix
        self.data = self.data.rename(columns={"EA":"EA_1","NEA":"NEA_1"})
        #sumstats1 with suffix _1, sumstats2 with no suffix
        molded_sumstats = _merge_mold_with_sumstats_by_chrpos(mold=self.data, 
                                                    sumstats=sumstatsObject2.data, 
                                                    log=self.log,
                                                    verbose=verbose,
                                                    merge_mode="outer",
                                                    stats_cols1 = self.other_cols[0],
                                                    stats_cols2 = self.other_cols[i],
                                                    suffixes=("_1",""),
                                                    return_not_matched_mold = False)


        molded_sumstats = _align_with_mold(molded_sumstats, log=self.log, verbose=verbose,suffixes=("_1",""))
        
        # flip sumstats2 statistics
        molded_sumstats = flipallelestats(molded_sumstats, log=self.log, verbose=verbose)
        
        # drop sumstats2 EA NEA
        molded_sumstats = molded_sumstats.drop(columns=["EA","NEA"])
        
        # rename sumstats1 EA NEA
        molded_sumstats = molded_sumstats.rename(columns={"EA_1":"EA","NEA_1":"NEA"})
        
        if not set(self.stats_cols[i]) == set(self.stats_cols[0]):
            cols_to_fill = set(self.stats_cols[0]).difference(set(self.stats_cols[i]))
            molded_sumstats = _fill_missing_columns(molded_sumstats, cols_to_fill, log=self.log, verbose=verbose)

        # rename sumstast2 with _2
        molded_sumstats = _renaming_cols(molded_sumstats, 
                                         self.stats_cols[0] + self.other_cols[i], 
                                         log=self.log, 
                                         verbose=verbose, 
                                         suffixes=("_1","_{}".format(i+1)))
        
        molded_sumstats = _sort_pair_cols(molded_sumstats, verbose=verbose, log=self.log, suffixes=["_{}".format(j) for j in range(1,i+2)])
        
        return molded_sumstats
    def run_meta_analysis(self,**kwargs):
        return meta_analyze_multi(self.data,nstudy = self.nstudy,**kwargs)