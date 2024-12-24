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

class SumstatsPair( ):
    def __init__(self, sumstatsObject1, sumstatsObject2, study=None, suffixes = ("_1","_2") ,verbose=True ):
        
        if not isinstance(sumstatsObject1, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats Object #1.")
        if not isinstance(sumstatsObject2, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats Object #2.")
        if sumstatsObject1.meta["gwaslab"]["study_name"]!=sumstatsObject2.meta["gwaslab"]["study_name"]:
            self.study_name = "{}_{}".format(sumstatsObject1.meta["gwaslab"]["study_name"], sumstatsObject2.meta["gwaslab"]["study_name"])
        else:
            self.study_name = "{}_{}".format("STUDY1", "STUDY2")
        
        self.snp_info_cols = []
        self.stats_cols =[]
        self.stats_cols2 =[]
        self.other_cols =[]
        self.other_cols2 =[]
        self.log = Log()
        self.suffixes = suffixes
        self.colocalization=pd.DataFrame()
        self.sumstats1 = pd.DataFrame()
        self.sumstats2 = pd.DataFrame()
        
        self.mr =dict()
        self.clumps =dict()
        self.ns = None
        self.finemapping = dict()
        #self.to_finemapping_file_path = ""
        #self.plink_log = ""

        self.log.write( "Start to create SumstatsPair object..." )
        
        self.log.write( " -Checking sumstats 1..." , verbose=verbose)
        check_datatype(sumstatsObject1.data, log=self.log, verbose=verbose)
        check_dataframe_shape(sumstats=sumstatsObject1.data, 
                        log=self.log, 
                        verbose=verbose)
        
        self.log.write( " -Checking sumstats 2..." , verbose=verbose)
        check_datatype(sumstatsObject2.data, log=self.log, verbose=verbose)
        check_dataframe_shape(sumstats=sumstatsObject2.data, 
                                log=self.log, 
                                verbose=verbose)

        for i in sumstatsObject1.data.columns:
            if i in _get_headers(mode="info"):
                # extract SNP info columns from sumstats1
                self.snp_info_cols.append(i)
            elif i in _get_headers(mode="stats"):
                self.stats_cols.append(i)
            else:
                self.other_cols.append(i)
        for i in sumstatsObject2.data.columns:
            if i in _get_headers(mode="info"):
                continue
            elif i in _get_headers(mode="stats"):
                self.stats_cols2.append(i)
            else:
                self.other_cols2.append(i)           
        
        self.log.write( " -Variant Info columns: {}".format(self.snp_info_cols) , verbose=verbose)
        self.log.write( " -Variant statistics columns: {}".format(self.stats_cols) , verbose=verbose)
        self.log.write( " -Sumstats1 other columns: {}".format(self.other_cols) , verbose=verbose)
        self.log.write( " -Sumstats2 other columns: {}".format(self.other_cols2) , verbose=verbose)
        
        # extract only info and stats cols
        self.data = sumstatsObject1.data
        
        #rename with _1
        self.data = self.data.rename(columns={"EA":"EA_1","NEA":"NEA_1"})
        self.data = self.data.rename(columns={i:i + suffixes[0] for i in self.stats_cols})
        self.data = self.data.rename(columns={i:i + suffixes[0] for i in self.other_cols})

        self.data, self.sumstats1 = self._merge_two_sumstats(sumstatsObject2, suffixes=suffixes)

        if "N{}".format(self.suffixes[0]) in self.data.columns and "N{}".format(self.suffixes[1]) in self.data.columns:
            n1 = int(floor(self.data["N{}".format(self.suffixes[0])].mean()))
            n2 = int(floor(self.data["N{}".format(self.suffixes[1])].mean()))
            self.ns=(n1, n2)
        else:
            self.ns = None

    def _merge_two_sumstats(self, sumstatsObject2, threshold=0.2, verbose=True,windowsizeb=10, ref_path=None,suffixes=("_1","_2")):

        # sumstats1 with suffix _1, sumstats2 with no suffix
        molded_sumstats, sumstats1 = _merge_mold_with_sumstats_by_chrpos(mold=self.data, 
                                                    sumstats=sumstatsObject2.data, 
                                                    log=self.log,
                                                    verbose=verbose,
                                                    suffixes=(suffixes[0],""),
                                                    return_not_matched_mold = True)

        molded_sumstats = _align_with_mold(molded_sumstats, log=self.log, verbose=verbose,suffixes=(suffixes[0],""))
        
        # flip sumstats2 statistics
        molded_sumstats = flipallelestats(molded_sumstats, log=self.log, verbose=verbose)
        
        # drop sumstats2 EA NEA
        molded_sumstats = molded_sumstats.drop(columns=["EA","NEA"])
        
        # rename sumstats1 EA NEA
        molded_sumstats = molded_sumstats.rename(columns={"EA_1":"EA","NEA_1":"NEA"})
        
        if not set(self.stats_cols2) == set(self.stats_cols):
            cols_to_fill = set(self.stats_cols).difference(set(self.stats_cols2))
            molded_sumstats = _fill_missing_columns(molded_sumstats, cols_to_fill, log=self.log, verbose=verbose)

        # rename sumstast2 with _2
        molded_sumstats = _renaming_cols(molded_sumstats, self.stats_cols + self.other_cols2, log=self.log, verbose=verbose, suffixes=suffixes)
        
        molded_sumstats = _sort_pair_cols(molded_sumstats, verbose=verbose, log=self.log)
        
        return molded_sumstats, sumstats1


    def clump(self,**kwargs):
        self.clumps["clumps"],self.clumps["clumps_raw"],self.clumps["plink_log"] = _clump(self.data, log=self.log, p="P_1",mlog10p="MLOG10P_1", study = self.study_name, **kwargs)

    def to_coloc(self,**kwargs):
        self.finemapping["path"],self.finemapping["file"],self.finemapping["plink_log"] = tofinemapping(self.data,study=self.study_name,suffixes=self.suffixes,log=self.log,**kwargs)

    def run_coloc_susie(self,**kwargs):
        self.colocalization = _run_coloc_susie(self.finemapping["path"],log=self.log,ncols=self.ns,**kwargs)

    def run_two_sample_mr(self, clump=False, **kwargs):
        exposure1 = self.study_name.split("_")[0]
        outcome2 = self.study_name.split("_")[1]
        _run_two_sample_mr(self,exposure1=exposure1,outcome2=outcome2, clump=clump,**kwargs)

    def extract_with_ld_proxy(self,**arg):
        return _extract_with_ld_proxy(common_sumstats = self.data, sumstats1=self.sumstats1,  **arg)

    def filter_value(self, expr, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filtervalues(new_Sumstats_object.data,expr,log=new_Sumstats_object.log, **kwargs)
            return new_Sumstats_object
        else:
            self.data = filtervalues(self.data, expr,log=self.log,**kwargs)
        gc.collect()

    ## Visualization #############################################################################################################################################
    def plot_miami(self,**kwargs):

        plot_miami2(merged_sumstats=self.data, 
                    suffixes=self.suffixes,
                    **kwargs)
    
    def compare_af(self, **kwargs):
        
        return plotdaf( self.data,
                     eaf="EAF_2",
                     raf="EAF_1",
                     xlabel="Effect Allele Frequency in Sumstats 1",
                     ylabel="Effect Allele Frequency in Sumstats 2",
                     **kwargs)
                     