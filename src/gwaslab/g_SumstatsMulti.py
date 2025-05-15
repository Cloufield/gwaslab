import pandas as pd
import numpy as np
import copy
import gc
from gwaslab.util_in_filter_value import filtervalues
from gwaslab.g_Log import Log
from math import floor
from gwaslab.g_Sumstats import Sumstats
from gwaslab.g_Sumstats_polars import Sumstatsp
from gwaslab.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm_casting import _align_with_mold
from gwaslab.hm_casting import _fill_missing_columns
from gwaslab.hm_casting import _check_daf
from gwaslab.hm_casting import _assign_warning_code
from gwaslab.hm_casting import _renaming_cols
from gwaslab.hm_casting import _sort_pair_cols

from gwaslab.hm_casting_polars import _merge_mold_with_sumstats_by_chrposp
from gwaslab.hm_casting_polars import _align_with_moldp
from gwaslab.hm_casting_polars import _fill_missing_columnsp
from gwaslab.hm_casting_polars import _renaming_colsp
from gwaslab.hm_casting_polars import _sort_pair_colsp

from gwaslab.qc_fix_sumstats import flipallelestats
from gwaslab.qc_fix_sumstats_polars import flipallelestatsp
from gwaslab.qc_check_datatype_polars import check_datatype
from gwaslab.qc_check_datatype_polars import check_dataframe_shape
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
from gwaslab.util_ex_run_hyprcoloc import _run_hyprcoloc
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_fill_data import _get_multi_min
from gwaslab.g_meta import _init_meta
from gwaslab.g_meta_update import _update_meta
from gwaslab.qc_fix_sumstats import _process_build
from gwaslab.util_ex_run_mtag import _run_mtag



class SumstatsMulti( ):
    def __init__(self, 
                 sumstatsObjects, 
                 group_name=None, 
                 build="99",
                 engine="pandas",
                 merge_mode="outer",
                 merge_by_id=False,
                 verbose=True ):
        
        for i,sumstatsObject in enumerate(sumstatsObjects):
            if not isinstance(sumstatsObject, Sumstats):
                if not isinstance(sumstatsObject, Sumstatsp):
                    raise ValueError("Please provide GWASLab Sumstats Object #{}.".format(i+1))
        
        self.log = Log()
        self.meta = _init_meta(object="SumstatsMulti") 
        
        if engine=="polars":
            import polars as pl
            merge_mode="full"


        self.engine=engine
            
        self.meta["gwaslab"]["number_of_studies"] = len(sumstatsObjects)
        self.meta["gwaslab"]["genome_build"] = _process_build(build, log=self.log, verbose=False)
        self.meta["gwaslab"]["objects"] =  dict()
        self.meta["gwaslab"]["study_index"] =  dict()

        if group_name is None:
            self.group_name = "Group1" 
            self.meta["gwaslab"]["group_name"] =  "Group1" 
        else:
            self.group_name = group_name
            self.meta["gwaslab"]["group_name"] =  group_name

        self.names=[]
        self.hyprcoloc = {}

        self.snp_info_cols = dict()
        self.stats_cols =  dict()
        self.other_cols= dict()

        self.log.write( "Start to create SumstatsMulti object..." )

        for i,sumstatsObject in enumerate(sumstatsObjects):
            self.log.write( " -Checking sumstats Object #{}...".format(i+1), verbose=verbose)
            check_datatype(sumstatsObject.data, log=self.log, verbose=verbose)
            check_dataframe_shape(sumstats=sumstatsObject.data, 
                            log=self.log, 
                            verbose=verbose)
            
            if sumstatsObject.meta["gwaslab"]["study_name"] in self.names:
                new_study_name = "{}_{}".format(sumstatsObject.meta["gwaslab"]["study_name"],i+1)
                self.log.write( "  -Sumstats Object #{} name: {}".format(i+1,new_study_name), verbose=verbose)
                self.names.append(new_study_name)
            else:
                self.log.write( "  -Sumstats Object #{} name: {}".format(i+1, sumstatsObject.meta["gwaslab"]["study_name"]), verbose=verbose)
                self.names.append(sumstatsObject.meta["gwaslab"]["study_name"])
            self.meta["gwaslab"]["objects"][i+1] = sumstatsObject.meta
            self.meta["gwaslab"]["study_index"][i+1] = self.names[-1]

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

        self.meta["gwaslab"]["study_names_in_group"] = ",".join(self.names)
        


        self.log.write( " -Variant Info columns: {}".format(self.snp_info_cols[0]) , verbose=verbose)
        for i in range(len(sumstatsObjects)):
            self.log.write( " -Sumstats #{} variant statistics columns: {}".format(i+1, self.stats_cols[i]) , verbose=verbose)
            self.log.write( " -Sumstats #{} other columns: {}".format(i+1, self.other_cols[i]) , verbose=verbose)
        
        #for i,sumstatsObject in enumerate(sumstatsObjects):
        #    sumstatsObject.data["_RAW_INDEX_{}".format(i+1)] = range(len(sumstatsObject.data))

        # extract only info and stats cols
        self.data = sumstatsObjects[0].data
        
        #rename with _1
        if engine=="polars":
            self.data = self.data.rename({"EA":"EA_1","NEA":"NEA_1","STATUS":"STATUS_1"})
            self.data = self.data.rename({i:i + "_1" for i in self.stats_cols[0]})
            self.data = self.data.rename({i:i + "_1" for i in self.other_cols[0]})
        else:
            self.data = self.data.rename(columns={"EA":"EA_1","NEA":"NEA_1","STATUS":"STATUS_1"})
            self.data = self.data.rename(columns={i:i + "_1" for i in self.stats_cols[0]})
            self.data = self.data.rename(columns={i:i + "_1" for i in self.other_cols[0]})        

        if engine=="polars":
            self.data = pl.DataFrame(self.data)
            for i, sumstatsObject in enumerate(sumstatsObjects):
                if i >0:
                    self.log.write("Merging Sumstats #{} to main DataFrame...".format(i+1))
                    self.data = self._merge_two_sumstats(sumstatsObject.data,i=i,merge_mode=merge_mode,engine=engine,merge_by_id=merge_by_id)
                    self.log.write("Finished merging Sumstats #{} to main DataFrame.".format(i+1))
        else:
            for i, sumstatsObject in enumerate(sumstatsObjects):
                if i >0:
                    self.log.write("Merging Sumstats #{} to main DataFrame...".format(i+1))
                    self.data = self._merge_two_sumstats(sumstatsObject.data,i=i,merge_mode=merge_mode,engine=engine)
                    self.log.write("Finished merging Sumstats #{} to main DataFrame.".format(i+1))
        

    def _merge_two_sumstats(self, 
                            sumstatsObject2, 
                            verbose=True,
                            merge_mode="outer",
                            engine="pandas",
                            merge_by_id=False,
                            i=0):

        # _1 _2 
        # add suffix
        if engine=="polars":
            if "EA" in self.data.columns:
                self.data = self.data.rename({"EA":"EA_1","NEA":"NEA_1"})
        else:
            self.data = self.data.rename(columns={"EA":"EA_1","NEA":"NEA_1"})

        #sumstats1 with suffix _1, sumstats2 with no suffix
        if engine=="polars":
            molded_sumstats = _merge_mold_with_sumstats_by_chrposp(mold=self.data, 
                                                        sumstats=sumstatsObject2, 
                                                        log=self.log,
                                                        verbose=verbose,
                                                        merge_mode=merge_mode,
                                                        merge_by_id=merge_by_id,
                                                        stats_cols1 = self.other_cols[0],
                                                        stats_cols2 = self.other_cols[i],
                                                        suffixes=("_1",""),
                                                        return_not_matched_mold = False)
            molded_sumstats = _align_with_moldp(molded_sumstats, log=self.log, verbose=verbose,suffixes=("_1",""))
            molded_sumstats = flipallelestatsp(molded_sumstats, log=self.log, verbose=verbose)
            molded_sumstats = molded_sumstats.drop(["EA","NEA"] )
            molded_sumstats = molded_sumstats.rename({"EA_1":"EA","NEA_1":"NEA"})
        else:
            molded_sumstats = _merge_mold_with_sumstats_by_chrpos(mold=self.data, 
                                                        sumstats=sumstatsObject2, 
                                                        log=self.log,
                                                        verbose=verbose,
                                                        merge_mode=merge_mode,
                                                        stats_cols1 = self.other_cols[0],
                                                        stats_cols2 = self.other_cols[i],
                                                        suffixes=("_1",""),
                                                        return_not_matched_mold = False)    
            molded_sumstats = _align_with_mold(molded_sumstats, log=self.log, verbose=verbose,suffixes=("_1",""))
            molded_sumstats = flipallelestats(molded_sumstats, log=self.log, verbose=verbose)        
            molded_sumstats = molded_sumstats.drop(columns=["EA","NEA"] )
            molded_sumstats = molded_sumstats.rename(columns={"EA_1":"EA","NEA_1":"NEA"})
            
        if not set(self.stats_cols[i]) == set(self.stats_cols[0]):
            cols_to_fill = set(self.stats_cols[0]).difference(set(self.stats_cols[i]))
            molded_sumstats = _fill_missing_columns(molded_sumstats, cols_to_fill, log=self.log, verbose=verbose)
        
        if engine=="polars":
            # rename sumstast2 with _2
            molded_sumstats = _renaming_colsp(molded_sumstats, 
                                            self.stats_cols[0] + self.other_cols[i], 
                                            log=self.log, 
                                            verbose=verbose, 
                                            suffixes=("_1","_{}".format(i+1)))
        else:
            molded_sumstats = _renaming_cols(molded_sumstats, 
                                            self.stats_cols[0] + self.other_cols[i], 
                                            log=self.log, 
                                            verbose=verbose, 
                                            suffixes=("_1","_{}".format(i+1)))   
                     
        molded_sumstats = _sort_pair_cols(molded_sumstats, verbose=verbose, log=self.log, suffixes=["_{}".format(j) for j in range(1,i+2)])
        return molded_sumstats
    
    def update_meta(self,**kwargs):
        self.meta = _update_meta(self.meta, self.data, log = self.log, **kwargs)

    def run_meta_analysis(self, **kwargs):
        if self.engine == "polars":
            from gwaslab.util_in_meta_polars import meta_analyze_polars
            return meta_analyze_polars(self.data,nstudy = self.meta["gwaslab"]["number_of_studies"] ,**kwargs)
        else:
            return meta_analyze_multi(self.data,nstudy = self.meta["gwaslab"]["number_of_studies"] ,**kwargs)
    
    def run_hyprcoloc(self,**kwargs):
        hyprcoloc_res_combined = _run_hyprcoloc(self.data,
                       nstudy = self.meta["gwaslab"]["number_of_studies"], 
                       study= self.meta["gwaslab"]["group_name"], 
                       traits=self.names, **kwargs)
        self.hyprcoloc = hyprcoloc_res_combined

    def run_mtag(self,**kwargs):
        _run_mtag(     self.data,
                       nstudy = self.meta["gwaslab"]["number_of_studies"], 
                       study= self.meta["gwaslab"]["group_name"], 
                       traits=self.names, 
                       **kwargs)

    def get_lead(self, build=None, gls=False, **kwargs):
        
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        
        # extract build information from meta data
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]

        self.data = _get_multi_min(self.data,
                                   col="P", 
                                   nstudy=self.meta["gwaslab"]["number_of_studies"])

        output = getsig(self.data,
                            id=id_to_use,
                            chrom="CHR",
                            pos="POS",
                            p="P_MIN",
                            log=self.log,
                            build=build,
                            **kwargs)
        # return sumstats object    

        if gls == True:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = output
            gc.collect()
            return new_Sumstats_object
        
        return output
    
