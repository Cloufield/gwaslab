import pandas as pd
import numpy as np
import copy
import gc
from math import floor
from typing import TYPE_CHECKING, Optional, Dict, Any, List, Tuple, Union, Callable

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
    from gwaslab.info.g_Log import Log

from gwaslab.info.g_meta import _init_meta
from gwaslab.qc.qc_reserved_headers import _get_headers
from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats

from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm.hm_casting import _align_with_mold
from gwaslab.hm.hm_casting import _fill_missing_columns
from gwaslab.hm.hm_casting import _renaming_cols
from gwaslab.hm.hm_casting import _sort_pair_cols

from gwaslab.qc.qc_fix_sumstats import _flip_allele_stats
from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_shape

from gwaslab.util.util_ex_calculate_ldmatrix import _to_finemapping
from gwaslab.util.rwrapper.util_ex_run_coloc import _run_coloc_susie
from gwaslab.util.rwrapper.util_ex_run_ccgwas import _run_ccgwas
from gwaslab.util.rwrapper.util_ex_run_2samplemr import _run_two_sample_mr
from gwaslab.util.util_ex_run_clumping import _clump
from gwaslab.util.util_ex_ldproxyfinder import _extract_with_ld_proxy
from gwaslab.util.util_ex_match_ldmatrix import tofinemapping_m
from gwaslab.util.rwrapper.util_ex_run_mesusie2 import _run_mesusie
from gwaslab.util.util_ex_run_multisusie import _run_multisusie_rss
from gwaslab.util.util_in_filter_value import _filter_values
from gwaslab.util.util_in_meta import meta_analyze_multi

from gwaslab.viz.viz_plot_stackedregional import plot_stacked_mqq
from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2
from gwaslab.viz.viz_plot_compare_af import  plotdaf
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

from gwaslab.io.io_to_pickle import _offload
from gwaslab.io.io_to_pickle import _reload
from gwaslab.io.io_to_pickle import dump_pickle_pair
from gwaslab.io.io_read_pipcs import _read_pipcs

class SumstatsPair( ):
    def __init__(
        self, 
        sumstatsObject1: 'Sumstats', 
        sumstatsObject2: 'Sumstats', 
        study: Optional[str] = None, 
        suffixes: Tuple[str, str] = ("_1","_2"), 
        keep_all_variants: bool = True,
        verbose: bool = True
    ) -> None:
        
        if not isinstance(sumstatsObject1, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats Object #1.")
        if not isinstance(sumstatsObject2, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats Object #2.")
        
        # Initialize metadata and basic attributes
        self.meta = _init_meta(object="SumstatsPair")
        self.id = id(self)
        self.tmp_path = f"./{self.id}"
        self.suffixes = suffixes
        self.log = Log()
        self.viz_params = VizParamsManager()
        load_viz_config(self.viz_params)

        # Set study names
        study_name1 = sumstatsObject1.meta["gwaslab"]["study_name"]
        study_name2 = sumstatsObject2.meta["gwaslab"]["study_name"]
        
        if study_name1 != study_name2:
            self.study_name = f"{study_name1}_{study_name2}"
            self.study_names = [study_name1, study_name2]
        else:
            self.study_name = f"{study_name1}1_{study_name2}2"
            self.study_names = [f"{study_name1}1", f"{study_name2}2"]
        
        self.meta["gwaslab"]["objects"] = {0: sumstatsObject1.meta, 1: sumstatsObject2.meta}
        self.meta["gwaslab"]["group_name"] = self.study_name
        
        # Initialize LDSC data
        self.ldsc = {0: sumstatsObject1.ldsc_h2, 1: sumstatsObject2.ldsc_h2}
        self.ldsc_rg = sumstatsObject1.ldsc_rg

        # Initialize column lists
        self.snp_info_cols = []
        self.stats_cols = []
        self.stats_cols2 = []
        self.other_cols = []
        self.other_cols2 = []
        
        # Initialize result containers
        self.colocalization = pd.DataFrame()
        self.sumstats1 = pd.DataFrame()
        self.sumstats2 = pd.DataFrame()
        self.ns = None
        self.mr = {}
        self.clumps = {}
        self.mesusie = {}
        self.mesusie_res = pd.DataFrame()
        self.multisusie = {}
        self.multisusie_res = pd.DataFrame()
        self.coloc = {}
        self.coloc_susie_res = pd.DataFrame()

        # Validate and check sumstats
        self.log.write("Start to create SumstatsPair object...")
        self.log.write(" -Checking sumstats 1...", verbose=verbose)
        check_datatype(sumstatsObject1.data, log=self.log, verbose=verbose)
        check_dataframe_shape(sumstats=sumstatsObject1.data, log=self.log, verbose=verbose)
        
        self.log.write(" -Checking sumstats 2...", verbose=verbose)
        check_datatype(sumstatsObject2.data, log=self.log, verbose=verbose)
        check_dataframe_shape(sumstats=sumstatsObject2.data, log=self.log, verbose=verbose)

        # Classify columns from sumstats1
        info_headers = _get_headers(mode="info")
        stats_headers = _get_headers(mode="stats")
        
        for col in sumstatsObject1.data.columns:
            if col in info_headers:
                self.snp_info_cols.append(col)
            elif col in stats_headers:
                self.stats_cols.append(col)
            else:
                self.other_cols.append(col)

        # Classify columns from sumstats2 (skip info columns)
        for col in sumstatsObject2.data.columns:
            if col in stats_headers:
                self.stats_cols2.append(col)
            elif col not in info_headers:
                self.other_cols2.append(col)
        
        self.log.write(f" -Variant Info columns: {self.snp_info_cols}", verbose=verbose)
        self.log.write(f" -Variant statistics columns: {self.stats_cols}", verbose=verbose)
        self.log.write(f" -Sumstats1 other columns: {self.other_cols}", verbose=verbose)
        self.log.write(f" -Sumstats2 other columns: {self.other_cols2}", verbose=verbose)
        
        # Add temporary index columns for tracking
        sumstatsObject1.data["_RAW_INDEX_1"] = range(len(sumstatsObject1.data))
        sumstatsObject2.data["_RAW_INDEX_2"] = range(len(sumstatsObject2.data))
        
        # Prepare data from sumstats1 with suffixes
        self.data = sumstatsObject1.data.copy()
        rename_dict = {
            "EA": "EA_1",
            "NEA": "NEA_1",
            **{col: col + suffixes[0] for col in self.stats_cols},
            **{col: col + suffixes[0] for col in self.other_cols}
        }
        self.data = self.data.rename(columns=rename_dict)

        # Merge with sumstats2
        self.data, self.sumstats1, self.sumstats2 = self._merge_two_sumstats(
            sumstatsObject2, suffixes=suffixes, keep_all_variants=keep_all_variants, verbose=verbose
        )

        # Calculate sample sizes
        n_col1 = f"N{suffixes[0]}"
        n_col2 = f"N{suffixes[1]}"
        if n_col1 in self.data.columns and n_col2 in self.data.columns:
            # Use skipna=True to handle NA values (from keep_all_variants=True)
            n1_mean = self.data[n_col1].mean(skipna=True)
            n2_mean = self.data[n_col2].mean(skipna=True)
            # Check if mean is valid (not NaN)
            if pd.notna(n1_mean) and pd.notna(n2_mean):
                n1 = int(floor(n1_mean))
                n2 = int(floor(n2_mean))
                self.ns = (n1, n2)
            else:
                self.ns = None
        else:
            self.ns = None
        
        # Clean up temporary index columns
        sumstatsObject1.data = sumstatsObject1.data.drop(columns=["_RAW_INDEX_1"])
        sumstatsObject2.data = sumstatsObject2.data.drop(columns=["_RAW_INDEX_2"])

    def _merge_two_sumstats(
        self, 
        sumstatsObject2: 'Sumstats', 
        suffixes: Tuple[str, str] = ("_1","_2"),
        keep_all_variants: bool = True,
        verbose: bool = True
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        suffix1 = suffixes[0]
        
        # Determine merge mode: use "outer" to keep all variants if requested
        merge_mode = "outer" if keep_all_variants else "inner"
        
        # Merge by chromosome and position (sumstats1 has suffix _1, sumstats2 has no suffix yet)
        molded_sumstats, sumstats1, sumstats2 = _merge_mold_with_sumstats_by_chrpos(
            mold=self.data,
            sumstats_or_dataframe=sumstatsObject2.data,
            log=self.log,
            verbose=verbose,
            stats_cols1=self.stats_cols,
            stats_cols2=self.stats_cols2,
            suffixes=(suffix1, ""),
            merge_mode=merge_mode,
            return_not_matched_mold=True,
            keep_all_variants=keep_all_variants
        )
        
        # Align alleles between the two sumstats (only for variants present in both)
        molded_sumstats = _align_with_mold(
            molded_sumstats, 
            log=self.log, 
            verbose=verbose,
            suffixes=(suffix1, ""),
            keep_all_variants=keep_all_variants
        )
        
        # Flip sumstats2 statistics when alleles are flipped (only for variants present in both)
        molded_sumstats = _flip_allele_stats(molded_sumstats, log=self.log, verbose=verbose)
        
        # Drop sumstats2 EA/NEA columns and rename sumstats1 EA/NEA to standard names
        molded_sumstats = molded_sumstats.drop(columns=["EA", "NEA"], errors='ignore')
        molded_sumstats = molded_sumstats.rename(columns={"EA_1": "EA", "NEA_1": "NEA"})
        
        # Fill missing columns if sumstats2 doesn't have all the same stats columns
        if set(self.stats_cols2) != set(self.stats_cols):
            cols_to_fill = set(self.stats_cols).difference(set(self.stats_cols2))
            molded_sumstats = _fill_missing_columns(
                molded_sumstats, cols_to_fill, log=self.log, verbose=verbose
            )

        # Apply suffix _2 to sumstats2 columns
        cols_to_rename = self.stats_cols + self.other_cols2
        molded_sumstats = _renaming_cols(
            molded_sumstats, 
            cols_to_rename, 
            log=self.log, 
            verbose=verbose, 
            suffixes=suffixes
        )
        
        # Sort columns in standard order
        molded_sumstats = _sort_pair_cols(molded_sumstats, verbose=verbose, log=self.log)
        
        return molded_sumstats, sumstats1, sumstats2


    def clump(self, **kwargs: Any) -> None:
        self.clumps["clumps"],self.clumps["clumps_raw"],self.clumps["plink_log"] = _clump(self, log=self.log, p="P_1",mlog10p="MLOG10P_1", study = self.meta["gwaslab"]["group_name"], **kwargs)

    def to_coloc(self, **kwargs: Any) -> None:
        self.coloc["path"],self.coloc["file"],self.coloc["plink_log"] = _to_finemapping(self,study=self.meta["gwaslab"]["group_name"],suffixes=self.suffixes,log=self.log,**kwargs)

    def to_mesusie(self, **kwargs: Any) -> None:
        self.mesusie["path"],self.mesusie["file"],self.mesusie["plink_log"] = tofinemapping_m(self.data,
                                                                                              studies = self.study_names,
                                                                                              group = self.meta["gwaslab"]["group_name"],
                                                                                              suffixes=self.suffixes,
                                                                                              log=self.log,
                                                                                              **kwargs)
        
    def run_mesusie(self, **kwargs: Any) -> None:
        # Run MESuSiE using new framework (returns DataFrame directly, processing done internally)
        self.mesusie_res = _run_mesusie(
            self.mesusie["path"],
            log=self.log,
            ncols=self.ns,
            reference_data=self.data,
            group_name=self.meta["gwaslab"]["group_name"],
            study_name=self.study_name,
            **kwargs
        )
    
    def run_multisusie_rss(self, **kwargs: Any) -> None:
        if "path" not in self.mesusie or self.mesusie["path"] is None:
            self.log.warning(" -MESuSiE filelist not found. Please run to_mesusie() first.", 
                           verbose=kwargs.get("verbose", True))
            self.multisusie_res = pd.DataFrame()
            return
        
        # Run MultiSuSiE using filepath from to_mesusie (processing done internally)
        self.multisusie_res = _run_multisusie_rss(
            filepath=self.mesusie["path"],
            log=self.log,
            population_sizes=self.ns,
            reference_data=self.data,
            group_name=self.meta["gwaslab"]["group_name"],
            study_name=self.study_name,
            **kwargs
        )
    
    def run_ccgwas(self, **kwargs: Any) -> None:
         _run_ccgwas(self.data, 
                      meta = self.meta,
                      ldsc = self.ldsc,
                      ldsc_rg = self.ldsc_rg,
                      group=self.meta["gwaslab"]["group_name"],
                      studies = self.study_names,
                      log=self.log,
                      **kwargs)

    def read_pipcs(self, prefix: str, **kwargs: Any) -> None:
        self.mesusie_res = _read_pipcs(self.data[["SNPID","CHR","POS"]], 
                                   prefix, 
                                   group=self.meta["gwaslab"]["group_name"], 
                                   studie_names = self.study_name,
                                   **kwargs)
         
    def run_coloc_susie(self, **kwargs: Any) -> None:
        self.coloc_susie_res = _run_coloc_susie(self, self.coloc["path"],log=self.log,ncols=self.ns,**kwargs)

    def run_two_sample_mr(self, clump: bool = False, **kwargs: Any) -> None:
        exposure1 = self.study_names[0]
        outcome2 = self.study_names[1]
        _run_two_sample_mr(self,exposure1=exposure1,outcome2=outcome2, clump=clump,**kwargs)

    def run_meta_analysis(self, **kwargs: Any) -> Any:
        return meta_analyze_multi(self.data, nstudy=2, log=self.log, **kwargs)

    def extract_with_ld_proxy(self, **arg: Any) -> pd.DataFrame:
        return _extract_with_ld_proxy(common_sumstats = self.data, sumstats1=self.sumstats1,  **arg)

    def filter_value(self, expr: str, inplace: bool = False, **kwargs: Any) -> Optional['SumstatsPair']:
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_values(new_Sumstats_object.data,expr,log=new_Sumstats_object.log, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_values(self.data, expr,log=self.log,**kwargs)
        gc.collect()
        return None

    def _apply_viz_params(self, func: Callable[..., Any], kwargs: Dict[str, Any], key: Optional[str] = None, mode: Optional[str] = None) -> Dict[str, Any]:
        params = self.viz_params.merge(key or func.__name__, kwargs, mode=mode)
        return self.viz_params.filter(func, params, key=key or func.__name__, mode=mode, log=self.log, verbose=kwargs.get("verbose", True))

    def stacked_mqq(self, **kwargs: Any) -> None:
        df1 = self.data[["SNPID","CHR","POS","EA","NEA","P_1"]].rename(columns={"P_1":"P"})
        df2 = self.data[["SNPID","CHR","POS","EA","NEA","P_2"]].rename(columns={"P_2":"P"})
        objects = [df1, df2]
        if isinstance(self.mesusie_res, pd.DataFrame) and "PIP" in self.mesusie_res.columns:
            objects.append(self.mesusie_res)
        if isinstance(self.multisusie_res, pd.DataFrame) and "PIP" in self.multisusie_res.columns:
            objects.append(self.multisusie_res)

        if "vcfs" not in kwargs or kwargs.get("vcfs") is None:
            kwargs["vcfs"] = [None]
        if "region_ld_threshold" not in kwargs or kwargs.get("region_ld_threshold") is None:
            kwargs["region_ld_threshold"] = [0.2,0.4,0.6,0.8]
        if "region_ld_colors" not in kwargs or kwargs.get("region_ld_colors") is None:
            kwargs["region_ld_colors"] = ["#E4E4E4","#020080","#86CEF9","#24FF02","#FDA400","#FF0000","#FF0000"]
        if "region_marker_shapes" not in kwargs or kwargs.get("region_marker_shapes") is None:
            kwargs["region_marker_shapes"] = ['o','^','s','D','*','P','X','h','8']
        if "titles" not in kwargs:
            kwargs["titles"] = self.study_names

        params = self._apply_viz_params(plot_stacked_mqq, kwargs, key="plot_stacked_mqq")
        plot_stacked_mqq(objects=objects, **params)

    ## Visualization #############################################################################################################################################
    def plot_miami(self, **kwargs: Any) -> None:
        params = self._apply_viz_params(plot_miami2, kwargs, key="plot_miami2")
        plot_miami2(merged_sumstats=self.data, 
                    suffixes=self.suffixes,
                    cols1=["CHR","POS","P_1"],
                    cols2=["CHR","POS","P_2"],
                    **params)
    
    def compare_af(self, **kwargs: Any) -> Any:
        params = self._apply_viz_params(plotdaf, kwargs, key="plot_compare_af")
        return plotdaf( self.data,
                     eaf="EAF_2",
                     raf="EAF_1",
                     xlabel="Effect Allele Frequency in Sumstats 1",
                     ylabel="Effect Allele Frequency in Sumstats 2",
                     **params)

    def to_pickle(self, path: str = "~/mysumpair.pickle", overwrite: bool = False) -> None:
        dump_pickle_pair(self, path=path, overwrite=overwrite)

    def offload(self) -> None:
        # Only offload if data exists (skip if already offloaded)
        if hasattr(self, 'data') and self.data is not None:
            _offload(self.data, self.tmp_path, self.log)
            del self.data
            gc.collect()

    def reload(self, delete_files: Optional[List[str]] = None) -> None:
        self.data = _reload(self.tmp_path, self.log, delete_files=delete_files)
