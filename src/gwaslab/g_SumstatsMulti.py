import pandas as pd
import numpy as np
import copy
import gc
from typing import TYPE_CHECKING, Optional, Dict, Any, List, Tuple, Union, Callable

if TYPE_CHECKING:
    from gwaslab.info.g_Log import Log

from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats
from gwaslab.io.io_input_type import _get_id_column
from gwaslab.g_Sumstats_polars import Sumstatsp

from gwaslab.util.general.util_path_manager import _path

from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm.hm_casting import _align_with_mold
from gwaslab.hm.hm_casting import _fill_missing_columns
from gwaslab.hm.hm_casting import _renaming_cols
from gwaslab.hm.hm_casting import _sort_pair_cols
from gwaslab.hm.hm_casting_polars import _merge_mold_with_sumstats_by_chrposp
from gwaslab.hm.hm_casting_polars import _align_with_moldp
from gwaslab.hm.hm_casting_polars import _fill_missing_columnsp
from gwaslab.hm.hm_casting_polars import _renaming_colsp
from gwaslab.hm.hm_casting_polars import _sort_pair_colsp

from gwaslab.qc.qc_fix_sumstats import _flip_allele_stats
from gwaslab.qc.qc_fix_sumstats_polars import flipallelestatsp

from gwaslab.qc.qc_check_datatype_polars import check_datatype_polars
from gwaslab.qc.qc_check_datatype_polars import check_dataframe_shape_polars

from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_fix_sumstats import _process_build

from gwaslab.util.util_ex_calculate_ldmatrix import _to_finemapping
from gwaslab.util.rwrapper.util_ex_run_coloc import _run_coloc_susie
from gwaslab.util.util_in_filter_value import _filter_values
from gwaslab.util.rwrapper.util_ex_run_2samplemr import _run_two_sample_mr
from gwaslab.util.util_ex_run_clumping import _clump
from gwaslab.util.util_ex_ldproxyfinder import _extract_with_ld_proxy
from gwaslab.util.util_ex_match_ldmatrix import tofinemapping_m
from gwaslab.util.rwrapper.util_ex_run_mesusie2 import _run_mesusie
from gwaslab.util.util_in_meta import meta_analyze_multi
from gwaslab.util.rwrapper.util_ex_run_hyprcoloc import _run_hyprcoloc
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.util.util_in_fill_data import _get_multi_min
from gwaslab.util.util_ex_run_mtag import _run_mtag
from gwaslab.util.util_ex_run_multisusie import _run_multisusie_rss

from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2
from gwaslab.viz.viz_plot_compare_af import  plotdaf
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

from gwaslab.qc.qc_reserved_headers import _get_headers
from gwaslab.info.g_meta import _init_meta
from gwaslab.info.g_meta import _update_meta

from gwaslab.io.io_to_pickle import _offload
from gwaslab.io.io_to_pickle import _reload
from gwaslab.io.io_to_pickle import dump_pickle_multi

class SumstatsMulti( ):
    def __init__(
        self, 
        sumstatsObjects: List[Union[Sumstats, Sumstatsp]], 
        group_name: Optional[str] = None, 
        build: str = "99",
        engine: str = "pandas",
        merge_by_id: bool = False,
        keep_all_variants: bool = True,
        verbose: bool = True
    ) -> None:
        
        for i,sumstatsObject in enumerate(sumstatsObjects):
            if not isinstance(sumstatsObject, Sumstats):
                if not isinstance(sumstatsObject, Sumstatsp):
                    raise ValueError("Please provide GWASLab Sumstats Object #{}.".format(i+1))
        
        self.log = Log()
        self.meta = _init_meta(object="SumstatsMulti") 
        self.id = id(self)
        self.tmp_path = _path(pid=self.id, 
                              log = self.log, 
                              verbose=verbose)
        self.viz_params = VizParamsManager()
        load_viz_config(self.viz_params)
        
        # Determine merge_mode from keep_all_variants (consistent with SumstatsPair)
        # keep_all_variants=True (default) is equivalent to merge_mode="outer" (pandas) or "full" (polars)
        # keep_all_variants=False is equivalent to merge_mode="inner"
        if engine == "polars":
            import polars as pl
            merge_mode = "full" if keep_all_variants else "inner"
        else:
            merge_mode = "outer" if keep_all_variants else "inner"

        self.engine=engine
        self.keep_all_variants = keep_all_variants
            
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

            if engine=="polars":
                check_datatype_polars(sumstatsObject.data, log=self.log, verbose=verbose)
                check_dataframe_shape_polars(sumstats=sumstatsObject.data, 
                                log=self.log, 
                                verbose=verbose)
            else:
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
                    self.data = self._merge_two_sumstats(
                        sumstatsObject.data,
                        i=i,
                        merge_mode=merge_mode,
                        engine=engine,
                        merge_by_id=merge_by_id,
                        keep_all_variants=keep_all_variants
                    )
                    self.log.write("Finished merging Sumstats #{} to main DataFrame.".format(i+1))
        else:
            for i, sumstatsObject in enumerate(sumstatsObjects):
                if i >0:
                    self.log.write("Merging Sumstats #{} to main DataFrame...".format(i+1))
                    self.data = self._merge_two_sumstats(
                        sumstatsObject.data,
                        i=i,
                        merge_mode=merge_mode,
                        engine=engine,
                        keep_all_variants=keep_all_variants
                    )
                    self.log.write("Finished merging Sumstats #{} to main DataFrame.".format(i+1))
        

    def _merge_two_sumstats(
        self, 
        sumstatsObject2: Union[pd.DataFrame, Any], 
        verbose: bool = True,
        merge_mode: str = "outer",
        engine: str = "pandas",
        merge_by_id: bool = False,
        keep_all_variants: bool = True,
        i: int = 0
    ) -> Union[pd.DataFrame, Any]:

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
                                                        return_not_matched_mold = False,
                                                        keep_all_variants=keep_all_variants)
            molded_sumstats = _align_with_moldp(molded_sumstats, log=self.log, verbose=verbose, suffixes=("_1",""), keep_all_variants=keep_all_variants)
            molded_sumstats = flipallelestatsp(molded_sumstats, log=self.log, verbose=verbose)
            molded_sumstats = molded_sumstats.drop(["EA","NEA"] )
            molded_sumstats = molded_sumstats.rename({"EA_1":"EA","NEA_1":"NEA"})
        else:
            molded_sumstats = _merge_mold_with_sumstats_by_chrpos(mold=self.data, 
                                                        sumstats_or_dataframe=sumstatsObject2, 
                                                        log=self.log,
                                                        verbose=verbose,
                                                        merge_mode=merge_mode,
                                                        stats_cols1 = self.other_cols[0],
                                                        stats_cols2 = self.other_cols[i],
                                                        suffixes=("_1",""),
                                                        return_not_matched_mold = False,
                                                        keep_all_variants=keep_all_variants)    
            molded_sumstats = _align_with_mold(molded_sumstats, log=self.log, verbose=verbose, suffixes=("_1",""), keep_all_variants=keep_all_variants)
            molded_sumstats = _flip_allele_stats(molded_sumstats, log=self.log, verbose=verbose)        
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
    
    def _apply_viz_params(self, func: Callable[..., Any], kwargs: Dict[str, Any], key: Optional[str] = None, mode: Optional[str] = None) -> Dict[str, Any]:
        params = self.viz_params.merge(key or func.__name__, kwargs, mode=mode)
        return self.viz_params.filter(func, params, key=key or func.__name__, mode=mode, log=self.log, verbose=kwargs.get("verbose", True))

    def update_meta(self, **kwargs: Any) -> None:
        self.meta = _update_meta(self.meta, self.data, log = self.log, **kwargs)

    def run_meta_analysis(self, **kwargs: Any) -> Any:
        if self.engine == "polars":
            from gwaslab.util.util_in_meta_polars import meta_analyze_polars
            # Filter out verbose and log from kwargs as meta_analyze_polars doesn't accept them
            filtered_kwargs = {k: v for k, v in kwargs.items() if k not in ["verbose", "log"]}
            return meta_analyze_polars(self.data, nstudy=self.meta["gwaslab"]["number_of_studies"], log=self.log, **filtered_kwargs)
        else:
            return meta_analyze_multi(self.data,nstudy = self.meta["gwaslab"]["number_of_studies"] ,**kwargs)
    
    def run_hyprcoloc(self, **kwargs: Any) -> None:
        hyprcoloc_res_combined = _run_hyprcoloc(self.data,
                       nstudy = self.meta["gwaslab"]["number_of_studies"], 
                       study= self.meta["gwaslab"]["group_name"], 
                       traits=self.names, **kwargs)
        self.hyprcoloc = hyprcoloc_res_combined

    def run_mtag(self, **kwargs: Any) -> None:
        _run_mtag(     self,
                       nstudy = self.meta["gwaslab"]["number_of_studies"], 
                       study= self.meta["gwaslab"]["group_name"], 
                       traits=self.names, 
                       **kwargs)

    def run_meta_regression(
        self,
        covariate: Union[List[float], np.ndarray, Dict[int, float]],
        covariate_name: str = "covariate",
        center_covariate: bool = True,
        min_studies: int = 4,
        verbose: bool = True,
        **kwargs: Any
    ) -> Sumstats:
        """
        Perform meta-regression on GWAS summary statistics.
        
        This method fits a meta-regression model where effect sizes vary as a
        function of a study-level covariate. The model uses weighted least squares
        with inverse-variance weights.
        
        Parameters
        ----------
        covariate : list, array, or dict
            Study-level covariate values. Can be:
            - List/array of length K (number of studies): [x_1, x_2, ..., x_K]
            - Dict mapping study index (1-based) to covariate: {1: x_1, 2: x_2, ...}
        covariate_name : str, default "covariate"
            Name of the covariate for output columns (e.g., "age", "year").
        center_covariate : bool, default True
            If True, center the covariate by subtracting its mean.
            This makes the intercept interpretable as the effect at the mean covariate value.
        min_studies : int, default 3
            Minimum number of studies required per variant to perform meta-regression.
        verbose : bool, default True
            Whether to print progress messages.
        **kwargs : Any
            Additional keyword arguments passed to meta_regress.
        
        Returns
        -------
        Sumstats
            Sumstats object containing meta-regression results with columns:
            - All original SNP information (SNPID, CHR, POS, EA, NEA)
            - BETA_INTERCEPT: Intercept (effect at mean covariate if centered)
            - SE_INTERCEPT: Standard error of intercept
            - Z_INTERCEPT: Z-score for intercept
            - P_INTERCEPT: P-value for intercept
            - BETA_SLOPE_{COVARIATE_NAME}: Slope (change in effect per unit change in covariate)
            - SE_SLOPE_{COVARIATE_NAME}: Standard error of slope
            - Z_SLOPE_{COVARIATE_NAME}: Z-score for slope
            - P_SLOPE_{COVARIATE_NAME}: P-value for slope
            - N_STUDIES: Number of studies included for each variant
            - Q: Cochran's Q statistic for residual heterogeneity
            - P_HET: P-value for heterogeneity test
            - I2: I-squared statistic for heterogeneity
        
        Examples
        --------
        >>> import gwaslab as gl
        >>> # Load multiple sumstats
        >>> sumstats_list = [
        ...     gl.Sumstats("study1.txt.gz", verbose=False),
        ...     gl.Sumstats("study2.txt.gz", verbose=False),
        ...     gl.Sumstats("study3.txt.gz", verbose=False)
        ... ]
        >>> multi = gl.SumstatsMulti(sumstats_list, verbose=False)
        >>> # Mean ages for each study
        >>> ages = [45.2, 52.1, 48.7]
        >>> result = multi.run_meta_regression(covariate=ages, covariate_name="age")
        """
        # Note: Covariate-based meta-regression has been removed.
        # MR-MEGA extension now only contains the MDS-based MR-MEGA algorithm.
        # For covariate-based meta-regression, use a different method or implement separately.
        raise NotImplementedError(
            "Covariate-based meta-regression is not available in the MR-MEGA extension. "
            "The MR-MEGA extension only implements the MDS-based MR-MEGA algorithm. "
            "Use meta_regress_mrmega() for MR-MEGA analysis."
        )

    def run_multisusie_rss(self, R_list: List[np.ndarray], **kwargs: Any) -> pd.DataFrame:
        """
        Run finemapping using MultiSuSiE RSS (summary statistics version).
        
        This method runs MultiSuSiE on multiple populations/studies using summary statistics.
        It extracts BETA, SE, and N from the SumstatsMulti object for each study and runs
        MultiSuSiE finemapping.
        
        Parameters
        ----------
        R_list : List[np.ndarray]
            List of LD correlation matrices, one for each population/study.
            Each matrix should be PxP where P is the number of variants.
            Variants should be in the same order across all matrices.
        **kwargs : Any
            Additional keyword arguments passed to _run_multisusie_rss.
            See _run_multisusie_rss documentation for full list of parameters.
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing finemapping results with columns:
            - SNPID: Variant identifier
            - CHR: Chromosome (if available)
            - POS: Position (if available)
            - PIP: Posterior inclusion probability
            - CREDIBLE_SET_INDEX: Credible set index (0 if not in any set)
            - CS_CATEGORY: Credible set category (if available)
            - Additional columns from MultiSuSiE results
        
        Examples
        --------
        >>> # Assuming you have a SumstatsMulti object and LD matrices
        >>> results = sumstats_multi.run_multisusie_rss(
        ...     R_list=[ld_matrix_1, ld_matrix_2],
        ...     L=10,
        ...     max_iter=100
        ... )
        """
        kwargs = {k: v for k, v in kwargs.items() if k != "log"}
        verbose = kwargs.pop("verbose", True)
        return _run_multisusie_rss(
            gls=self,
            R_list=R_list,
            log=self.log,
            verbose=verbose,
            **kwargs
        )

    def get_lead(self, build: Optional[str] = None, gls: bool = False, **kwargs: Any) -> Union[pd.DataFrame, 'SumstatsMulti']:
        
        id_to_use = _get_id_column(self.data)
        
        # extract build information from meta data
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]

        self.data = _get_multi_min(self.data,
                                   col="P", 
                                   nstudy=self.meta["gwaslab"]["number_of_studies"])

        output = _get_sig(self.data,
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
    
    def to_pickle(self, path: str = "~/mysummulti.pickle", overwrite: bool = False) -> None:
        """
        Save SumstatsMulti object to a pickle file.
        
        Parameters
        ----------
        path : str, optional
            File path for the pickle file. Supports `~` for home directory expansion.
            Defaults to "~/mysummulti.pickle"
        overwrite : bool, optional
            If True, overwrite the file if it already exists. If False, skip saving if file exists.
            Defaults to False
        
        Returns
        -------
        None
            Saves the object to disk
        """
        dump_pickle_multi(self, path=path, overwrite=overwrite)

    def offload(self) -> None:
        _offload(self.data, self.tmp_path, self.log)
        del self.data
        gc.collect()

    def reload(self, delete_files: Optional[List[str]] = None) -> None:
        """
        Reload data from temporary pickle file.
        
        Parameters
        ----------
        delete_files : list of str, optional
            Additional files to delete after successful reload
        """
        self.data = _reload(self.tmp_path, self.log, delete_files=delete_files)