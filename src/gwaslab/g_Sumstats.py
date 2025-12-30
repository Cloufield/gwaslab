import copy
import gc
import time
from typing import TYPE_CHECKING, Optional, Union, Dict, Any, List, Tuple, Callable
import numpy as np
import pandas as pd
from gwaslab.info.g_object_helper import add_doc, suppress_display

if TYPE_CHECKING:
    from gwaslab.info.g_Log import Log

# ----- Version / Metadata / Logging -----
from gwaslab.info.g_Log import Log
from gwaslab.info.g_Sumstats_summary import lookupstatus, summarize
from gwaslab.info.g_version import _show_version, gwaslab_info
from gwaslab.info.g_meta import (
    _append_meta_record,
    _check_sumstats_qc_status,
    _init_meta,
    _set_harmonization_status,
    _set_qc_status,
    _update_step_status,
    _update_qc_step,
    _update_harmonize_step,
)
from gwaslab.info.g_meta import _update_meta
from gwaslab.qc.qc_build import _set_build
# ----- QC: Fix Sumstats -----
from gwaslab.qc.qc_fix_sumstats import (
    _fix_allele,
    _fix_chr,
    _fix_ID,
    _fix_pos,
    _flip_allele_stats,
    _flip_SNPID,
    _parallelize_normalize_allele,
    _remove_dup,
    _sort_column,
    _sort_coordinate,
    _strip_SNPID,
    _process_build
)
from gwaslab.hm.hm_liftover_v2 import _liftover_variant
from gwaslab.qc.qc_sanity_check import (
    _sanity_check_stats,
    _check_data_consistency,
)
# ----- Harmonization -----
from gwaslab.hm.hm_harmonize_sumstats import (
    _check_ref,
    _parallelize_check_af,
    _parallelize_infer_af,
    _parallelize_infer_strand,
    _parallelize_assign_rsid,
    _parallelize_rsid_to_chrpos,
    _parallele_infer_af_with_maf,
)

from gwaslab.hm.hm_casting import (
    _align_with_mold,
    _merge_mold_with_sumstats_by_chrpos,
)

# ----- Utility: Filtering -----
from gwaslab.util.util_in_filter_value import (
    _filter_in,
    _filter_out,
    _filter_region_in,
    _filter_region_out,
    _filter_values,
    _infer_build,
    _sampling,
    _exclude_hla,
    _filter_indel,
    _filter_palindromic,
    _filter_region,
    _filter_snp,
    _get_flanking,
    _get_flanking_by_chrpos,
    _get_flanking_by_id,
    _get_region_start_and_end,
    _search_variants,
)

# ----- Utility: Calculation -----
from gwaslab.util.util_in_calculate_gc import _lambda_GC
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2
from gwaslab.util.util_in_estimate_ess import _get_ess
from gwaslab.util.util_in_get_density import (
    _get_signal_density2,
)
from gwaslab.util.util_in_get_sig import (
    _anno_gene,
    _get_novel,
    _get_sig,
    _get_top,
    _check_cis,
    _check_novel_set,
)
from gwaslab.util.util_in_fill_data import _fill_data
from gwaslab.view.view_sumstats import _view_sumstats
from gwaslab.view.view_report import generate_qc_report
from gwaslab.util.util_ex_phewwas import _extract_associations
from gwaslab.bd.bd_sex_chromosomes import Chromosomes
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

# ----- Utility: LDSC / PRS / LD -----
from gwaslab.util.util_ex_calculate_ldmatrix import _to_finemapping
from gwaslab.util.util_ex_calculate_prs import _calculate_prs
from gwaslab.util.util_ex_infer_ancestry import _infer_ancestry
from gwaslab.util.util_ex_ldproxyfinder import _extract_ld_proxy
from gwaslab.util.util_ex_ldsc import (
    _estimate_h2_by_ldsc,
    _estimate_h2_cts_by_ldsc,
    _estimate_partitioned_h2_by_ldsc,
    _estimate_rg_by_ldsc,
)
from gwaslab.util.util_ex_run_clumping import _clump
from gwaslab.util.util_ex_run_magma import _run_magma
from gwaslab.util.util_ex_run_prscs import _run_prscs
from gwaslab.util.util_ex_run_scdrs import _run_scdrs
from gwaslab.util.util_ex_run_susie import _get_cs_lead, _run_susie_rss

# ----- Utility: Fine-mapping ABF -----
from gwaslab.util.util_abf_finemapping import _abf_finemapping, _make_cs

# ----- Base Data (Reference Panels / Common Data) -----
from gwaslab.bd.bd_common_data import (
    get_chr_list,
    get_chr_to_number,
    get_format_dict,
    get_formats_list,
    get_high_ld,
    get_number_to_chr,
)
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap3
from gwaslab.bd.bd_path_manager import _path

# ----- Visualization -----
from gwaslab.viz.viz_plot_associations import _plot_associations
from gwaslab.viz.viz_plot_compare_af import plotdaf
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.viz.viz_plot_density import _process_density
from gwaslab.viz.viz_plot_effect import _plot_effect
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.viz.viz_plot_phe_heatmap import _gwheatmap
from gwaslab.viz.viz_plot_qqplot import _plot_qq
from gwaslab.viz.viz_plot_regional2 import _plot_regional
from gwaslab.viz.viz_plot_trumpetplot import _plot_trumpet
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

# ----- Downstream Analysis -----
from gwaslab.downstream.ds_result_manager import DownstreamResultManager

# ----- IO: Reading / Writing -----
from gwaslab.io.io_load_ld import _to_finemapping_using_ld
from gwaslab.io.io_preformat_input import _preformat
from gwaslab.io.io_read_pipcs import _read_pipcs
from gwaslab.io.io_to_formats import _to_format
from gwaslab.io.io_to_pickle import _offload, _reload, dump_pickle
from gwaslab.io.io_gwaslab_standard import write_gsf
from gwaslab.io.io_process_kwargs import remove_overlapping_kwargs
from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
from gwaslab.io.io_input_type import _get_id_column
from gwaslab.hm.hm_assign_rsid import _assign_rsid
from gwaslab.hm.hm_infer_with_af  import _infer_strand_with_annotation
from gwaslab.hm.hm_check_af import _check_af_with_annotation, _infer_af_with_annotation, _infer_af_with_maf_annotation
from datetime import datetime

# ----- Internal Helper Functions -----

#20220309
class Sumstats():
    
    @add_doc(_preformat)
    def __init__(self,
             sumstats: Optional[Union[str, pd.DataFrame]] = None,
             fmt: Optional[str] = None,
             tab_fmt: str = "tsv",
             snpid: Optional[str] = None,
             rsid: Optional[str] = None,
             chrom: Optional[str] = None,
             pos: Optional[str] = None,
             ea: Optional[str] = None,
             nea: Optional[str] = None,
             ref: Optional[str] = None,
             alt: Optional[str] = None,
             eaf: Optional[str] = None,
             neaf: Optional[str] = None,
             maf: Optional[str] = None,
             n: Optional[str] = None,
             beta: Optional[str] = None,
             se: Optional[str] = None,
             chisq: Optional[str] = None,
             z: Optional[str] = None,
             f: Optional[str] = None,
             t: Optional[str] = None,
             p: Optional[str] = None,
             q: Optional[str] = None,
             mlog10p: Optional[str] = None,
             test: Optional[str] = None,
             info: Optional[str] = None,
             OR: Optional[str] = None,
             OR_95L: Optional[str] = None,
             OR_95U: Optional[str] = None,
             beta_95L: Optional[str] = None,
             beta_95U: Optional[str] = None,
             HR: Optional[str] = None,
             HR_95L: Optional[str] = None,
             HR_95U: Optional[str] = None,
             ncase: Optional[str] = None,
             ncontrol: Optional[str] = None,
             neff: Optional[str] = None,
             i2: Optional[str] = None,
             phet: Optional[str] = None,
             dof: Optional[str] = None,
             snpr2: Optional[str] = None,
             status: Optional[str] = None,
             other: Optional[List[str]] = None,
             exclude: Optional[List[str]] = None,
             include: Optional[List[str]] = None,
             chrom_pat: Optional[str] = None,
             snpid_pat: Optional[str] = None,
             direction: Optional[str] = None,
             verbose: bool = True,
             study: str = "Study_1",
             trait: str = "Trait_1",
             build: str = "99",
             species: str = "homo sapiens",
             readargs: Optional[Dict[str, Any]] = None,
             **kwreadargs: Any) -> None:
        
        self.log = Log()
        # Store reference to Sumstats object in log for shape tracking
        self.log._sumstats_obj = self
        # Store verbose as attribute for use in other methods
        self.verbose = verbose
        # print gwaslab version information
        _show_version(self.log, verbose=verbose)

        # basic attributes
        self.data = pd.DataFrame()
        
        
        # Track shape for logging optimization (skip logging if shape unchanged)
        self._last_shape = None
        
        #preformat the data
        self.data  = _preformat(
          sumstats=sumstats,
          fmt=fmt,
          tab_fmt = tab_fmt,
          snpid=snpid,
          rsid=rsid,
          chrom=chrom,
          pos=pos,
          ea=ea,
          nea=nea,
          ref=ref,
          alt=alt,
          eaf=eaf,
          neaf=neaf,
          maf=maf,
          n=n,
          beta=beta,
          se=se,
          chisq=chisq,
          z=z,
          f=f,
          t=t,
          p=p,
          q=q,
          mlog10p=mlog10p,
          test=test,
          info=info,
          OR=OR,
          OR_95L=OR_95L,
          OR_95U=OR_95U,
          beta_95L=beta_95L,
          beta_95U=beta_95U,
          HR=HR,
          HR_95L=HR_95L,
          HR_95U=HR_95U,
          i2=i2,
          phet=phet,
          dof=dof,
          snpr2=snpr2,
          ncase=ncase,
          ncontrol=ncontrol,
          neff=neff,
          direction=direction,
          study=study,
          build=build,
          trait=trait,
          status=status,
          other=other,
          exclude=exclude,
          include=include,
          chrom_pat=chrom_pat,
          snpid_pat=snpid_pat,
          verbose=verbose,
          log=self.log, 
          readargs=readargs,
          **kwreadargs)

        # meta information
        self.meta = _init_meta() 
        self.meta["gwaslab"]["study_name"] = study
        self.meta["gwaslab"]["species"] = species

        # set build (using setter to ensure consistency)
        self.build = build
        self.viz_params = VizParamsManager()
        load_viz_config(self.viz_params)

        self.id = id(self)
        self.tmp_path = _path(pid=self.id, log = self.log, verbose=False)
        
        # Initialize downstream analysis result manager
        self.downstream = DownstreamResultManager()
        
        # Initialize Chromosomes instance based on species
        self.chromosomes = Chromosomes(species=species)
        
        # Initialize ChromosomeMapper instance
        # Build is already set above via the setter, so use the processed build value
        self.mapper = ChromosomeMapper(
            species=species,
            build=self._build,
            log=self.log,
            verbose=verbose
        )
        
        # Auto-detect and build sumstats layer if data is available
        self._update_mapper_from_data()
        
        # initialize attributes for clumping and finmapping
        #self.to_finemapping_file_path = ""
        #self.to_finemapping_file  = pd.DataFrame()
        #self.plink_log = ""

        # path / file / plink_log        
    
    @property
    def build(self) -> str:
        if self._build =="Unknown" or self._build =="99":
            if self.meta["gwaslab"]["species"] == "homo sapiens":
                self.log.warning("Build is unknown. .infer_build first.")
                try:
                    self.infer_build()
                except:
                    pass
        return self._build

    @build.setter
    def build(self, value: Union[str, int]) -> None:
        # Process the build value and update both _build and meta
        from gwaslab.qc.qc_build import _process_build
        # Ensure meta["gwaslab"] exists
        if "gwaslab" not in self.meta:
            self.meta["gwaslab"] = {}
        species = self.meta.get("gwaslab", {}).get("species", "homo sapiens")
        processed_build = _process_build(value, log=self.log, verbose=False, species=species)
        self._build = processed_build
        self.meta["gwaslab"]["genome_build"] = processed_build
        
        # Update mapper with new build if mapper exists
        if hasattr(self, 'mapper'):
            self.mapper.build = processed_build
            self.mapper._build_nc_mappings()  # Rebuild NC mappings for new build
    
    def _update_mapper_from_data(self, chrom_col: str = "CHR"):
        """
        Auto-detect chromosome format from data and rebuild sumstats layer in mapper.
        
        This method should be called whenever the chromosome column is modified
        (e.g., after fix_chr, after loading new data, etc.).
        
        Parameters
        ----------
        chrom_col : str, default="CHR"
            Column name for chromosome data in self.data.
        """
        if not hasattr(self, 'data') or self.data.empty:
            return
        
        if not hasattr(self, 'mapper'):
            return
        
        if chrom_col not in self.data.columns:
            return
        
        try:
            # Detect format and build sumstats layer
            self.mapper.detect_sumstats_format(self.data[chrom_col])
        except Exception as e:
            if self.verbose:
                self.log.warning(f"Could not auto-detect chromosome format: {e}")
    
    # ============================================================================
    # Downstream Analysis Result Properties (backward compatibility)
    # ============================================================================
    
    @property
    def ldsc_h2(self) -> Optional[float]:
        """LDSC heritability estimate (backward compatibility)."""
        return self.downstream.ldsc_h2
    
    @ldsc_h2.setter
    def ldsc_h2(self, value: Optional[float]) -> None:
        """Set LDSC heritability estimate (backward compatibility)."""
        self.downstream.ldsc_h2 = value
    
    @property
    def ldsc_h2_results(self) -> Optional[Any]:
        """Detailed LDSC heritability results (backward compatibility)."""
        return self.downstream.ldsc_h2_results
    
    @ldsc_h2_results.setter
    def ldsc_h2_results(self, value: Optional[Any]) -> None:
        """Set detailed LDSC heritability results (backward compatibility)."""
        self.downstream.ldsc_h2_results = value
    
    @property
    def ldsc_rg(self) -> pd.DataFrame:
        """LDSC genetic correlation results (backward compatibility)."""
        return self.downstream.ldsc_rg
    
    @ldsc_rg.setter
    def ldsc_rg(self, value: pd.DataFrame) -> None:
        """Set LDSC genetic correlation results (backward compatibility)."""
        self.downstream.ldsc_rg = value
    
    @property
    def ldsc_h2_cts(self) -> Optional[Any]:
        """LDSC cell-type-specific heritability results (backward compatibility)."""
        return self.downstream.ldsc_h2_cts
    
    @ldsc_h2_cts.setter
    def ldsc_h2_cts(self, value: Optional[Any]) -> None:
        """Set LDSC cell-type-specific heritability results (backward compatibility)."""
        self.downstream.ldsc_h2_cts = value
    
    @property
    def ldsc_partitioned_h2_summary(self) -> Optional[Any]:
        """LDSC partitioned heritability summary (backward compatibility)."""
        return self.downstream.ldsc_partitioned_h2_summary
    
    @ldsc_partitioned_h2_summary.setter
    def ldsc_partitioned_h2_summary(self, value: Optional[Any]) -> None:
        """Set LDSC partitioned heritability summary (backward compatibility)."""
        self.downstream.ldsc_partitioned_h2_summary = value
    
    @property
    def ldsc_partitioned_h2_results(self) -> Optional[Any]:
        """Detailed LDSC partitioned heritability results (backward compatibility)."""
        return self.downstream.ldsc_partitioned_h2_results
    
    @ldsc_partitioned_h2_results.setter
    def ldsc_partitioned_h2_results(self, value: Optional[Any]) -> None:
        """Set detailed LDSC partitioned heritability results (backward compatibility)."""
        self.downstream.ldsc_partitioned_h2_results = value
    
    @property
    def finemapping(self) -> Dict[str, Any]:
        """Finemapping results dictionary (backward compatibility)."""
        return self.downstream.finemapping
    
    @finemapping.setter
    def finemapping(self, value: Dict[str, Any]) -> None:
        """Set finemapping results dictionary (backward compatibility)."""
        self.downstream.finemapping = value
    
    @property
    def clumps(self) -> Dict[str, Any]:
        """Clumping results dictionary (backward compatibility)."""
        return self.downstream.clumps
    
    @clumps.setter
    def clumps(self, value: Dict[str, Any]) -> None:
        """Set clumping results dictionary (backward compatibility)."""
        self.downstream.clumps = value
    
    @property
    def pipcs(self) -> pd.DataFrame:
        """PIPCS results (backward compatibility)."""
        return self.downstream.pipcs
    
    @pipcs.setter
    def pipcs(self, value: pd.DataFrame) -> None:
        """Set PIPCS results (backward compatibility)."""
        self.downstream.pipcs = value

    def _apply_viz_params(self, func: Callable[..., Any], kwargs: Dict[str, Any], key: Optional[str] = None, mode: Optional[str] = None) -> Dict[str, Any]:
        params = self.viz_params.merge(key or func.__name__, kwargs, mode=mode)
        return self.viz_params.filter(func, params, key=key or func.__name__, mode=mode, log=self.log, verbose=kwargs.get("verbose", True))

    def __getitem__(self, index: Any) -> Any:
        return self.data[index]
        
    def __setitem__(self, index: Any, value: Any) -> None:
        self.data[index] = value

    def __len__(self) -> int:
        return len(self.data)
        
    def __getattr__(self, name: str) -> Any:
        # Use object.__getattribute__ to access self.data directly without triggering __getattr__ again
        # This prevents infinite recursion during pickle loading when self.data might not be initialized yet
        try:
            data = object.__getattribute__(self, 'data')
            return getattr(data, name)
        except AttributeError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")
    
    def __deepcopy__(self, memo: Dict[int, Any]) -> 'Sumstats':
        """Custom deepcopy implementation to properly copy Sumstats objects."""
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))
        return result
    
    def copy(self) -> 'Sumstats':
        """Return a deep copy of the Sumstats object."""
        return copy.deepcopy(self)
#### healper #################################################################################
    #@add_doc(_update_meta)
    #def update_meta(self, **kwargs):
    #    self.meta = _update_meta(self.meta, self.data,log = self.log, **kwargs)
    
    @add_doc(summarize)
    def summary(self) -> Any:
        return summarize(self)
    
    @add_doc(lookupstatus)
    def lookup_status(self, status: str = "STATUS") -> Any:
        return lookupstatus(self.data[status])
    
    @add_doc(_set_build)
    def set_build(self, build: Union[str, int], verbose: bool = True) -> None:
        self.data, processed_build = _set_build(self, build=build, log=self.log,verbose=verbose)
        # Update self.build to ensure consistency (setter will also update meta)
        self.build = processed_build
    
    @add_doc(_infer_build)
    def infer_build(self, verbose: bool = True, **kwargs: Any) -> None:
        kwargs = remove_overlapping_kwargs(kwargs,{"log","verbose"})
        # _infer_build now handles updating self.build and self.meta internally
        self.data = _infer_build(self,log=self.log,verbose=verbose,**kwargs)

    @add_doc(_liftover_variant)
    def liftover(self, to_build: Optional[str] = None, from_build: Optional[str] = None, chain_path: Optional[str] = None, **kwargs: Any) -> None:
        if from_build is None:
            from_build = self.build
        kwargs = remove_overlapping_kwargs(kwargs, {"from_build", "to_build", "chain_path", "log"})
        self.data = _liftover_variant(self, from_build=from_build, to_build=to_build, chain_path=chain_path, log=self.log, **kwargs)


# QC ######################################################################################
    #clean the sumstats with one line
    def basic_check(self,
                    remove=False,
                    remove_dup=False,
                    threads=1,
                    fix_id_kwargs={},
                    remove_dup_kwargs={},
                    fix_chr_kwargs={},
                    fix_pos_kwargs={},
                    fix_allele_kwargs={},
                    sanity_check_stats_kwargs={},
                    consistency_check_kwargs={},
                    normalize=True,
                    normalize_allele_kwargs={},
                    verbose=True):
        """
        All-in-one function for Sumstats quality control (QC), which is a wrapper of separate functions including:
        `fix_id` for SNPID and rsID check, `fix_chr` for chromosome notation (CHR) check,
        `fix_pos` for basepair position (POS) check, `fix_allele` for allele notation (EA and NEA) check,
        `check_sanity` for statistics sanity check and datatype check (BETA, SE, P and so forth),
        `check_data_consistency` for checking if convitable data are consistent (e.g., if the calculated BETA/SE is close to the original Z),
        `normalize_allele` for indel normalization, `remove_dup` for removal of multi-allelic variants, indels, and duplicated variants,
        `sort_coordinate` for sorting the genomic coordinates, and `sort_column` for sorting the order of columns in the dataframe.

        For a detailed description of all checks performed, see the documentation in docs/QC&Filtering.md.

        Parameters
        ----------
        remove : bool, optional
            Whether to remove bad quality variants detected in _fix_chr, _fix_pos, and _fix_allele.
        remove_dup : bool, optional
            Whether to remove duplicated or multi-allelic variants using remove_dup.
        threads : int, optional
            Number of threads to use for parallel processing.
        fix_id_kwargs : dict, optional
            Keyword arguments passed to `fix_id`.
        remove_dup_kwargs : dict, optional
            Keyword arguments passed to `remove_dup`.
        fix_chr_kwargs : dict, optional
            Keyword arguments passed to `fix_chr`.
        fix_pos_kwargs : dict, optional
            Keyword arguments passed to `fix_pos`.
        fix_allele_kwargs : dict, optional
            Keyword arguments passed to `fix_allele`.
        sanity_check_stats_kwargs : dict, optional
            Keyword arguments passed to `check_sanity`.
        consistency_check_kwargs : dict, optional
            Keyword arguments passed to `check_data_consistency`.
        normalize : bool, optional
            Whether to perform indel normalization.
        normalize_allele_kwargs : dict, optional
            Keyword arguments passed to `normalize_allele`.
        verbose : bool, optional
            Whether to print progress information.
        """
        saved_kwargs = {**locals()}
        # prepare status section
        ###############################################
        # Sequential version (original implementation)
        # try to fix data without dropping any information 
        
        fix_id_kwargs = remove_overlapping_kwargs(fix_id_kwargs,{"log", "remove", "verbose"})
        self.data = _fix_ID(self,log=self.log,verbose=verbose, **fix_id_kwargs)
        
        fix_chr_kwargs = remove_overlapping_kwargs(fix_chr_kwargs,{"log", "remove", "verbose"})
        self.data = _fix_chr(self,log=self.log,remove=remove,verbose=verbose,**fix_chr_kwargs)
        # Auto-detect and rebuild mapper after chromosome column is modified
        self._update_mapper_from_data()
        
        fix_pos_kwargs = remove_overlapping_kwargs(fix_pos_kwargs,{"log", "remove", "verbose"})
        self.data = _fix_pos(self,log=self.log,remove=remove,verbose=verbose,**fix_pos_kwargs)
        
        fix_allele_kwargs = remove_overlapping_kwargs(fix_allele_kwargs,{"log", "remove", "verbose"})
        self.data = _fix_allele(self,log=self.log, remove=remove, verbose=verbose,**fix_allele_kwargs)
        
        sanity_check_stats_kwargs = remove_overlapping_kwargs(sanity_check_stats_kwargs,{"log", "remove", "verbose"})
        self.data = _sanity_check_stats(self,log=self.log,verbose=verbose,**sanity_check_stats_kwargs)
        
        consistency_check_kwargs = remove_overlapping_kwargs(consistency_check_kwargs,{"log", "remove", "verbose"})
        _check_data_consistency(self,log=self.log,verbose=verbose,**consistency_check_kwargs)
        
        if normalize is True:
            normalize_allele_kwargs = remove_overlapping_kwargs(normalize_allele_kwargs,{"log", "remove", "verbose","threads"})
            self.data = _parallelize_normalize_allele(self,threads=threads,verbose=verbose,log=self.log,**normalize_allele_kwargs)
        
        if remove_dup is True:
            remove_dup_kwargs = remove_overlapping_kwargs(remove_dup_kwargs,{"log", "remove", "verbose"})
            self.data = _remove_dup(self,log=self.log,verbose=verbose,**remove_dup_kwargs)
        
        self.data = _sort_coordinate(self,verbose=verbose,log=self.log)
        self.data = _sort_column(self,verbose=verbose,log=self.log)
        
        # is_sorted is already set by _sort_column
        
        _set_qc_status(self, saved_kwargs)
        ###############################################
        return self
        
    
    def harmonize(self,
              basic_check=True,
              ref_seq=None,
              ref_rsid_tsv=None,
              ref_rsid_vcf=None,
              ref_infer=None,
              ref_alt_freq=None,
              ref_maf_threshold = 0.4,
              maf_threshold=0.40,
              threads=1,
              remove=False,
              check_ref_kwargs={},
              remove_dup_kwargs={},
              assign_rsid_kwargs={},
              infer_strand_kwargs={},
              flip_allele_stats_kwargs={},
              liftover_kwargs={},
              fix_id_kwargs={},
              fix_chr_kwargs={},
              fix_pos_kwargs={},
              fix_allele_kwargs={},
              sanity_check_stats_kwargs={},
              normalize_allele_kwargs={},
              verbose=True,
              sweep_mode=False
              ):
        """
        Standard pipeline for harmonizing sumstats including:
        1. Basic check and standardization using fix_id, fix_chr, fix_pos, fix_allele, check_sanity, normalize_allele, sort_coordinate, and sort_column.
        2. Reference-based annotation and flipping using check_ref, flip_allele_stats, infer_strand, and assign_rsid.
        3. Optional duplicate removal with removedup and final sorting via sort_coordinate and sort_column.
        For infer_strand, and assign_rsid, check_ref with ref_seq is required.

        Parameters
        ----------
        basic_check : bool, optional
            Whether to run basic QC pipeline (_fix_ID, _fix_chr, _fix_pos, etc.)
        ref_seq : str or None
            Full path to reference sequence file in fasta format for allele flipping
        ref_rsid_tsv : str or None
            Full path to rsID TSV reference file
        ref_rsid_vcf : str or None
            Full path to rsID VCF/BCF reference file
        ref_infer : str or None
            Full path to Reference VCF/BCF file for strand inference
        ref_alt_freq : str or None
            Allele frequency field name in VCF/BCF INFO for strand inference. Default is "AF".
        ref_maf_threshold : float
            MAF threshold (applied to reference VCF/BCF) for strand inference
        maf_threshold : float
            MAF threshold (applied to sumstats) for strand inference
        threads : int, optional
            Number of threads for parallel processing
        remove : bool, optional
            Whether to remove bad variants during QC
        check_ref_kwargs : dict, optional
            Arguments passed to check_ref
        remove_dup_kwargs : dict, optional
            Arguments passed to remove_dup
        assign_rsid_kwargs : dict, optional
            Arguments passed to assign_rsid
        infer_strand_kwargs : dict, optional
            Arguments passed to infer_strand
        flip_allele_stats_kwargs : dict, optional
            Arguments passed to flip_allele_stats
        fix_id_kwargs : dict, optional
            Arguments passed to fix_ID
        fix_chr_kwargs : dict, optional
            Arguments passed to fix_chr
        fix_pos_kwargs : dict, optional
            Arguments passed to fix_pos
        fix_allele_kwargs : dict, optional
            Arguments passed to fix_allele
        sanity_check_stats_kwargs : dict, optional
            Arguments passed to check_sanity
        normalize_allele_kwargs : dict, optional
            Arguments passed to normalize_allele
        verbose : bool, optional
            Whether to print progress information.
        sweep_mode:bool, optional
            If false, use lookup (per-variant) mode. If true, use sweep model (fast for large dataset).
        """
        saved_kwargs = {**locals()}
        # prepare status helpers

        if basic_check is True:
            fix_id_kwargs = remove_overlapping_kwargs(fix_id_kwargs,{"log", "remove", "verbose"})
            self.data = _fix_ID(self,log=self.log,verbose=verbose,**fix_id_kwargs)
            fix_chr_kwargs = remove_overlapping_kwargs(fix_chr_kwargs,{"log", "remove", "verbose"})
            self.data = _fix_chr(self,remove=remove,log=self.log,verbose=verbose,**fix_chr_kwargs)
            fix_pos_kwargs = remove_overlapping_kwargs(fix_pos_kwargs,{"log", "remove", "verbose"})
            self.data = _fix_pos(self,remove=remove,log=self.log,verbose=verbose,**fix_pos_kwargs)
            fix_allele_kwargs = remove_overlapping_kwargs(fix_allele_kwargs,{"log", "remove", "verbose"})
            self.data = _fix_allele(self,log=self.log,verbose=verbose,**fix_allele_kwargs)
            sanity_check_stats_kwargs = remove_overlapping_kwargs(sanity_check_stats_kwargs,{"log", "remove", "verbose"})
            self.data = _sanity_check_stats(self,log=self.log,verbose=verbose,**sanity_check_stats_kwargs)
            normalize_allele_kwargs = remove_overlapping_kwargs(normalize_allele_kwargs,{"log", "remove", "verbose","threads"})
            self.data = _parallelize_normalize_allele(self,log=self.log,threads=threads,verbose=verbose,**normalize_allele_kwargs)
            
            self.data = _sort_column(self,log=self.log,verbose=verbose)

            self.data = _sort_coordinate(self,log=self.log,verbose=verbose)

            _set_qc_status(self, saved_kwargs)
            gc.collect()
        
        #####################################################
        #part 2 : annotating and flipping
        #    2.1  ref check -> flip allele and allel-specific stats
        #    2.2  assign rsid
        #    2.3 infer strand for palindromic SNP
        #
        ########## liftover ###############
        #    3 : liftover by chr and pos to target build  -> reset status
        ###################################
        #   3.1 ref check (target build) -> flip allele and allel-specific stats  
        #   3.2  assign rsid (target build)
        #   3.2 infer strand for palindromic SNP (target build)
        #####################################################
        if ref_seq is not None:
            check_ref_kwargs = remove_overlapping_kwargs(check_ref_kwargs,{"log", "ref_seq"})
            self.data = _check_ref(self,ref_seq,log=self.log,**check_ref_kwargs)
            # Metadata and harmonization status are already set by _check_ref
            flip_allele_stats_kwargs = remove_overlapping_kwargs(flip_allele_stats_kwargs,{"log","verbose"})
            self.data = _flip_allele_stats(self,log=self.log,verbose=verbose,**flip_allele_stats_kwargs)
            # Harmonization status is already set by _flip_allele_stats
            
            gc.collect()
            
        if ref_infer is not None: 
            infer_strand_kwargs = remove_overlapping_kwargs(infer_strand_kwargs,{"log","verbose","ref_infer","ref_alt_freq","maf_threshold","ref_maf_threshold","threads","path","assign_cols"})
            if sweep_mode:
                self.data = _infer_strand_with_annotation(self, 
                                                        path = ref_infer, 
                                                        assign_cols = ref_alt_freq, 
                                                        threads=threads,
                                                        log=self.log,
                                                        verbose=verbose,
                                                        **infer_strand_kwargs)
            else:
                self.data= _parallelize_infer_strand(self,ref_infer = ref_infer,ref_alt_freq=ref_alt_freq,
                                                maf_threshold=maf_threshold, ref_maf_threshold=ref_maf_threshold,
                                                  threads=threads,log=self.log,verbose=verbose,**infer_strand_kwargs)

            # Metadata and harmonization status are already set by _parallelize_infer_strand or _infer_strand_with_annotation
            
            flip_allele_stats_kwargs = remove_overlapping_kwargs(flip_allele_stats_kwargs,{"log","verbose"})
            self.data =_flip_allele_stats(self,log=self.log,verbose=verbose,**flip_allele_stats_kwargs)
            # Harmonization status is already set by _flip_allele_stats
            
            gc.collect()
            
        if (ref_seq is not None or ref_infer is not None) and (ref_rsid_tsv is not None or ref_rsid_vcf is not None):

            self.data = _fix_ID(self, log=self.log, verbose=verbose, **{"fixid":True, "fixsep":True, "overwrite":True})

            gc.collect()
        
        #####################################################
        
        if ref_rsid_tsv is not None:
            assign_rsid_kwargs = remove_overlapping_kwargs(assign_rsid_kwargs,{"ref_mode","path","threads","log","verbose"})
            self.data = _parallelize_assign_rsid(self,path=ref_rsid_tsv,ref_mode="tsv",
                                                 threads=threads,log=self.log,verbose=verbose,**assign_rsid_kwargs)

            # Metadata and harmonization status are already set by _parallelize_assign_rsid
            gc.collect()

        if ref_rsid_vcf is not None:
            assign_rsid_kwargs = remove_overlapping_kwargs(assign_rsid_kwargs,{"ref_mode","path","threads","log","verbose"})

            if sweep_mode:
                self.data = _assign_rsid(self, path = ref_rsid_vcf, threads=threads,log=self.log,verbose=verbose,**assign_rsid_kwargs)
            else:
                self.data = _parallelize_assign_rsid(self,path=ref_rsid_vcf,ref_mode="vcf",
                                                     threads=threads,log=self.log,verbose=verbose,**assign_rsid_kwargs)   

            # Metadata and harmonization status are already set by _assign_rsid or _parallelize_assign_rsid
            gc.collect()
        
        ######################################################    
        if remove is True:
            remove_dup_kwargs = remove_overlapping_kwargs(remove_dup_kwargs,{"log","verbose"})
            self.data = _remove_dup(self,log=self.log,verbose=verbose,**remove_dup_kwargs)
        ################################################ 
        
        self.data = _sort_coordinate(self,log=self.log,verbose=verbose)
        
        self.data = _sort_column(self,log=self.log,verbose=verbose)
        gc.collect()
        # is_sorted is already set by _sort_column, just set is_harmonised
        self.meta["is_harmonised"] = True

        _set_harmonization_status(self, saved_kwargs)
        return self
    
    def align_with_template(self, template, **kwargs):
        ## merge
        molded_sumstats, sumstats1 = _merge_mold_with_sumstats_by_chrpos(mold=template, 
                                            sumstats_or_dataframe=self,
                                            log=self.log,
                                            suffixes=("_MOLD",""),
                                            return_not_matched_mold = True)
        ## align
        aligned_data = _align_with_mold(molded_sumstats)
        # Assign aligned data back to self.data
        self.data = aligned_data

        ## flip
        self.data =_flip_allele_stats(self, log=self.log)
    
    @add_doc(_check_sumstats_qc_status)
    def check_sumstats_qc_status(self) -> Any:
        return _check_sumstats_qc_status(self)
    ############################################################################################################
    #customizable API to build your own QC pipeline
    @add_doc(_fix_ID)
    def fix_id(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _fix_ID(self,log=self.log,**kwargs)
        return self
    @add_doc(_flip_SNPID)
    def flip_snpid(self, **kwargs: Any) -> None:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _flip_SNPID(self,log=self.log,**kwargs)
    @add_doc(_strip_SNPID)
    def strip_snpid(self, **kwargs: Any) -> None:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _strip_SNPID(self,log=self.log,**kwargs)
    @add_doc(_fix_chr)
    def fix_chr(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _fix_chr(self,log=self.log,**kwargs)
        # Auto-detect and rebuild mapper after chromosome column is modified
        self._update_mapper_from_data()
        return self
    @add_doc(_fix_pos)
    def fix_pos(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _fix_pos(self,log=self.log,**kwargs)
        return self
    @add_doc(_fix_allele)
    def fix_allele(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _fix_allele(self,log=self.log,**kwargs)
        return self
    @add_doc(_remove_dup)
    def remove_dup(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _remove_dup(self,log=self.log,**kwargs)
        return self
    @add_doc(_sanity_check_stats)
    def check_sanity(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _sanity_check_stats(self,log=self.log,**kwargs)
        return self
    @add_doc(_check_data_consistency)
    def check_data_consistency(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        _check_data_consistency(self,log=self.log,**kwargs)
        return self
    def check_id(self, **kwargs: Any) -> None:
        pass
    @add_doc(_check_ref)
    def check_ref(self, ref_seq: Any, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log","ref_seq"})
        self.data = _check_ref(self,ref_seq,log=self.log,**kwargs)
        # Harmonization status is already set by _check_ref
        return self
    @add_doc(_parallelize_infer_strand)
    def infer_strand(self, ref_infer: Any, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log","ref_infer"})
        self.data = _parallelize_infer_strand(self,ref_infer=ref_infer,log=self.log,**kwargs)
        # Harmonization status is already set by _parallelize_infer_strand
        return self
    
    @add_doc(_flip_allele_stats)
    def flip_allele_stats(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _flip_allele_stats(self,log=self.log,**kwargs)
        # Harmonization status is already set by _flip_allele_stats
        return self
    @add_doc(_parallelize_normalize_allele)
    def normalize_allele(self, **kwargs: Any) -> 'Sumstats':
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _parallelize_normalize_allele(self,log=self.log,**kwargs)
        return self

    def annotate_sumstats(self, **kwargs: Any) -> None:
        from gwaslab.hm.hm_assign_rsid import _annotate_sumstats
        self.data = _annotate_sumstats(self,**kwargs)
    
    @add_doc(_assign_rsid)
    def assign_rsid2(self,**kwargs):
        self.data = _assign_rsid(self,**kwargs)
        # Harmonization status is already set by _assign_rsid
        return self

    @add_doc(_infer_strand_with_annotation)
    def infer_strand2(self,**kwargs):
        self.data = _infer_strand_with_annotation(self,**kwargs)

    @add_doc(_parallelize_assign_rsid)
    def assign_rsid(self,
                    ref_rsid_tsv=None,
                    ref_rsid_vcf=None,
                    **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"ref_mode","path","log"})
        if ref_rsid_tsv is not None:
            self.data = _parallelize_assign_rsid(self,path=ref_rsid_tsv,ref_mode="tsv",log=self.log,**kwargs)
            # Metadata and harmonization status are already set by _parallelize_assign_rsid
        if ref_rsid_vcf is not None:
            self.data = _parallelize_assign_rsid(self,path=ref_rsid_vcf,ref_mode="vcf",log=self.log,**kwargs)   
            # Metadata and harmonization status are already set by _parallelize_assign_rsid
        return self
    @add_doc(_parallelize_rsid_to_chrpos)
    def rsid_to_chrpos(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _parallelize_rsid_to_chrpos(self,log=self.log,**kwargs)
        return self
    @add_doc(_parallelize_rsid_to_chrpos)
    def rsid_to_chrpos2(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _parallelize_rsid_to_chrpos(self,log=self.log,**kwargs)
        return self

    ############################################################################################################
    @add_doc(_sort_coordinate)
    def sort_coordinate(self,**sort_kwargs):
        sort_kwargs = remove_overlapping_kwargs(sort_kwargs,{"log"})
        self.data = _sort_coordinate(self,log=self.log,**sort_kwargs)
        # is_sorted is already set by _sort_coordinate
        return self

    @add_doc(_sort_column)
    def sort_column(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _sort_column(self,log=self.log,**kwargs)
        return self
    
    ############################################################################################################
    @add_doc(_fill_data)
    def fill_data(self, verbose: bool = True, **kwargs: Any) -> None:
        kwargs = remove_overlapping_kwargs(kwargs,{"log","verbose"})
        self.data = _fill_data(self, verbose=verbose, log=self.log, **kwargs)
        self.data = _sort_column(self, verbose=verbose, log=self.log)

# utilities ############################################################################################################
    # filter series ######################################################################
    @add_doc(_filter_region)
    def filter_region(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_region(new_Sumstats_object, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_region(self, **kwargs)
            return None
    
    @add_doc(_get_flanking)
    def filter_flanking(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking(new_Sumstats_object, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _get_flanking(self, **kwargs)
            return None
    
    @add_doc(_get_flanking_by_chrpos)
    def filter_flanking_by_chrpos(self, chrpos: Tuple[int, int], inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking_by_chrpos(new_Sumstats_object, chrpos, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _get_flanking_by_chrpos(self, chrpos,**kwargs)
            return None
    
    @add_doc(_get_flanking_by_id)
    def filter_flanking_by_id(self, snpid: str, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking_by_id(new_Sumstats_object, snpid, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _get_flanking_by_id(self, snpid, **kwargs)
            return None
    
    @add_doc(_filter_values)
    def filter_value(self, expr: str, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log","expr"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            _filter_values(new_Sumstats_object,expr,log=new_Sumstats_object.log, **kwargs)
            return new_Sumstats_object
        else:
            _filter_values(self, expr,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_out)
    def filter_out(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            _filter_out(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            _filter_out(self,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_in)
    def filter_in(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            _filter_in(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            _filter_in(self,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_region_in)
    def filter_region_in(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_region_in(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_region_in(self,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_region_out)
    def filter_region_out(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_region_out(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_region_out(self,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_palindromic)
    def filter_palindromic(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_palindromic(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_palindromic(self,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_snp)
    def filter_snp(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_snp(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_snp(self,log=self.log,**kwargs)
            return None
    
    @add_doc(_filter_indel)
    def filter_indel(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_indel(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_indel(self,log=self.log,**kwargs)
            return None
    
    
    @add_doc(_exclude_hla)
    def exclude_hla(self, inplace: bool = False, **kwargs: Any) -> Optional['Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _exclude_hla(new_Sumstats_object,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _exclude_hla(self,log=self.log,**kwargs)
            return None

    @add_doc(_search_variants)
    def search(self, inplace=False, **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _search_variants(new_Sumstats_object,
                                                        log=new_Sumstats_object.log,
                                                        **kwargs)
            return new_Sumstats_object
        else:
            self.data = _search_variants(self,log=self.log,**kwargs)
    
    @add_doc(_extract_ld_proxy)
    def get_proxy(self,snplist, **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","snplist"})
        proxy_variants = _extract_ld_proxy(common_sumstats = self,
                                                            snplist=snplist,
                                                            log=self.log,
                                                            **kwargs)
        return proxy_variants
    
    @add_doc(_sampling)
    def random_variants(self,inplace=False,n=1,p=None,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","n","p"})
        if inplace is True:
            self.data = _sampling(self,n=n,p=p,log=self.log,**kwargs)
        else:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _sampling(new_Sumstats_object,n=n,p=p,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
    
    @add_doc(_get_hapmap3)
    def filter_hapmap3(self, inplace=False, **kwargs ):
        kwargs = remove_overlapping_kwargs(kwargs,{"build","log"})
        if inplace is True:
            self.data = _get_hapmap3(self, build=self.build,log=self.log, **kwargs)
        else:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_hapmap3(new_Sumstats_object, build=self.build,log=self.log, **kwargs)
            return new_Sumstats_object
    
    @add_doc(_get_region_start_and_end)
    def get_region_start_and_end(self,**kwargs):
        return _get_region_start_and_end(**kwargs)

    ######################################################################
    
    def check_af(self,ref_infer,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","ref_infer"})
        self.data = _parallelize_check_af(self,ref_infer=ref_infer,log=self.log,**kwargs)
        # Metadata is already set by _parallelize_check_af
    
    def check_af2(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _check_af_with_annotation(self,**kwargs)
        # Metadata is already set by _check_af_with_annotation
    
    def infer_af(self,ref_infer,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","ref_infer"})
        self.data = _parallelize_infer_af(self,ref_infer=ref_infer,log=self.log,**kwargs)
        # Metadata is already set by _parallelize_infer_af
    
    def infer_af2(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _infer_af_with_annotation(self,**kwargs)
        # Metadata is already set by _infer_af_with_annotation
    
    def infer_eaf_from_maf(self,ref_infer,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","ref_infer"})
        self.data = _parallele_infer_af_with_maf(self,ref_infer=ref_infer,log=self.log,**kwargs)
        # Metadata is already set by _parallele_infer_af_with_maf
    
    def infer_eaf_from_maf2(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        self.data = _infer_af_with_maf_annotation(self,**kwargs)
        # Metadata is already set by _infer_af_with_maf_annotation    
    
    @suppress_display
    def plot_daf(self, **kwargs: Any) -> Tuple[Any, Any]:
        fig, outliers = plotdaf(self, **self._apply_viz_params(plotdaf, kwargs, key="plot_daf"))
        return fig, outliers
    
    @add_doc(_infer_ancestry)
    def infer_ancestry(self, **kwargs: Any) -> None:
        kwargs = remove_overlapping_kwargs(kwargs,{"build","log"})
        self.meta["gwaslab"]["inferred_ancestry"] = _infer_ancestry(self,
                                                                    build=self.build, 
                                                                    log=self.log,
                                                                    **kwargs)

    @suppress_display
    def plot_gwheatmap(self, **kwargs: Any) -> Any:
        params = self._apply_viz_params(_gwheatmap, kwargs, key="plot_gwheatmap")
        fig = _gwheatmap(self, **params)
        return fig
    
    @add_doc(_mqqplot)
    @suppress_display
    def plot_mqq(self, build: Optional[str] = None, **kwargs: Any) -> Any:
        mode = kwargs.get("mode", None)
        params = self._apply_viz_params(_mqqplot, kwargs, key="plot_mqq", mode=mode)
        if "build" not in params:
            params["build"] = self.build
        plot, log = _mqqplot(self, **params)
        return plot
    
    @add_doc(_mqqplot)
    @suppress_display
    def plot_manhattan(self, build: Optional[str] = None, **kwargs: Any) -> Any:
        params = self._apply_viz_params(_mqqplot, kwargs, key="plot_manhattan")
        if "build" not in params:
            params["build"] = self.build
        plot, log = _mqqplot(self, **params)
        return plot
    
    @add_doc(_process_density)
    @suppress_display
    def plot_snp_density(self, build: Optional[str] = None, **kwargs: Any) -> Any:
        plot, log = _mqqplot(self, **{**self._apply_viz_params(_mqqplot, kwargs, key="plot_snp_density", mode="b"), "mode": "b", "build": kwargs.get("build", self.build)})
        return plot
    
    @add_doc(_plot_qq)
    @suppress_display
    def plot_qq(self, build: Optional[str] = None, **kwargs: Any) -> Any:
        plot, log = _mqqplot(self, **{**self._apply_viz_params(_mqqplot, kwargs, key="plot_qq", mode="qq"), "mode": "qq", "build": kwargs.get("build", self.build)})
        return plot
    
    @add_doc(_plot_regional)
    @suppress_display
    def plot_region(self, build: Optional[str] = None, **kwargs: Any) -> Any:
        plot, log = _mqqplot(self, **{**self._apply_viz_params(_mqqplot, kwargs, key="plot_region", mode="r"), "mode": "r", "build": kwargs.get("build", self.build)})
        return plot

    @add_doc(_plot_trumpet)
    @suppress_display
    def plot_trumpet(self, build: Optional[str] = None, **kwargs: Any) -> Any:
        fig = _plot_trumpet(
            self,
            **{**self._apply_viz_params(_plot_trumpet, kwargs, key="plot_trumpet", mode=("b" if kwargs.get("mode") == "b" else "q")),
               "build": kwargs.get("build", self.build)}
        )
        return fig

    @add_doc(_plot_effect)
    def plot_effect(self, **kwargs: Any) -> None:
        _plot_effect(self, **self._apply_viz_params(_plot_effect, kwargs, key="plot_effect"))

    @add_doc(_get_sig)
    def get_lead(self, gls: bool = False, build: Optional[str] = None, **kwargs: Any) -> Union[pd.DataFrame, 'Sumstats']:
        kwargs = remove_overlapping_kwargs(kwargs,{"log","build"})
        output = _get_sig(self,
                        log=self.log,
                        build=self.build,
                        **kwargs)
        
        # return sumstats object    
        if gls == True:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = output
            gc.collect()
            return new_Sumstats_object
        return output
    
    @add_doc(_get_top)
    def get_top(self, gls=False, build=None, **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","build"})
        output = _get_top(self,
                        log=self.log,
                        build=self.build,
                        **kwargs)
        if gls == True:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = output
            gc.collect()
            return new_Sumstats_object
        return output
    
    @add_doc(generate_qc_report)
    def report(self, output_path="gwas_qc_report.html", **kwargs):
        # Extract specific kwargs for generate_qc_report
        basic_check_kwargs = kwargs.pop("basic_check_kwargs", None)
        harmonize_kwargs = kwargs.pop("harmonize_kwargs", None)
        get_lead_kwargs = kwargs.pop("get_lead_kwargs", None)
        mqq_plot_kwargs = kwargs.pop("mqq_plot_kwargs", None)
        regional_plot_kwargs = kwargs.pop("regional_plot_kwargs", None)
        output_kwargs = kwargs.pop("output_kwargs", None)
        report_title = kwargs.pop("report_title", "GWAS Quality Control Report")
        verbose = kwargs.pop("verbose", True)
        
        # Pass remaining kwargs to appropriate functions if needed
        if basic_check_kwargs is None:
            basic_check_kwargs = {}
        if get_lead_kwargs is None:
            get_lead_kwargs = {}
        if mqq_plot_kwargs is None:
            mqq_plot_kwargs = {}
        if regional_plot_kwargs is None:
            regional_plot_kwargs = {}
        
        # Generate the report
        return generate_qc_report(
            sumstats=self,
            output_path=output_path,
            basic_check_kwargs=basic_check_kwargs,
            harmonize_kwargs=harmonize_kwargs,
            get_lead_kwargs=get_lead_kwargs,
            mqq_plot_kwargs=mqq_plot_kwargs,
            regional_plot_kwargs=regional_plot_kwargs,
            output_kwargs=output_kwargs,
            report_title=report_title,
            verbose=verbose
        )

    @add_doc(_get_signal_density2)
    def get_density(self, sig_list=None, windowsizekb=100, **kwargs):
        self.data = _get_signal_density2(
            self,
            bwindowsizekb=windowsizekb,
            sig_sumstats=sig_list if isinstance(sig_list, pd.DataFrame) else None,
            log=self.log
        )

    @add_doc(_get_novel)
    def get_novel(self, **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","build"})
        output = _get_novel(self,
                           build=self.build,
                           log=self.log,
                           **kwargs)
        return output
    
    def get_associations(self, **kwargs):
        """
        Extract GWAS Catalog associations for variants in sumstats.
        
        Returns
        -------
        pandas.DataFrame
            Summary DataFrame containing GWAS Catalog associations for variants in sumstats.
            The full associations are stored in self.associations.
        """
        associations_full,associations_summary  = _extract_associations(self,**kwargs)
        self.associations = associations_full
        return associations_summary
    
    @add_doc(_plot_associations)
    def plot_associations(self,**kwargs):
        _plot_associations(self.associations, **self._apply_viz_params(_plot_associations, kwargs, key="plot_associations"))

    def check_cis(self, gls=False, **kwargs):
        output = _check_cis(self,
                           log=self.log,
                           **kwargs)
        
        # return sumstats object   
        if gls == True:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = output
            gc.collect()
            return new_Sumstats_object
        return output
    
    def check_novel_set(self, **kwargs):
        output = _check_novel_set(self,
                           log=self.log,
                           **kwargs)
        # return sumstats object    
        return output

    def check_cs_overlap(self, **kwargs):
        output = _check_novel_set(self.pipcs,
                           log=self.log,
                           **kwargs)
        # return sumstats object    
        return output

    def anno_gene(self, **kwargs):
        output = _anno_gene(self,
                           log=self.log,
                           **kwargs)
        return output
    
    @add_doc(_get_per_snp_r2)
    def get_per_snp_r2(self,**kwargs):
        self.data = _get_per_snp_r2(self, beta="BETA", af="EAF", n="N", log=self.log, **kwargs)
        #add data inplace
    
    @add_doc(_get_ess)
    def get_ess(self, **kwargs):
        self.data = _get_ess(self, log=self.log, **kwargs)
    
    @add_doc(_lambda_GC)
    def get_gc(self, mode=None, **kwargs):
        output = _lambda_GC(self, mode=mode, log=self.log, **kwargs)
        self.meta["Genomic inflation factor"] = output
        return output
        
    def abf_finemapping(self, region=None, chrpos=None, snpid=None,**kwargs):        
        region_data = _abf_finemapping(self.data.copy(),region=region,chrpos=chrpos,snpid=snpid,log=self.log, **kwargs)
        credible_sets = _make_cs(region_data,threshold=0.95,log=self.log)
        return region_data, credible_sets
        
    def get_ld_matrix_from_vcf(self, **kwargs):
        return _get_ld_matrix_from_vcf(self.data.copy(), log=self.log, **kwargs)
######################################################################################################
    def run_prscs(self, build=None, verbose=True, match_allele=True, how="inner", **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log"})
        if build is None:
            build = self.build
        # Create a temporary copy for hapmap3 filtering
        temp_sumstats = copy.deepcopy(self)
        temp_sumstats.data = self.data.copy()
        insumstats = _get_hapmap3(temp_sumstats, build=build, verbose=verbose , match_allele=match_allele, how=how )
        _run_prscs(sst_file = insumstats[["rsID","CHR","POS","EA","NEA","BETA","SE"]], 
                   log=self.log, 
                   **kwargs)
    
    def run_magma(self, build=None, verbose=True, **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"study","build","verbose","log"})
        _run_magma(self, 
                   study=self.meta["gwaslab"]["study_name"], 
                   build=build, verbose=verbose, log=self.log, **kwargs)
        
    def run_scdrs(self, build=None, verbose=True, **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"study","build","verbose","log"})
        _run_scdrs(self,
                   study=self.meta["gwaslab"]["study_name"], 
                   build=build, 
                   verbose=verbose, log=self.log, **kwargs)
## LDSC ##############################################################################################
    @add_doc(_estimate_h2_by_ldsc)
    def estimate_h2_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right", **kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"insumstats","meta","log","verbose"})
        if build is None:
            build = self.build
        # Create a temporary copy for hapmap3 filtering
        temp_sumstats = copy.deepcopy(self)
        temp_sumstats.data = self.data.copy()
        insumstats = _get_hapmap3(temp_sumstats, build=build, verbose=verbose , match_allele=match_allele, how=how )
        self.ldsc_h2, self.ldsc_h2_results = _estimate_h2_by_ldsc(insumstats=insumstats, 
                                                                  meta=self.meta,
                                                                  log=self.log, 
                                                                  verbose=verbose, 
                                                                  **kwargs)
        return self.ldsc_h2_results
    
    def estimate_rg_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right", get_hm3=True,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"insumstats","meta","log","verbose"})
        if build is None:
            build = self.build
        if get_hm3==True:
            # Create a temporary copy for hapmap3 filtering
            temp_sumstats = copy.deepcopy(self)
            temp_sumstats.data = self.data.copy()
            insumstats = _get_hapmap3(temp_sumstats, build=build, verbose=verbose , match_allele=match_allele, how=how )
        else:
            insumstats = self.data
        ldsc_rg = _estimate_rg_by_ldsc(insumstats=insumstats,
                                             meta=self.meta,
                                             log=self.log, 
                                             verbose=verbose, 
                                             **kwargs)
        self.ldsc_rg = pd.concat([self.ldsc_rg, ldsc_rg],ignore_index=True)
        return self.ldsc_rg
    
    @add_doc(_estimate_h2_cts_by_ldsc)
    def estimate_h2_cts_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right",**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"insumstats","log","verbose"})
        if build is None:
            build = self.build
        # Create a temporary copy for hapmap3 filtering
        temp_sumstats = copy.deepcopy(self)
        temp_sumstats.data = self.data.copy()
        insumstats = _get_hapmap3(temp_sumstats, build=build, verbose=verbose , match_allele=match_allele, how=how )
        self.ldsc_h2_cts  = _estimate_h2_cts_by_ldsc(insumstats=insumstats, 
                                                     log=self.log,
                                                       verbose=verbose, 
                                                       **kwargs)
        return self.ldsc_h2_cts
    
    @add_doc(_estimate_partitioned_h2_by_ldsc)
    def estimate_partitioned_h2_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right",**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"insumstats","meta","log","verbose"})
        if build is None:
            build = self.build
        # Create a temporary copy for hapmap3 filtering
        temp_sumstats = copy.deepcopy(self)
        temp_sumstats.data = self.data.copy()
        insumstats = _get_hapmap3(temp_sumstats, build=build, verbose=verbose , match_allele=match_allele, how=how )
        self.ldsc_partitioned_h2_summary, self.ldsc_partitioned_h2_results  = _estimate_partitioned_h2_by_ldsc(insumstats=insumstats, 
                                                                                                               meta=self.meta,
                                                                                                               log=self.log, 
                                                                                                               verbose=verbose, 
                                                                                                               **kwargs)
        return self.ldsc_partitioned_h2_results
# external ################################################################################################
    
    def calculate_ld_matrix(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"study"})
        self.finemapping["path"],self.finemapping["file"],self.finemapping["plink_log"]= _to_finemapping(self, study = self.meta["gwaslab"]["study_name"],**kwargs)
        #self.to_finemapping_file_path, self.to_finemapping_file, self.plink_log  = tofinemapping(self.data,study = self.meta["gwaslab"]["study_name"],**kwargs)
    
    def extract_ld_matrix(self,**kwargs):
        self.finemapping["path"],self.finemapping["file"],self.finemapping["plink_log"]= _to_finemapping_using_ld(self,study = self.meta["gwaslab"]["study_name"],**kwargs)

    def run_susie_rss(self,**kwargs):
        self.pipcs=_run_susie_rss(self, self.finemapping["path"], **kwargs)
        self.finemapping["pipcs"] = self.pipcs
        #self.pipcs=_run_susie_rss(self.to_finemapping_file_path,**kwargs)
    
    def get_cs_lead(self,**kwargs):
        return _get_cs_lead(self.pipcs,**kwargs)
    
    @add_doc(_clump)
    def clump(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","study"})
        self.clumps["clumps"], self.clumps["clumps_raw"], self.clumps["plink_log"] = _clump(self, 
                                                                                            log=self.log, 
                                                                                            study = self.meta["gwaslab"]["study_name"], 
                                                                                            **kwargs)

    def calculate_prs(self,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"log","study"})
        combined_results_summary = _calculate_prs(self, log=self.log, study = self.meta["gwaslab"]["study_name"], **kwargs)
        return combined_results_summary



# loading aux data
    def read_pipcs(self,prefix,**kwargs):
        kwargs = remove_overlapping_kwargs(kwargs,{"study"})
        # Create a temporary Sumstats object with just the needed columns
        temp_data = self.data[["SNPID","CHR","POS"]].copy()
        self.pipcs = _read_pipcs(temp_data,prefix, study= self.meta["gwaslab"]["study_name"], **kwargs)
        self.finemapping["pipcs"] = self.pipcs

    def plot_pipcs(self, region=None, locus=None, **kwargs):
        _plot_cs(self.pipcs, region=region, locus=locus, **self._apply_viz_params(_plot_cs, kwargs, key="plot_pipcs"))

# to_format ###############################################################################################       
    @add_doc(_to_format)
    def to_format(self, path, build=None, verbose=True, **kwargs):
        if build is None:
            build = self.build
        kwargs = remove_overlapping_kwargs(kwargs,{"log","verbose","meta","build","path"})
        _to_format(self, path, log=self.log, verbose=verbose, meta=self.meta, build=build, **kwargs)

    def to_pickle(self, path="~/mysumstats.pickle", overwrite=False):
        dump_pickle(self, path=path, overwrite=overwrite)

    @add_doc(write_gsf)
    def to_gsf(self, path, partition_cols=None, compression="zstd", verbose=True):
        write_gsf(
            self,
            path=path,
            meta=self.meta,
            partition_cols=partition_cols,
            compression=compression,
            log=self.log,
            verbose=verbose
        )

######################################################################################
    def offload(self):
        _offload(self.data, self.tmp_path, self.log)
        del self.data
        gc.collect()

    @add_doc(_view_sumstats)
    def view_sumstats(self, expr=None):
        return _view_sumstats(self, expr=expr)
    
    
    @add_doc(_reload)
    def reload(self, delete_files=None):
        self.data = _reload(self.tmp_path, self.log, delete_files=delete_files)
