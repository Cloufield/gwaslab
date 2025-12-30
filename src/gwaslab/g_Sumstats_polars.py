import pandas as pd
import numpy as np
import time
import copy
from typing import TYPE_CHECKING, Optional, List, Any, Union

if TYPE_CHECKING:
    from gwaslab.info.g_Log import Log

from gwaslab.info.g_Sumstats_summary import summarize
from gwaslab.info.g_Sumstats_summary import lookupstatus
from gwaslab.io.io_preformat_input_polars import preformatp
from gwaslab.io.io_to_formats import _to_format
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_fix_sumstats import _fix_ID
from gwaslab.qc.qc_fix_sumstats import _flip_SNPID
from gwaslab.qc.qc_fix_sumstats import _strip_SNPID
from gwaslab.qc.qc_fix_sumstats import _remove_dup
from gwaslab.qc.qc_fix_sumstats import _fix_chr
from gwaslab.qc.qc_fix_sumstats import _fix_pos
from gwaslab.qc.qc_fix_sumstats import _fix_allele
from gwaslab.qc.qc_fix_sumstats import _parallelize_normalize_allele
from gwaslab.qc.qc_sanity_check import _sanity_check_stats

from gwaslab.hm.hm_liftover_v2 import _liftover_variant
from gwaslab.qc.qc_fix_sumstats import _flip_allele_stats
from gwaslab.qc.qc_fix_sumstats import _sort_coordinate
from gwaslab.qc.qc_fix_sumstats import _sort_column
from gwaslab.qc.qc_build import _set_build
from gwaslab.qc.qc_fix_sumstats import _process_build
from gwaslab.hm.hm_harmonize_sumstats import _parallelize_check_af
from gwaslab.hm.hm_harmonize_sumstats import _parallelize_infer_af
from gwaslab.hm.hm_harmonize_sumstats import _check_ref
from gwaslab.hm.hm_harmonize_sumstats import _parallelize_assign_rsid
from gwaslab.hm.hm_harmonize_sumstats import _parallelize_infer_strand
from gwaslab.hm.hm_harmonize_sumstats import _parallelize_rsid_to_chrpos
from gwaslab.hm.hm_harmonize_sumstats import _parallele_infer_af_with_maf
from gwaslab.util.util_in_filter_value import _filter_values
from gwaslab.util.util_in_filter_value import _filter_out
from gwaslab.util.util_in_filter_value import _filter_in
from gwaslab.util.util_in_filter_value import _filter_region_in
from gwaslab.util.util_in_filter_value import _filter_region_out
from gwaslab.util.util_in_filter_value import _filter_indel
from gwaslab.util.util_in_filter_value import _filter_palindromic
from gwaslab.util.util_in_filter_value import _filter_snp
from gwaslab.util.util_in_filter_value import _filter_region
from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.util.util_in_filter_value import _search_variants
from gwaslab.util.util_in_filter_value import _infer_build
from gwaslab.util.util_in_filter_value import _sampling
from gwaslab.util.util_in_filter_value import _get_flanking
from gwaslab.util.util_in_filter_value import _get_flanking_by_chrpos
from gwaslab.util.util_in_filter_value import _get_flanking_by_id
from gwaslab.util.util_in_calculate_gc import _lambda_GC
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.util.util_in_get_density import _get_signal_density2
from gwaslab.util.util_in_get_sig import _anno_gene
from gwaslab.util.util_in_get_sig import _get_novel
from gwaslab.util.util_in_get_sig import _check_cis
from gwaslab.util.util_in_get_sig import _check_novel_set
from gwaslab.util.util_in_fill_data import _fill_data
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap3
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_high_ld
from gwaslab.bd.bd_common_data import get_format_dict
from gwaslab.bd.bd_common_data import get_formats_list
from gwaslab.info.g_version import _show_version
from gwaslab.info.g_version import gwaslab_info
from gwaslab.info.g_meta import _init_meta
from gwaslab.info.g_meta import _append_meta_record
from gwaslab.info.g_meta import _update_meta
from gwaslab.util.util_ex_run_clumping import _clump
from gwaslab.util.util_ex_calculate_ldmatrix import _to_finemapping
from gwaslab.io.io_load_ld import _to_finemapping_using_ld
from gwaslab.util.util_ex_calculate_prs import _calculate_prs
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.viz.viz_plot_trumpetplot import _plot_trumpet
from gwaslab.viz.viz_plot_compare_af import plotdaf
from gwaslab.util.util_ex_run_susie import _run_susie_rss
from gwaslab.util.util_ex_run_susie import _get_cs_lead
from gwaslab.qc.qc_sanity_check import _check_data_consistency
from gwaslab.util.util_ex_ldsc import _estimate_h2_by_ldsc
from gwaslab.util.util_ex_ldsc import _estimate_rg_by_ldsc
from gwaslab.util.util_ex_ldsc import _estimate_h2_cts_by_ldsc
from gwaslab.util.util_ex_ldsc import _estimate_partitioned_h2_by_ldsc 
from gwaslab.util.util_ex_ldproxyfinder import _extract_ld_proxy
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap3
from gwaslab.util.util_abf_finemapping import _abf_finemapping
from gwaslab.util.util_abf_finemapping import _make_cs
from gwaslab.io.io_read_pipcs import _read_pipcs
from gwaslab.util.util_in_estimate_ess import _get_ess
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.hm.hm_casting import _align_with_mold
from gwaslab.hm.hm_casting  import _merge_mold_with_sumstats_by_chrpos
import gc
from gwaslab.viz.viz_plot_phe_heatmap import _gwheatmap
from gwaslab.util.util_ex_run_prscs import _run_prscs

#20220309
class Sumstatsp():
    def __init__(
        self,
        sumstats: Union[str, pd.DataFrame, Any],
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
        other: List[str] = [],
        chrom_pat: Optional[str] = None,
        snpid_pat: Optional[str] = None,
        usekeys: Optional[List[str]] = None,
        direction: Optional[str] = None,
        verbose: bool = True,
        study: str = "Study_1",
        trait: str = "Trait_1",
        build: str = "99",
        species: str = "homo sapiens",
        build_infer: bool = False,
        **readargs: Any
    ) -> None:

        # basic attributes
        self.data = pd.DataFrame()
        self.log = Log()
        self.ldsc_h2 = None
        self.ldsc_h2_results = None
        self.ldsc_rg = pd.DataFrame()
        self.ldsc_h2_cts = None
        self.ldsc_partitioned_h2_summary = None
        self.ldsc_partitioned_h2_results = None
        # meta information
        self.meta = _init_meta() 
        self.build = build
        self.meta["gwaslab"]["study_name"] =  study
        self.meta["gwaslab"]["species"] = species
        
        # initialize attributes for clumping and finmapping
        #self.to_finemapping_file_path = ""
        #self.to_finemapping_file  = pd.DataFrame()
        #self.plink_log = ""

        # path / file / plink_log
        self.finemapping = dict()

        # clumps / clumps_raw / plink_log
        self.clumps = dict()
        
        #
        self.pipcs = pd.DataFrame()

        # print gwaslab version information
        _show_version(self.log, verbose=verbose)

        #preformat the data
        self.data  = preformatp(
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
          usekeys=usekeys,
          chrom_pat=chrom_pat,
          snpid_pat=snpid_pat,
          verbose=verbose,
          readargs=readargs,
          log=self.log)

        gc.collect()   