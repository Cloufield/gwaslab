import pandas as pd
import numpy as np
import time
import copy
from gwaslab.g_Sumstats_summary import summarize
from gwaslab.g_Sumstats_summary import lookupstatus
from gwaslab.io_preformat_input_polars import preformatp
from gwaslab.io_to_formats import _to_format
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import fixID
from gwaslab.qc_fix_sumstats import flipSNPID
from gwaslab.qc_fix_sumstats import stripSNPID
from gwaslab.qc_fix_sumstats import removedup
from gwaslab.qc_fix_sumstats import fixchr
from gwaslab.qc_fix_sumstats import fixpos
from gwaslab.qc_fix_sumstats import fixallele
from gwaslab.qc_fix_sumstats import parallelnormalizeallele
from gwaslab.qc_fix_sumstats import sanitycheckstats
from gwaslab.qc_fix_sumstats import parallelizeliftovervariant
from gwaslab.qc_fix_sumstats import flipallelestats
from gwaslab.qc_fix_sumstats import sortcoordinate
from gwaslab.qc_fix_sumstats import sortcolumn
from gwaslab.qc_fix_sumstats import _set_build
from gwaslab.qc_fix_sumstats import _process_build
from gwaslab.hm_harmonize_sumstats import parallelecheckaf
from gwaslab.hm_harmonize_sumstats import paralleleinferaf
from gwaslab.hm_harmonize_sumstats import checkref
from gwaslab.hm_harmonize_sumstats import oldcheckref
from gwaslab.hm_harmonize_sumstats import rsidtochrpos
from gwaslab.hm_harmonize_sumstats import parallelizeassignrsid
from gwaslab.hm_harmonize_sumstats import parallelinferstrand
from gwaslab.hm_harmonize_sumstats import parallelrsidtochrpos
from gwaslab.hm_harmonize_sumstats import _paralleleinferafwithmaf
from gwaslab.util_in_filter_value import filtervalues
from gwaslab.util_in_filter_value import filterout
from gwaslab.util_in_filter_value import filterin
from gwaslab.util_in_filter_value import filterregionin
from gwaslab.util_in_filter_value import filterregionout
from gwaslab.util_in_filter_value import _filter_indel
from gwaslab.util_in_filter_value import _filter_palindromic
from gwaslab.util_in_filter_value import _filter_snp
from gwaslab.util_in_filter_value import _filter_region
from gwaslab.util_in_filter_value import _exclude_hla
from gwaslab.util_in_filter_value import _search_variants
from gwaslab.util_in_filter_value import inferbuild
from gwaslab.util_in_filter_value import sampling
from gwaslab.util_in_filter_value import _get_flanking
from gwaslab.util_in_filter_value import _get_flanking_by_chrpos
from gwaslab.util_in_filter_value import _get_flanking_by_id
from gwaslab.util_in_calculate_gc import lambdaGC
from gwaslab.util_in_convert_h2 import _get_per_snp_r2
from gwaslab.util_in_get_sig import getsig
from gwaslab.util_in_get_density import getsignaldensity
from gwaslab.util_in_get_density import assigndensity
from gwaslab.util_in_get_sig import annogene
from gwaslab.util_in_get_sig import getnovel
from gwaslab.util_in_get_sig import _check_cis
from gwaslab.util_in_get_sig import _check_novel_set
from gwaslab.util_in_fill_data import filldata
from gwaslab.bd_get_hapmap3 import gethapmap3
from gwaslab.bd_common_data import get_chr_list
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_high_ld
from gwaslab.bd_common_data import get_format_dict
from gwaslab.bd_common_data import get_formats_list
from gwaslab.g_version import _show_version
from gwaslab.g_version import gwaslab_info
from gwaslab.g_meta import _init_meta
from gwaslab.g_meta import _append_meta_record
from gwaslab.g_meta_update import _update_meta
from gwaslab.util_ex_run_clumping import _clump
from gwaslab.util_ex_calculate_ldmatrix import tofinemapping
from gwaslab.io_load_ld import tofinemapping_using_ld
from gwaslab.util_ex_calculate_prs import _calculate_prs
from gwaslab.viz_plot_mqqplot import mqqplot
from gwaslab.viz_plot_trumpetplot import plottrumpet
from gwaslab.viz_plot_compare_af import plotdaf
from gwaslab.util_ex_run_susie import _run_susie_rss
from gwaslab.util_ex_run_susie import _get_cs_lead
from gwaslab.qc_fix_sumstats import _check_data_consistency
from gwaslab.util_ex_ldsc import _estimate_h2_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_rg_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_h2_cts_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_partitioned_h2_by_ldsc 
from gwaslab.util_ex_ldproxyfinder import _extract_ld_proxy
from gwaslab.bd_get_hapmap3 import gethapmap3
from gwaslab.util_abf_finemapping import abf_finemapping
from gwaslab.util_abf_finemapping import make_cs
from gwaslab.io_read_pipcs import _read_pipcs
from gwaslab.util_in_estimate_ess import _get_ess
from gwaslab.viz_plot_credible_sets import _plot_cs
from gwaslab.hm_casting import _align_with_mold
from gwaslab.hm_casting  import _merge_mold_with_sumstats_by_chrpos
import gc
from gwaslab.viz_plot_phe_heatmap import _gwheatmap
from gwaslab.util_ex_run_prscs import _run_prscs

#20220309
class Sumstatsp():
    def __init__(self,
             sumstats,
             fmt=None,
             tab_fmt="tsv",
             snpid=None,
             rsid=None,
             chrom=None,
             pos=None,
             ea=None,
             nea=None,
             ref=None,
             alt=None,
             eaf=None,
             neaf=None,
             maf=None,
             n=None,
             beta=None,
             se=None,
             chisq=None,
             z=None,
             f=None,
             t=None,
             p=None,
             q=None,
             mlog10p=None,
             test=None,
             info=None,
             OR=None,
             OR_95L=None,
             OR_95U=None,
             beta_95L=None,
             beta_95U=None,
             HR=None,
             HR_95L=None,
             HR_95U=None,
             ncase=None,
             ncontrol=None,
             neff=None,
             i2=None,
             phet=None,
             dof=None,
             snpr2=None,
             status=None,
             other=[],
             chrom_pat=None,
             snpid_pat=None,
             usekeys=None,
             direction=None,
             verbose=True,
             study="Study_1",
             trait="Trait_1",
             build="99",
             species="homo sapiens",
             build_infer=False,
             **readargs):

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