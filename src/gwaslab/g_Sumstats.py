import pandas as pd
import numpy as np
import time
import copy
import gc
from gwaslab.g_Sumstats_summary import summarize
from gwaslab.g_Sumstats_summary import lookupstatus
from gwaslab.io_preformat_input import preformat
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
from gwaslab.util_ex_phewwas import _extract_associations
from gwaslab.util_in_fill_data import filldata
from gwaslab.bd_get_hapmap3 import gethapmap3
from gwaslab.bd_common_data import get_chr_list
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_high_ld
from gwaslab.bd_common_data import get_format_dict
from gwaslab.bd_common_data import get_formats_list
from gwaslab.bd_path_manager import _path
from gwaslab.g_version import _show_version
from gwaslab.g_version import gwaslab_info
from gwaslab.g_meta import _init_meta
from gwaslab.g_meta import _append_meta_record
from gwaslab.g_meta_update import _update_meta
from gwaslab.io_load_ld import tofinemapping_using_ld
from gwaslab.qc_fix_sumstats import _check_data_consistency
from gwaslab.bd_get_hapmap3 import gethapmap3
from gwaslab.util_abf_finemapping import abf_finemapping
from gwaslab.util_abf_finemapping import make_cs
from gwaslab.io_read_pipcs import _read_pipcs
from gwaslab.util_in_estimate_ess import _get_ess
from gwaslab.hm_casting import _align_with_mold
from gwaslab.hm_casting  import _merge_mold_with_sumstats_by_chrpos
from gwaslab.viz_plot_phe_heatmap import _gwheatmap
from gwaslab.viz_plot_mqqplot import mqqplot
from gwaslab.viz_plot_trumpetplot import plottrumpet
from gwaslab.viz_plot_compare_af import plotdaf
from gwaslab.viz_plot_credible_sets import _plot_cs
from gwaslab.util_ex_run_clumping import _clump
from gwaslab.util_ex_calculate_ldmatrix import tofinemapping
from gwaslab.util_ex_run_susie import _run_susie_rss
from gwaslab.util_ex_run_susie import _get_cs_lead
from gwaslab.util_ex_ldsc import _estimate_h2_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_rg_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_h2_cts_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_partitioned_h2_by_ldsc 
from gwaslab.util_ex_ldproxyfinder import _extract_ld_proxy
from gwaslab.util_ex_run_magma import _run_magma
from gwaslab.util_ex_infer_ancestry import _infer_ancestry
from gwaslab.util_ex_run_prscs import _run_prscs
from gwaslab.util_ex_calculate_prs import _calculate_prs
from gwaslab.util_ex_run_scdrs import _run_scdrs

from gwaslab.io_to_pickle import _offload
from gwaslab.io_to_pickle import _reload

#20220309
class Sumstats():
    def __init__(self,
             sumstats=None,
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
             other=None,
             exclude=None,
             include=None,
             chrom_pat=None,
             snpid_pat=None,
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
        
        self.id = id(self)
        self.tmp_path = _path(pid=self.id)
        
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
        self.data  = preformat(
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
          readargs=readargs,
          log=self.log)

        # checking genome build
        self.meta["gwaslab"]["genome_build"] = _process_build(build, log=self.log, verbose=False)

        # if build is unknown and build_infer is True, infer the build
        if species=="homo sapiens" and self.meta["gwaslab"]["genome_build"]=="99" and build_infer is True:
            try:
                self.infer_build()
            except:
                pass
        gc.collect()   

#### healper #################################################################################
    def update_meta(self, **kwargs):
        self.meta = _update_meta(self.meta, self.data,log = self.log, **kwargs)

    def summary(self):
        return summarize(self.data)

    def lookup_status(self,status="STATUS"):
        return lookupstatus(self.data[status])
    
    def set_build(self, build, verbose=True):
        self.data, self.meta["gwaslab"]["genome_build"] = _set_build(self.data, build=build, log=self.log,verbose=verbose)
        gc.collect()

    def infer_build(self,verbose=True,**kwargs):
        self.data, self.meta["gwaslab"]["genome_build"] = inferbuild(self.data,log=self.log,verbose=verbose,**kwargs)
    
    def liftover(self,to_build, from_build=None,**kwargs):
        if from_build is None:
            if self.meta["gwaslab"]["genome_build"]=="99":
                self.data, self.meta["gwaslab"]["genome_build"] = inferbuild(self.data,**kwargs)
            from_build = self.meta["gwaslab"]["genome_build"]
        self.data = parallelizeliftovervariant(self.data,from_build=from_build, to_build=to_build, log=self.log,**kwargs)
        self.meta["is_sorted"] = False
        self.meta["is_harmonised"] = False
        self.meta["gwaslab"]["genome_build"]=to_build

# QC ######################################################################################
    #clean the sumstats with one line
    def basic_check(self,
                    remove=False,
                    remove_dup=False,
                    n_cores=1,
                    fixid_args={},
                    removedup_args={},
                    fixchr_args={},
                    fixpos_args={},
                    fixallele_args={},
                    sanitycheckstats_args={},
                    consistencycheck_args={},
                    normalize=True,
                    normalizeallele_args={},
                    verbose=True):
        ###############################################
        # try to fix data without dropping any information
        self.data = fixID(self.data,log=self.log,verbose=verbose, **fixid_args)
        self.data = fixchr(self.data,log=self.log,remove=remove,verbose=verbose,**fixchr_args)
        self.data = fixpos(self.data,log=self.log,remove=remove,verbose=verbose,**fixpos_args)
        self.data = fixallele(self.data,log=self.log,remove=remove,verbose=verbose,**fixallele_args)
        self.data = sanitycheckstats(self.data,log=self.log,verbose=verbose,**sanitycheckstats_args)
        _check_data_consistency(self.data,log=self.log,verbose=verbose,**consistencycheck_args)
        
        if normalize is True:
            self.data = parallelnormalizeallele(self.data,n_cores=n_cores,verbose=verbose,log=self.log,**normalizeallele_args)
        if remove_dup is True:
            self.data = removedup(self.data,log=self.log,verbose=verbose,**removedup_args)
        self.data = sortcoordinate(self.data,verbose=verbose,log=self.log)
        self.data = sortcolumn(self.data,verbose=verbose,log=self.log)
        self.meta["is_sorted"] = True
        ###############################################
        
    
    def harmonize(self,
              basic_check=True,
              ref_seq=None,
              ref_rsid_tsv=None,
              ref_rsid_vcf=None,
              ref_infer=None,
              ref_alt_freq=None,
              maf_threshold=0.40,
              ref_seq_mode="v",
              n_cores=1,
              remove=False,
              checkref_args={},
              removedup_args={},
              assignrsid_args={},
              inferstrand_args={},
              flipallelestats_args={},
              liftover_args={},
              fixid_args={},
              fixchr_args={},
              fixpos_args={},
              fixallele_args={},
              sanitycheckstats_args={},
              normalizeallele_args={}
              ):
        
        #Standard pipeline
        ####################################################
        #part 1 : basic_check
        #    1.1 fix ID
        #    1.2 remove duplication
        #    1.3 standardization : CHR POS EA NEA
        #    1.4 normalization : EA NEA
        #    1.5 sanity check : BETA SE OR EAF N OR_95L OR_95H
        #    1.6 sorting genomic coordinates and column order 
        if basic_check is True:
            
            self.data = fixID(self.data,log=self.log,**fixid_args)
            
            self.data = fixchr(self.data,remove=remove,log=self.log,**fixchr_args)
            
            self.data = fixpos(self.data,remove=remove,log=self.log,**fixpos_args)
            
            self.data = fixallele(self.data,log=self.log,**fixallele_args)
            
            self.data = sanitycheckstats(self.data,log=self.log,**sanitycheckstats_args)
            
            self.data = parallelnormalizeallele(self.data,log=self.log,n_cores=n_cores,**normalizeallele_args)
            
            self.data = sortcolumn(self.data,log=self.log)

            self.data = sortcoordinate(self.data,log=self.log)
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

            if ref_seq_mode=="v":
                self.data = checkref(self.data,ref_seq,log=self.log,**checkref_args)
            elif ref_seq_mode=="s":
                self.data = oldcheckref(self.data,ref_seq,log=self.log,**checkref_args)
            else:
                raise ValueError("ref_seq_mode should be 'v' (vectorized, faster) or 's' (sequential, slower)")
            
            self.meta["gwaslab"]["references"]["ref_seq"] = ref_seq
            
            self.data = flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
            gc.collect()
            
        if ref_infer is not None: 
            
            self.data= parallelinferstrand(self.data,ref_infer = ref_infer,ref_alt_freq=ref_alt_freq,maf_threshold=maf_threshold,
                                              n_cores=n_cores,log=self.log,**inferstrand_args)

            self.meta["gwaslab"]["references"]["ref_infer"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer"] , ref_infer)
            
            self.data =flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
            gc.collect()
            
        if (ref_seq is not None or ref_infer is not None) and (ref_rsid_tsv is not None or ref_rsid_vcf is not None):

            self.data = fixID(self.data, log=self.log, **{"fixid":True, "fixsep":True, "overwrite":True})

            gc.collect()
        
        #####################################################
        if ref_rsid_tsv is not None:
            
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_tsv,ref_mode="tsv",
                                                 n_cores=n_cores,log=self.log,**assignrsid_args)
            
            

            self.meta["gwaslab"]["references"]["ref_rsid_tsv"] = ref_rsid_tsv
            gc.collect()

        if ref_rsid_vcf is not None:
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_vcf,ref_mode="vcf",
                                                 n_cores=n_cores,log=self.log,**assignrsid_args)   

            self.meta["gwaslab"]["references"]["ref_rsid_vcf"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_rsid_vcf"] , ref_rsid_vcf)
            
            gc.collect()
        ######################################################    
        if remove is True:
            
            self.data = removedup(self.data,log=self.log,**removedup_args)
        ################################################ 
        
        self.data = sortcoordinate(self.data,log=self.log)
        
        self.data = sortcolumn(self.data,log=self.log)
        gc.collect()
        self.meta["is_sorted"] = True
        self.meta["is_harmonised"] = True
        return self
    
    def align_with_template(self, template, **kwargs):
        ## merge
        molded_sumstats, sumstats1 = _merge_mold_with_sumstats_by_chrpos(mold=template, 
                                            sumstats=self.data,
                                            log=self.log,
                                            suffixes=("_MOLD",""),
                                            return_not_matched_mold = True)
        ## align
        aligned_data = _align_with_mold(molded_sumstats)

        ## flip
        self.data =flipallelestats(aligned_data, log=self.log)
        
    ############################################################################################################
    #customizable API to build your own QC pipeline
    def fix_id(self,**kwargs):
        self.data = fixID(self.data,log=self.log,**kwargs)
    def flip_snpid(self,**kwargs):
        self.data = flipSNPID(self.data,log=self.log,**kwargs)
    def strip_snpid(self,**kwargs):
        self.data = stripSNPID(self.data,log=self.log,**kwargs)
    def fix_chr(self,**kwargs):
        self.data = fixchr(self.data,log=self.log,**kwargs)
    def fix_pos(self,**kwargs):
        self.data = fixpos(self.data,log=self.log,**kwargs)
    def fix_allele(self,**kwargs):
        self.data = fixallele(self.data,log=self.log,**kwargs)
    def remove_dup(self,**kwargs):
        self.data = removedup(self.data,log=self.log,**kwargs)
    def check_sanity(self,**kwargs):
        self.data = sanitycheckstats(self.data,log=self.log,**kwargs)
    def check_data_consistency(self, **kwargs):
        _check_data_consistency(self.data,log=self.log,**kwargs)
    def check_id(self,**kwargs):
        pass
    def check_ref(self,ref_seq,ref_seq_mode="v",**kwargs):
        if ref_seq_mode=="v":
            self.meta["gwaslab"]["references"]["ref_seq"] = ref_seq
            self.data = checkref(self.data,ref_seq,log=self.log,**kwargs)
        elif ref_seq_mode=="s":
            self.meta["gwaslab"]["references"]["ref_seq"] = ref_seq
            self.data = oldcheckref(self.data,ref_seq,log=self.log,**kwargs)            
    def infer_strand(self,ref_infer,**kwargs):
        self.meta["gwaslab"]["references"]["ref_infer"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer"] , ref_infer)
        self.data = parallelinferstrand(self.data,ref_infer=ref_infer,log=self.log,**kwargs)
    def flip_allele_stats(self,**kwargs):
        self.data = flipallelestats(self.data,log=self.log,**kwargs)
    def normalize_allele(self,**kwargs):
        self.data = parallelnormalizeallele(self.data,log=self.log,**kwargs)
    def assign_rsid(self,
                    ref_rsid_tsv=None,
                    ref_rsid_vcf=None,
                    **kwargs):
        if ref_rsid_tsv is not None:
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_tsv,ref_mode="tsv",log=self.log,**kwargs)
            self.meta["gwaslab"]["references"]["ref_rsid_tsv"] = ref_rsid_tsv
        if ref_rsid_vcf is not None:
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_vcf,ref_mode="vcf",log=self.log,**kwargs)   
            self.meta["gwaslab"]["references"]["ref_rsid_vcf"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_rsid_vcf"] , ref_rsid_vcf)
    def rsid_to_chrpos(self,**kwargs):
        self.data = rsidtochrpos(self.data,log=self.log,**kwargs)
    def rsid_to_chrpos2(self,**kwargs):
        self.data = parallelrsidtochrpos(self.data,log=self.log,**kwargs)

    ############################################################################################################
    
    def sort_coordinate(self,**sort_args):
        self.data = sortcoordinate(self.data,log=self.log,**sort_args)
        self.meta["is_sorted"] = True
    def sort_column(self,**kwargs):
        self.data = sortcolumn(self.data,log=self.log,**kwargs)
    
    ############################################################################################################
    def fill_data(self, verbose=True, **kwargs):
        self.data = filldata(self.data, verbose=verbose, log=self.log, **kwargs)
        self.data = sortcolumn(self.data, verbose=verbose, log=self.log)

# utilities ############################################################################################################
    # filter series ######################################################################
    
    def filter_region(self, inplace=False,**kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_region(new_Sumstats_object.data, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_region(self.data, **kwargs)

    def filter_flanking(self, inplace=False,**kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking(new_Sumstats_object.data, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _get_flanking(self.data, **kwargs)
    def filter_flanking_by_chrpos(self, chrpos,  inplace=False,**kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking_by_chrpos(new_Sumstats_object.data, chrpos, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _get_flanking_by_chrpos(self.data, chrpos,**kwargs)
    def filter_flanking_by_id(self, snpid, inplace=False,**kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking_by_id(new_Sumstats_object.data, snpid, **kwargs)
            return new_Sumstats_object
        else:
            self.data = _get_flanking_by_id(self.data, snpid, **kwargs)
    def filter_value(self, expr, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filtervalues(new_Sumstats_object.data,expr,log=new_Sumstats_object.log, **kwargs)
            return new_Sumstats_object
        else:
            self.data = filtervalues(self.data, expr,log=self.log,**kwargs)
    def filter_out(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterout(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = filterout(self.data,log=self.log,**kwargs)
    def filter_in(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterin(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = filterin(self.data,log=self.log,**kwargs)
    def filter_region_in(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterregionin(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = filterregionin(self.data,log=self.log,**kwargs)
    def filter_region_out(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterregionout(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = filterregionout(self.data,log=self.log,**kwargs)
    def filter_palindromic(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_palindromic(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_palindromic(self.data,log=self.log,**kwargs)
    def filter_snp(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_snp(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_snp(self.data,log=self.log,**kwargs)
    def filter_indel(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_indel(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _filter_indel(self.data,log=self.log,**kwargs)
    
    def exclude_hla(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _exclude_hla(new_Sumstats_object.data,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        else:
            self.data = _exclude_hla(self.data,log=self.log,**kwargs)

    def search(self, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _search_variants(new_Sumstats_object.data,
                                                        log=new_Sumstats_object.log,
                                                        **kwargs)
            return new_Sumstats_object
        else:
            self.data = _search_variants(self.data,log=self.log,**kwargs)
    
    def get_proxy(self,snplist, inplace=False, **kwargs):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _extract_ld_proxy(common_sumstats = new_Sumstats_object.data,
                                                         snplist=snplist,
                                                        log=new_Sumstats_object.log,
                                                        **kwargs)
            return new_Sumstats_object

    def random_variants(self,inplace=False,n=1,p=None,**kwargs):
        if inplace is True:
            self.data = sampling(self.data,n=n,p=p,log=self.log,**kwargs)
        else:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = sampling(new_Sumstats_object.data,n=n,p=p,log=new_Sumstats_object.log,**kwargs)
            return new_Sumstats_object
        
    def filter_hapmap3(self, inplace=False, build=None, **kwargs ):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        if inplace is True:
            self.data = gethapmap3(self.data, build=build,log=self.log, **kwargs)
        else:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = gethapmap3(new_Sumstats_object.data, build=build,log=self.log, **kwargs)
            return new_Sumstats_object
    ######################################################################
    
    def check_af(self,ref_infer,**kwargs):
        self.data = parallelecheckaf(self.data,ref_infer=ref_infer,log=self.log,**kwargs)
        self.meta["gwaslab"]["references"]["ref_infer_daf"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer_daf"] , ref_infer)
    def infer_af(self,ref_infer,**kwargs):
        self.data = paralleleinferaf(self.data,ref_infer=ref_infer,log=self.log,**kwargs)
        self.meta["gwaslab"]["references"]["ref_infer_af"] = ref_infer
        self.meta["gwaslab"]["references"]["ref_infer_af"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer_af"] , ref_infer)
    def maf_to_eaf(self,ref_infer,**kwargs):
        self.data = _paralleleinferafwithmaf(self.data,ref_infer=ref_infer,log=self.log,**kwargs)
        self.meta["gwaslab"]["references"]["ref_infer_maf"] = ref_infer
        self.meta["gwaslab"]["references"]["ref_infer_maf"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer_af"] , ref_infer)    
    
    def plot_daf(self, **kwargs):
        fig,outliers = plotdaf(self.data, **kwargs)
        return fig, outliers
    
    def infer_ancestry(self, **kwargs):
        self.meta["gwaslab"]["inferred_ancestry"] = _infer_ancestry(self.data, **kwargs)

    def plot_gwheatmap(self, **kwargs):
        fig = _gwheatmap(self.data, **kwargs)
        return fig
    
    def plot_mqq(self, build=None, **kwargs):

        chrom="CHR"
        pos="POS"
        p="P"
        
        if "SNPID" in self.data.columns:
            snpid="SNPID"
        elif "rsID" in self.data.columns:
            snpid="rsID"
        
        if "EAF" in self.data.columns:
            eaf="EAF"
        else:
            eaf=None

        # extract build information from meta data
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]

        plot = mqqplot(self.data,
                       snpid=snpid, 
                       chrom=chrom, 
                       pos=pos, 
                       p=p, 
                       eaf=eaf,
                       build = build, 
                       **kwargs)
        
        return plot
    
    def plot_trumpet(self, build=None, **kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        fig = plottrumpet(self.data,build = build,  **kwargs)
        return fig

    def get_lead(self, build=None, gls=False, **kwargs):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        
        # extract build information from meta data
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        
        output = getsig(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
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

    def get_density(self, sig_list=None, windowsizekb=100,**kwargs):
        
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        
        if sig_list is None:
            self.data["DENSITY"] = getsignaldensity(self.data,
                                                    id=id_to_use,
                                                    chrom="CHR",
                                                    pos="POS",
                                                    bwindowsizekb=windowsizekb,
                                                    log=self.log)
        else:
            if isinstance(sig_list, pd.DataFrame):
                self.data["DENSITY"] = assigndensity(self.data,
                                                    sig_list,
                                                    id=id_to_use, 
                                                    chrom="CHR", 
                                                    pos="POS", 
                                                    bwindowsizekb=windowsizekb,
                                                    log=self.log)

        
    def get_novel(self, **kwargs):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = getnovel(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           build=self.meta["gwaslab"]["genome_build"],
                           log=self.log,
                           **kwargs)
        # return sumstats object    
        return output
    
    def get_associations(self, **kwargs):
        associations_full,associations_summary  = _extract_associations(self.data,**kwargs)
        
        self.associations = associations_full
        
        return associations_summary

    def check_cis(self, gls=False, **kwargs):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = _check_cis(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
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
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = _check_novel_set(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           log=self.log,
                           **kwargs)
        # return sumstats object    
        return output

    def check_cs_overlap(self, **kwargs):
        if "SNPID" in self.pipcs.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = _check_novel_set(self.pipcs,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           log=self.log,
                           **kwargs)
        # return sumstats object    
        return output

    def anno_gene(self, **kwargs):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = annogene(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           log=self.log,
                           **kwargs)
        return output
        
    def get_per_snp_r2(self,**kwargs):
        self.data = _get_per_snp_r2(self.data, beta="BETA", af="EAF", n="N", log=self.log, **kwargs)
        #add data inplace
        
    def get_ess(self, **kwargs):
        self.data = _get_ess(self.data, log=self.log, **kwargs)

    def get_gc(self, mode=None, **kwargs):
        if mode is None:
            if "P" in self.data.columns:
                output = lambdaGC(self.data[["CHR","P"]],mode="P",**kwargs)
            elif "Z" in self.data.columns:
                output = lambdaGC(self.data[["CHR","Z"]],mode="Z",**kwargs)
            elif "CHISQ" in self.data.columns:
                output = lambdaGC(self.data[["CHR","CHISQ"]],mode="CHISQ",**kwargs)
            elif "MLOG10P" in self.data.columns:
                output = lambdaGC(self.data[["CHR","MLOG10P"]],mode="MLOG10P",**kwargs)
            
            #return scalar
            self.meta["Genomic inflation factor"] = output
            return output    
        else:
            output = lambdaGC(self.data[["CHR",mode]],mode=mode,**kwargs)
            self.meta["Genomic inflation factor"] = output
            return output
        
    def abf_finemapping(self, region=None, chrpos=None, snpid=None,**kwargs):        
        region_data = abf_finemapping(self.data.copy(),region=region,chrpos=chrpos,snpid=snpid,log=self.log, **kwargs)
        credible_sets = make_cs(region_data,threshold=0.95,log=self.log)
        return region_data, credible_sets
    
######################################################################################################
    def run_prscs(self, build=None, verbose=True, match_allele=True, how="inner", **kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=match_allele, how=how )
        _run_prscs(sst_file = insumstats[["rsID","CHR","POS","EA","NEA","BETA","SE"]], 
                   log=self.log, 
                   **kwargs)
    
    def run_magma(self, build=None, verbose=True, **kwargs):
        _run_magma(self, 
                   study=self.meta["gwaslab"]["study_name"], 
                   build=build, verbose=verbose, log=self.log, **kwargs)
        
    def run_scdrs(self, build=None, verbose=True, **kwargs):
        _run_scdrs(self,
                   study=self.meta["gwaslab"]["study_name"], 
                   build=build, 
                   verbose=verbose, log=self.log, **kwargs)
## LDSC ##############################################################################################
    def estimate_h2_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right", **kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=match_allele, how=how )
        self.ldsc_h2, self.ldsc_h2_results = _estimate_h2_by_ldsc(insumstats=insumstats, 
                                                                  meta=self.meta,
                                                                  log=self.log, 
                                                                  verbose=verbose, 
                                                                  **kwargs)
    
    def estimate_rg_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right", get_hm3=True,**kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        if get_hm3==True:
            insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=match_allele, how=how )
        else:
            insumstats = self.data
        ldsc_rg = _estimate_rg_by_ldsc(insumstats=insumstats,
                                             meta=self.meta,
                                             log=self.log, 
                                             verbose=verbose, 
                                             **kwargs)
        self.ldsc_rg = pd.concat([self.ldsc_rg, ldsc_rg],ignore_index=True)

    def estimate_h2_cts_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right",**kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=match_allele, how=how )
        self.ldsc_h2_cts  = _estimate_h2_cts_by_ldsc(insumstats=insumstats, 
                                                     log=self.log,
                                                       verbose=verbose, 
                                                       **kwargs)
    
    def estimate_partitioned_h2_by_ldsc(self, build=None, verbose=True, match_allele=True, how="right",**kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=match_allele, how=how )
        self.ldsc_partitioned_h2_summary, self.ldsc_partitioned_h2_results  = _estimate_partitioned_h2_by_ldsc(insumstats=insumstats, 
                                                                                                               meta=self.meta,
                                                                                                               log=self.log, 
                                                                                                               verbose=verbose, 
                                                                                                               **kwargs)
# external ################################################################################################
    
    def calculate_ld_matrix(self,**kwargs):
        self.finemapping["path"],self.finemapping["file"],self.finemapping["plink_log"]= tofinemapping(self, study = self.meta["gwaslab"]["study_name"],**kwargs)
        #self.to_finemapping_file_path, self.to_finemapping_file, self.plink_log  = tofinemapping(self.data,study = self.meta["gwaslab"]["study_name"],**kwargs)
    def extract_ld_matrix(self,**kwargs):
        self.finemapping["path"],self.finemapping["file"],self.finemapping["plink_log"]= tofinemapping_using_ld(self.data,study = self.meta["gwaslab"]["study_name"],**kwargs)

    def run_susie_rss(self,**kwargs):
        self.pipcs=_run_susie_rss(self, self.finemapping["path"], **kwargs)
        self.finemapping["pipcs"] = self.pipcs
        #self.pipcs=_run_susie_rss(self.to_finemapping_file_path,**kwargs)
    
    def get_cs_lead(self,**kwargs):
        return _get_cs_lead(self.pipcs,**kwargs)
    
    def clump(self,**kwargs):
        self.clumps["clumps"], self.clumps["clumps_raw"], self.clumps["plink_log"] = _clump(self, 
                                                                                            log=self.log, 
                                                                                            study = self.meta["gwaslab"]["study_name"], 
                                                                                            **kwargs)

    def calculate_prs(self,**kwargs):
        combined_results_summary = _calculate_prs(self.data, log=self.log, study = self.meta["gwaslab"]["study_name"], **kwargs)
        return combined_results_summary

# loading aux data
    def read_pipcs(self,prefix,**kwargs):
        self.pipcs = _read_pipcs(self.data[["SNPID","CHR","POS"]],prefix, study= self.meta["gwaslab"]["study_name"], **kwargs)
        self.finemapping["pipcs"] = self.pipcs

    def plot_pipcs(self, region=None, locus=None, **kwargs):
        _plot_cs(self.pipcs, region=region,locus=locus, **kwargs)
# to_format ###############################################################################################       

    def to_format(self, path, build=None, verbose=True, **kwargs):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        _to_format(self.data, path, log=self.log, verbose=verbose, meta=self.meta, build=build, **kwargs)

######################################################################################
    def offload(self):
        _offload(self.data, self.tmp_path, self.log)
        del self.data
        gc.collect()

    def reload(self):
         self.data = _reload(self.tmp_path, self.log)