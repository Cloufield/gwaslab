import pandas as pd
import numpy as np
import time
import copy
from gwaslab.g_Sumstats_summary import summarize
from gwaslab.g_Sumstats_summary import lookupstatus
from gwaslab.io_preformat_input import preformat
from gwaslab.io_to_formats import _to_format
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import fixID
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
from gwaslab.hm_harmonize_sumstats import rsidtochrpos
from gwaslab.hm_harmonize_sumstats import parallelizeassignrsid
from gwaslab.hm_harmonize_sumstats import parallelinferstrand
from gwaslab.hm_harmonize_sumstats import parallelrsidtochrpos
from gwaslab.util_in_filter_value import filtervalues
from gwaslab.util_in_filter_value import filterout
from gwaslab.util_in_filter_value import filterin
from gwaslab.util_in_filter_value import filterregionin
from gwaslab.util_in_filter_value import filterregionout
from gwaslab.util_in_filter_value import _filter_indel
from gwaslab.util_in_filter_value import _filter_palindromic
from gwaslab.util_in_filter_value import _filter_snp
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
from gwaslab.util_ex_run_clumping import _clump
from gwaslab.util_ex_calculate_ldmatrix import tofinemapping
from gwaslab.util_ex_calculate_prs import _calculate_prs
from gwaslab.viz_plot_mqqplot import mqqplot
from gwaslab.viz_plot_trumpetplot import plottrumpet
from gwaslab.viz_plot_compare_af import plotdaf
from gwaslab.util_ex_run_susie import _run_susie_rss
from gwaslab.qc_fix_sumstats import _check_data_consistency
from gwaslab.util_ex_ldsc import _estimate_h2_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_rg_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_h2_cts_by_ldsc
from gwaslab.util_ex_ldsc import _estimate_partitioned_h2_by_ldsc 
from gwaslab.bd_get_hapmap3 import gethapmap3
import gc

#20220309
class Sumstats():
    def __init__(self,
             sumstats,
             fmt=None,
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
             i2=None,
             phet=None,
             dof=None,
             snpr2=None,
             status=None,
             other=[],
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
        self.ldsc_rg = None
        self.ldsc_h2_cts = None
        self.ldsc_partitioned_h2_summary = None
        self.ldsc_partitioned_h2_results = None
        # meta information
        self.meta = _init_meta() 
        self.build = build
        self.meta["gwaslab"]["study_name"] =  study
        self.meta["gwaslab"]["species"] = species
        
        # initialize attributes for clumping and finmapping
        self.to_finemapping_file_path = ""
        self.to_finemapping_file  = pd.DataFrame()
        self.plink_log = ""
        self.clumps = pd.DataFrame()
        self.pipcs = pd.DataFrame()

        # print gwaslab version information
        _show_version(self.log, verbose=verbose)

        #preformat the data
        self.data  = preformat(
          sumstats=sumstats,
          fmt=fmt,
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
          direction=direction,
          study=study,
          build=build,
          trait=trait,
          status=status,
          other=other,
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
    def update_meta(self):
        self.meta["gwaslab"]["variants"]["variant_number"]=len(self.data)
        if "CHR" in self.data.columns:
            self.meta["gwaslab"]["variants"]["number_of_chromosomes"]=len(self.data["CHR"].unique())
        if "P" in self.data.columns:
            self.meta["gwaslab"]["variants"]["min_P"]=np.min(self.data["P"])
        if "EAF" in self.data.columns:
            self.meta["gwaslab"]["variants"]["min_minor_allele_freq"]=min (np.min(self.data["EAF"]) , 1- np.max(self.data["EAF"]))
        if "N" in self.data.columns:
            self.meta["gwaslab"]["samples"]["sample_size"] = int(self.data["N"].max())
            self.meta["gwaslab"]["samples"]["sample_size_median"] = self.data["N"].median()
            self.meta["gwaslab"]["samples"]["sample_size_min"] = int(self.data["N"].min())

    def summary(self):
        return summarize(self.data)

    def lookup_status(self,status="STATUS"):
        return lookupstatus(self.data[status])
    
    def set_build(self, build, verbose=True):
        self.data, self.meta["gwaslab"]["genome_build"] = _set_build(self.data, build=build, log=self.log,verbose=verbose)
        gc.collect()

    def infer_build(self,**args):
        self.data, self.meta["gwaslab"]["genome_build"] = inferbuild(self.data,**args)
    
    def liftover(self,to_build, from_build=None,**args):
        if from_build is None:
            if self.meta["gwaslab"]["genome_build"]=="99":
                self.data, self.meta["gwaslab"]["genome_build"] = inferbuild(self.data,**args)
            from_build = self.meta["gwaslab"]["genome_build"]
        self.data = parallelizeliftovervariant(self.data,from_build=from_build, to_build=to_build, log=self.log,**args)
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
                    fixchr_agrs={},
                    fixpos_args={},
                    fixallele_args={},
                    sanitycheckstats_args={},
                    consistencycheck_args={},
                    normalize=True,
                    normalizeallele_args={},
                    verbose=True):
        ###############################################
        # try to fix data without dropping any information
        self.data = fixID(self.data,verbose=verbose, **fixid_args)
        self.data = fixchr(self.data,log=self.log,remove=remove,verbose=verbose,**fixchr_agrs)
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
              n_cores=1,
              remove=False,
              checkref_args={},
              removedup_args={},
              assignrsid_args={},
              inferstrand_args={},
              flipallelestats_args={},
              liftover_args={},
              fixid_args={},
              fixchr_agrs={},
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
            
            self.data = fixID(self.data,**fixid_args)
            
            self.data = fixchr(self.data,remove=remove,log=self.log,**fixchr_agrs)
            
            self.data = fixpos(self.data,remove=remove,log=self.log,**fixpos_args)
            
            self.data = fixallele(self.data,log=self.log,**fixallele_args)
            
            self.data = sanitycheckstats(self.data,log=self.log,**sanitycheckstats_args)
            
            self.data = parallelnormalizeallele(self.data,log=self.log,n_cores=n_cores,**normalizeallele_args)
            
            self.data = sortcolumn(self.data,log=self.log)
            
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
            
            self.data = checkref(self.data,ref_seq,log=self.log,**checkref_args)
            
            self.meta["gwaslab"]["references"]["ref_seq"] = ref_seq
            
            self.data = flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
            gc.collect()
            
        if ref_infer is not None: 
            
            self.data= parallelinferstrand(self.data,ref_infer = ref_infer,ref_alt_freq=ref_alt_freq,maf_threshold=maf_threshold,
                                              n_cores=n_cores,log=self.log,**inferstrand_args)

            self.meta["gwaslab"]["references"]["ref_infer"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer"] , ref_infer)
            
            self.data =flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
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
    ############################################################################################################
    #customizable API to build your own QC pipeline
    def fix_id(self,**args):
        self.data = fixID(self.data,log=self.log,**args)
    def fix_chr(self,**args):
        self.data = fixchr(self.data,log=self.log,**args)
    def fix_pos(self,**args):
        self.data = fixpos(self.data,log=self.log,**args)
    def fix_allele(self,**args):
        self.data = fixallele(self.data,log=self.log,**args)
    def remove_dup(self,**args):
        self.data = removedup(self.data,log=self.log,**args)
    def check_sanity(self,**args):
        self.data = sanitycheckstats(self.data,log=self.log,**args)
    def check_data_consistency(self, **args):
        _check_data_consistency(self.data,log=self.log,**args)
    def check_id(self,**args):
        pass
    def check_ref(self,ref_seq,**args):
        self.meta["gwaslab"]["references"]["ref_seq"] = ref_seq
        self.data = checkref(self.data,ref_seq,log=self.log,**args)
    def infer_strand(self,ref_infer,**args):
        self.meta["gwaslab"]["references"]["ref_infer"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer"] , ref_infer)
        self.data = parallelinferstrand(self.data,ref_infer=ref_infer,log=self.log,**args)
    def flip_allele_stats(self,**args):
        self.data = flipallelestats(self.data,log=self.log,**args)
    def normalize_allele(self,**args):
        self.data = parallelnormalizeallele(self.data,log=self.log,**args)
    def assign_rsid(self,
                    ref_rsid_tsv=None,
                    ref_rsid_vcf=None,
                    **args):
        if ref_rsid_tsv is not None:
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_tsv,ref_mode="tsv",log=self.log,**args)
            self.meta["gwaslab"]["references"]["ref_rsid_tsv"] = ref_rsid_tsv
        if ref_rsid_vcf is not None:
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_vcf,ref_mode="vcf",log=self.log,**args)   
            self.meta["gwaslab"]["references"]["ref_rsid_vcf"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_rsid_vcf"] , ref_rsid_vcf)
    def rsid_to_chrpos(self,**args):
        self.data = rsidtochrpos(self.data,log=self.log,**args)
    def rsid_to_chrpos2(self,**args):
        self.data = parallelrsidtochrpos(self.data,log=self.log,**args)

    ############################################################################################################
    
    def sort_coordinate(self,**sort_args):
        self.data = sortcoordinate(self.data,log=self.log,**sort_args)
        self.meta["is_sorted"] = True
    def sort_column(self,**args):
        self.data = sortcolumn(self.data,log=self.log,**args)
    
    ############################################################################################################
    def fill_data(self, verbose=True, **args):
        self.data = filldata(self.data, verbose=verbose, **args)
        self.data = sortcolumn(self.data, verbose=verbose, log=self.log)

# utilities ############################################################################################################
    # filter series ######################################################################
    def filter_flanking(self, inplace=False,**args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking(new_Sumstats_object.data, **args)
            return new_Sumstats_object
        else:
            self.data = _get_flanking(self.data, **args)
    def filter_flanking_by_chrpos(self, chrpos,  inplace=False,**args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking_by_chrpos(new_Sumstats_object.data, chrpos, **args)
            return new_Sumstats_object
        else:
            self.data = _get_flanking_by_chrpos(self.data, chrpos,**args)
    def filter_flanking_by_id(self, snpid, inplace=False,**args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _get_flanking_by_id(new_Sumstats_object.data, snpid, **args)
            return new_Sumstats_object
        else:
            self.data = _get_flanking_by_id(self.data, snpid, **args)
    def filter_value(self, expr, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filtervalues(new_Sumstats_object.data,expr,log=new_Sumstats_object.log, **args)
            return new_Sumstats_object
        else:
            self.data = filtervalues(self.data, expr,log=self.log,**args)
    def filter_out(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterout(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = filterout(self.data,log=self.log,**args)
    def filter_in(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterin(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = filterin(self.data,log=self.log,**args)
    def filter_region_in(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterregionin(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = filterregionin(self.data,log=self.log,**args)
    def filter_region_out(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = filterregionout(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = filterregionout(self.data,log=self.log,**args)
    def filter_palindromic(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_palindromic(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = _filter_palindromic(self.data,log=self.log,**args)
    def filter_snp(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_snp(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = _filter_snp(self.data,log=self.log,**args)
    def filter_indel(self, inplace=False, **args):
        if inplace is False:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = _filter_indel(new_Sumstats_object.data,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        else:
            self.data = _filter_indel(self.data,log=self.log,**args)

    def random_variants(self,inplace=False,n=1,p=None,**args):
        if inplace is True:
            self.data = sampling(self.data,n=n,p=p,log=self.log,**args)
        else:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = sampling(new_Sumstats_object.data,n=n,p=p,log=new_Sumstats_object.log,**args)
            return new_Sumstats_object
        
    def filter_hapmap3(self, inplace=False, build=None, **args ):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        if inplace is True:
            self.data = gethapmap3(self.data, build=build,log=self.log, **args)
        else:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = gethapmap3(new_Sumstats_object.data, build=build,log=self.log, **args)
            return new_Sumstats_object
    ######################################################################
    
    def check_af(self,ref_infer,**args):
        self.data = parallelecheckaf(self.data,ref_infer=ref_infer,log=self.log,**args)
        self.meta["gwaslab"]["references"]["ref_infer_daf"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer_daf"] , ref_infer)
    def infer_af(self,ref_infer,**args):
        self.data = paralleleinferaf(self.data,ref_infer=ref_infer,log=self.log,**args)
        self.meta["gwaslab"]["references"]["ref_infer_af"] = ref_infer
        self.meta["gwaslab"]["references"]["ref_infer_af"] = _append_meta_record(self.meta["gwaslab"]["references"]["ref_infer_af"] , ref_infer)
    def plot_daf(self, **args):
        fig,outliers = plotdaf(self.data, **args)
        return fig, outliers
    def plot_mqq(self, build=None, **args):

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
                       **args)
        
        return plot
    
    def plot_trumpet(self, build=None, **args):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        fig = plottrumpet(self.data,build = build,  **args)
        return fig

    def get_lead(self, build=None, gls=False, **args):
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
                           **args)
        # return sumstats object    
        if gls == True:
            new_Sumstats_object = copy.deepcopy(self)
            new_Sumstats_object.data = output
            gc.collect()
            return new_Sumstats_object
        return output

    def get_density(self, sig_list=None, windowsizekb=100,**args):
        
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

        
    def get_novel(self, **args):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = getnovel(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           log=self.log,
                           **args)
        # return sumstats object    
        return output
    
    def check_cis(self, **args):
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
                           **args)
        # return sumstats object    
        return output
    
    def check_novel_set(self, **args):
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
                           **args)
        # return sumstats object    
        return output
    
    def anno_gene(self, **args):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = annogene(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           log=self.log,
                           **args)
        return output
        
    def get_per_snp_r2(self,**args):
        self.data = _get_per_snp_r2(self.data, beta="BETA", af="EAF", n="N", log=self.log, **args)
        #add data inplace

    def get_gc(self, mode=None, **args):
        if mode is None:
            if "P" in self.data.columns:
                output = lambdaGC(self.data[["CHR","P"]],mode="P",**args)
            elif "Z" in self.data.columns:
                output = lambdaGC(self.data[["CHR","Z"]],mode="Z",**args)
            elif "CHISQ" in self.data.columns:
                output = lambdaGC(self.data[["CHR","CHISQ"]],mode="CHISQ",**args)
            elif "MLOG10P" in self.data.columns:
                output = lambdaGC(self.data[["CHR","MLOG10P"]],mode="MLOG10P",**args)
            
            #return scalar
            self.meta["Genomic inflation factor"] = output
            return output    
        else:
            output = lambdaGC(self.data[["CHR",mode]],mode=mode,**args)
            self.meta["Genomic inflation factor"] = output
            return output 
## LDSC ##############################################################################################
    def estimate_h2_by_ldsc(self, build=None, verbose=True, match_allele=True, **args):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=True )
        self.ldsc_h2 = _estimate_h2_by_ldsc(insumstats=insumstats, log=self.log, verbose=verbose, **args)
    
    def estimate_rg_by_ldsc(self, build=None, verbose=True, match_allele=True, **args):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=True )
        self.ldsc_rg = _estimate_rg_by_ldsc(insumstats=insumstats, log=self.log, verbose=verbose, **args)

    def estimate_h2_cts_by_ldsc(self, build=None, verbose=True, match_allele=True, **args):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=True )
        self.ldsc_h2_cts  = _estimate_h2_cts_by_ldsc(insumstats=insumstats, log=self.log, verbose=verbose, **args)
    
    def estimate_partitioned_h2_by_ldsc(self, build=None, verbose=True, match_allele=True, **args):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        insumstats = gethapmap3(self.data.copy(), build=build, verbose=verbose , match_allele=True )
        self.ldsc_partitioned_h2_summary, self.ldsc_partitioned_h2_results  = _estimate_partitioned_h2_by_ldsc(insumstats=insumstats, log=self.log, verbose=verbose, **args)
# external ################################################################################################
    
    def calculate_ld_matrix(self,**args):
        self.to_finemapping_file_path, self.to_finemapping_file, self.plink_log  = tofinemapping(self.data,study = self.meta["gwaslab"]["study_name"],**args)
    
    def run_susie_rss(self,**args):
        self.pipcs=_run_susie_rss(self.to_finemapping_file_path,**args)
    
    def clump(self,**args):
        self.clumps,self.plink_log = _clump(self.data, log=self.log, study = self.meta["gwaslab"]["study_name"], **args)

    def calculate_prs(self,**args):
        combined_results_summary = _calculate_prs(self.data, log=self.log, study = self.meta["gwaslab"]["study_name"], **args)
        return combined_results_summary
    
# to_format ###############################################################################################       

    def to_format(self, path, build=None, **args):
        if build is None:
            build = self.meta["gwaslab"]["genome_build"]
        _to_format(self.data, path, log=self.log, meta=self.meta, build=build, **args)