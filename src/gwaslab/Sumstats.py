import pandas as pd
import numpy as np
import time
from gwaslab.Sumstatssummary import summarize
from gwaslab.Sumstatssummary import lookupstatus
from gwaslab.preformat_input import preformat
from gwaslab.Log import Log
from gwaslab.fixdata import fixID
from gwaslab.fixdata import removedup
from gwaslab.fixdata import fixchr
from gwaslab.fixdata import fixpos
from gwaslab.fixdata import fixallele
from gwaslab.fixdata import parallelnormalizeallele
from gwaslab.fixdata import sanitycheckstats
from gwaslab.fixdata import parallelizeliftovervariant
from gwaslab.fixdata import flipallelestats
from gwaslab.fixdata import sortcoordinate
from gwaslab.fixdata import sortcolumn
from gwaslab.retrievedata import parallelecheckaf
from gwaslab.retrievedata import checkref
from gwaslab.retrievedata import rsidtochrpos
from gwaslab.retrievedata import parallelizeassignrsid
from gwaslab.retrievedata import parallelinferstrand
from gwaslab.retrievedata import parallelrsidtochrpos
from gwaslab.filtervalue import filterout
from gwaslab.filtervalue import filterin
from gwaslab.filtervalue import filterregionin
from gwaslab.filtervalue import filterregionout
from gwaslab.filtervalue import inferbuild
from gwaslab.filtervalue import sampling
from gwaslab.mqqplot import mqqplot
from gwaslab.calculate_gc import gc
from gwaslab.getsig import getsig
from gwaslab.getsig import annogene
from gwaslab.getsig import getnovel
from gwaslab.fill import filldata
from gwaslab.to_formats import tofmt
from gwaslab.get_hapmap3 import gethapmap3
from gwaslab.compare_effect import plotdaf
from gwaslab.CommonData import get_chr_list
from gwaslab.CommonData import get_number_to_chr
from gwaslab.CommonData import get_chr_to_number
from gwaslab.CommonData import get_high_ld
from gwaslab.CommonData import get_format_dict
from gwaslab.CommonData import get_formats_list
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
             n=None,
             beta=None,
             se=None,
             chisq=None,
             z=None,
             p=None,
             mlog10p=None,
             info=None,
             OR=None,
             OR_95L=None,
             OR_95U=None,
             status=None,
             other=[],
             direction=None,
             verbose=True,
             build="99",
             **readargs):
        
        self.build = build
        self.meta = {"Genome build":build,
                 "Name":"Sumstats_1","StudyType":"Unknown"}
        self.data = pd.DataFrame()
        self.log = Log()
        
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
          n=n,
          beta=beta,
          se=se,
          chisq=chisq,
          z=z,
          p=p,
          mlog10p=mlog10p,
          info=info,
          OR=OR,
          OR_95L=OR_95L,
          OR_95U=OR_95U,
          direction=direction,
          build=build,
          status=status,
          other=other,
          verbose=verbose,
          readargs=readargs,
          log=self.log)
        gc.collect()     
#### healper #################################################################################
    #value to arguement name 
    def get_columns(self, columns=[]):
        dic = {
            "SNPID": "snpid",
            "rsID": "rsid",
            "CHR": "chrom",
            "POS": "pos",
            "NEA": "nea",
            "EA": "ea",
            "EAF": "eaf",
            "N": "n",
            "BETA": "beta",
            "SE": "se",
            "Z": "z",
            "CHISQ": "chisq",
            "P": "p",
            "MLOG10P": "mlog10p",
            "INFO": "info",
            "OR": "OR",
            "OR_95L": "OR_95L",
            "OR_95U": "OR_95U",
        }
        other=[]
        col_subset={}
        for key in columns:
            if key not in dic.keys():
                other.append(key)
            else:
                col_subset[dic[key]] = key
        if len(other)>0: col_subset["other"] = other      
        return col_subset

##########################################################################################
    def update_meta(self):
        self.meta["Number_of_variants"]=len(self.data)
        if "CHR" in self.data.columns:
            self.meta["Number_of_chromosomes"]=len(self.data["CHR"].unique())
        if "P" in self.data.columns:
            self.meta["Min_P"]=np.min(self.data["P"])
    def summary(self):
        return summarize(self.data)
    def lookup_status(self,status="STATUS"):
        return lookupstatus(self.data[status])
        
# QC ######################################################################################
    #clean the sumstats with one line
    def basic_check(self,remove=False,
                    n_cores=1,
                    fixid_args={},
                    removedup_args={},
                    fixchr_agrs={},
                    fixpos_args={},
                    fixallele_args={},
                    sanitycheckstats_args={},
                    normalize=True,
                    normalizeallele_args={}):
        ###############################################
        # try to fix data without dropping any information
        self.data = fixID(self.data,**fixid_args)
        if remove is True:
            self.data = removedup(self.data,log=self.log,**removedup_args)
        self.data = fixchr(self.data,log=self.log,remove=remove,**fixchr_agrs)
        self.data = fixpos(self.data,log=self.log,remove=remove,**fixpos_args)
        self.data = fixallele(self.data,log=self.log,**fixallele_args)
        self.data = sanitycheckstats(self.data,log=self.log,**sanitycheckstats_args)
        if normalize is True:
            self.data = parallelnormalizeallele(self.data,n_cores=n_cores,log=self.log,**normalizeallele_args)
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
            
            self.data = flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
            gc.collect()
            
        if ref_infer is not None: 
            
            self.data= parallelinferstrand(self.data,ref_infer = ref_infer,ref_alt_freq=ref_alt_freq,
                                              n_cores=n_cores,log=self.log,**inferstrand_args)
            
            self.data =flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
            gc.collect()
        
        #####################################################
        if ref_rsid_tsv is not None:
            
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_tsv,ref_mode="tsv",
                                                 n_cores=n_cores,log=self.log,**assignrsid_args)
            gc.collect()
        if ref_rsid_vcf is not None:
            
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_vcf,ref_mode="vcf",
                                                 n_cores=n_cores,log=self.log,**assignrsid_args)   
            gc.collect()
        ######################################################    
        if remove is True:
            
            self.data = removedup(self.data,log=self.log,**removedup_args)
        ################################################ 
        
        self.data = sortcoordinate(self.data,log=self.log)
        
        self.data = sortcolumn(self.data,log=self.log)
        gc.collect()
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
    def check_id(self,**args):
        pass
    def check_ref(self,**args):
        self.data = checkref(self.data,log=self.log,**args)
    def infer_strand(self,**args):
        self.data =inferstrand(self.data,maf_threshold=0.43,log=self.log,**args)
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
        if ref_rsid_vcf is not None:
            self.data = parallelizeassignrsid(self.data,path=ref_rsid_vcf,ref_mode="vcf",log=self.log,**args)   
        
    def rsid_to_chrpos(self,**args):
        self.data = rsidtochrpos(self.data,log=self.log,**args)
    def rsid_to_chrpos2(self,**args):
        self.data = parallelrsidtochrpos(self.data,log=self.log,**args)
    def liftover(self,**args):
        self.data = parallelizeliftovervariant(self.data,log=self.log,**args)
    def sort_coordinate(self,**sort_args):
        self.data = sortcoordinate(self.data,log=self.log,**sort_args)
    def sort_column(self,**args):
        self.data = sortcolumn(self.data,log=self.log,**args)
    ############################################################################################################
    def fill_data(self, **args):
        self.data = filldata(self.data,**args)
    def infer_build(self,**args):
        self.data = inferbuild(self.data,**args)

# utilities ############################################################################################################
    def filter_out(self, inplace=True, **args):
        if inplace is False:
            return filterout(self.data,**args)
        else:
            self.data = filterout(self.data,log=self.log,**args)
            
    def filter_in(self, inplace=True, **args):
        if inplace is False:
            return filterin(self.data,**args)
        else:
            self.data = filterin(self.data,log=self.log,**args)
    def filter_region_in(self, inplace=True, **args):
        if inplace is False:
            return filterregionin(self.data,**args)
        else:
            self.data = filterregionin(self.data,log=self.log,**args)
    def filter_region_out(self, inplace=True, **args):
        if inplace is False:
            return filterregionout(self.data,**args)
        else:
            self.data = filterregionout(self.data,log=self.log,**args)
    
    def random_variants(self,in_place=False,n=1,**args):
        if in_place is True:
            self.data = sampling(self.data,n=n,log=self.log,**args)
        else:
            return sampling(self.data,n=n,log=self.log,**args)
    
###########################################################################################################################            
    def check_af(self,**args):
        self.data = parallelecheckaf(self.data,log=self.log,**args)
      
    def plot_daf(self, **args):
        plot = plotdaf(self.data, **args)
        
    def plot_mqq(self, **args):

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
        
        plot = mqqplot(self.data,snpid=snpid, chrom=chrom, pos=pos, p=p, eaf=eaf,**args)
        
        return plot

    def get_lead(self, **args):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = getsig(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           log=self.log,
                           **args)
        # return sumstats object    
        return output
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
        
    def get_per_snp_h2(self,**args):
        self.data = getpersnph2(self.data,beta="BETA",af="EAF",**args)
        #add data inplace

    
    def get_gc(self, **args):
        if "P" in self.data.columns:
            output = gc(self.data["P"],mode="p",**args)
        elif "Z" in self.data.columns:
            output = gc(self.data["Z"],mode="z",**args)
        elif "CHISQ" in self.data.columns:
            output = gc(self.data["CHISQ"],mode="chi2",**args)
        #return scalar
        self.meta["Genomic inflation factor"] = output
        return output    
    
    def exclude_hla(self, **args):
        is_hla = (self.data["CHR"] == "6") & (self.data["POS"] > 25000000) & (
            self.data["POS"] < 34000000)
        output = self.data.loc[~is_hla]  #not in hla
        # return sumstats object    
        return Sumstats(output, **self.get_columns(), verbose=False)

    def extract_hla(self, **args):
        is_hla = (self.data["CHR"] == "6") & (self.data["POS"] > 25000000) & (
            self.data["POS"] < 34000000)
        output = self.data.loc[is_hla]  #in hla
        # return sumstats object    
        return Sumstats(output, **self.get_columns(self.data.columns), verbose=False)
    
    def extract_chr(self,chrom="1", **args):
        is_chr = (self.data["CHR"] == str(chrom))
        output = self.data.loc[is_chr]  #in hla
        # return sumstats object    
        return Sumstats(output, **self.get_columns(self.data.columns), verbose=False)
    
    

        
# to_format ###############################################################################################       
    def to_format(self,
              path="./sumstats",
              fmt="ldsc",   
              extract=None,
              exclude=None,
              id_use="rsID",
              hapmap3=False,
              exclude_hla=False,  
              build="19", 
              n=None,
              verbose=True,
              output_log=True,
              to_csvargs={},
              float_formats={},
              xymt_number=False,
              xymt=["X","Y","MT"],
              chr_prefix=None,
              bgzip=False,
              tabix=False):
        
        if fmt in get_formats_list() + ["vep","bed","annovar","vcf"]:
            if verbose: self.log.write("Start to format the output sumstats in: ",fmt, " format")
        else:
            raise ValueError("Please select a format to output")
        
        
        #######################################################################################################
        # filter
        output = self.data.copy()
        if extract is not None:
            output = output.loc[output[id_use].isin(extract),:]

        if exclude is not None:
            output = output.loc[~output[id_use].isin(exclude),:]
        
        #hla and hapmap3 #######################################################################################
        suffix=fmt
        
        #exclude hla
        if exclude_hla is True:
            if verbose: self.log.write(" -Excluding variants in HLA region ...")
            before = len(output)
            is_hla = (output["CHR"]== 6) & (output["POS"].astype("Int64") >25000000) & (output["POS"].astype("Int64") < 34000000)
            output = output.loc[~is_hla,:]
            after = len(output)
            if verbose: self.log.write(" -Exclude "+ str(before - after) + " variants in HLA region.")
            suffix = "noMHC."+suffix
        
        #extract hapmap3 SNPs
        if hapmap3 is True:
            output = gethapmap3(output,build=build,verbose=True)
            after = len(output)
            if verbose: self.log.write(" -Extract "+ str(after) + " variants in Hapmap3 datasets for build "+build+".")
            suffix = "hapmap3."+suffix
        
        # add a n column
        if n is not None:
            output["N"] = n
                
        #######################################################################################################
        #formatting statistics
        if verbose: self.log.write(" -Formatting statistics ...")
        
        formats = {'EAF': '{:.4g}', 
                'BETA': '{:.4f}', 
                'Z': '{:.4f}',
                'CHISQ': '{:.4f}',
                'SE': '{:.4f}',
                'OR': '{:.4f}',
                'OR_95U': '{:.4f}',
                'OR_95L': '{:.4f}',
                'INFO': '{:.4f}',
                'P': '{:.4e}',
                'MLOG10P': '{:.4f}',
                'DAF': '{:.4f}'
                  }
        
        for col, f in float_formats.items():
            if col in output.columns: 
                formats[col]=f
        for col, f in formats.items():
            if col in output.columns: 
                if output[col].dtype in ["float64","float32","float16","float"]:
                    output[col] = output[col].map(f.format)
        if verbose: 
            self.log.write(" - Float statistics formats:")  
            keys=[]
            values=[]
            for key,value in formats.items():
                if key in output.columns: 
                    keys.append(key)
                    values.append(value)
            self.log.write("  - Columns:",keys) 
            self.log.write("  - Output formats:",values) 
            
        ##########################################################################################################          
        # output, mapping column names
        
        if fmt in get_formats_list() + ["vep","bed","annovar","vcf"]:
            tofmt(output,
                  path=path,
                  fmt=fmt,
                  suffix=suffix,
                  build=build,
                  verbose=True,
                  log=self.log,
                  to_csvargs=to_csvargs,
                  chr_prefix=chr_prefix,
                  meta = self.meta,
                  bgzip=bgzip,
                  tabix=tabix,
                  xymt_number=xymt_number,
                  xymt=xymt)
        if output_log is True:
            if verbose: self.log.write(" -Saving log file:",path + "."+ suffix +".log")
        if verbose: self.log.write("Finished outputting successfully!")
