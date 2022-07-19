import pandas as pd
import numpy as np
import gwaslab as gl
import time
#20220309
class Sumstats():
    def __init__(self,
             sumstats,
             snpid=None,
             rsid=None,
             chrom=None,
             pos=None,
             ea=None,
             nea=None,
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
             OR_se=None,
             OR_95L=None,
             OR_95U=None,
             status=None,
             other=[],
             direction=None,
             verbose=True,
             build="99",
             **readargs):
        
        self.build = build
        self.meta = {"Genome build":build}
        self.data = pd.DataFrame()
        self.log = gl.Log()
        
        #preformat the data
        self.data  = gl.preformat(sumstats=sumstats,
          snpid=snpid,
          rsid=rsid,
          chrom=chrom,
          pos=pos,
          ea=ea,
          nea=nea,
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
          OR_se=OR_se,
          OR_95L=OR_95L,
          OR_95U=OR_95U,
          direction=direction,
          build=build,
          status=status,
          other=other,
          verbose=verbose,
          readargs=readargs,
          log=self.log)
                   
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
            "OR_SE": "OR_se",
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
    
    def update_meta(self):
        self.meta["Number_of_variants"]=len(self.data)
        if "CHR" in self.data.columns:
            self.meta["Number_of_chromosomes"]=len(self.data["CHR"].unique())
        if "P" in self.data.columns:
            self.meta["Min_P"]=np.min(self.data["P"])
        
        
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
        self.data = gl.fixID(self.data,**fixid_args)
        if remove is True:
            self.data = gl.removedup(self.data,log=self.log,**removedup_args)
        self.data = gl.fixchr(self.data,log=self.log,remove=remove,**fixchr_agrs)
        self.data = gl.fixpos(self.data,log=self.log,remove=remove,**fixpos_args)
        self.data = gl.fixallele(self.data,log=self.log,**fixallele_args)
        self.data = gl.sanitycheckstats(self.data,log=self.log,**sanitycheckstats_args)
        if normalize is True:
            self.data = gl.parallelnormalizeallele(self.data,n_cores=n_cores,log=self.log,**normalizeallele_args)
        ###############################################
        
    
    def harmonize(self,
              basic_check=True,
              ref_seq=None,
              ref_rsid=None,
              ref_rsid_chr_dict=None,
              ref_infer=None,
              ref_infer_chr_dict=None,
              ref_alt_freq=None,
              from_build=None,
              to_build=None,
              to_ref_seq=None,
              to_ref_infer=None,
              to_ref_rsid=None,
              to_ref_alt_freq=None,
              maf_threshold=0.43,
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
            self.data = gl.fixID(self.data,**fixid_args)
            self.data = gl.fixchr(self.data,remove=remove,log=self.log,**fixchr_agrs)
            self.data = gl.fixpos(self.data,remove=remove,log=self.log,**fixpos_args)
            self.data = gl.fixallele(self.data,log=self.log,**fixallele_args)
            self.data = gl.sanitycheckstats(self.data,log=self.log,**sanitycheckstats_args)
            self.data = gl.parallelnormalizeallele(self.data,log=self.log,**normalizeallele_args)
            self.data = gl.sortcolumn(self.data,log=self.log)
        
        
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
            self.data = gl.checkref(self.data,ref_seq,log=self.log,**checkref_args)
            self.data = gl.flipallelestats(self.data,log=self.log,**flipallelestats_args)
        
        if ref_rsid is not None:
            self.data = gl.parallelizeassignrsid(self.data,ref_rsid,n_cores=n_cores,log=self.log,chr_dict=ref_rsid_chr_dict,**assignrsid_args)
            
        if ref_infer is not None: 
            self.data= gl.inferstrand(self.data,ref_infer = ref_infer,ref_alt_freq=ref_alt_freq,maf_threshold=0.43,chr_dict=ref_infer_chr_dict,log=self.log,**inferstrand_args)
            self.data =gl.flipallelestats(self.data,log=self.log,**flipallelestats_args)
            
        #if liftover           
        if to_build is not None :
            self.data = gl.parallelizeliftovervariant(self.data,n_cores=n_cores,from_build=from_build,to_build=to_build,log=self.log,**liftover_args)
            if to_ref_seq is not None:
                self.data = gl.checkref(self.data,to_ref_seq,log=self.log,**checkref_args)
                self.data = gl.flipallelestats(self.data,log=self.log,**flipallelestats_args)
                if to_ref_rsid is not None  :
                    self.data = gl.parallelizeassignrsid(self.data,to_ref_rsid,n_cores=n_cores,log=self.log,**assignrsid_args)
                if (to_ref_infer is not None) and (to_ref_alt_freq is not None): 
                    self.data = gl.inferstrand(self.data,ref_infer = to_ref_infer,ref_alt_freq=to_ref_alt_freq,maf_threshold=0.43,log=self.log,**inferstrand_args)
                    self.data = gl.flipallelestats(self.data,log=self.log,**flipallelestats_args)
        self.data = gl.fixchr(self.data, remove=remove)  
        if remove is True:
            self.data = gl.removedup(self.data,log=self.log,**removedup_args)
        ################################################ 
        self.data = gl.sortcoordinate(self.data,log=self.log)
        self.data = gl.sortcolumn(self.data,log=self.log)
        return self
    
    ############################################################################################################
    #customizable API to build your own QC pipeline
    def fix_ID(self,**args):
        self.data = gl.fixID(self.data,log=self.log,**args)
    def fix_chr(self,**args):
        self.data = gl.fixchr(self.data,log=self.log,**args)
    def fix_pos(self,**args):
        self.data = gl.fixpos(self.data,log=self.log,**args)
    def fix_allele(self,**args):
        self.data = gl.fixallele(self.data,log=self.log,**args)
    def remove_dup(self,**args):
        self.data = gl.removedup(self.data,log=self.log,**args)
    def check_sanity(self,**args):
        self.data = gl.sanitycheckstats(self.data,log=self.log,**args)
    def check_ID(self,**args):
        pass
    def check_ref(self,**args):
        self.data = gl.checkref(self.data,log=self.log,**args)
    def infer_strand(self,**args):
        self.data =gl.inferstrand(self.data,maf_threshold=0.43,log=self.log,**args)
    def flip_allele_stats(self,**args):
        self.data = gl.flipallelestats(self.data,log=self.log,**args)
    def normalize_allele(self,**args):
        self.data = gl.parallelnormalizeallele(self.data,log=self.log,**args)
    def assign_rsid(self,**args):
        self.data = gl.parallelizeassignrsid(self.data,log=self.log,**args)
    def rsid_to_chrpos(self,**args):
        self.data = gl.rsidtochrpos(self.data,log=self.log,**args)
    def rsid_to_chrpos2(self,**args):
        self.data = gl.parallelrsidtochrpos(self.data,log=self.log,**args)
    def liftover(self,**args):
        self.data = gl.parallelizeliftovervariant(self.data,log=self.log,**args)
    def sort_coordinate(self,**sort_args):
        self.data = gl.sortcoordinate(self.data,log=self.log,**sort_args)
    def sort_column(self,**args):
        self.data = gl.sortcolumn(self.data,log=self.log,**args)
    ############################################################################################################
    def fill_data(self, **args):
        self.data = gl.filldata(self.data,**args)

# utilities ############################################################################################################
    def filter_out(self, inplace=True, **args):
        if inplace is False:
            return gl.filterout(self.data,**args)
        else:
            self.data = gl.filterout(self.data,log=self.log,**args)
            
    def filter_in(self, inplace=True, **args):
        if inplace is False:
            return gl.filterin(self.data,**args)
        else:
            self.data = gl.filterin(self.data,log=self.log,**args)
###########################################################################################################################            
    def check_af(self,**args):
        self.data = gl.checkaf(self.data,log=self.log,**args)
      
    def plot_daf(self, **args):
        plot = gl.plotdaf(self.data, **args)
        
    def plot_mqq(self, **args):
        plot = gl.mqqplot(self.data, chrom="CHR", pos="POS", p="P", **args)
        return plot

    def get_lead(self, **args):
        if "SNPID" in self.data.columns:
            id_to_use = "SNPID"
        else:
            id_to_use = "rsID"
        output = gl.getsig(self.data,
                           id=id_to_use,
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           log=self.log,
                           **args)
        # return sumstats object    
        return output

    
    def get_per_snp_h2(self,**args):
        self.data = gl.getpersnph2(self.data,beta="BETA",af="EAF",**args)
        #add data inplace

    
    def get_gc(self, **args):
        if "P" in self.data.columns:
            output = gl.gc(self.data["P"],mode="p",**args)
        elif "Z" in self.data.columns:
            output = gl.gc(self.data["Z"],mode="z",**args)
        elif "CHISQ" in self.data.columns:
            output = gl.gc(self.data["CHISQ"],mode="chi2",**args)
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
              format="ldsc",   
              extract=None,
              exclude=None,
              id_use="rsID",
              hapmap3=False,
              exclude_hla=False,        
              verbose=True,
              output_log=True,
              to_csvargs={}):
        
        if format in ["ldsc","fuma","metal","bed","vcf"]:
            if verbose: self.log.write("Output sumstats to format:",format)
        else:
            raise ValueError("Please select a format to output")
        
        #######################################################################################################
        # filter
        output = self.data.copy()
        if extract is not None:
            output = output.loc[output[id_use].isin(extract),:]

        if exclude is not None:
            output = output.loc[~output[id_use].isin(exclude),:]

        if exclude_hla is True:
            if verbose: self.log.write(" -Excluding variants in HLA region ...")
            is_hla = (output["CHR"]== "6") & (output["POS"] >25000000) & (output["POS"] < 34000000)
            output = output.loc[~is_hla,:]

        if hapmap3 is True:
            output = gl.gethapmap3(output)
        #######################################################################################################

        # output
        if format=="ldsc":
            gl.toldsc(output, path=path+"."+format,verbose=True,log=self.log,to_csvargs=to_csvargs)
        if format=="fuma":
            gl.toldsc(output, path=path+"."+format,verbose=True,log=self.log,to_csvargs=to_csvargs)
        if format=="bed":
            gl.toldsc(output, path=path+"."+format,verbose=True,log=self.log,to_csvargs=to_csvargs)
        if format=="metal":
            gl.toldsc(output, path=path+"."+format,verbose=True,log=self.log,to_csvargs=to_csvargs)
        if format=="vcf":
            gl.toldsc(output, path=path+"."+format,verbose=True,log=self.log,to_csvargs=to_csvargs)
        if output_log is True:
            if verbose: self.log.write("Saving log file...")
            self.log.save(path + "."+ format +".log")
            
