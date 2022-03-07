import pandas as pd
import gwaslab as gl

#20220307
class Sumstats():
    # Sumstats.data will store a pandas df
    data = pd.DataFrame()
    build = "Unknown"
    def __init__(self,
                 sumstats,
                 snpid=None,
                 rsid=None,
                 chrom=None,
                 pos=None,
                 ea=None,
                 nea=None,
                 eaf=None,
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
                 other=[],
                 verbose=True,
                 build="hg19",
                 **arguments):
        
        self.data  = gl.preformat(sumstats,snpid,rsid,chrom,pos,ea,nea,eaf,n,beta,se,chisq,z,p,mlog10p,info,OR,OR_se,OR_95L,OR_95U,other,verbose,**arguments)
        self.build = build

    def __getattr__(self, name):
        if name == 'gc':
            return self.get_gc(verbose=False)
        else:
            raise AttributeError        
            
#### healper #################################################################################
    #value to arguement name 
    def get_columns(self, columns=[]):
        dic = {
            "MARKERNAME": "snpid",
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
        if len(columns)>0:
            col_subset = {dic[key]: key for key in columns}
            return col_subset
        else:
            other=[]
            for key in self.data.columns:
                if key not in dic.keys():
                    other.append(key)
                else:
                    columns.append(key)
              
            col_subset = {dic[key]: key for key in columns}
            if len(other)>0: col_subset["other"] = other  
            return col_subset
    
# QC ######################################################################################
    
    def fix_chr(self,chrom="CHR",add_prefix=""):
        self.data[chr] = self.data[chr].astype("string")
        self.data[chr] = self.data[chr].str.lstrip("chr|CHR|Chr")
        self.data[chr] = self.data[chr].str.lstrip("0")
        self.data[chr] = self.data[chr].str.upper()
        if add_prefix:
            self.data[chr] = add_prefix + sumstats[chr]
        #return sumstats

# utilities ###############################################################################

    def plot_mqq(self, **args):
        plot = gl.mqqplot(self.data, chrom="CHR", pos="POS", p="P", **args)
        return plot

    def get_lead(self, **args):
        output = gl.getsig(self.data,
                           id="MARKERNAME",
                           chrom="CHR",
                           pos="POS",
                           p="P",
                           **args)
        return Sumstats(output, **self.get_columns(), verbose=False)

    def get_gc(self, **args):
        output = gl.gc(self.data["P"], **args)
        return output
    
    def exclude_hla(self, **args):
        is_hla = (self.data["CHR"] == "6") & (self.data["POS"] > 25000000) & (
            self.data["POS"] < 34000000)
        output = self.data.loc[~is_hla]  #not in hla
        return Sumstats(output, **self.get_columns(), verbose=False)

    def extract_hla(self, **args):
        is_hla = (self.data["CHR"] == "6") & (self.data["POS"] > 25000000) & (
            self.data["POS"] < 34000000)
        output = self.data.loc[is_hla]  #in hla
        return Sumstats(output, **self.get_columns(), verbose=False)
    
    def extract_chr(self,chrom="1", **args):
        is_chr = (self.data["CHR"] == str(chrom))
        output = self.data.loc[is_chr]  #in hla
        return Sumstats(output, **self.get_columns(), verbose=False)
    
    def lift_over(self,from_build="hg19" ,to_build="hg38",remove_unmapped=True):
        output = gl.liftover(self.data,chrom="CHR",pos="POS",
                                from_build=from_build,to_build=to_build,
                                remove_unmapped=remove_unmapped)
        return Sumstats(output, other=["CHR_"+to_build,"POS_"+to_build], **self.get_columns(), verbose=False)
        
# to_format ###############################################################################################       
    ## ldsc ##################################################################################################
    def to_ldsc(self, path="", exclude_hla=True, **args):
        if "rsID" in list(self.data.columns.values):
            self.hapmap3 = gl.ldsc(self.data, rsid="rsID", **args)
        else:
            self.hapmap3 = gl.ldsc(self.data, chrom="CHR", pos="POS", **args)
        if path:
            self.hapmap3.to_csv(path, sep="\t", index=False)
            print(
                "  - Saving sumstats in ldsc format at Hapmap3 variants(excluding hla region) to : "
            )
            print("  - " + path)
        if exclude_hla is True:
            print("  - Excluding variants in HLA region ...")
            is_hla = (self.hapmap3["CHR"]
                      == "6") & (self.hapmap3["POS"] >
                               25000000) & (self.hapmap3["POS"] < 34000000)
            self.hapmap3 = self.hapmap3.loc[~is_hla]
            print("  - Hapmap3 variants in sumstas : ", len(self.hapmap3))
        return Sumstats(self.hapmap3,
                        **self.get_columns(self.hapmap3.columns),
                        verbose=False)
    
    ## bed ##################################################################################################
    def to_bed(self, path="",snpid="MARKERNAME",chrom="CHR",pos="POS",ref="NEA",alt="EA",other=[],flank=0):
        output = pd.DataFrame()
        flank = int(flank)
        output["Chrom"]= self.data[chrom]
        output["Start"]= self.data[pos].astype("int") - flank
        output["End"]=   self.data[pos].astype("int") + flank
        if ref and alt:
            output["Ref"]=   self.data[ref]
            output["Alt"]=   self.data[alt]
        output["Name"]=  self.data[snpid]
        if path:
            output.to_csv(path,"\t",header=None)
        else:
            return output
    #fuma ############################################################################################
    def to_fuma(self,path=None):
        sumstats_fuma = gl.tofuma(self.data,path)
        return sumstats_fuma
