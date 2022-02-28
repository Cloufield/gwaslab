import pandas as pd
import gwaslab as gl

class Sumstats():
    
    data = pd.DataFrame()
    columns = ["rsID","CHR","POS","EA","NEA","N","P"]
    
    def __getattr__(self, name):
        if name == 'gc':
            return gl.gc(self.data["P"])
        raise AttributeError
    
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
                 verbose=True,**arguments):
            
            usecols=[]
            
            if snpid: usecols.append(snpid)
            if rsid: usecols.append(rsid)
            if chrom: usecols.append(chrom)
            if pos: usecols.append(pos)
            if ea: usecols.append(ea)
            if nea: usecols.append(nea)
            if eaf: usecols.append(eaf)
            if n and (type(n) is str): usecols.append(n)
            if beta: usecols.append(beta)
            if se: usecols.append(se)
            if chisq: usecols.append(chisq)
            if z:usecols.append(z)
            if p: usecols.append(p)
            if mlog10p: usecols.append(mlog10p)
            if info: usecols.append(info)
            if OR: usecols.append(OR)
            if OR_se:usecols.append(OR_se)
            if OR_95L: usecols.append(OR_95L)
            if OR_95U: usecols.append(OR_95U)
            if other: usecols = usecols + other
                
            rename_dictionary={snpid:   "MARKERNAME",
                               rsid:    "rsID",
                               chrom:   "CHR",
                               pos:     "POS",
                               nea:     "NEA",
                               ea:      "EA",
                               eaf:     "EAF",
                               n:       "N",
                               beta:    "BETA",
                               se:      "SE",
                               z:       "Z",   
                               chisq:   "CHISQ",
                               p:       "P",
                               mlog10p: "MLOG10P",
                               info:    "INFO",
                               OR:      "OR",
                               OR_se:   "OR_SE",
                               OR_95L:  "OR_95L",
                               OR_95U:  "OR_95U"}
                        
            if verbose: print("Reading columns :", usecols)
            if verbose: print("Renaming columns to :", list(map(lambda x:rename_dictionary[x], usecols)))
            
            if type(sumstats) is str:
                inpath = sumstats
                if verbose:print("Initiating Sumstats Object form file :"+ inpath)
                self.data = pd.read_table(inpath,usecols=usecols,sep="\s+")
                
            elif type(sumstats) is pd.DataFrame:
                if verbose:print("Initiating Sumstats Object form pandas dataframe object ...")
                self.data = sumstats.loc[:,usecols]
                
            self.data = self.data.rename(columns=rename_dictionary)
            
            if type(n) is int:
                self.data["N"]=n
            
            if not snpid:
                self.data["MARKERNAME"] = self.data["CHR"].astype("string") + ":" +self.data["POS"].astype("string")
    
    def get_columns(self,columns=None):
        dic={
         "MARKERNAME":"snpid", 
         "rsID":"rsid",
         "CHR":"chrom",
         "POS":"pos",
         "NEA":"nea",
         "EA":"ea",
         "EAF":"eaf",
         "N":"n",
         "BETA":"beta",
         "SE":"se",
         "Z":"z",
         "CHISQ":"chisq",
         "P":"p",
         "MLOG10P":"mlog10p",
         "INFO":"info",
         "OR":"OR",
         "OR_SE":"OR_se",
         "OR_95L":"OR_95L",
         "OR_95U":"OR_95U",
        }
        if columns is not None: 
            col_subset = { dic[key]:key for key in columns}
        else: 
            col_subset = { dic[key]:key for key in self.data.columns}
        return col_subset
    
    def plot_mqq(self, **args):
        plot = gl.mqqplot(self.data,chrom="CHR",pos="POS",p="P",**args)
        return plot
    
    def get_sig(self,**args):
        output = gl.getsig(self.data,id="MARKERNAME",chrom="CHR",pos="POS",p="P",**args)
        return Sumstats(output,**self.get_columns())
    
    def exclude_hla(self,**args):
        is_hla=(self.data["CHR"]==6)&(self.data["POS"]>25000000)&(self.data["POS"]<34000000)
        output=self.data.loc[~is_hla]
        return Sumstats(output,**self.get_columns(),verbose=False)
    
    def extract_hla(self,**args):
        is_hla=(self.data["CHR"]==6)&(self.data["POS"]>25000000)&(self.data["POS"]<34000000)
        output = self.data.loc[is_hla]
        return Sumstats(output,**self.get_columns(),verbose=False)
    
    def get_gc(self,**args):
        output = gl.gc(self.data["P"],**args)
        return output
    
    def to_ldsc(self,path="",exclude_hla=True,**args):
        if "rsID" in list(self.data.columns.values):
            self.hapmap3 = gl.ldsc(self.data,rsid="rsID",**args)
        else:
            self.hapmap3 = gl.ldsc(self.data,chrom="CHR",pos="POS",**args)
        if path:
            self.hapmap3.to_csv(path,sep="\t",index=False)
            print("  - Saving sumstats in ldsc format at Hapmap3 variants(excluding hla region) to : ")
            print("  - "+path)
        if exclude_hla is True:
            print("  - Excluding variants in HLA region ...")
            is_hla=(self.hapmap3["CHR"]==6) & (self.hapmap3["POS"]>25000000)&(self.hapmap3["POS"]<34000000)
            self.hapmap3 = self.hapmap3.loc[~is_hla]
            print("  - Hapmap3 variants in sumstas : ",len(self.hapmap3 ))
        return Sumstats(self.hapmap3,**self.get_columns(self.hapmap3.columns),verbose=False)