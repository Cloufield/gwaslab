import pandas as pd
import gwaslab as gl

#20220306
class Sumstats():

    # Sumstats.data will store a pandas df
    data = pd.DataFrame()

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
                 **arguments):

        #renaming dictionary
        rename_dictionary = {
            snpid: "MARKERNAME",
            rsid: "rsID",
            chrom: "CHR",
            pos: "POS",
            nea: "NEA",
            ea: "EA",
            eaf: "EAF",
            n: "N",
            beta: "BETA",
            se: "SE",
            z: "Z",
            chisq: "CHISQ",
            p: "P",
            mlog10p: "MLOG10P",
            info: "INFO",
            OR: "OR",
            OR_se: "OR_SE",
            OR_95L: "OR_95L",
            OR_95U: "OR_95U"
        }

        usecols = []
        datatype_dictionary = {}
        if snpid:
            usecols.append(snpid)
            datatype_dictionary[snpid] = "string"
        if rsid:
            usecols.append(rsid)
            datatype_dictionary[rsid] = "string"
        if chrom:
            usecols.append(chrom)
            datatype_dictionary[chrom] = "string"
        if pos:
            usecols.append(pos)
            datatype_dictionary[pos] = "int"
        if ea:
            usecols.append(ea)
            datatype_dictionary[ea] = "string"
        if nea:
            usecols.append(nea)
            datatype_dictionary[nea] = "string"
        if eaf:
            usecols.append(eaf)
            datatype_dictionary[eaf] = "float"
        if n and (type(n) is str):
            usecols.append(n)
            datatype_dictionary[n] = "int"
        if beta:
            usecols.append(beta)
            datatype_dictionary[beta] = "float"
        if se:
            usecols.append(se)
            datatype_dictionary[se] = "float"
        if chisq:
            usecols.append(chisq)
            datatype_dictionary[chisq] = "float"
        if z:
            usecols.append(z)
            datatype_dictionary[z] = "float"
        if p:
            usecols.append(p)
            datatype_dictionary[p] = "float"
        if mlog10p:
            usecols.append(mlog10p)
            datatype_dictionary[mlog10p] = "float"
        if info:
            usecols.append(info)
            datatype_dictionary[info] = "float"
        if OR:
            usecols.append(OR)
            datatype_dictionary[OR] = "float"
        if OR_se:
            usecols.append(OR_se)
            datatype_dictionary[OR_se] = "float"
        if OR_95L:
            usecols.append(OR_95L)
            datatype_dictionary[OR_95L] = "float"
        if OR_95U:
            usecols.append(OR_95U)
            datatype_dictionary[OR_95U] = "float"
        if other:
            usecols = usecols + other

#loading data ##################################################################################
        if type(sumstats) is str:
            ## loading data from path
            inpath = sumstats
            if verbose:
                print("Initiating Sumstats Object form file :" + inpath)
            self.data = pd.read_table(inpath,
                                      usecols=usecols,
                                      sep="\s+",
                                      dtype=datatype_dictionary)
        elif type(sumstats) is pd.DataFrame:
            ## loading data from dataframe
            if verbose:
                print(
                    "Initiating Sumstats Object from pandas dataframe object ..."
                )
            self.data = sumstats.loc[:, usecols].astype(datatype_dictionary)

        ##renaming columns
        conveted_columns = list(map(lambda x: rename_dictionary[x], usecols))
        datatype_columns = list(map(lambda x: datatype_dictionary[x], usecols))
        
        if verbose: print("Reading columns          :", ",".join(usecols))
        if verbose: print("Renaming columns to      :", ",".join(conveted_columns))
        if verbose: print("Datatype of each columns :", ",".join(datatype_columns))
        if verbose: print("Current dataframe shape  : Rows ", len(self.data)," x ",len(self.data.columns)," Columns")
        self.data = self.data.rename(columns=rename_dictionary)

        if type(n) is int:
            self.data["N"] = n

        if not snpid:
            self.data["MARKERNAME"] = self.data["CHR"].astype(
                "string") + ":" + self.data["POS"].astype("string")

    def __getattr__(self, name):
        if name == 'gc':
            return self.get_gc(verbose=False)
        else:
            raise AttributeError        
            
#### healper #################################################################################

    def get_columns(self, columns=None):
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
        if columns is not None:
            col_subset = {dic[key]: key for key in columns}
        else:
            col_subset = {dic[key]: key for key in self.data.columns}
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
        is_hla = (self.data["CHR"] == 6) & (self.data["POS"] > 25000000) & (
            self.data["POS"] < 34000000)
        output = self.data.loc[~is_hla]  #not in hla
        return Sumstats(output, **self.get_columns(), verbose=False)

    def extract_hla(self, **args):
        is_hla = (self.data["CHR"] == 6) & (self.data["POS"] > 25000000) & (
            self.data["POS"] < 34000000)
        output = self.data.loc[is_hla]  #in hla
        return Sumstats(output, **self.get_columns(), verbose=False)
    
    def lift_over(self,from_build="hg19" ,to_build="hg38",remove_unmapped=True):
        output = gl.liftover(self.data,chrom="CHR",pos="POS",
                                from_build=from_build,to_build=to_build,
                                remove_unmapped=remove_unmapped)
        return Sumstats(output, **self.get_columns(), verbose=False)
        
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
                      == 6) & (self.hapmap3["POS"] >
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
