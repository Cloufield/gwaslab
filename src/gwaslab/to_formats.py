import pandas as pd
import gwaslab as gl

def toldsc(sumstats, path, verbose=True,log=gl.Log(),to_csvargs={}):
    '''
    Source: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
    
    The ldsc .sumstats format requires six pieces of information for each SNP:
    A unique identifier (e.g., the rs number)
    Allele 1 (effect allele)
    Allele 2 (non-effect allele)
    Sample size (which often varies from SNP to SNP)
    A P-value
    A signed summary statistic (beta, OR, log odds, Z-score, etc)
    
    gwaslab will output a minimum set of columns
    '''
    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2"})
    # select signed 
    for col in ["BETA","OR","Z"]:
        if col in sumstats.columns:
            signed = col
            break
    
    output_cols=["rsID","A1","A2","P",signed]
    #if output N
    if "N" in sumstats.columns:
        output_cols.append("N")
        
    sumstats.loc[:,output_cols].to_csv(path, sep="\t", index=False,**to_csvargs)
    if verbose: log.write(" -Saving sumstats in ldsc format at Hapmap3 variants to : ")
    if verbose: log.write(" -" + path)   



###################################################################

def tobed(sumstats,snpid,chrom,pos,ref=None,alt=None,flank=500, to_path="",verbose=True):
    if verbose: print("Start fornmatting to bed...")
    output = pd.DataFrame()
    flank = int(flank)
    output["Chrom"]= sumstats[chrom]
    output["Start"]= sumstats[pos].astype("int") - flank
    output["End"]=   sumstats[pos].astype("int") + flank
    if ref and alt:
        output["Ref"]=   sumstats[ref]
        output["Alt"]=   sumstats[alt]
    output["Name"]=  sumstats[snpid]
    return output

def tofuma(sumstats,path=None,snpid="MARKERNAME", chrom="CHR", pos="POS", ea="EA", nea="NEA", beta="BETA", se="SE" ,p="P" ,n="N",verbose=True):
        rename_dictionary = {
        snpid: "MARKERNAME",
        chrom: "CHR",
        pos: "POS",
        nea: "A2",
        ea: "A1",
        n: "N",
        beta: "BETA",
        se: "SE",
        p: "P"
        }
        if verbose: print("Start fornmatting to fuma input format...")
        sumstats = sumstats.rename(columns=rename_dictionary)
        sumstats = sumstats.loc[:,list(rename_dictionary.values())]
        if path is not None: sumstats.to_csv(path," ",index=None)
        return sumstats