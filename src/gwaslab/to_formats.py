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



################################################################################################################################################
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
#################################################################################################################################################
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


###################################################################################################################################################
def tossf(sumstats,path=None,snpid="MARKERNAME", chrom="CHR", pos="POS", ea="EA", nea="NEA", beta="BETA",OR="OR", se="SE" ,p="P" ,n="N",OR_95L="OR_95L",OR_95U="OR_95U",verbose=True):
    '''Hayhurst, J., Buniello, A., Harris, L., Mosaku, A., Chang, C., Gignoux, C. R., ... & Barroso, I. (2022). A community driven GWAS summary statistics standard. bioRxiv.
    '''
    rename_dictionary = {
    snpid: "variant_id",
    rsid: "variant_id",
    chrom: "chromosome",
    pos: "bas_pair_location",
    nea: "other_allele",
    ea: "effect_allele",
    n: "n",
    beta: "beta",
    se: "standard_error",
    p: "p_value",
    info:"info",
    OR:"odds_ratio",
    OR_95L:"ci_lower",
    OR_95U:"ci_upper",
    status:"hm_code"
    }
    if verbose: print("Start fornmatting to GWAS-SSF format...")
    #check mandatory columns
    if beta in sumstats.columns: 
        mandatory=["chromosome","bas_pair_location","effect_allele","other_allele","beta","standard_error","effect_allele_frequency","p_value"]
    elif OR in sumstats.columns: 
        mandatory=["chromosome","bas_pair_location","effect_allele","other_allele","odds_ratio","standard_error","effect_allele_frequency","p_value"]
    encouraged=["variant_id","rsid","info","ci_upper","ci_lower","n"]      
    sumstats = sumstats.rename(columns=rename_dictionary)

    other=[]
    ouput_order=mandatory

    for i in mandatory:
        if i not in sumstats.columns:
            if verbose: print(i+" column is missing... please check again...")
    for i in encoraged:
        if i in sumstats.columns:
            output_order.append(i)
    for i in sumstats.columns:
        if (i not in mandatory) and( i not in encoraged):
            other.append(i)
    if len(other)>1:
        other = other.sort()
    ouput_order = ouput_order + other

    sumstats = sumstats.loc[:,ouput_order]
    if verbose: print("Start fornmatting to GWAS-SSF format...")
    if path is not None: sumstats.to_csv(path,"\t",index=None)
    #if verbose: print("Start outputing metadata...")
    ##    
            
###################################################################################################################################################
def tovcf(sumstats,path=None,snpid="MARKERNAME", chrom="CHR", pos="POS", ea="EA", nea="NEA", beta="BETA", se="SE" ,p="P" ,n="N",verbose=True):
    pass
    
###################################################################################################################################################
def tofmt(sumstats,path=None,suffix=None,fmt=None,cols=[],verbose=True,log=gl.Log(),to_csvargs={}):
    if fmt is not None:
        if verbose: log.write("Start outputting sumstats in "+fmt+" format...")
        if verbose: log.write(" -"+fmt+" format will be loaded...")
        meta_data,rename_dictionary = gl.CommonData.get_format_dict(fmt,inverse=True)
        if verbose:             
            log.write(" -"+fmt+" format meta info:")   
            for key,value in meta_data.items():
                log.write("  -",key," : ",value)
        if verbose: 
            log.write(" -gwaslab to "+fmt+" format dictionary:",)  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            log.write("  - gwaslab keys:",keys) 
            log.write("  - "+fmt+" values:",values) 
            
        
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        ouput_cols = ouput_cols + cols
        sumstats = sumstats.loc[:,ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        path = path + "."+suffix+".tsv.gz"
        if verbose: log.write(" -Output columns:",sumstats.columns)
        if verbose: log.write(" -Output path:",path) 
        if path is not None: sumstats.to_csv(path,sep="\t",index=None,**to_csvargs)
        return sumstats
    else:
        if path is not None: sumstats.to_csv(path,sep="\t",index=None,**to_csvargs)
    
    
    
    
    
    
    
    
    
    
    
    
    