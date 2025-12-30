from typing import Optional, List, Dict, Any
import pandas as pd
import numpy as np
import gc
import os, psutil

#20220306
def usage() -> float:
    process = psutil.Process(os.getpid())
    return process.memory_info()[0] / float(2 ** 20)

#flip direction
def flip_direction(string: str) -> str:
    flipped_string=""
    for char in string:
        if char=="?":flipped_string+="?"
        elif char=="+":flipped_string+="-"
        elif char=="-":flipped_string+="+"
    return flipped_string

def rsID2chrpos(
    path: str,
    rsid_col: str,
    readargs: Dict[str, Any] = {},
    chr_to_anno: List[int] = [],
    ea: Optional[str] = None,
    nea: Optional[str] = None,
    eaf: Optional[str] = None,
    beta: Optional[str] = None,
    direction: Optional[str] = None,
    chrpos_path_12: List[str] = ["/home/heyunye/mydata/reference/snp/organisms/human_9606_b151_GRCh37p13/VCF/chr/dbsnp_b151_rsid_chrposrefalt_uniq_hg19_chr",".feather"]
) -> pd.DataFrame:
    
    print("Loading into memory :" + path)
    sumstats = pd.read_csv(path,"\s+",**readargs)
    print("Total number of raw input variants:" +str(len(sumstats)))
    
    print("Removing duplicated rsIDs...")
    sumstats.drop_duplicates(rsid_col, keep="first", inplace=True)
    sumstats.set_index(rsid_col,inplace=True)
    print('Total number of variants after removing duplicates(keep first) :' +str(len(sumstats)))
    
    # initiate data columns
    sumstats['#CHROM']=0
    sumstats['#POS']  =0
    sumstats['#POS']  = sumstats['#POS'].astype('int')
    sumstats['#Ref']  =""
    sumstats['#Alt']  =""

    if chr_to_anno:
        # annotate specified chromosomes
        for chrom in chr_to_anno:
            print("Annotating SNPs on chromosome "+str(chrom))
            chrpos_path=chrpos_path_12[0]+str(chrom)+chrpos_path_12[1] # create path to annotation file for each chromosome
            chunks = pd.read_feather(chrpos_path)
            chunks.rename(columns={"Ref":"#Ref","Alt":"#Alt","POS":"#POS"},inplace=True)
            chunks.set_index("ID",inplace=True)
            variants_count=len(chunks)
            print("  - Loading :"+ str(variants_count) +" variants...")
            chunks["#CHROM"] = int(chrom)
            sumstats.update(chunks)
            del chunks
            gc.collect()

    
    else:    
        # annotate all chromosomes
        print("Annotating 1-22 chromosomes...")
        for chrom in range(1,23):
            print("Annotating SNPs on chromosome "+str(chrom))
            chrpos_path=chrpos_path_12[0]+str(chrom)+chrpos_path_12[1] # create path to annotation file for each chromosome
            chunks = pd.read_feather(chrpos_path)
            chunks.rename(columns={"Ref":"#Ref","Alt":"#Alt","POS":"#POS"},inplace=True)
            chunks.set_index("ID",inplace=True)
            variants_count = len(chunks)
            print("  - Loading :"+ str(variants_count) +" variants...")
            chunks["#CHROM"] = str(chrom)
            sumstats.update(chunks)
            del chunks
            gc.collect()
    
    print("Total number of annotated variants:" +str(len(sumstats.loc[sumstats["#POS"]!=0])) )
    print("Total number of unmapped variants:"  +str(len(sumstats.loc[sumstats["#POS"]==0])) )
    
    sumstats["#CHROM"]=sumstats["#CHROM"].astype(int)
    sumstats["#POS"]=sumstats["#POS"].astype(int)
    
#  check alleles
    if ea and nea: 
        sumstats[ea]=sumstats[ea].str.upper()
        sumstats[nea]=sumstats[nea].str.upper()
        sumstats["Allele_match"]=9
        sumstats["Aligned_NEA"]=sumstats[nea]
        sumstats["Aligned_EA"]=sumstats[ea]
        
    if beta:       sumstats["Beta_aligned"]=np.nan
    if eaf:        sumstats["EAF_aligned"] =np.nan
    if direction:  sumstats["Direction_aligned"]=""
    #allele match code:
        #0 match
        #1 need to be flipped
    #match ref nea 
    if ea and nea: 
        sumstats.loc[sumstats[nea]==sumstats["#Ref"],"Allele_match"]=0
        sumstats.loc[sumstats[nea]==sumstats["#Ref"],"Aligned_NEA"] = sumstats.loc[sumstats[nea]==sumstats["#Ref"],nea]
        sumstats.loc[sumstats[nea]==sumstats["#Ref"],"Aligned_EA"]  = sumstats.loc[sumstats[nea]==sumstats["#Ref"],ea]
        print("#Ref-NEA match:"     +str(len(sumstats.loc[sumstats[nea]==sumstats["#Ref"]])))
    #if ea -> ref , need to flip
        sumstats.loc[(sumstats[ea]==sumstats["#Ref"]),"Allele_match"]=1
        sumstats.loc[sumstats[ea]==sumstats["#Ref"],"Aligned_NEA"] = sumstats.loc[sumstats[ea]==sumstats["#Ref"],ea]
        sumstats.loc[sumstats[ea]==sumstats["#Ref"],"Aligned_EA"]  = sumstats.loc[sumstats[ea]==sumstats["#Ref"],nea]
        print("#Ref-NEA not match but EA can be filpped:" +str(len(sumstats.loc[sumstats[ea]==sumstats["#Ref"]])))
    
    print("#Ref-NEA and #Alt-EA not on same strand or not match:" +str(len(sumstats.loc[(sumstats[nea]!=sumstats["#Ref"])&(sumstats[ea]!=sumstats["#Ref"])])))  
    print("Flipping signed statistics for "+str(len(sumstats.loc[sumstats[ea]==sumstats["#Ref"]])))
    
    to_round = []
    if beta:
        print("Flipping beta...")
        sumstats.loc[sumstats["Allele_match"]==0,"Beta_aligned"] =  sumstats.loc[sumstats["Allele_match"]==0,beta]
        sumstats.loc[sumstats["Allele_match"]==1,"Beta_aligned"] = -sumstats.loc[sumstats["Allele_match"]==1,beta]
        sumstats.loc[sumstats["Allele_match"]>1,"Beta_aligned"]  =  sumstats.loc[sumstats["Allele_match"]>1,beta]
        to_round += [beta,"Beta_aligned"] 
    if eaf:
        print("Flipping eaf...")
        sumstats.loc[sumstats["Allele_match"]==0,"EAF_aligned"]  =   sumstats.loc[sumstats["Allele_match"]==0,eaf]
        sumstats.loc[sumstats["Allele_match"]==1,"EAF_aligned"]  = 1-sumstats.loc[sumstats["Allele_match"]==1,eaf]
        sumstats.loc[sumstats["Allele_match"]>1,"EAF_aligned"]   =   sumstats.loc[sumstats["Allele_match"]>1,eaf]
        to_round += [eaf,"EAF_aligned"] 
    if direction:
        print("Flipping direction...")
        sumstats.loc[sumstats["Allele_match"]==0,"Direction_aligned"] =   sumstats.loc[sumstats["Allele_match"]==0,direction]
        sumstats.loc[sumstats["Allele_match"]==1,"Direction_aligned"] =   sumstats.loc[sumstats["Allele_match"]==1,direction].apply(flip_direction)
        sumstats.loc[sumstats["Allele_match"]>1,"Direction_aligned"]  =   sumstats.loc[sumstats["Allele_match"]>1,direction]
    
    print("Rounding beta and frequency to 4 digits...")
    print("Writing to "+path+".rsid")
    sumstats = sumstats.loc[sumstats["#POS"]!=0]
    
    sumstats[to_round] = sumstats[to_round].round(4)
    sumstats.sort_values(["#CHROM","#POS"]).to_csv(path+".rsid","\t",na_rep='\.')
    
    return sumstats