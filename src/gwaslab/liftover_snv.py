import pandas as pd
from liftover import get_lifter

def liftover_snv(row,converter):
    chrom = row[0]  
    pos_0_based = int(row[1]) - 1
    results = converter[chrom][pos_0_based]
    if converter[chrom][pos_0_based]:
        # return chrom, pos_1_based
        return results[0][0].strip("chr"),results[0][1]+1
    else:
        return "unmapped","unmapped"
    
def liftover(insumstats, 
             chrom, 
             pos, 
             from_build="hg19", 
             to_build="hg38",
             remove_unmapped=True):
    
    sumstats = insumstats.copy()
    print("Creating converter : " + from_build +" to "+ to_build)
    converter = get_lifter(from_build, to_build)
    
    print("Converting variants: "+str(len(sumstats)))
    sumstats.loc[:,"CHR_POS_"+to_build] = sumstats.loc[:,["CHR","POS"]].apply(lambda x: liftover_snv(x[["CHR","POS"]],converter),axis=1)
    sumstats.loc[:,"CHR_"+to_build] = sumstats.loc[:,"CHR_POS_"+to_build].apply(lambda x :str(x[0]))
    sumstats.loc[:,"POS_"+to_build] = sumstats.loc[:,"CHR_POS_"+to_build].apply(lambda x :str(x[1]))
    
    map_num= len(sumstats.loc[sumstats["POS_"+to_build]!="unmapped",:])
    unmap_num = len(sumstats.loc[sumstats["POS_"+to_build]=="unmapped",:])
    
    print("Converted variants: "+str(map_num))
    if remove_unmapped:
        print("Remove unmapped variants: "+str(unmap_num))
        output = sumstats.loc[sumstats["POS_"+to_build]!="unmapped",:].copy()
        output.loc[:,["POS_"+to_build]] = output.loc[:,"POS_"+to_build].astype("int")
    else:
        output = sumstats.copy()
    return output.drop("CHR_POS_"+to_build,axis=1)
