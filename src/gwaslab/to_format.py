import pandas as pd

def tobed(sumstats,snpid,chrom,pos,ref=None,alt=None,flank=500, to_path=""):
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

