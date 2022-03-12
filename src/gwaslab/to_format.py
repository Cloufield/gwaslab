import pandas as pd

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