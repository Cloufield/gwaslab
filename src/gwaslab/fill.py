import pandas as pd
import numpy as np
import scipy.stats as ss

def preformat(
    inpath,
    outpath=None,
    snpid=None,
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
    verbose=False,
    readargs={"sep":"\s+"},
    writeargs={"sep":"\t","index":None}
    ):
    
    usecols=[]
    dtype={}
    if snpid: 
        usecols.append(snpid)
    if chrom: 
        usecols.append(chrom)
    if pos: 
        usecols.append(pos)
    if ea: 
        usecols.append(ea)
    if nea: 
        usecols.append(nea)
    if eaf: 
        usecols.append(eaf)
    if n: 
        usecols.append(n)
    if beta: 
        usecols.append(beta)
    if se: 
        usecols.append(se)
    if chisq: 
        usecols.append(chisq)
    if z:
        usecols.append(z)
    if p: 
        usecols.append(p)
    if mlog10p: 
        usecols.append(mlog10p)
    if info: 
        usecols.append(info)
    if OR: 
        usecols.append(OR)
    if OR_se:
        usecols.append(OR_se)
    if OR_95L: 
        usecols.append(OR_95L)
    if OR_95U: 
        usecols.append(OR_95U)
    if other: usecols = usecols + other
    
    
    rename_dictionary={
                                             snpid:   "MARKERNAME",
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
                                             OR_95U:  "OR_95U"
                                            }
    usecols
    print("Loading file :"+inpath)
    print("Reading columns :", usecols)
    print("Renaming columns to :", list(map(lambda x:rename_dictionary[x], usecols)))
    print("Preformatted sumstats will be written to :" + outpath)
    if verbose: print("Verbose mode: filling convertable columns...")
    
    
    if outpath:
        #read file in chunks
        chunks = pd.read_table(inpath,usecols=usecols,chunksize=1000000,**readargs)

        #variants counter
        variants_count=0

        for chunk_num,sumstats in enumerate(chunks):
            variants_count+=len(sumstats)
            print("Loading chunk number " + str(chunk_num+1) +":"+ str(variants_count) +" variants...")

            #rename
            sumstats = sumstats.rename(columns=rename_dictionary)

            if verbose:    
                if beta and se:
                    sumstats["OR"]     = np.exp(sumstats["BETA"])
                    sumstats["OR_95L"] = np.exp(sumstats["BETA"]-ss.norm.ppf(0.975)*sumstats["SE"])
                    sumstats["OR_95U"] = np.exp(sumstats["BETA"]+ss.norm.ppf(0.975)*sumstats["SE"])

                elif OR and OR_se:
                    if not beta: sumstats["BETA"] = np.log(sumstats["OR"])
                    sumstats["SE"]   = np.log(sumstats["OR"]+sumstats["OR_SE"]) - np.log(sumstats["OR"])

                elif OR and OR_95U: 
                    if not beta: sumstats["BETA"] = np.log(sumstats["OR"])
                    tl,tu = ss.norm.interval(0.95,0)
                    sumstats["SE"]=(sumstats["OR_95U"] - sumstats["OR"])/tu

                elif OR and OR_95L:
                    if not beta: sumstats["BETA"] = np.log(sumstats["OR"])
                    tl,tu = ss.norm.interval(0.95,0)
                    sumstats["SE"]=(sumstats["OR_95L"] -  sumstats["OR"])/tl

                if not chisq:
                    if not z:
                        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]
                    sumstats["CHISQ"] = (sumstats["Z"])**2

                if p:
                    sumstats["MLOG10P"] = -np.log10(sumstats["P"])
                elif mlog10p:
                    sumstats["P"] = np.power(10,-sumstats["MLOG10P"])

            order=["MARKERNAME","CHR","POS","EA","NEA","EAF","N","BETA","SE","Z","CHISQ","P","MLOG10P","INFO", "OR","OR_SE","OR_95L","OR_95U"]+other
            output_columns=[]

            if verbose:
                for i in order:
                    if i in sumstats.columns:
                        output_columns.append(i)
            else:
                output_columns=usecols
            if chunk_num<1: 
                sumstats.to_csv(outpath,**writeargs,header=True)
            else:
                sumstats.to_csv(outpath,**writeargs,mode='a',header=False)
    else:
        
        sumstats = pd.read_table(inpath,usecols=usecols,**readargs)
        
        
        
        
    
    print("Preformated "+ str(variants_count) +" successfully!")
    return 
    