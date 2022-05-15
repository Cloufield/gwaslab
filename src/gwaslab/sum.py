import pandas as pd

def summarize(sumstats):
    output = {}
    output["Current data shape"] = AF.data.shape
    output["Current columns"] = AF.data.columns.values
    output["Total Variants"] = len(sumstats)
    if "EAF" in sumstats.columns:
        maf = sumstats["EAF"].copy()
        maf[maf>0.5] = 1 - maf[maf>0.5]
        output["Common Variants"] =  sum(maf>0.05)
        output["Low-frequency Variants"] = sum((maf>0.01 )& (maf<0.05))
        output["Rare Variants"] = sum(maf<0.01)
        if ("N" in sumstats.columns) and ("P" in sumstats.columns):
            output["Min MAC"] = round(min(output["Rare Variants"] * maf))
    if "P" in sumstats.columns:
        output["Minimum P"] = sumstats["P"].min()
    return output

def sum_status():
    pass