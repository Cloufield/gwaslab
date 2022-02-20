import scipy.stats as ss
import numpy as np

def or2beta(insumstats,OR,OR_lower=None,OR_upper=None,OR_se=None,level=0.95):
    output = insumstats
    output["BETA"]=np.log(sumstats[OR])
    if OR_se:
        output["SE"]=np.log(sumstats[OR]+sumstats[OR_se])-np.log(sumstats[OR])
    elif OR_upper:
        tl,tu = ss.norm.interval(0.95,0)
        output["SE"]=(sumstats[upper] - sumstats[OR])/tu
    elif OR_lower:
        tl,tu = ss.norm.interval(0.95,0)
        output["SE"]=(sumstats[lower] - sumstats[OR])/tl
    return output
    
