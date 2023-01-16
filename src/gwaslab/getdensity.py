
import pandas as pd
from gwaslab.Log import Log
import gc

def getsignaldensity(insumstats, id="SNPID", chrom="CHR",pos="POS", bwindowsizekb=100, large_number=10000000000,log=Log(),verbose=True):    
    sumstats = insumstats.loc[:,[id,chrom,pos]].copy()
    if verbose:log.write(" -Calculating DENSITY with windowsize of ",bwindowsizekb ," kb")
    stack=[]
    sumstats["TCHR+POS"] = sumstats[chrom]*large_number +  sumstats[pos]
    sumstats = sumstats.sort_values(by=["TCHR+POS"])
    
    for index,row in sumstats.iterrows():
        stack.append([row[id],row["TCHR+POS"],0])  
        for i in range(2,len(stack)+1):
        # closest signals in range	
            if stack[-i][1]>= (row["TCHR+POS"]- 1000*bwindowsizekb):
                #add 1 to both point
                stack[-i][2]+=1
                stack[-1][2]+=1
            else:
            	# closest signals still out out range
                break
    df = pd.DataFrame(stack,columns=["SNPID","TCHR+POS","DENSITY"])
    sumstats["DENSITY"] = df["DENSITY"].values

    # mean and median
    bmean=sumstats["DENSITY"].mean()
    bmedian=sumstats["DENSITY"].median()
    if verbose:log.write(" -Mean : {} signals per {} kb".format(bmean,bwindowsizekb))
    if verbose:log.write(" -Median : {} signals per {} kb".format(bmedian,bwindowsizekb))
    
    sumstats = sumstats.drop("TCHR+POS",axis=1)
    return sumstats["DENSITY"]

def assigndensity(insumstats,
				sig_sumstats,
				id="SNPID", 
				chrom="CHR", 
				pos="POS", 
				bwindowsizekb=100,
				large_number=10000000000,
				log=Log(),verbose=True):
	sumstats = insumstats.loc[:,[id,chrom,pos]].copy()
	sumstats["DENSITY"] = 0
	sumstats["TCHR+POS"] = sumstats[chrom]*large_number +  sumstats[pos]
	sig_sumstats["TCHR+POS"] = sig_sumstats[chrom]*large_number +  sig_sumstats[pos]
	counter = 0
	
	for index,row in sig_sumstats.iterrows():
		counter+=1
		to_add =(sumstats["TCHR+POS"]>=(row["TCHR+POS"]- 1000*bwindowsizekb)) & (sumstats["TCHR+POS"]<=(row["TCHR+POS"]+ 1000*bwindowsizekb))
		sumstats.loc[to_add,"DENSITY"] += 1
		if counter%1000==0:
			if verbose:log.write(" -Processed {} signals".format(counter//1000))
			gc.collect()
	
	return sumstats["DENSITY"]
