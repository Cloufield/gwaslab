
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
import gc

def getsignaldensity(insumstats, id="SNPID", chrom="CHR",pos="POS", bwindowsizekb=100,log=Log(),verbose=True):    
    log.write("Start to calculate signal DENSITY..." ,verbose=verbose)
    sumstats = insumstats[[id,chrom,pos]].copy()
    log.write(" -Calculating DENSITY with windowsize of ",bwindowsizekb ," kb",verbose=verbose)
    #stack=[]

    large_number = 1000000000
    for i in range(7):
        if insumstats[pos].max()*10 >  large_number:
            large_number = int(large_number * 10)
        else:
            break

    sumstats["TCHR+POS"] = sumstats[chrom]*large_number +  sumstats[pos]
    sumstats = sumstats.sort_values(by=["TCHR+POS"])
    positions = sumstats["TCHR+POS"].values
    #sumstats = sumstats.reset_index()
    density = np.zeros(len(sumstats))

    last_left_pos = 0
    for current_pos, position in enumerate(positions):
        for current_left_pos in range(last_left_pos,current_pos):
            if  positions[current_left_pos] >= position-1000 * bwindowsizekb:
                
                # from p_left to p_current -1, plus 1 density (update the rightside numbers)
                density[current_left_pos:current_pos]+=1
                
                # for p_current, plus p_current - p_left (update the leftside numbers)
                density[current_pos]+= current_pos - current_left_pos
                
                # update left pointer
                last_left_pos = current_left_pos
                break
        
        #stack.append([row[id],row["TCHR+POS"],0])  
        #for i in range(2,len(stack)+1):
        ## closest signals in range	
        #    if stack[-i][1]>= (row["TCHR+POS"]- 1000*bwindowsizekb):
        #        #add 1 to both point
        #        stack[-i][2]+=1
        #        stack[-1][2]+=1
        #    else:
        #    	# closest signals still out out range
        #        break
    
    sumstats["DENSITY"] = density
    sumstats["DENSITY"]  = sumstats["DENSITY"].astype("Int32")
    # mean and median
    bmean = sumstats["DENSITY"].mean()
    bmedian = sumstats["DENSITY"].median()
    bsd = sumstats["DENSITY"].std()
    bmax = sumstats["DENSITY"].max()
    bmaxid = sumstats["DENSITY"].idxmax()

    log.write(" -Mean : {} signals per {} kb".format(bmean,bwindowsizekb),verbose=verbose)
    log.write(" -SD : {}".format(bsd),verbose=verbose)
    log.write(" -Median : {} signals per {} kb".format(bmedian,bwindowsizekb),verbose=verbose)
    log.write(" -Max : {} signals per {} kb at variant(s) {}".format(bmax,bwindowsizekb,sumstats.loc[bmaxid,id]),verbose=verbose)
    
    sumstats = sumstats.drop("TCHR+POS",axis=1)
    log.write("Finished calculating signal DENSITY successfully!",verbose=verbose)
    return sumstats["DENSITY"]

def assigndensity(insumstats,
				sig_sumstats,
				id="SNPID", 
				chrom="CHR", 
				pos="POS", 
				bwindowsizekb=100,
				log=Log(),verbose=True):

    large_number = 1000000000
    for i in range(7):
        if insumstats[pos].max()*10 >  large_number:
            large_number = int(large_number * 10)
        else:
            break
    sumstats = insumstats[[id,chrom,pos]].copy()
    sumstats["DENSITY"] = 0
    sumstats["TCHR+POS"] = sumstats[chrom]*large_number +  sumstats[pos]
    sig_sumstats["TCHR+POS"] = sig_sumstats[chrom]*large_number +  sig_sumstats[pos]
    counter = 0

    for index,row in sig_sumstats.iterrows():
        counter+=1
        to_add =(sumstats["TCHR+POS"]>=(row["TCHR+POS"]- 1000*bwindowsizekb)) & (sumstats["TCHR+POS"]<=(row["TCHR+POS"]+ 1000*bwindowsizekb))
        sumstats.loc[to_add,"DENSITY"] += 1
        if counter%1000==0:
            log.write(" -Processed {} signals".format(counter//1000),verbose=verbose)
            gc.collect()
	
    return sumstats["DENSITY"]
