import pandas as pd

def fixchr(sumstats,chrom="CHR",add_prefix="",remove=True, verbose=True):
        # convert to string datatype
        sumstats[chrom] = sumstats[chrom].astype("string") 
        
        # strip prefix
        if verbose: print("Stripping chr prefix...") 
        sumstats[chrom] = sumstats[chrom].str.lstrip("chrCHR_-")
        
        # strip leading zeros
        if verbose: print("Removing leading zeroes...") 
        sumstats[chrom] = sumstats[chrom].str.lstrip("0")
        
        # x,y,mt to X,Y,MT
        if verbose: print("Converting to UPPERCASE...") 
        sumstats[chrom] = sumstats[chrom].str.upper()
        
        chrom_list={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"}
        
        if remove:
            # remove variants with unrecognized chrom
            unrecognized_num = len(sumstats.loc[~sumstats[chrom].isin(chrom_list),:])
            if unrecognized_num>0:
                if verbose: print("Removing variants with unrecognized chrom : "+str(unrecognized_num)) 
            sumstats = sumstats.loc[sumstats[chrom].isin(chrom_list),:]
        
        if add_prefix:
            if verbose: print("Adding prefix : "+ add_prefix+"...")
            sumstats[chrom] = add_prefix + sumstats[chrom]
            
        return sumstats

def fixallele(sumstats,EA="EA", NEA="NEA",verbose=True):
        
        # remove variants with alleles other than actgACTG
        all_var_num = len(sumstats)
        sumstats = sumstats.loc[~sumstats[EA].str.contains("[^actgACTG]"),:]
        sumstats = sumstats.loc[~sumstats[NEA].str.contains("[^actgACTG]"),:]
        remain_var_num = len(sumstats)
        if verbose: print("Removed "+str(all_var_num - remain_var_num)+" variants with alleles that contain bases other than A/C/T/G .")        
        return sumstats
    
def normalizeallele(sumstats,pos="POS",ea="EA" ,nea="NEA",verbose=True):
    if verbose: print("Start to normalize variants...")
    normalized = sumstats.loc[:,[pos,nea,ea]].apply(lambda x: normalizevariant(x[0],x[1],x[2]),axis=1)
    normalized_pd = pd.DataFrame(normalized.to_list(), columns=[pos,nea,ea])
    if verbose:
        changed_num = len(normalized_pd.loc[(sumstats[ea]!=normalized_pd[ea]) | (sumstats[nea]!=normalized_pd[nea]),:])
        print("Number of variants that were changed : ",changed_num )
    sumstats.loc[:,[pos,nea,ea]] = normalized_pd
    return sumstats
    

##############################################################################
def normalizevariant(pos,a,b):
    # https://genome.sph.umich.edu/wiki/Variant_Normalization
    if len(a)==1 or len(b)==1:
        return pos,a,b
    
    pos_change=0
    pointer_a_l, pointer_a_r = 0,len(a)-1
    pointer_b_l, pointer_b_r = 0,len(b)-1

    #remove from right
    for i in range(max(len(a),len(b))):
        if a[pointer_a_r] == b[pointer_b_r]:
            pointer_a_r-=1
            pointer_b_r-=1
        else:
            break
        if pointer_a_r-pointer_a_l==0 or pointer_b_r-pointer_b_l==0:
            return pos,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1]

    # remove from left
    for i in range(max(len(a),len(b))):
        if a[pointer_a_l] == b[pointer_b_l]:
            pointer_a_l+=1
            pointer_b_l+=1
            pos_change+=1
        else:
            break
        if pointer_a_r-pointer_a_l==0 or pointer_b_r-pointer_b_l==0:
            break

    return pos+pos_change,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1]