import re
import gwaslab as gl
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from liftover import get_lifter
from functools import partial


#fixID
#rsidtochrpos
#remove_dup
#fixchr
#fixpos
#fixallele
#normalizeallele
#normalizevariant
#checkref
#sanitycheckstats
#flipallelestats
#parallelizeassignrsid
#sortcoordinate
#sortcolumn

###############################################################################################################
#20220428 
def fixID(sumstats,
       snpid="SNPID",rsid="rsID",chrom="CHR",pos="POS",nea="NEA",ea="EA",status="STATUS",
       fixchrpos=True,fixid=False,fixeanea=False,fixeanea_flip=False,overwrite=False,remove=True,verbose=True,log=gl.Log()):  
    '''
    1. Fix SNPid
    2. fix chr and pos using snpid
    3. checking rsid and chr:pos:nea:ea
    '''
    if verbose: log.write("Start to check IDs...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    

    
    ############################  checking ###################################################  
    if snpid in sumstats.columns: 
        if verbose: log.write(" -Checking if SNPID is chr:pos:ref:alt...(separator: - ,: , _)")
        is_chrposrefalt = sumstats[snpid].str.match(r'(chr)?([0-9XYMT]+)[:_-]([0-9]+)[:_-]([ATCG]+)[:_-]([ATCG]+)', case=False, flags=0, na=False)
    
    if rsid in sumstats.columns: 
        if verbose: log.write(" -Checking if rsID is rsxxxxxx or RSxxxxxxx...")
        is_rsid = sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)
        
        if verbose: log.write(" -Checking if chr:pos:ref:alt is mixed in rsID column ...")
        is_rs_chrpos = sumstats[rsid].str.match(r'(chr)?([0-9XYMT]+)[:_-]([0-9]+)[:_-]([ATCG]+)[:_-]([ATCG]+)', case=False, flags=0, na=False)
        
        if verbose: log.write(" -Number of chr:pos:ref:alt mixed in rsID column :",sum(is_rs_chrpos))
        if verbose: log.write(" -Number of Unrecognized rsID :",len(sumstats) - sum(is_rs_chrpos) - sum(is_rsid) ) 
        if verbose: log.write(" -A look at the unrecognized rsID :",set(sumstats.loc[(~is_rsid)&(~is_rs_chrpos),rsid].head()),"...") 
      
    ############################  fixing chr pos###################################################  
    if fixchrpos is True:
        if snpid in sumstats.columns: 
            if verbose: log.write(" -Fixing CHR and POS...")
            if overwrite is True: 
                if verbose: log.write(" -Overwrite is applied...")
                to_fix = is_chrposrefalt
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                to_fix = is_chrposrefalt&sumstats[chrom].isna()&sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants. ...")
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                if verbose: log.write(" -Initiating CHR columns...")
                sumstats.loc[:,chrom]=pd.Series(dtype="string")   
                to_fix = is_chrposrefalt&sumstats[chrom].isna()&sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants. ...")
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                if verbose: log.write(" -Initiating CHR and POS column...")
                sumstats.loc[:,pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt&sumstats[chrom].isna()&sumstats[pos].isna() 
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants. ...")
            else:
                if verbose: log.write(" -Initiating CHR and POS columns...")
                sumstats.loc[:,chrom]=pd.Series(dtype="string")   
                sumstats.loc[:,pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants. ...")   
            if sum(to_fix)>0:
                if verbose: log.write(" -Filling CHR and POS columns using valid SNPID's chr:pos...")
                sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[0].strip("chrCHR")).astype("string")
                sumstats.loc[to_fix,pos] =np.floor(pd.to_numeric(sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]), errors='coerce')).astype('Int64')
                                            
        if rsid in sumstats.columns:
            if verbose: log.write(" -Fixing CHR and POS using chr:pos:ref:alt format variants in rsID column...")
            if overwrite is True: 
                if verbose: log.write(" -Overwrite is applied...")
                to_fix = is_rs_chrpos
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                to_fix = is_rs_chrpos&sumstats[chrom].isna()&sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants ...")
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                if verbose: log.write(" -Initiating CHR columns...")
                sumstats.loc[:,chrom]=pd.Series(dtype="string")   
                to_fix = is_rs_chrpos&sumstats[chrom].isna()&sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants ...")
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                if verbose: log.write(" -Initiating CHR and POS column...")
                sumstats.loc[:,pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos&sumstats[chrom].isna()&sumstats[pos].isna() 
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants ...")
            else:
                if verbose: log.write(" -Initiating CHR and POS columns...")
                sumstats.loc[:,chrom]=pd.Series(dtype="string")   
                sumstats.loc[:,pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants ...")   
            
            if sum(to_fix)>0:    
                if verbose: log.write(" -Filling CHR and POS columns using chr:pos:ref:alt format variants in rsID column...")
                sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[0]).astype("string")
                sumstats.loc[to_fix,pos] = np.floor(pd.to_numeric(sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[1]), errors='coerce')).astype('Int64')
            
    ############################  fixing chr pos###################################################   
    if fixeanea is True:
        if verbose: log.write(" -Warning: Please make sure a1 is ref or not in Chr:pos:a1:a2")
        if overwrite is True:
            if verbose: log.write(" -Overwrite is applied...")
            to_fix = is_chrposrefalt
        elif (nea in sumstats.columns) and (nea in sumstats.columns):
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
        elif (nea in sumstats.columns) and (ea not in sumstats.columns):
            if verbose: log.write(" -Initiating EA columns...")
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
        elif (nea not in sumstats.columns) and (ea in sumstats.columns):
            if verbose: log.write(" -Initiating NEA columns...")
            sumstats[nea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
        else:
            if verbose: log.write(" -Initiating EA and NEA columns...")
            sumstats[nea]=pd.Series(dtype="string")
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt
            if sum(to_fix)>0: 
                if verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                    
        if sum(to_fix)>0:    
            if verbose: log.write(" -Filling "+str(sum(to_fix))+" EA and NEA columns using SNPID's chr:pos:nea:ea...")
            
            if fixeanea_flip is True:
                if verbose: log.write(" -Flipped : chr:pos:a1:a2...a1->EA , a2->NEA ")
                sumstats.loc[to_fix,ea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[2]).astype("string")  
                sumstats.loc[to_fix,nea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[3]).astype("string")
            else:
                if verbose: log.write(" -Chr:pos:a1:a2...a1->EA , a2->NEA ")
                sumstats.loc[to_fix,ea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[3]).astype("string")  
                sumstats.loc[to_fix,nea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[2]).astype("string")
                
            #sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]).astype("string")
            #sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[1]).astype("Int64")
    ############################  fixing id ################################################### 
    if fixid is True:
        if snpid not in sumstats.columns: 
        # initiate a SNPID column
            sumstats.loc[:,snpid]=pd.Series(dtype="string")
        if (chrom in sumstats.columns) and (pos in sumstats.columns):
            pre_number=sum(sumstats[snpid].isna())                                
            if overwrite is False:
                to_fix = sumstats[snpid].isna()
            else:
                to_fix = sumstats[snpid].isna() | sumstats[snpid].notnull()
            if (ea in sumstats.columns) and (nea in sumstats.columns):
                pattern = r"\d\d1[123]\d"  
                matched_index = sumstats[status].str.match(pattern) 
                to_full_fix = matched_index & to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull() & sumstats[ea].notnull()& sumstats[nea].notnull()
                to_part_fix = (~sumstats[status].str.match(pattern)) & to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull()
                if sum(to_full_fix)>0:
                    sumstats.loc[to_full_fix,snpid] = sumstats.loc[to_full_fix,chrom].astype("string") + ":"+sumstats.loc[to_full_fix,pos].astype("string") +":"+ sumstats.loc[to_full_fix,nea].astype("string") +":"+ sumstats.loc[to_full_fix,ea].astype("string")
                if verbose: log.write(" -Filling "+str(sum(to_part_fix)) +" SNPID using CHR POS...")
                if sum(to_part_fix)>0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")     
            else:
                to_part_fix = to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull()
                if verbose: log.write(" -Filling "+str(sum(to_part_fix)) +" SNPID using CHR POS...")
                if sum(to_part_fix)>0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string") 
            after_number=sum(sumstats[snpid].isna())
            if verbose: log.write(" -Fixed "+ str(pre_number - after_number) +" variants ID...")
        elif verbose: log.write(" -ID Unfixable: no CHR and POS columns or no SNPID. ")
    return sumstats

###############################################################################################################

    
###############################################################################################################
#20220428 
def removedup(sumstats,snpid="SNPID",ea="EA",nea="NEA",rsid="rsID",keep='first',verbose=True,log=gl.Log()):
    '''
    remove duplicate SNPs based on  1. SNPID, EA, and NEA or 2. rsID for non-NA variants
    '''
    if verbose: log.write("Start to remove duplicate...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    total_number = len(sumstats)   
    pre_number =len(sumstats)   
    if snpid in sumstats.columns:
        sumstats = sumstats.loc[sumstats[snpid].isna() | (~sumstats.duplicated(subset=[snpid,ea,nea], keep=keep)),:]
        
        after_number=len(sumstats)   
        if verbose:  log.write(" -Removed ",pre_number -after_number ," based on SNPID, EA, and NEA...")
    
    pre_number =len(sumstats) 
    if rsid in sumstats.columns:
        sumstats = sumstats.loc[sumstats[rsid].isna() | (~sumstats.duplicated(subset=rsid, keep=keep)),:]
        after_number=len(sumstats)   
        if verbose:  log.write(" -Removed ",pre_number -after_number ," based on rsID...")
    
    after_number=len(sumstats)   
    if verbose:  log.write(" -Removed ",total_number -after_number," duplicates in total.")
    return sumstats

###############################################################################################################
#20220423
def fixchr(sumstats,chrom="CHR",add_prefix="",remove=False, verbose=True,log=gl.Log()):
        if verbose: log.write("Start to fix chromosome notation...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
        if chrom not in sumstats.columns:
            raise ValueError("Chromosome column was not detected.")
        
        is_chr_fixable = sumstats[chrom].str.match(r'(chr)?([012][0-9]|[0-9]|X|Y|MT)', case=False, flags=0, na=False)
        if sum(~is_chr_fixable)>0 and verbose: 
            log.write(" -Unrecognized chromosome notations :",set(sumstats.loc[~is_chr_fixable,chrom].head()),"...") 
        if verbose: log.write(" -Number of fixable chr notation :",sum(is_chr_fixable) )

        # convert to string datatype
        sumstats.loc[:,chrom] = sumstats.loc[:,chrom].astype("string") 
        
        # strip prefix
        if verbose: log.write(" -Stripping chr prefix if exists...") 
        sumstats.loc[:,chrom] = sumstats.loc[:,chrom].str.lstrip("chrCHR_-.")
        
        # strip leading zeros
        if verbose: log.write(" -Removing leading zeroes if exists...") 
        sumstats.loc[:,chrom] = sumstats.loc[:,chrom].str.lstrip("0")
        
        # x,y,mt to X,Y,MT
        if verbose: log.write(" -Converting to string datatype and UPPERCASE...") 
        sumstats.loc[:,chrom] =sumstats.loc[:,chrom].str.upper()
        
        chrom_list = gl.get_chr_list() #bottom 
        unrecognized_num = sum(~sumstats[chrom].isin(chrom_list))
        
        if unrecognized_num>0 and verbose is True: 
            log.write(" -Unrecognized chromosome notations :" , set(sumstats.loc[~sumstats[chrom].isin(chrom_list),chrom].head()) )
        elif verbose:
            log.write(" -No unrecognized chromosome notations...")
            
        if (remove is True) and unrecognized_num>0:
            # remove variants with unrecognized chrom 
            if verbose: log.write(" -Removed "+ str(unrecognized_num)+ " variants with unrecognized chromosome notations.") 
            #sumstats = sumstats.loc[sumstats.index[sumstats[chrom].isin(chrom_list)],:]
            good_chr = sumstats[chrom].isin(chrom_list)
            sumstats = sumstats.loc[good_chr, :].copy()
        if add_prefix:
            if verbose: log.write(" -Adding prefix : "+ add_prefix+"...")
            sumstats.loc[:,chrom] = add_prefix + sumstats[chrom]
            
        return sumstats
    
    
###############################################################################################################    
#20220426
def fixpos(sumstats,pos="POS",remove=False, verbose=True,log=gl.Log()):
        if verbose: log.write("Start to fix basepair positions...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
        
        all_var_num = len(sumstats)
       
        #convert to numeric
        sumstats.loc[:,pos] = np.floor(pd.to_numeric(sumstats.loc[:,pos], errors='coerce')).astype('Int64')
        
        #remove na
        if remove is True: 
            sumstats = sumstats.loc[~sumstats[pos].isna(),:]
            remain_var_num = len(sumstats)
            if verbose: log.write(" -Removed "+str(all_var_num - remain_var_num)+" variants with bad positions.")        
        
        if verbose: log.write(" -Converted all position to datatype Int64.")
        if verbose: log.write(" -Fixed basepair position successfully.")
        return sumstats
    
###############################################################################################################    
#20220426
def fixallele(sumstats,ea="EA", nea="NEA",verbose=True,log=gl.Log()):
        # remove variants with alleles other than actgACTG
        if verbose: log.write("Start to fix alleles...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
        
        all_var_num = len(sumstats)

        good_ea  = ~sumstats[ea].str.contains("[^actgACTG]")
        good_nea = ~sumstats[nea].str.contains("[^actgACTG]")
        exclude  = sum(sumstats[ea].str.contains("[^actgACTG]") | sumstats[nea].str.contains("[^actgACTG]"))
        
        sumstats = sumstats.loc[good_ea & good_nea,:]
        
        sumstats[ea] = sumstats[ea].str.upper()
        sumstats[nea]= sumstats[nea].str.upper()

        remain_var_num = len(sumstats)
        if verbose: log.write(" -Removed "+str(exclude)+" variants with alleles that contain bases other than A/C/T/G .")        
        if verbose: log.write(" -Converted all bases to string datatype and UPPERCASE.")
        if verbose: log.write(" -Fixed allele successfully.")
        return sumstats

###############################################################################################################   
#20220426

def parallelnormalizeallele(sumstats,pos="POS",ea="EA" ,nea="NEA",status="STATUS",n_cores=1,verbose=True,log=gl.Log()):
    if verbose: log.write("Start to normalize variants...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    
    df_split = np.array_split(sumstats.loc[:,[pos,nea,ea,status]], n_cores)
    pool = Pool(n_cores)
    #df = pd.concat(pool.starmap(func, df_split))
    func=normalizeallele
    normalized_pd = pd.concat(pool.map(partial(func,pos=pos,ea=ea,nea=nea,status=status),df_split))
    pool.close()
    pool.join()
    
    if verbose:
        changed_num = len(normalized_pd.loc[(sumstats[ea]!=normalized_pd[ea]) | (sumstats[nea]!=normalized_pd[nea]),:])
        if changed_num>0:
            log.write(" -Not normalized allele:",end="")
            for i in sumstats.loc[(sumstats[ea]!=normalized_pd[ea]) | (sumstats[nea]!=normalized_pd[nea]),[ea,nea]].head().values:
                log.write(i,end="",show_time=False)
            log.write("... \n",end="",show_time=False)         
            log.write(" -Modified "+str(changed_num) +" variants according to parsimony and left alignment principal.")
        else:
            log.write(" -All variants are already normalized..")
    
    sumstats.update(normalized_pd,overwrite=True)
    sumstats.loc[:,pos] = np.floor(pd.to_numeric(sumstats.loc[:,pos], errors='coerce')).astype('Int64')
    sumstats.loc[:,[nea,ea,status]] = sumstats.loc[:,[nea,ea,status]].astype("string")
    return sumstats

def normalizeallele(sumstats,pos="POS",ea="EA" ,nea="NEA",status="STATUS"):
       
    normalized = sumstats.loc[:,[pos,nea,ea,status]].apply(lambda x: normalizevariant(x[0],x[1],x[2],x[3]),axis=1)
    normalized_pd = pd.DataFrame(normalized.to_list(), columns=[pos,nea,ea,status],index=sumstats.index)
    return normalized_pd
def normalizevariant(pos,a,b,status):
    # https://genome.sph.umich.edu/wiki/Variant_Normalization
    # a - ref - nea    starting -> pos
    # b - alt - ea
    # status xx1xx /xx2xx /xx9xx
    
    status_pre=status[:2]
    status_end=status[3:]
    if len(a)==1 or len(b)==1:
        return pos,a,b,status_pre+"1"+status_end
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
            return pos,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"1"+status_end

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
    return pos+pos_change,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"1"+status_end
###############################################################################################################

###############################################################################################################
#20220426
def sanitycheckstats(sumstats,coltocheck=["P","BETA","SE","EAF","CHISQ","N","OR","OR_95L","OR_95U"],verbose=True,log=gl.Log()):
    '''
    Sanity check:
        N:      Int64    , N>0 , 
        EAF:    float64  , 0<= EAF <=1, 
        P:      float64  , 0< P <5e-300, 
        BETA:   float64  , abs(BETA) <10
        SE:     float64  , SE >0
        OR:     float64  , OR>0
        OR_95L: float64  , OR_95L>0
        OR_95U: float64  , OR_95L>0
        INFO:   float64  , INFO>0
    '''
    ## add direction
    
    if verbose: log.write("Start sanity check for statistics ...") 
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    cols_to_check=[]
    oringinal_number=len(sumstats)
    
    pre_number=len(sumstats)
    if "N" in coltocheck and "N" in sumstats.columns:
        cols_to_check.append("N")
        if verbose: log.write(" -Checking if N is Int64 and N>0 ...") 
        sumstats.loc[:,"N"] = np.floor(pd.to_numeric(sumstats.loc[:,"N"], errors='coerce')).astype("Int64")
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad N.") 
            
    pre_number=len(sumstats)    
    if "EAF" in coltocheck and "EAF" in sumstats.columns:
        cols_to_check.append("EAF")
        if verbose: log.write(" -Checking if 0<= EAF <=1 ...") 
        sumstats.loc[:,"EAF"] = pd.to_numeric(sumstats.loc[:,"EAF"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["EAF"]>=0) & (sumstats["EAF"]<=1),:]
        sumstats.loc[:,"EAF"] = sumstats.loc[:,"EAF"].round(4)
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad EAF.") 
    
    pre_number=len(sumstats)    
    if "CHISQ" in coltocheck and "CHISQ" in sumstats.columns:
        cols_to_check.append("CHISQ")
        if verbose: log.write(" -Checking if CHISQ>0 ...") 
        sumstats.loc[:,"CHISQ"] = np.floor(pd.to_numeric(sumstats.loc[:,"CHISQ"], errors='coerce'))
        sumstats = sumstats.loc[(sumstats["CHISQ"]>=0),:]
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad CHISQ.") 
            
    pre_number=len(sumstats) 
    if "P" in coltocheck and "P" in sumstats.columns:
        cols_to_check.append("P")
        if verbose: log.write(" -Checking if 0< P <5e-300 ...") 
        sumstats.loc[:,"P"] = pd.to_numeric(sumstats.loc[:,"P"], errors='coerce')
        sumstats = sumstats.loc[sumstats["P"]>=0,:]
        sumstats.loc[sumstats["P"]<5e-300,"P"] = 5e-300
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad P.") 
            
    pre_number=len(sumstats)    
    if "BETA" in coltocheck and "BETA" in sumstats.columns:
        cols_to_check.append("BETA")
        if verbose: log.write(" -Checking if abs(BETA)<10 ...") 
        sumstats.loc[:,"BETA"] = pd.to_numeric(sumstats.loc[:,"BETA"], errors='coerce')
        sumstats = sumstats.loc[np.abs(sumstats["BETA"])<10,:]
        sumstats.loc[:,"BETA"] = sumstats.loc[:,"BETA"].round(4)
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad BETA.") 
            
    pre_number=len(sumstats)    
    if "SE" in coltocheck and "SE" in sumstats.columns:
        cols_to_check.append("SE")
        if verbose: log.write(" -Checking if SE >0 ...") 
        sumstats.loc[:,"SE"] = pd.to_numeric(sumstats.loc[:,"SE"], errors='coerce')
        sumstats = sumstats.loc[sumstats["SE"]>0,:]
        sumstats.loc[:,"SE"] = sumstats.loc[:,"SE"].round(4)
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad SE.") 
            
    pre_number=len(sumstats)    
    if "OR" in coltocheck and "OR" in sumstats.columns:
        cols_to_check.append("OR")
        if verbose: log.write(" -Checking if 0<OR<10 ...") 
        sumstats.loc[:,"OR"] = pd.to_numeric(sumstats.loc[:,"OR"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["OR"]>0) &(sumstats["OR"]<10),:]
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad OR.") 
            
    pre_number=len(sumstats)   
    if "OR_95L" in coltocheck and "OR_95L" in sumstats.columns:
        cols_to_check.append("OR_95L")
        if verbose: log.write(" -Checking if OR_95L>0 ...") 
        sumstats.loc[:,"OR_95L"] = pd.to_numeric(sumstats.loc[:,"OR_95L"], errors='coerce')
        sumstats = sumstats.loc[sumstats["OR_95L"]>0,:]
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad OR_95L.") 
            
    pre_number=len(sumstats)    
    if "OR_95U" in coltocheck and "OR_95U" in sumstats.columns:
        cols_to_check.append("OR_95U")
        if verbose: log.write(" -Checking if OR_95L>0 ...") 
        sumstats.loc[:,"OR_95U"] = pd.to_numeric(sumstats.loc[:,"OR_95U"], errors='coerce')
        sumstats = sumstats.loc[sumstats["OR_95U"]>0,:]
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad OR_95U.") 
    
    sumstats = sumstats.dropna(subset=cols_to_check)
    after_number=len(sumstats)
    if verbose: log.write(" -Removed "+str(oringinal_number - after_number)+" variants with bad statistics in total.") 
    return sumstats

###############################################################################################################
#20220426
def get_reverse_complementary_allele(a):
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a[::-1].translate(dic)

def flip_direction(string):
    flipped_string=""
    for char in string:
        if char=="?":flipped_string+="?"
        elif char=="+":flipped_string+="-"
        elif char=="-":flipped_string+="+"
    return flipped_string
    
def flipallelestats(sumstats,status="STATUS",verbose=True,log=gl.Log()):
    #
    
    
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))
    
    ###################get reverse complementary####################
    pattern = r"\w\w\w[67]\w"  
    matched_index = sumstats[status].str.match(pattern)
    if sum(matched_index)>0:
        if verbose: log.write("Start to convert alleles to reverse complement for SNPs with status xxx6x or xxx7x ... ") 
        if verbose: log.write(" -Flipping "+ str(sum(matched_index)) +" variants...") 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            if verbose: log.write(" -Converting to reverse complement : EA and NEA...") 
            reverse_complement_ea = sumstats.loc[matched_index,'EA'].apply(lambda x :get_reverse_complementary_allele(x))  
            reverse_complement_nea = sumstats.loc[matched_index,'NEA'].apply(lambda x :get_reverse_complementary_allele(x)) 
            sumstats.loc[matched_index,['NEA']] = reverse_complement_nea.astype("string")
            sumstats.loc[matched_index,['EA']] = reverse_complement_ea.astype("string")
    

    ###################flip ref####################
    pattern = r"\w\w\w[78]\w"  
    matched_index = sumstats[status].str.match(pattern)
    if sum(matched_index)>0:
        if verbose: log.write("Start to flip allele-specific stats for SNPs with status xxx8x: alt->ea , ref->nea ... ") 
        if verbose: log.write(" -Flipping "+ str(sum(matched_index)) +" variants...") 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            if verbose: log.write(" -Swapping column: NEA <=> EA...") 
            sumstats.loc[matched_index,['NEA','EA']] = sumstats.loc[matched_index,['EA','NEA']].values
        if "BETA" in sumstats.columns:
            if verbose: log.write(" -Flipping column: BETA = - BETA...") 
            sumstats.loc[matched_index,"BETA"] =     - sumstats.loc[matched_index,"BETA"].values
        if "EAF" in sumstats.columns:
            if verbose: log.write(" -Flipping column: EAF = 1 - EAF...") 
            sumstats.loc[matched_index,"EAF"] =   1 - sumstats.loc[matched_index,"EAF"].values
        if "OR" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR = 1 / OR...") 
            sumstats.loc[matched_index,"OR"] =   1 / sumstats.loc[matched_index,"freq"].values
        if "OR_95L" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR_95L = 1 / OR_95L...") 
            sumstats.loc[matched_index,"OR_95L"] =   1 / sumstats.loc[matched_index,"OR_95L"].values
        if "OR_95U" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR_95U = 1 / OR_95U...") 
            sumstats.loc[matched_index,"OR_95U"] =   1 / sumstats.loc[matched_index,"OR_95U"].values
        if "DIRECTION" in sumstats.columns:
            if verbose: log.write(" -Flipping column: DIRECTION +-? <=> -+? ...") 
            sumstats.loc[matched_index,"DIRECTION"] =   sumstats.loc[matched_index,"DIRECTION"].apply(flip_direction)
        #change status    
        if verbose: log.write(" -Changed the status for flipped variants: xxx8x -> xxx2x") 
        sumstats.loc[matched_index,status] = sumstats.loc[matched_index,status].apply(lambda x:change_status(x,4,"2")).astype("string")
         # flip ref
        # flip ref
    
    ###################flip ref for undistingushable indels####################
    pattern = r"\w\w[1][45][6]"  
    matched_index = sumstats[status].str.match(pattern)
    if sum(matched_index)>0:
        if verbose: log.write("Start to flip allele-specific stats for SNPs with status xxx8x: alt->ea , ref->nea ... ") 
        if verbose: log.write(" -Flipping "+ str(sum(matched_index)) +" variants...") 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            if verbose: log.write(" -Swapping column: NEA <=> EA...") 
            sumstats.loc[matched_index,['NEA','EA']] = sumstats.loc[matched_index,['EA','NEA']].values
        if "BETA" in sumstats.columns:
            if verbose: log.write(" -Flipping column: BETA = - BETA...") 
            sumstats.loc[matched_index,"BETA"] =     - sumstats.loc[matched_index,"BETA"].values
        if "EAF" in sumstats.columns:
            if verbose: log.write(" -Flipping column: EAF = 1 - EAF...") 
            sumstats.loc[matched_index,"EAF"] =   1 - sumstats.loc[matched_index,"EAF"].values
        if "OR" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR = 1 / OR...") 
            sumstats.loc[matched_index,"OR"] =   1 / sumstats.loc[matched_index,"freq"].values
        if "OR_95L" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR_95L = 1 / OR_95L...") 
            sumstats.loc[matched_index,"OR_95L"] =   1 / sumstats.loc[matched_index,"OR_95L"].values
        if "OR_95U" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR_95U = 1 / OR_95U...") 
            sumstats.loc[matched_index,"OR_95U"] =   1 / sumstats.loc[matched_index,"OR_95U"].values
        if "DIRECTION" in sumstats.columns:
            if verbose: log.write(" -Flipping column: DIRECTION +-? <=> -+? ...") 
            sumstats.loc[matched_index,"DIRECTION"] =   sumstats.loc[matched_index,"DIRECTION"].apply(flip_direction)
        #change status    
        if verbose: log.write(" -Changed the status for flipped variants: xxx8x -> xxx2x") 
        sumstats.loc[matched_index,status] = sumstats.loc[matched_index,status].apply(lambda x:change_status(x,5,"5")).astype("string")
         # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    pattern = r"\w\w\w[12][8]"  
    matched_index = sumstats[status].str.match(pattern)
    if sum(matched_index)>0:
        if verbose: log.write("Start to flip allele-specific stats for palindromic SNPs with status xxx[12]8: (-)strand <=> (+)strand ... ") 
        if verbose: log.write(" -Flipping "+ str(sum(matched_index)) +" variants...") 
        if "BETA" in sumstats.columns:
            if verbose: log.write(" -Flipping column: BETA = - BETA...") 
            sumstats.loc[matched_index,"BETA"] =     - sumstats.loc[matched_index,"BETA"].values
        if "EAF" in sumstats.columns:
            if verbose: log.write(" -Flipping column: EAF = 1 - EAF...") 
            sumstats.loc[matched_index,"EAF"] =   1 - sumstats.loc[matched_index,"EAF"].values
        if "OR" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR = 1 / OR...") 
            sumstats.loc[matched_index,"OR"] =   1 / sumstats.loc[matched_index,"freq"].values
        if "OR_95L" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR_95L = 1 / OR_95L...") 
            sumstats.loc[matched_index,"OR_95L"] =   1 / sumstats.loc[matched_index,"OR_95L"].values
        if "OR_95U" in sumstats.columns:
            if verbose: log.write(" -Flipping column: OR_95U = 1 / OR_95U...") 
            sumstats.loc[matched_index,"OR_95U"] =   1 / sumstats.loc[matched_index,"OR_95U"].values
        if "DIRECTION" in sumstats.columns:
            if verbose: log.write(" -Flipping column: DIRECTION +-? <=> -+? ...") 
            sumstats.loc[matched_index,"DIRECTION"] =   sumstats.loc[matched_index,"DIRECTION"].apply(flip_direction)
        #change status    
        if verbose: log.write(" -Changed the status for flipped variants:  xxx[12]8: ->  xxx[12]3") 
        sumstats.loc[matched_index,status] = sumstats.loc[matched_index,status].apply(lambda x:change_status(x,5,"3")).astype("string")


    return sumstats
###############################################################################################################

###############################################################################################################
#20220426
def liftover_snv(row,converter,to_build):
    # status 00xxx /19xxx /38xxx/99xxx
    status_pre=""
    status_end=row[2][2]+"00"
    chrom = row[0]  
    pos_0_based = int(row[1]) - 1
    results = converter[chrom][pos_0_based]
    if converter[chrom][pos_0_based]:
        # return chrom, pos_1_based
        return results[0][0].strip("chr"),results[0][1]+1,to_build+status_end
    else:
        return pd.NA,pd.NA,"99"+status_end

def liftover_variant(sumstats, 
             chrom="CHR", 
             pos="POS", 
             status="STATUS",
             from_build="19", 
             to_build="38",
             remove=True):
#sumstats.loc[:,pos] = np.floor(pd.to_numeric(sumstats.loc[:,pos], errors='coerce')).astype('Int64')
    converter = get_lifter("hg"+from_build,"hg"+to_build)
    not_na = (~sumstats[pos].isna()) & (~sumstats[chrom].isna())
    na=~not_na
    sumstats.loc[na,status] = sumstats.loc[na,status].apply(lambda x:change_status(x,1,to_build[0])).astype("string")
    sumstats.loc[na,status] = sumstats.loc[na,status].apply(lambda x:change_status(x,2,to_build[1])).astype("string")
    lifted = sumstats.loc[not_na,[chrom,pos,status]].apply(lambda x: liftover_snv(x[[chrom,pos,status]],converter,to_build),axis=1)
    sumstats.loc[not_na,chrom] = lifted.apply(lambda x :pd.NA if x[0] is pd.NA else str(x[0])).astype("string")
    sumstats.loc[not_na,pos] = np.floor(pd.to_numeric(lifted.apply(lambda x :pd.NA if x[0] is pd.NA else str(x[1])), errors='coerce')).astype('Int64')
    sumstats.loc[not_na,status] = lifted.apply(lambda x :str(x[2])).astype("string")
    return sumstats.loc[:,:]

def parallelizeliftovervariant(sumstats,n_cores=1,chrom="CHR", pos="POS", from_build="19", to_build="38",status="STATUS",remove=True, verbose=True,log=gl.Log()):
    if verbose: log.write("Start to perform liftover...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose: log.write(" -CPU Cores to use :",n_cores)
    if verbose: log.write(" -Performing liftover ...")
    if verbose: log.write(" -Creating converter : hg" + from_build +" to hg"+ to_build)
    if verbose: log.write(" -Converting variants : "+str(len(sumstats)))
    func=liftover_variant
    df_split = np.array_split(sumstats, n_cores)
    pool = Pool(n_cores)
    #df = pd.concat(pool.starmap(func, df_split))
    sumstats = pd.concat(pool.map(partial(func,chrom=chrom,pos=pos,from_build=from_build,to_build=to_build,remove=remove,status=status),df_split))
    pool.close()
    pool.join()
    
    map_num   = len(sumstats.loc[~sumstats["POS"].isna(),:])
    unmap_num = len(sumstats.loc[sumstats["POS"].isna(),:])
    if verbose:log.write(" -Converted variants: "+str(map_num))
    if remove:
        if verbose: log.write(" -Remove unmapped variants: "+str(unmap_num))
        sumstats = sumstats.loc[~sumstats["POS"].isna(),:]
    sumstats = fixchr(sumstats,chrom=chrom,add_prefix="",remove=remove, verbose=verbose,log=log)
    if verbose: log.write(" -Liftover is performed successfully!")
    return sumstats

###############################################################################################################
#20220426
def sortcoordinate(sumstats,chrom="CHR",pos="POS",reindex=True,verbose=True,log=gl.Log()):
    
    chromosome_number_to_string = {i:str(i) for i in range(1,23)}
    chromosome_number_to_string[23] = "X"
    chromosome_number_to_string[24] = "Y"
    chromosome_number_to_string[25] = "MT"
    
    chromosome_string_to_number = {str(i):i for i in range(1,23)}
    chromosome_string_to_number["X"] = 23
    chromosome_string_to_number["Y"] = 24
    chromosome_string_to_number["MT"] = 25
    
    for i in sumstats[chrom]:
        if i not in chromosome_string_to_number.keys():
            chromosome_string_to_number[i]=i
            chromosome_number_to_string[i]=i

    good_chrpos = sumstats[chrom].isin(gl.get_chr_list())&sumstats[pos].notnull()
    sumstats_good = sumstats.loc[good_chrpos,:]
    sumstats_bad = sumstats.loc[~good_chrpos,:]
    
    sumstats = sumstats.sort_values(by=[chrom,pos],ascending=True,ignore_index=True)
                                    #####################????????????)
    
    
    if verbose: log.write("Start to sort the genome coordinates...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    
    if verbose: log.write(" -Mapping chr1-22,X,Y,MT to numbers 1-25, other chromosome notations will be sorted as strings...")
    
    sumstats[chrom] = sumstats[chrom].map(chromosome_string_to_number)
    if verbose: log.write(" -Force converting POS to integers...")
    sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    if verbose: log.write(" -Sorting genome coordinates...")
    
    sumstats = sumstats.sort_values(by=[chrom,pos],ascending=True,ignore_index=True)
    
    if verbose: log.write(" -Converting chromosome column back to string type...")
    
    sumstats[chrom] = sumstats[chrom].map(chromosome_number_to_string)
    sumstats[chrom] = sumstats[chrom].astype("string")
    
    return sumstats
###############################################################################################################
#20220426
def sortcolumn(sumstats,verbose=True,log=gl.Log()):
    if verbose: log.write("Start to reorder the columns...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    order = [
        "SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "Z",
        "CHISQ", "P", "MLOG10P", "OR", "OR_SE", "OR_95L", "OR_95U", "INFO", "N","DIRECTION","STATUS"
           ]
    
    output_columns = []
    for i in order:
        if i in sumstats.columns: output_columns.append(i)
    for i in sumstats.columns:
        if i not in order: output_columns.append(i)
    if verbose: log.write(" -Reordering columns to    :", ",".join(output_columns))
    sumstats = sumstats.loc[:, output_columns]
    return sumstats

################################################################################################################
def change_status(status_value,digit,to):
    if digit>1:
        status_pre=status_value[:digit-1]
    else:
        status_pre=""
        
    status_end=status_value[digit:]
    return status_pre+str(to)+status_end

################################################################################################################
