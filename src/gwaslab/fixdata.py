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
#20220514 

def fixID(sumstats,
       snpid="SNPID",rsid="rsID",chrom="CHR",pos="POS",nea="NEA",ea="EA",status="STATUS",
       fixchrpos=False,fixid=False,fixeanea=False,fixeanea_flip=False,overwrite=False,verbose=True,log=gl.Log()):  
    
    '''
    1. fx SNPid
    2. fix chr and pos using snpid
    3. checking rsid and chr:pos:nea:ea
    '''
    if verbose: log.write("Start to check IDs...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    
    check_col(sumstats,[snpid,rsid],status)
    
    ############################  checking ###################################################  
    if snpid in sumstats.columns:  
        if verbose: log.write(" -Checking if SNPID is chr:pos:ref:alt...(separator: - ,: , _)")
        is_chrposrefalt = sumstats[snpid].str.match(r'(chr)?([0-9XYMT]+)[:_-]([0-9]+)[:_-]([ATCG]+)[:_-]([ATCG]+)', case=False, flags=0, na=False)
        is_snpid_na = sumstats[snpid].isna()
        sumstats.loc[ is_chrposrefalt,status] = vchange_status(sumstats.loc[ is_chrposrefalt,status],3 ,"975" ,"630")
        sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status] = vchange_status(sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status],3 ,"975" ,"842")
        
    if rsid in sumstats.columns: 
        if verbose: log.write(" -Checking if rsID is rsxxxxxx or RSxxxxxxx...")
        is_rsid = sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)
        
        sumstats.loc[ is_rsid,status] = vchange_status(sumstats.loc[ is_rsid,status], 3, "986","520")
        sumstats.loc[~is_rsid,status] = vchange_status(sumstats.loc[~is_rsid,status], 3, "986","743")
        
        if verbose: log.write(" -Checking if chr:pos:ref:alt is mixed in rsID column ...")
        is_rs_chrpos = sumstats[rsid].str.match(r'(chr)?([0-9XYMT]+)[:_-]([0-9]+)[:_-]([ATCG]+)[:_-]([ATCG]+)', case=False, flags=0, na=False)
        
        if verbose: log.write(" -Number of chr:pos:ref:alt mixed in rsID column :",sum(is_rs_chrpos))
        if verbose: log.write(" -Number of Unrecognized rsID :",len(sumstats) - sum(is_rs_chrpos) - sum(is_rsid) ) 
        if verbose: log.write(" -A look at the unrecognized rsID :",set(sumstats.loc[(~is_rsid)&(~is_rs_chrpos),rsid].head()),"...") 
      
    ############################  fixing chr pos###################################################  
    if fixchrpos is True:
    # from snpid or rsid, extract chr:pos to fix CHR and POS    
        if snpid in sumstats.columns: 
            if verbose: log.write(" -Fixing CHR and POS...")
            if overwrite is True: 
                if verbose: log.write(" -Overwrite is applied...")
                # fix all
                to_fix = is_chrposrefalt
                
                #fix variants with chr and pos being empty
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants. ...")
            
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                if verbose: log.write(" -Initiating CHR columns...")
                sumstats.loc[:,chrom]=pd.Series(dtype="string")   
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants. ...")
            
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                if verbose: log.write(" -Initiating CHR and POS column...")
                sumstats.loc[:,pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna() 
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
                # format and qc filled chr and pos
                sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[0].strip("chrCHR")).astype("string")
                sumstats.loc[to_fix,pos] =np.floor(pd.to_numeric(sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]), errors='coerce')).astype('Int64')
                # change status
                sumstats.loc[to_fix,status] = vchange_status(sumstats.loc[to_fix,status], 4, "98765432","00000000")    
        
        if rsid in sumstats.columns:
            if verbose: log.write(" -Fixing CHR and POS using chr:pos:ref:alt format variants in rsID column...")
            if overwrite is True: 
                if verbose: log.write(" -Overwrite is applied...")
                to_fix = is_rs_chrpos
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants ...")
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                if verbose: log.write(" -Initiating CHR columns...")
                sumstats.loc[:,chrom]=pd.Series(dtype="string")   
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                elif verbose: log.write(" -No fixable vairants ...")
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                if verbose: log.write(" -Initiating CHR and POS column...")
                sumstats.loc[:,pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna() 
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
                sumstats.loc[to_fix,status] = vchange_status(sumstats.loc[to_fix,status], 4, "98765432","00000000") 
                
    ############################  fixing chr pos###################################################   
    #if fixeanea is True:
    #    if verbose: log.write(" -Warning: Please make sure a1 is ref or not in Chr:pos:a1:a2")
    #    if overwrite is True:
    #        if verbose: log.write(" -Overwrite is applied...")
    #        to_fix = is_chrposrefalt
    #    elif (nea in sumstats.columns) and (nea in sumstats.columns):
    #        to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
    #        if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
    #    elif (nea in sumstats.columns) and (ea not in sumstats.columns):
    #        if verbose: log.write(" -Initiating EA columns...")
    #        sumstats[ea]=pd.Series(dtype="string")
    #        to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
    #        if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
    #    elif (nea not in sumstats.columns) and (ea in sumstats.columns):
    #        if verbose: log.write(" -Initiating NEA columns...")
    #        sumstats[nea]=pd.Series(dtype="string")
    #        to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
    #        if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
    #    else:
    #        if verbose: log.write(" -Initiating EA and NEA columns...")
    #        sumstats[nea]=pd.Series(dtype="string")
    #        sumstats[ea]=pd.Series(dtype="string")
    #        to_fix = is_chrposrefalt
    #        if sum(to_fix)>0: 
    #            if verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
    #                
    #    if sum(to_fix)>0:    
    #        if verbose: log.write(" -Filling "+str(sum(to_fix))+" EA and NEA columns using SNPID's chr:pos:nea:ea...")
    #        
    #        if fixeanea_flip is True:
    #            if verbose: log.write(" -Flipped : chr:pos:a1:a2...a1->EA , a2->NEA ")
    #            sumstats.loc[to_fix,ea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[2]).astype("string")  
    #            sumstats.loc[to_fix,nea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[3]).astype("string")
    #        else:
    #            if verbose: log.write(" -Chr:pos:a1:a2...a1->EA , a2->NEA ")
    #            sumstats.loc[to_fix,ea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[3]).astype("string")  
    #            sumstats.loc[to_fix,nea] = sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[2]).astype("string")
    #        #to_change_status = sumstats[status].str.match(r"\w\w\w[45]\w\w\w")
    #        #sumstats.loc[to_fix&to_change_status,status] = vchange_status(sumstats.loc[to_fix&to_change_status,status],4,"2")  
    #        #sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]).astype("string")
    #        #sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[1]).astype("Int64")
    
    ############################  fixing id ################################################### 
    if fixid is True:
        if snpid not in sumstats.columns: 
        # initiate a SNPID column
            sumstats.loc[:,snpid]=pd.Series(dtype="string")
            
        if (chrom in sumstats.columns) and (pos in sumstats.columns):
            #only fix when CHR and POS is available
            pre_number=sum(sumstats[snpid].isna())                                
            
            if overwrite is False:
                #fix empty 
                to_fix = sumstats[snpid].isna() & sumstats[status].str.match( r"\w\w\w[0]\w\w\w", case=False, flags=0, na=False ) 
            else:
                #fix all
                to_fix = sumstats[status].str.match( r"\w\w\w[0]\w\w\w", case=False, flags=0, na=False ) 
            
            if (ea in sumstats.columns) and (nea in sumstats.columns):
            # when ea and nea is available  -> check status -> fix to chr:pos:nea:ea 
                
                pattern = r"\w\w\w[0]\w\w\w"  
                matched_index = sumstats[status].str.match(pattern) 
                to_part_fix = matched_index & to_fix 
                
                pattern = r"\w\w\w[0][01267][01234]\w"  
                matched_index = sumstats[status].str.match(pattern) 
                to_full_fix = matched_index & to_fix 
                
                if sum(to_full_fix)>0:
                    sumstats.loc[to_full_fix,snpid] = sumstats.loc[to_full_fix,chrom].astype("string") + ":"+sumstats.loc[to_full_fix,pos].astype("string") +":"+ sumstats.loc[to_full_fix,nea].astype("string") +":"+ sumstats.loc[to_full_fix,ea].astype("string")
                if verbose: log.write(" -Filling "+str(sum(to_part_fix)) +" SNPID using CHR POS...")
                if sum(to_part_fix)>0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")    
                sumstats.loc[(to_full_fix),status] = vchange_status(sumstats.loc[(to_full_fix),status],3,"975","630")    
                sumstats.loc[(to_part_fix),status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")
            else:
            #when these is no ea or ena, just fix to chr:pos
                to_part_fix = to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull()
                if verbose: log.write(" -Filling "+str(sum(to_part_fix)) +" SNPID using CHR POS...")
                if sum(to_part_fix)>0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                    sumstats.loc[to_part_fix,status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")
                    
            after_number=sum(sumstats[snpid].isna())
            if verbose: log.write(" -Fixed "+ str(pre_number - after_number) +" variants ID...")
        elif verbose: log.write(" -ID unfixable: no CHR and POS columns or no SNPID. ")
    return sumstats

###############################################################################################################

    
###############################################################################################################
#20220514 
def removedup(sumstats,mode="dm",chrom="CHR",pos="POS",snpid="SNPID",ea="EA",nea="NEA",rsid="rsID",keep='first',verbose=True,log=gl.Log()):
    '''
    remove duplicate SNPs based on  1. SNPID, EA, and NEA
    remove duplicate SNPs based on  2. rsID for non-NA variants
    remove multiallelic SNPs based on  3. chr pos 
    '''
    total_number = len(sumstats)  
    if (snpid in sumstats.columns) and (nea in sumstats.columns) and (ea in sumstats.columns) and "d" in mode:
        if verbose: log.write("Start to remove duplicated variants based on snpid,ea and nea...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns)) 
        if verbose: log.write(" -Which variant to keep: ",  keep )   
        pre_number =len(sumstats)   
        if snpid in sumstats.columns:
            sumstats = sumstats.loc[sumstats[snpid].isna() | (~sumstats.duplicated(subset=[snpid,ea,nea], keep=keep)),:]
            after_number=len(sumstats)   
            if verbose:  log.write(" -Removed ",pre_number -after_number ," based on SNPID, EA, and NEA...")
     
    
    if (rsid in sumstats.columns) and ("d" in mode):
        pre_number =len(sumstats)
        if verbose: log.write("Start to remove duplicated variants based on rsID...")
        sumstats = sumstats.loc[sumstats[rsid].isna() | (~sumstats.duplicated(subset=rsid, keep=keep)),:]
        after_number=len(sumstats)   
        if verbose:  log.write(" -Removed ",pre_number -after_number ," based on rsID...")

      
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and "m" in mode:
        pre_number =len(sumstats) 
        if verbose: log.write("Start to remove multiallelic variants based on chr:pos...")    
        if verbose: log.write(" -Which variant to keep: ",  keep ) 
        sumstats = sumstats.loc[(~sumstats.loc[:,[snpid,chrom,pos]].all(axis=1)) | (~sumstats.duplicated(subset=[chrom,pos], keep=keep)),:]
        after_number=len(sumstats)  
        if verbose:  log.write(" -Removed ",total_number -after_number," multiallelic variants in total.")   
    after_number=len(sumstats)   
    if verbose:  log.write(" -Removed ",total_number -after_number," duplicates in total.")
    return sumstats

###############################################################################################################
#20220514
def fixchr(sumstats,chrom="CHR",status="STATUS",add_prefix="",x=23,y=24,mt=25,remove=False, verbose=True,log=gl.Log()):
        
        if verbose: log.write("Start to fix chromosome notation...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
        
        check_col(sumstats,chrom,status)
        
        is_chr_fixable = sumstats[chrom].str.match(r'(chr)?([012][0-9]|[0-9]|X|Y|MT)', case=False, flags=0, na=False)
        if verbose: log.write(" -Vairants with fixable chromosome notations:",sum(is_chr_fixable))   
        
        is_chr_na    = sumstats[chrom].isna()
        if sum(is_chr_na)>0 and verbose: 
            log.write(" -Vairants with NA chromosome notations:",sum(is_chr_na))  
        
        is_chr_invalid = (~is_chr_fixable)&(~is_chr_na)
        if sum(is_chr_invalid)>0 and verbose: 
            log.write(" -Vairants with invalid chromosome notations:",sum(is_chr_invalid))  
        if sum(~is_chr_fixable)>0 and verbose: 
            log.write(" -Unrecognized chromosome notations :",set(sumstats.loc[~is_chr_fixable,chrom].head()),"...") 

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
        
        # 23,24,25 to X, Y, MT
        convert_num_to_xymt={str(x):"X",str(y):"Y",str(mt):"MT"}
        sex_chr = sumstats[chrom].isin([str(x),str(y),str(mt)])
        sumstats.loc[sex_chr,chrom] =sumstats.loc[sex_chr,chrom].map(convert_num_to_xymt)
        if verbose: log.write(" -Standardizing sex chromosome notations...") 
        
        chrom_list = gl.get_chr_list() #bottom 
        unrecognized_num = sum(~sumstats[chrom].isin(chrom_list))
        
        sumstats.loc[is_chr_fixable,status] = vchange_status(sumstats.loc[is_chr_fixable,status],4,"986","520")
        sumstats.loc[is_chr_invalid,status] = vchange_status(sumstats.loc[is_chr_invalid,status],4,"986","743")
        
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
#20220514
def fixpos(sumstats,pos="POS",status="STATUS",remove=False, verbose=True,log=gl.Log()):
        if verbose: log.write("Start to fix basepair positions...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
        
        check_col(sumstats,pos,status)
        
        all_var_num = len(sumstats)
        #convert to numeric
        is_pos_na = sumstats.loc[:,pos].isna()
        
        sumstats.loc[:,pos] = np.floor(pd.to_numeric(sumstats.loc[:,pos], errors='coerce')).astype('Int64')
        
        is_pos_fixed = ~sumstats.loc[:,pos].isna()
        is_pos_invalid = (~is_pos_na)&(~is_pos_fixed)
        
        sumstats.loc[is_pos_fixed,status] = vchange_status(sumstats.loc[is_pos_fixed,status],4,"975","630")
        sumstats.loc[is_pos_invalid,status] = vchange_status(sumstats.loc[is_pos_invalid,status],4,"975","842")
        
        # remove 300,000,000
        
        #remove na
        if remove is True: 
            sumstats = sumstats.loc[~sumstats[pos].isna(),:]
            remain_var_num = len(sumstats)
            if verbose: log.write(" -Removed "+str(all_var_num - remain_var_num)+" variants with bad positions.")        
        
        if verbose: log.write(" -Converted all position to datatype Int64.")
        if verbose: log.write(" -Fixed basepair position successfully.")
        return sumstats
    
###############################################################################################################    
#20220514
def fixallele(sumstats,ea="EA", nea="NEA",status="STATUS",remove=False,verbose=True,log=gl.Log()):
        # remove variants with alleles other than actgACTG
        if verbose: log.write("Start to fix alleles...")
        if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
        
        check_col(sumstats,ea,nea,status)
        
        #if (ea not in sumstats.columns) or (nea not in sumstats.columns):

        all_var_num = len(sumstats)
        
        good_ea  = ~sumstats[ea].str.contains("[^actgACTG]")
        good_nea = ~sumstats[nea].str.contains("[^actgACTG]")
        is_eanea_na = sumstats[ea].isna() |  sumstats[nea].isna()
        is_invalid = sumstats[ea].str.contains("[^actgACTG]") | sumstats[nea].str.contains("[^actgACTG]") | (sumstats[nea] == sumstats[ea])
        exclude  = sum(sumstats[ea].str.contains("[^actgACTG]") | sumstats[nea].str.contains("[^actgACTG]"))
        if remove is True:
            sumstats = sumstats.loc[good_ea & good_nea,:].copy()
            if verbose: log.write(" -Removed "+str(exclude)+" variants with alleles that contain bases other than A/C/T/G .")  
                
        sumstats.loc[:,ea] = sumstats[ea].str.upper()
        sumstats.loc[:,nea]= sumstats[nea].str.upper()
        if verbose: log.write(" -Converted all bases to string datatype and UPPERCASE.")
        
        is_eanea_fixed = good_ea | good_nea
        is_snp = (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()==1)
        is_indel = (sumstats[ea].str.len()!=sumstats[nea].str.len())
        is_not_normalized = (sumstats[ea].str.len()>1) &(sumstats[nea].str.len()>1)
        is_normalized = is_indel &( (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()>1) | (sumstats[ea].str.len()>1) &(sumstats[nea].str.len()==1) )
        
        sumstats.loc[is_invalid, status]                 = vchange_status(sumstats.loc[is_invalid,status],                5,"9","6") 
        sumstats.loc[is_eanea_na,status]                 = vchange_status(sumstats.loc[is_eanea_na, status],              5,"9","7")
        sumstats.loc[is_eanea_fixed&is_not_normalized,status]   = vchange_status(sumstats.loc[is_eanea_fixed&is_not_normalized,status], 5,"9","5")
        sumstats.loc[is_eanea_fixed&is_snp, status]          = vchange_status(sumstats.loc[is_eanea_fixed&is_snp,status],        5,"9","0")
        sumstats.loc[is_eanea_fixed&is_indel,status]         = vchange_status(sumstats.loc[is_eanea_fixed&is_indel, status],      5,"9","4")   
        sumstats.loc[is_eanea_fixed&is_normalized,status]      = vchange_status(sumstats.loc[is_eanea_fixed&is_normalized, status],  5,"4","3")
        remain_var_num = len(sumstats)
              
        
        if verbose: log.write(" -Fixed allele successfully.")
        return sumstats

###############################################################################################################   
#20220514

def parallelnormalizeallele(sumstats,pos="POS",ea="EA" ,nea="NEA",status="STATUS",n_cores=1,verbose=True,log=gl.Log()):
    
    check_col(sumstats,pos,ea,nea,status)
    
    if verbose: log.write("Start to normalize variants...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    variants_to_check = sumstats[status].str.match(r'\w\w\w\w[45]\w\w', case=False, flags=0, na=False)
    if sum(variants_to_check)==0:
        log.write(" -No available variants to normalize..")
        return sumstats

    
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
    to_normalize = r"\w\w\w\w[45]\w\w"
    variants_to_check = sumstats[status].str.match(to_normalize)
    if sum(variants_to_check)>0:
        normalized = sumstats.loc[variants_to_check,[pos,nea,ea,status]].apply(lambda x: normalizevariant(x[0],x[1],x[2],x[3]),axis=1)
        sumstats.loc[variants_to_check,[pos,nea,ea,status]] = pd.DataFrame(normalized.to_list(), columns=[pos,nea,ea,status],index=sumstats.loc[variants_to_check,:].index)
    return sumstats

def normalizevariant(pos,a,b,status):
    # https://genome.sph.umich.edu/wiki/Variant_Normalization
    # a - ref - nea    starting -> pos
    # b - alt - ea
    # status xx1xx /xx2xx /xx9xx
    
    status_pre=status[:3]
    status_end=status[4:]
    
    if len(a)==1 or len(b)==1:
        return pos,a,b,status_pre+"3"+status_end
    
    pos_change=0
    pointer_a_l, pointer_a_r = 0, len(a)-1
    pointer_b_l, pointer_b_r = 0, len(b)-1

    #remove from right
    for i in range(max(len(a),len(b))):
        if a[pointer_a_r] == b[pointer_b_r]:
            pointer_a_r-=1
            pointer_b_r-=1
        else:
            break
        if pointer_a_r-pointer_a_l==0 or pointer_b_r-pointer_b_l==0:
            return pos,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"3"+status_end

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
    return pos+pos_change,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"3"+status_end
###############################################################################################################

###############################################################################################################
#20220426
def sanitycheckstats(sumstats,coltocheck=["P","MLOG10P","BETA","SE","EAF","CHISQ","N","OR","OR_95L","OR_95U"],verbose=True,log=gl.Log()):
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
    if "MLOG10P" in coltocheck and "MLOG10P" in sumstats.columns:
        cols_to_check.append("MLOG10P")
        if verbose: log.write(" -Checking if MLOG10P>=0 ...") 
        sumstats.loc[:,"MLOG10P"] = pd.to_numeric(sumstats.loc[:,"MLOG10P"], errors='coerce')
        sumstats = sumstats.loc[(sumstats["MLOG10P"]>=0),:]
        after_number=len(sumstats)
        if verbose: log.write(" -Removed "+str(pre_number - after_number)+" variants with bad MLOG10P.") 
    
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
        if char=="?":
            flipped_string+="?"
        elif char=="+":
            flipped_string+="-"
        elif char=="-":
            flipped_string+="+"
        else: #sometime it is 0
            flipped_string+=char
    return flipped_string
    
def flipallelestats(sumstats,status="STATUS",verbose=True,log=gl.Log()):
    #
    
    
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))
    
    ###################get reverse complementary####################
    pattern = r"\w\w\w\w\w[45]\w"  
    matched_index = sumstats[status].str.match(pattern)
    if sum(matched_index)>0:
        if verbose: log.write("Start to convert alleles to reverse complement ... ") 
        if verbose: log.write(" -Flipping "+ str(sum(matched_index)) +" variants...") 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            if verbose: log.write(" -Converting to reverse complement : EA and NEA...") 
            reverse_complement_ea = sumstats.loc[matched_index,'EA'].apply(lambda x :get_reverse_complementary_allele(x))  
            reverse_complement_nea = sumstats.loc[matched_index,'NEA'].apply(lambda x :get_reverse_complementary_allele(x)) 
            sumstats.loc[matched_index,['NEA']] = reverse_complement_nea.astype("string")
            sumstats.loc[matched_index,['EA']] = reverse_complement_ea.astype("string")
            sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 6, "4","2")

    ###################flip ref####################
    pattern = r"\w\w\w\w\w[35]\w"  
    matched_index = sumstats[status].str.match(pattern)
    if sum(matched_index)>0:
        if verbose: log.write("Start to flip allele-specific stats for SNPs with status: alt->ea , ref->nea ... ") 
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
        if verbose: log.write(" -Changed the status for flipped variants") 
        sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 6, "35","12")
        
    ###################flip ref for undistingushable indels####################
    pattern = r"\w\w\w\w[123][67][6]"  
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
        if verbose: log.write(" -Changed the status for flipped variants") 
        sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 7, "6","4")
         # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    pattern = r"\w\w\w[012][5]"  
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
        sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 7, "5","2")

    return sumstats
###############################################################################################################

###############################################################################################################
#20220426
def liftover_snv(row,converter,to_build):
    status_pre=""
    status_end=row[2][2]+"9"+row[2][4]+"99"
    
    chrom = row[0]  
    pos_0_based = int(row[1]) - 1
    results = converter[chrom][pos_0_based]
    if converter[chrom][pos_0_based]:
        # return chrom, pos_1_based
        return results[0][0].strip("chr"),results[0][1]+1,to_build+status_end
    else:
        return pd.NA,pd.NA,"97"+status_end

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
    na = sumstats[pos].isna() | sumstats[chrom].isna()
    sumstats.loc[na,status] = vchange_status(sumstats.loc[na, status],  1,"13","99")
    sumstats.loc[na,status] = vchange_status(sumstats.loc[na, status],  2,"39","88")  
    
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
    
    check_col(sumstats,chrom,pos)
    
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

################################################################################################################
def vchange_status(status,digit,before,after):
    for i in range(len(before)):
        pattern= (digit-1) * r'\w' + before[i] + (7 - digit)* r'\w'     
        to_change = status.str.match(pattern, case=False, flags=0, na=False)  
        if sum(to_change)>0:
            if digit>1:
                status_pre = status[to_change].str[:digit-1]
            else:
                status_pre = ""

            status_end=status[to_change].str[digit:]

            status[to_change] = status_pre+after[i]+status_end
    return status

def check_col(df,*args):
    not_in_df=[]
    for i in args:
        if type(i) is str:
            if i in df.columns:
                continue
            else:
                not_in_df.append(i)
        else:
            count=0
            for j in i:
                if j in df.columns:
                    count+=1
            if count==0:
                raise ValueError("Specified columns names was not detected. Please check:"+",".join(i))
    
    if len(not_in_df)>0:
        raise ValueError("Specified columns names was not detected. Please check:"+",".join(not_in_df))
    
            
    
    