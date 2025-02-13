import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from liftover import get_lifter
from liftover import ChainFile
from functools import partial
from gwaslab.g_vchange_status import vchange_status
from gwaslab.g_vchange_status import status_match
from gwaslab.g_vchange_status import change_status
from gwaslab.g_Log import Log
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_chr_list
from gwaslab.qc_check_datatype import check_datatype
from gwaslab.qc_check_datatype import check_dataframe_shape
from gwaslab.g_version import _get_version
from gwaslab.util_in_fill_data import _convert_betase_to_mlog10p
from gwaslab.util_in_fill_data import _convert_betase_to_p
from gwaslab.util_in_fill_data import _convert_mlog10p_to_p
from gwaslab.bd_common_data import get_chain
#process build
#setbuild
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
#_check_data_consistency
#flipallelestats
#parallelizeassignrsid
#sortcoordinate
#sortcolumn

###############################################################################################################
# 20220514 
def _process_build(build,log,verbose):
    if str(build).lower() in ["hg19","19","37","b37","grch37"]:
        log.write(" -Genomic coordinates are based on GRCh37/hg19...", verbose=verbose)
        final_build = "19"
    elif str(build).lower() in ["hg18","18","36","b36","grch36"]:
        log.write(" -Genomic coordinates are based on GRCh36/hg18...", verbose=verbose)
        final_build = "18"
    elif str(build).lower() in ["hg38","38","b38","grch38"]:
        log.write(" -Genomic coordinates are based on GRCh38/hg38...", verbose=verbose)
        final_build = "38"
    elif str(build).lower() in ["t2t","hs1","chm13","13"]:
        log.write(" -Genomic coordinates are based on T2T-CHM13...", verbose=verbose)
        final_build = "13"
    else:
        log.warning("Version of genomic coordinates is unknown...", verbose=verbose)
        final_build = "99"
    return final_build

def _set_build(sumstats, build="99", status="STATUS",verbose=True,log=Log()):
    build = _process_build(build,log=log,verbose=verbose)
    sumstats[status] = vchange_status(sumstats[status], 1, "139",build[0]*3)
    sumstats[status] = vchange_status(sumstats[status], 2, "89",build[1]*3)
    return sumstats, build

def fixID(sumstats,
       snpid="SNPID",rsid="rsID",chrom="CHR",pos="POS",nea="NEA",ea="EA",status="STATUS",fixprefix=False,
       fixchrpos=False,fixid=False,fixeanea=False,fixeanea_flip=False,fixsep=False,
       overwrite=False,verbose=True,forcefixid=False,log=Log()):  
    '''
    1. fx SNPid
    2. fix chr and pos using snpid
    3. checking rsid and chr:pos:nea:ea
    '''
    ##start function with col checking##########################################################
    _start_line = "check SNPID/rsID"
    _end_line = "checking SNPID/rsID"
    _start_cols =[]
    _start_function = ".fix_id()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    ############################  checking datatype ###################################################  
    if rsid in sumstats.columns:
        # convert to string datatype
        try:
            log.write(" -Checking rsID data type...",verbose=verbose)
            if sumstats[rsid].dtype == "string":
                pass
            else:
                log.write(" -Converting rsID to pd.string data type...",verbose=verbose)
                sumstats[rsid] = sumstats[rsid].astype("string")
        except:
            log.write(" -Force converting rsID to pd.string data type...",verbose=verbose)
            sumstats[rsid] = sumstats[rsid].astype("string")
    if snpid in sumstats.columns: 
        # convert to string datatype
        try:
            log.write(" -Checking SNPID data type...",verbose=verbose)
            if sumstats[snpid].dtype == "string":
                pass
            else:
                log.write(" -Converting SNPID to pd.string data type...",verbose=verbose)
                sumstats[snpid] = sumstats[snpid].astype("string")
        except:
            log.write(" -Force converting SNPID to pd.string data type...",verbose=verbose)
            sumstats[snpid] = sumstats[snpid].astype("string")

    ############################  checking ###################################################  
    if snpid in sumstats.columns:  
        log.write(" -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)",verbose=verbose)
        # check if SNPID is CHR:POS:EA:NEA
        is_chrposrefalt = sumstats[snpid].str.match(r'^\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+$', case=False, flags=0, na=False)
        # check if SNPID is NA
        is_snpid_na = sumstats[snpid].isna()

        # change STATUS code
        sumstats.loc[ is_chrposrefalt,status] = vchange_status(sumstats.loc[ is_chrposrefalt,status],3 ,"975" ,"630")
        sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status] = vchange_status(sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status],3 ,"975" ,"842")
        
    if rsid in sumstats.columns: 
        log.write(" -Checking if rsID is rsxxxxxx...", verbose=verbose)
        is_rsid = sumstats[rsid].str.match(r'^rs\d+$', case=False, flags=0, na=False)
        
        sumstats.loc[ is_rsid,status] = vchange_status(sumstats.loc[ is_rsid,status], 3, "986","520")
        sumstats.loc[~is_rsid,status] = vchange_status(sumstats.loc[~is_rsid,status], 3, "986","743")
        
        log.write(" -Checking if CHR:POS:NEA:EA is mixed in rsID column ...", verbose=verbose)
        is_rs_chrpos = sumstats[rsid].str.match(r'^\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+$', case=False, flags=0, na=False)
        
        log.write(" -Number of CHR:POS:NEA:EA mixed in rsID column :",sum(is_rs_chrpos), verbose=verbose)
        log.write(" -Number of Unrecognized rsID :",len(sumstats) - sum(is_rs_chrpos) - sum(is_rsid) , verbose=verbose) 
        log.write(" -A look at the unrecognized rsID :",set(sumstats.loc[(~is_rsid)&(~is_rs_chrpos),rsid].head()),"...", verbose=verbose) 
      
    ############################  fixing chr pos###################################################  
    if fixchrpos == True:
    # from snpid or rsid, extract CHR:POS to fix CHR and POS    
        if snpid in sumstats.columns: 
            log.write(" -Fixing CHR and POS...", verbose=verbose)
            if overwrite is True: 
                log.write(" -Overwrite is applied...", verbose=verbose)
                # fix all
                to_fix = is_chrposrefalt
                
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                #fix variants with chr and pos being NA
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                to_fix_num = sum(to_fix)
                if to_fix_num and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)
            
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                log.write(" -Initiating CHR columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                to_fix_num = sum(to_fix)
                if to_fix_num>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)
            
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                log.write(" -Initiating CHR and POS column...", verbose=verbose)
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna() 
                to_fix_num = sum(to_fix)
                if to_fix_num>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)     
                
            else:
                log.write(" -Initiating CHR and POS columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt
                to_fix_num = sum(to_fix)
                if to_fix_num>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)   
                    
            if sum(to_fix)>0:
                log.write(" -Filling CHR and POS columns using valid SNPID's chr:pos...", verbose=verbose)
                # format and qc filled chr and pos
               
                sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[1]
                sumstats.loc[to_fix,pos] = sumstats.loc[to_fix,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[2]
                
                #sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,snpid].str.split(':|_|-').str[0].str.strip("chrCHR").astype("string")
                #sumstats.loc[to_fix,pos] =np.floor(pd.to_numeric(sumstats.loc[to_fix,snpid].str.split(':|_|-').str[1], errors='coerce')).astype('Int64')
                # change status
                #sumstats.loc[to_fix,status] = vchange_status(sumstats.loc[to_fix,status], 4, "98765432","00000000")    
        
        if rsid in sumstats.columns:
            log.write(" -Fixing CHR and POS using chr:pos:ref:alt format variants in rsID column...", verbose=verbose)
            if overwrite is True: 
                log.write(" -Overwrite is applied...", verbose=verbose)
                to_fix = is_rs_chrpos
            elif (chrom in sumstats.columns) and (pos in sumstats.columns) :
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                log.write(" -Initiating CHR columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                log.write(" -Initiating CHR and POS column...", verbose=verbose)
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna() 
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)
            else:
                log.write(" -Initiating CHR and POS columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos
                if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)   
            
            if sum(to_fix)>0:    
                log.write(" -Filling CHR and POS columns using chr:pos:ref:alt format variants in rsID column...", verbose=verbose)
                sumstats.loc[to_fix,chrom] = sumstats.loc[to_fix,rsid].str.split(':|_|-',n=2).str[0]
                sumstats.loc[to_fix,pos] = sumstats.loc[to_fix,rsid].str.split(':|_|-',n=2).str[1]
                #sumstats.loc[to_fix,pos] = np.floor(pd.to_numeric(sumstats.loc[to_fix,rsid].str.split(':|_|-',x).get(1), errors='coerce')).astype('Int64')
                #sumstats.loc[to_fix,status] = vchange_status(sumstats.loc[to_fix,status], 4, "98765432","00000000").astype("string")  
                
    ############################  fixing chr pos###################################################   
    if fixeanea == True:
        log.warning("gwaslab assumes SNPID is in the format of CHR:POS:NEA:EA / CHR:POS:REF:ALT", verbose=verbose)
        if overwrite is True:
            log.write(" -Overwrite mode is applied...", verbose=verbose)
            to_fix = is_chrposrefalt
        elif (nea in sumstats.columns) and (nea in sumstats.columns):
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
        elif (nea in sumstats.columns) and (ea not in sumstats.columns):
            log.write(" -Initiating EA columns...", verbose=verbose)
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
        elif (nea not in sumstats.columns) and (ea in sumstats.columns):
            log.write(" -Initiating NEA columns...", verbose=verbose)
            sumstats[nea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if sum(to_fix)>0 and verbose: log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...")
        else:
            log.write(" -Initiating EA and NEA columns...", verbose=verbose)
            sumstats[nea]=pd.Series(dtype="string")
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt
            if sum(to_fix)>0: 
                log.write(" -Number of variants could be fixed: "+str(sum(to_fix))+" ...", verbose=verbose)
    #                
        if sum(to_fix)>0:    
            log.write(" -Filling "+str(sum(to_fix))+" EA and NEA columns using SNPID's CHR:POS:NEA:EA...", verbose=verbose)
    #        
            if fixeanea_flip == True:
                log.write(" -Flipped : CHR:POS:NEA:EA -> CHR:POS:EA:NEA ", verbose=verbose)
                sumstats.loc[to_fix,ea] = sumstats.loc[to_fix,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[3]
                sumstats.loc[to_fix,nea] = sumstats.loc[to_fix,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[4]
            else:
                log.write(" -Chr:pos:a1:a2...a1->EA , a2->NEA ", verbose=verbose)
                sumstats.loc[to_fix,ea] = sumstats.loc[to_fix,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[4]
                sumstats.loc[to_fix,nea] = sumstats.loc[to_fix,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[3]
    
    #        #to_change_status = sumstats[status].str.match(r"\w\w\w[45]\w\w\w")
    #        #sumstats.loc[to_fix&to_change_status,status] = vchange_status(sumstats.loc[to_fix&to_change_status,status],4,"2")  
    #        #sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]).astype("string")
    #        #sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[1]).astype("Int64")
    
    ############################  fixing id ################################################### 
    if fixsep == True:
        if snpid in sumstats.columns: 
            log.write(' -Replacing [_-] in SNPID with ":" ...', verbose=verbose)
            sumstats[snpid] = sumstats[snpid].str.replace(r"[_-]",":",regex=True)
    
    if fixprefix == True:
        if snpid in sumstats.columns: 
            log.write(' -Removing /^chr/ in SNPID ...', verbose=verbose)
            prefix_removed = sumstats[snpid].str.extract(r'^(chr)?(\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[1]
            sumstats.loc[~prefix_removed.isna(),snpid] = prefix_removed[~prefix_removed.isna()]

    if fixid == True:
        if snpid not in sumstats.columns: 
        # initiate a SNPID column
            sumstats[snpid]=pd.Series(dtype="string")
        
        if (rsid in sumstats.columns) and (sum(is_rs_chrpos)>0) :
            sumstats[snpid]= sumstats.loc[is_rs_chrpos,rsid]
            
        if (chrom in sumstats.columns) and (pos in sumstats.columns):
            #only fix when CHR and POS is available
            pre_number=sum(sumstats[snpid].isna())                                
            
            if overwrite is False:
                #fix empty 
                to_fix = sumstats[snpid].isna() & sumstats[status].str.match( r"\w\w\w[0]\w\w\w", case=False, flags=0, na=False ) 
                #status_match(sumstats[status],4,0)
                #
            else:
                #fix all
                to_fix = sumstats[status].str.match( r"\w\w\w[0]\w\w\w", case=False, flags=0, na=False ) 
                # status_match(sumstats[status],4,0)
                #
            
            if (ea in sumstats.columns) and (nea in sumstats.columns):
            # when ea and nea is available  -> check status -> fix to chr:pos:nea:ea 
                
                pattern = r"\w\w\w[0]\w\w\w"  
                matched_index = sumstats[status].str.match(pattern)
                #matched_index = status_match(sumstats[status],4,0) #
                to_part_fix = matched_index & to_fix 
                
                #pattern = r"\w\w\w[0][01267][01234]\w"  
                pattern = r"\w\w\w\w[0123][01267][01234]"
                matched_index = sumstats[status].str.match(pattern)
                #status_match(sumstats[status],5,[0,1,2,3]) | status_match(sumstats[status],6,[0,1,2,6,7])| status_match(sumstats[status],7,[0,1,2,3,4])
                #
                if forcefixid is True:
                    matched_index = to_fix
                to_full_fix = matched_index & to_fix 
                
                
                if sum(to_part_fix)>0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                if sum(to_full_fix)>0:
                    sumstats.loc[to_full_fix,snpid] = sumstats.loc[to_full_fix,chrom].astype("string") + ":"+sumstats.loc[to_full_fix,pos].astype("string") +":"+ sumstats.loc[to_full_fix,nea].astype("string") +":"+ sumstats.loc[to_full_fix,ea].astype("string")
                log.write(" -Filling "+str(sum(to_part_fix)-sum(to_full_fix)) +" SNPID using CHR:POS...", verbose=verbose)
                log.write(" -Filling "+str(sum(to_full_fix)) +" SNPID using CHR:POS:NEA:EA...", verbose=verbose)
                sumstats.loc[(to_full_fix),status] = vchange_status(sumstats.loc[(to_full_fix),status],3,"975","630") 
                sumstats.loc[(to_part_fix),status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")  
                
            else:
            #when these is no ea or ena, just fix to chr:pos
                to_part_fix = to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull()
                log.write(" -Filling "+str(sum(to_part_fix)) +" SNPID using CHR POS...", verbose=verbose)
                if sum(to_part_fix)>0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                    sumstats.loc[to_part_fix,status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")
                    
            after_number=sum(sumstats[snpid].isna())
            log.write(" -Fixed "+ str(pre_number - after_number) +" variants ID...", verbose=verbose)
        else:
            log.write(" -ID unfixable: no CHR and POS columns or no SNPID. ", verbose=verbose)

    finished(log,verbose,_end_line)
    return sumstats

""

def stripSNPID(sumstats,snpid="SNPID",overwrite=False,verbose=True,log=Log()):  
    '''
    flip EA and NEA SNPid   CHR:POS:EA:NEA -> CHR:POS:NEA:EA
    '''
    ##start function with col checking##########################################################
    _start_line = "strip SNPID"
    _end_line = "stripping SNPID"
    _start_cols =["SNPID"]
    _start_function = ".strip_snpid()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    log.write(" -Checking if SNPID is (xxx:)CHR:POS:ATCG_Allele:ATCG_Allele(:xxx)...(separator: - ,: , _)",verbose=verbose)
    is_chrposrefalt = sumstats[snpid].str.contains(r'[:_-]?\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+[:_-]?', case=False, flags=0, na=False)
    # check if SNPID is NA
    is_snpid_na = sumstats[snpid].isna()
    
    log.write(" -Stripping {} non-NA fixable SNPIDs...".format(sum(is_chrposrefalt)),verbose=verbose)

    # flip 
    sumstats.loc[is_chrposrefalt,snpid] = \
        sumstats.loc[is_chrposrefalt,snpid].str.extract(r'[:_-]?(chr)?(\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+)[:_-]?',flags=re.IGNORECASE|re.ASCII)[1].astype("string")  

    finished(log,verbose,_end_line)
    return sumstats

def flipSNPID(sumstats,snpid="SNPID",overwrite=False,verbose=True,log=Log()):  
    '''
    flip EA and NEA SNPid   CHR:POS:EA:NEA -> CHR:POS:NEA:EA
    '''
    ##start function with col checking##########################################################
    _start_line = "flip SNPID"
    _end_line = "flipping SNPID"
    _start_cols =["SNPID"]
    _start_function = ".flip_snpid()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    log.warning("This function only flips alleles in SNPID without changing EA, NEA, STATUS or any statistics.")
    log.write(" -Checking if SNPID is CHR:POS:ATCG_Allele:ATCG_Allele...(separator: - ,: , _)",verbose=verbose)
    is_chrposrefalt = sumstats[snpid].str.match(r'^\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+$', case=False, flags=0, na=False)
    # check if SNPID is NA
    is_snpid_na = sumstats[snpid].isna()
    
    log.write(" -Flipping {} non-NA fixable SNPIDs...".format(sum(is_chrposrefalt)),verbose=verbose)

    # flip 
    sumstats.loc[is_chrposrefalt,snpid] = \
        sumstats.loc[is_chrposrefalt,snpid].str.extract(r'^(chr)?(\w+[:_-]\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[1].astype("string")  \
        + ":"+sumstats.loc[is_chrposrefalt,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[4].astype("string") \
        + ":"+sumstats.loc[is_chrposrefalt,snpid].str.extract(r'^(chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$',flags=re.IGNORECASE|re.ASCII)[3].astype("string")

    finished(log,verbose,_end_line)
    return sumstats

###############################################################################################################
# 20230128 
def removedup(sumstats,mode="dm",chrom="CHR",pos="POS",snpid="SNPID",ea="EA",nea="NEA",rsid="rsID",keep='first',keep_col="P",remove=False,keep_ascend=True,verbose=True,log=Log()):
    '''
    remove duplicate SNPs based on 1. SNPID,
    remove duplicate SNPs based on 2. CHR, POS, EA, and NEA
    remove duplicate SNPs based on 3. rsID
    remove multiallelic SNPs based on 4. CHR, POS
    '''

    ##start function with col checking##########################################################
    _start_line = "remove duplicated/multiallelic variants"
    _end_line = "removing duplicated/multiallelic variants"
    _start_cols =[]
    _start_function = ".remove_dup()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    log.write(" -Removing mode:{}".format(mode), verbose=verbose)
    # sort the variants using the specified column before removing
    if keep_col is not None : 
        if keep_col in sumstats.columns:
            log.write("Start to sort the sumstats using {}...".format(keep_col), verbose=verbose)
            sumstats = sumstats.sort_values(by=keep_col,ascending=keep_ascend)
        else:
            log.write("Column" + keep_col +" was not detected... skipping... ", verbose=verbose)
    total_number = len(sumstats)  
    
    # remove by duplicated SNPID
    if (snpid in sumstats.columns) and ("d" in mode or "s" in mode):
        log.write("Start to remove duplicated variants based on snpid...{}".format(_get_version()), verbose=verbose)
        check_dataframe_shape(sumstats, log, verbose)
        log.write(" -Which variant to keep: ",  keep , verbose=verbose)   
        pre_number =len(sumstats)   
        if snpid in sumstats.columns:
            # keep na and remove duplicated
            sumstats = sumstats.loc[sumstats[snpid].isna() | (~sumstats.duplicated(subset=[snpid], keep=keep)),:]
            after_number=len(sumstats)   
            log.write(" -Removed ",pre_number -after_number ," based on SNPID...", verbose=verbose)
    
    # remove by duplicated rsID
    if (rsid in sumstats.columns) and ("d" in mode or "r" in mode):
        # keep na and remove duplicated
        pre_number =len(sumstats)
        log.write("Start to remove duplicated variants based on rsID...", verbose=verbose)
        check_dataframe_shape(sumstats, log, verbose)
        sumstats = sumstats.loc[sumstats[rsid].isna() | (~sumstats.duplicated(subset=rsid, keep=keep)),:]
        after_number=len(sumstats)   
        log.write(" -Removed ",pre_number -after_number ," based on rsID...", verbose=verbose)
    
    # remove by duplicated variants by CHR:POS:NEA:EA
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and (nea in sumstats.columns) and (ea in sumstats.columns) and ("d" in mode or "c" in mode):
        log.write("Start to remove duplicated variants based on CHR,POS,EA and NEA...", verbose=verbose)
        check_dataframe_shape(sumstats, log, verbose)
        log.write(" -Which variant to keep: ",  keep , verbose=verbose)   
        pre_number =len(sumstats)   
        if snpid in sumstats.columns:
            # keep na and remove duplicated
            sumstats = sumstats.loc[(~sumstats[[chrom,pos,ea,nea]].all(axis=1)) | (~sumstats.duplicated(subset=[chrom,pos,ea,nea], keep=keep)),:]
            after_number=len(sumstats)   
            log.write(" -Removed ",pre_number -after_number ," based on CHR,POS,EA and NEA...", verbose=verbose) 
    
    # remove by multiallelic variants by CHR:POS
    if (chrom in sumstats.columns) and (pos in sumstats.columns) and "m" in mode:
        # keep na and remove duplicated
        pre_number =len(sumstats) 
        log.write("Start to remove multiallelic variants based on chr:pos...", verbose=verbose)    
        check_dataframe_shape(sumstats, log, verbose)
        log.write(" -Which variant to keep: ",  keep , verbose=verbose) 
        sumstats = sumstats.loc[(~sumstats[[chrom,pos]].all(axis=1)) | (~sumstats.duplicated(subset=[chrom,pos], keep=keep)),:]
        after_number=len(sumstats)  
        log.write(" -Removed ",pre_number -after_number," multiallelic variants...", verbose=verbose)   
    after_number=len(sumstats)   
    
    # resort the coordinates
    log.write(" -Removed ",total_number -after_number," variants in total.", verbose=verbose)
    if keep_col is not None : 
        log.write(" -Sort the coordinates based on CHR and POS...", verbose=verbose)
        sumstats = sortcoordinate(sumstats,verbose=False)
    
    if "n" in mode or remove==True:
        # if remove==True, remove NAs
        log.write(" -Removing NAs...", verbose=verbose)
        pre_number =len(sumstats) 
        specified_columns = []
        if "d" in mode:
            specified_columns.append(rsid)
            specified_columns.append(snpid)
            specified_columns.append(chrom)
            specified_columns.append(pos)
            specified_columns.append(ea)
            specified_columns.append(nea)
        if "r" in mode:
            specified_columns.append(rsid)
        if "s" in mode:
            specified_columns.append(snpid)
        if "m" in mode:
            specified_columns.append(chrom)
            specified_columns.append(pos)
        if "c" in mode:
            specified_columns.append(chrom)
            specified_columns.append(pos)
            specified_columns.append(ea)
            specified_columns.append(nea)
        sumstats = sumstats.loc[~sumstats[specified_columns].isna().any(axis=1),:]
        after_number=len(sumstats) 
        log.write(" -Removed ",pre_number -after_number," variants with NA values in {} .".format(set(specified_columns)), verbose=verbose)  

    finished(log,verbose,_end_line)
    return sumstats

###############################################################################################################
# 20230128
def fixchr(sumstats,chrom="CHR",status="STATUS",add_prefix="",x=("X",23),y=("Y",24),mt=("MT",25), remove=False, verbose=True, chrom_list = None, minchr=1,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "fix chromosome notation (CHR)"
    _end_line = "fixing chromosome notation (CHR)"
    _start_cols =[chrom,status]
    _start_function = ".fix_chr()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    #chrom_list = get_chr_list() #bottom 
    if chrom_list is None:
        chrom_list = get_chr_list()
    
    # convert to string datatype
    try:
        log.write(" -Checking CHR data type...", verbose=verbose)
        if sumstats[chrom].dtype == "string":
            pass
        else:
            sumstats[chrom] = sumstats[chrom].astype("string")
    except:
        log.write(" -Force converting to pd string data type...", verbose=verbose)
        sumstats[chrom] = sumstats[chrom].astype("string")
    
    # check if CHR is numeric
    is_chr_fixed = sumstats[chrom].str.isnumeric()
    # fill NAs with False
    is_chr_fixed = is_chr_fixed.fillna(False)
    log.write(" -Variants with standardized chromosome notation:",sum(is_chr_fixed), verbose=verbose)  
    
    # if there are variants whose CHR need to be fixed
    if sum(is_chr_fixed)<len(sumstats):
        
        #extract the CHR number or X Y M MT
        chr_extracted = sumstats.loc[~is_chr_fixed,chrom].str.extract(r'^(chr)?(\d{1,3}|[XYM]|MT)$',flags=re.IGNORECASE|re.ASCII)[1]

        is_chr_fixable = ~chr_extracted.isna()
        log.write(" -Variants with fixable chromosome notations:",sum(is_chr_fixable), verbose=verbose)  

        # For not fixed variants, check if na
        is_chr_na  = sumstats.loc[~is_chr_fixed, chrom].isna()
        if sum(is_chr_na)>0 and verbose: 
            log.write(" -Variants with NA chromosome notations:",sum(is_chr_na))  
        
        # Check variants with CHR being not NA and not fixable
        is_chr_invalid = (~is_chr_fixable)&(~is_chr_na)
        if sum(is_chr_invalid)>0 and verbose: 
            log.write(" -Variants with invalid chromosome notations:",sum(is_chr_invalid), verbose=verbose) 
            try:
                log.write(" -A look at invalid chromosome notations:" , set(sumstats.loc[~is_chr_fixed,chrom][is_chr_invalid].head()), verbose=verbose)
            except:
                pass
        else:
            log.write(" -No unrecognized chromosome notations...", verbose=verbose)
        
        # Assign good chr back to sumstats
        sumstats.loc[is_chr_fixable.index,chrom] = chr_extracted[is_chr_fixable.index]

        # X, Y, MT to 23,24,25
        xymt_list = [x[0].lower(),y[0].lower(),mt[0].lower(),x[0].upper(),y[0].upper(),mt[0].upper()]
        
        # check if sumstats contain sex CHR
        sex_chr = sumstats[chrom].isin(xymt_list)
        
        # if sumstats contain sex CHR
        if sum(sex_chr)>0:
            log.write(" -Identifying non-autosomal chromosomes : {}, {}, and {} ...".format(x[0],y[0],mt[0]), verbose=verbose)
            log.write(" -Identified ",str(sum(sex_chr))," variants on sex chromosomes...", verbose=verbose)
            
            # convert "X, Y, MT" to numbers
            convert_num_to_xymt={}
            if x[0].lower() in sumstats[chrom].values or x[0].upper() in sumstats[chrom].values:
                convert_num_to_xymt[x[0].lower()] = str(x[1])
                convert_num_to_xymt[x[0].upper()] = str(x[1])
                log.write(" -Standardizing sex chromosome notations: {} to {}...".format(x[0], x[1]), verbose=verbose)
            if y[0].lower() in sumstats[chrom].values or y[0].upper() in sumstats[chrom].values:
                convert_num_to_xymt[y[0].lower()] = str(y[1])
                convert_num_to_xymt[y[0].upper()] = str(y[1])
                log.write(" -Standardizing sex chromosome notations: {} to {}...".format(y[0], y[1]), verbose=verbose)
            if mt[0].lower() in sumstats[chrom].values or mt[0].upper() in sumstats[chrom].values:
                convert_num_to_xymt[mt[0].lower()] = str(mt[1])
                convert_num_to_xymt[mt[0].upper()] = str(mt[1])
                log.write(" -Standardizing sex chromosome notations: {} to {}...".format(mt[0], mt[1]), verbose=verbose)
            sumstats.loc[sex_chr,chrom] =sumstats.loc[sex_chr,chrom].map(convert_num_to_xymt)
        
        # change status code
        sumstats.loc[is_chr_fixed,status] = vchange_status(sumstats.loc[is_chr_fixed,status],4,"986","520")
        if len(is_chr_fixable.index)>0:
            sumstats.loc[is_chr_fixable.index,status] = vchange_status(sumstats.loc[is_chr_fixable.index,status],4,"986","520")
        if len(is_chr_fixable.index)>0:
            sumstats.loc[is_chr_invalid.index,status] = vchange_status(sumstats.loc[is_chr_invalid.index,status],4,"986","743")
        
        # check variants with unrecognized CHR
        unrecognized_num = sum(~sumstats[chrom].isin(chrom_list))
        if (remove is True) and unrecognized_num>0:
            # remove variants with unrecognized CHR 
            try:
                log.write(" -Valid CHR list: {} - {}".format(min([int(x) for x in chrom_list if x.isnumeric()]),max([int(x) for x in chrom_list if x.isnumeric()])), verbose=verbose)
            except:
                pass
            log.write(" -Removed "+ str(unrecognized_num)+ " variants with chromosome notations not in CHR list.", verbose=verbose) 
            try:
                log.write(" -A look at chromosome notations not in CHR list:" , set(sumstats.loc[~sumstats[chrom].isin(chrom_list),chrom].head()), verbose=verbose)
            except:
                pass
            #sumstats = sumstats.loc[sumstats.index[sumstats[chrom].isin(chrom_list)],:]
            good_chr = sumstats[chrom].isin(chrom_list)
            sumstats = sumstats.loc[good_chr, :].copy()
    else:
        log.write(" -All CHR are already fixed...", verbose=verbose) 
        sumstats.loc[is_chr_fixed,status] = vchange_status(sumstats.loc[is_chr_fixed,status],4,"986","520")
    
    # Convert string to int
    try:
        sumstats[chrom] = sumstats[chrom].astype('Int64')
    except:
    #    # force convert
        sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    
    # filter out variants with CHR <=0
    out_of_range_chr = sumstats[chrom] < minchr
    out_of_range_chr = out_of_range_chr.fillna(False)
    if sum(out_of_range_chr)>0:
        log.write(" -Sanity check for CHR...", verbose=verbose) 
        log.write(" -Removed {} variants with CHR < {}...".format(sum(out_of_range_chr),minchr), verbose=verbose)
        sumstats = sumstats.loc[~out_of_range_chr,:]

    finished(log,verbose,_end_line)
    return sumstats

###############################################################################################################    
# 20230128
def fixpos(sumstats,pos="POS",status="STATUS",remove=False, verbose=True, lower_limit=0 , upper_limit=None , limit=250000000, log=Log()):
    ##start function with col checking##########################################################
    _start_line = "fix basepair positions (POS)"
    _end_line = "fixing basepair positions (POS)"
    _start_cols =[pos,status]
    _start_function = ".fix_pos()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    if upper_limit is None:
        upper_limit = limit
        
    all_var_num = len(sumstats)
    #convert to numeric
    is_pos_na = sumstats[pos].isna()
    
    try:
        if str(sumstats[pos].dtype) == "string" or str(sumstats[pos].dtype) == "object":
            sumstats[pos] = sumstats[pos].astype('string')
            # if so, remove thousands separator
            log.write(' -Removing thousands separator "," or underbar "_" ...', verbose=verbose)
            sumstats.loc[~is_pos_na, pos] = sumstats.loc[~is_pos_na, pos].str.replace(r'[,_]', '' ,regex=True)
    except:
        pass

    # convert POS to integer
    try:
        log.write(' -Converting to Int64 data type ...', verbose=verbose)
        sumstats[pos] = sumstats[pos].astype('Int64')
    except:
        log.write(' -Force converting to Int64 data type ...', verbose=verbose)
        sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    is_pos_fixed = ~sumstats[pos].isna()
    is_pos_invalid = (~is_pos_na)&(~is_pos_fixed)
    
    sumstats.loc[is_pos_fixed,status]   = vchange_status(sumstats.loc[is_pos_fixed,status]  ,4,"975","630")
    sumstats.loc[is_pos_invalid,status] = vchange_status(sumstats.loc[is_pos_invalid,status],4,"975","842")
    
    # remove outlier, limit:250,000,000
    log.write(" -Position bound:({} , {:,})".format(lower_limit, upper_limit), verbose=verbose)
    is_pos_na = sumstats[pos].isna()
    out_lier= ((sumstats[pos]<=lower_limit) | (sumstats[pos]>=upper_limit)) & (~is_pos_na)
    log.write(" -Removed outliers:",sum(out_lier), verbose=verbose)
    sumstats = sumstats.loc[~out_lier,:]
    #remove na
    if remove is True: 
        sumstats = sumstats.loc[~sumstats[pos].isna(),:]
        remain_var_num = len(sumstats)
        log.write(" -Removed "+str(all_var_num - remain_var_num)+" variants with bad positions.", verbose=verbose)        
 
    finished(log,verbose,_end_line)
    return sumstats

###############################################################################################################    
# 20220514
def fixallele(sumstats,ea="EA", nea="NEA",status="STATUS",remove=False,verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "fix alleles (EA and NEA)"
    _end_line = "fixing alleles (EA and NEA)"
    _start_cols =[ea, nea,status]
    _start_function = ".fix_allele()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################
    #try:
    #    ea_missing = sum(sumstats[ea].isna())
    #    nea_missing = sum(sumstats[nea].isna())
    #    if sum(ea_missing)>0:
    #        log.write(" -Converting {} missing EA to letter N.".format(ea_missing))
    #        sumstats[ea] = sumstats[ea].add_categories("N").fillna("N")
    #    if sum(sumstats[nea].isna())>0:
    #        log.write(" -Converting {} missing NEA to letter N.".format(nea_missing))
    #        sumstats[nea] = sumstats[nea].add_categories("N").fillna("N")
    #except:
    #    pass

    log.write(" -Converted all bases to string datatype and UPPERCASE.", verbose=verbose)
    categories = set(sumstats[ea].str.upper())|set(sumstats[nea].str.upper())|set("N")
    categories = {x for x in categories if pd.notna(x)}
    sumstats[ea]=pd.Categorical(sumstats[ea].str.upper(),categories = categories) 
    sumstats[nea]=pd.Categorical(sumstats[nea].str.upper(),categories = categories) 
    all_var_num = len(sumstats)
    
    ## check ATCG
    bad_ea = sumstats[ea].str.contains("[^actgACTG]",na=True)
    bad_nea = sumstats[nea].str.contains("[^actgACTG]",na=True)
    good_ea  = ~bad_ea
    good_nea = ~bad_nea

    log.write(" -Variants with bad EA  : {}".format(sum(bad_ea)), verbose=verbose)
    log.write(" -Variants with bad NEA : {}".format(sum(bad_nea)), verbose=verbose)
    
    ## check NA
    is_eanea_na = sumstats[ea].isna() |  sumstats[nea].isna()
    log.write(" -Variants with NA for EA or NEA: {}".format(sum(is_eanea_na)), verbose=verbose)
    
    ## check same alleles
    not_variant = sumstats[nea] == sumstats[ea]
    log.write(" -Variants with same EA and NEA: {}".format(sum(not_variant)), verbose=verbose)

    ## sum up invalid variants
    is_invalid = bad_ea | bad_nea | not_variant
    
    exclude  = bad_nea | bad_ea
    
    if len(set(sumstats.loc[bad_ea,ea].head())) >0:
        log.write(" -A look at the non-ATCG EA:",set(sumstats.loc[bad_ea,ea].head()),"...", verbose=verbose) 
    if len(set(sumstats.loc[bad_nea,nea].head())) >0:
        log.write(" -A look at the non-ATCG NEA:",set(sumstats.loc[bad_nea,nea].head()),"...", verbose=verbose) 
    
    if remove == True:
        sumstats = sumstats.loc[(good_ea & good_nea),:].copy()
        good_eanea_num = len(sumstats)
        log.write(" -Removed "+str(all_var_num - good_eanea_num)+" variants with NA alleles or alleles that contain bases other than A/C/T/G.", verbose=verbose)  
        sumstats = sumstats.loc[(good_ea & good_nea & (~not_variant)),:].copy()
        good_eanea_notsame_num = len(sumstats)
        log.write(" -Removed "+str(good_eanea_num - good_eanea_notsame_num)+" variants with same allele for EA and NEA.", verbose=verbose) 
    else:
        sumstats[[ea,nea]] = sumstats[[ea,nea]].fillna("N")
        log.write(" -Detected "+str(sum(exclude))+" variants with alleles that contain bases other than A/C/T/G .", verbose=verbose) 
    categories = set(sumstats[ea].str.upper())|set(sumstats[nea].str.upper())|set("N")
    sumstats[ea]=pd.Categorical(sumstats[ea].str.upper(),categories = categories) 
    sumstats[nea]=pd.Categorical(sumstats[nea].str.upper(),categories = categories) 
    
    is_eanea_fixed = good_ea | good_nea
    is_snp = (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()==1)
    is_indel = (sumstats[ea].str.len()!=sumstats[nea].str.len())
    is_not_normalized = (sumstats[ea].str.len()>1) &(sumstats[nea].str.len()>1)
    is_normalized = is_indel &( (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()>1) | (sumstats[ea].str.len()>1) &(sumstats[nea].str.len()==1) )
    
    if sum(is_invalid)>0:
        sumstats.loc[is_invalid, status]                      = vchange_status(sumstats.loc[is_invalid,status],                5,"9","6") 
    if sum(is_eanea_na)>0:
        sumstats.loc[is_eanea_na,status]                      = vchange_status(sumstats.loc[is_eanea_na, status],              5,"9","7")
    if sum(is_eanea_fixed&is_not_normalized)>0:
        sumstats.loc[is_eanea_fixed&is_not_normalized,status] = vchange_status(sumstats.loc[is_eanea_fixed&is_not_normalized,status], 5,"9","5")
    if sum(is_eanea_fixed&is_snp)>0:
        sumstats.loc[is_eanea_fixed&is_snp, status]           = vchange_status(sumstats.loc[is_eanea_fixed&is_snp,status],        5,"9","0")
    if sum(is_eanea_fixed&is_indel)>0:
        sumstats.loc[is_eanea_fixed&is_indel,status]          = vchange_status(sumstats.loc[is_eanea_fixed&is_indel, status],      5,"9","4")
    if sum(is_eanea_fixed&is_normalized)>0:
        sumstats.loc[is_eanea_fixed&is_normalized,status]     = vchange_status(sumstats.loc[is_eanea_fixed&is_normalized, status],  5,"4","3")

    finished(log,verbose,_end_line)
    return sumstats

###############################################################################################################   
# 20220721

def parallelnormalizeallele(sumstats,mode="s",snpid="SNPID",rsid="rsID",pos="POS",nea="NEA",ea="EA" ,status="STATUS",chunk=3000000,n_cores=1,verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "normalize indels"
    _end_line = "normalizing indels"
    _start_cols =[ea, nea,status]
    _start_function = ".normalize()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    #variants_to_check = status_match(sumstats[status],5,[4,5]) #
    #r'\w\w\w\w[45]\w\w'
    variants_to_check = sumstats[status].str[4].str.match(r'4|5', case=False, flags=0, na=False)
    if sum(variants_to_check)==0:
        log.write(" -No available variants to normalize..", verbose=verbose)
        log.write("Finished normalizing variants successfully!", verbose=verbose)
        return sumstats
    ###############################################################################################################
    if mode=="v":
        if sum(variants_to_check)<100000: 
            n_cores=1  
        if n_cores==1:
            normalized_pd, changed_index = fastnormalizeallele(sumstats.loc[variants_to_check,[pos,nea,ea,status]],pos=pos ,nea=nea,ea=ea,status=status,chunk=chunk, log=log, verbose=verbose)
        else:
            pool = Pool(n_cores)
            map_func = partial(fastnormalizeallele,pos=pos,nea=nea,ea=ea,status=status)
            df_split = _df_split(sumstats.loc[variants_to_check,[pos,nea,ea,status]], n_cores)
            results = pool.map(map_func,df_split)
            normalized_pd = pd.concat([i[0] for i in results])
            changed_index = np.concatenate([i[1] for i in results])
            del results
            pool.close()
            pool.join()
            gc.collect()
        ###############################################################################################################
        try:
            example_sumstats = sumstats.loc[changed_index,:].head()
            changed_num = len(changed_index)
            if changed_num>0:
                if snpid in example_sumstats.columns:
                    before_normalize_id = example_sumstats.loc[variants_to_check,snpid]
                elif rsid in example_sumstats.columns:
                    before_normalize_id = example_sumstats.loc[variants_to_check,rsid]
                else:
                    before_normalize_id = example_sumstats.index
                
                log.write(" -Not normalized allele IDs:",end="", verbose=verbose)
                for i in before_normalize_id.values:
                    log.write(i,end=" ",show_time=False)
                log.write("... \n",end="",show_time=False, verbose=verbose) 
            
                log.write(" -Not normalized allele:",end="", verbose=verbose)
                for i in example_sumstats[[ea,nea]].values:
                    log.write(i,end="",show_time=False, verbose=verbose)
                log.write("... \n",end="",show_time=False, verbose=verbose)     
                log.write(" -Modified "+str(changed_num) +" variants according to parsimony and left alignment principal.", verbose=verbose)
            else:
                log.write(" -All variants are already normalized..", verbose=verbose)
        except:
            pass

    ##########################################################################################################################################################
    elif mode=="s":
        if sum(variants_to_check)<10000: 
            n_cores=1  
        pool = Pool(n_cores)
        map_func = partial(normalizeallele,pos=pos,nea=nea,ea=ea,status=status)
        #df_split = np.array_split(sumstats.loc[variants_to_check,[pos,nea,ea,status]], n_cores)
        df_split = _df_split(sumstats.loc[variants_to_check,[pos,nea,ea,status]], n_cores)
        normalized_pd = pd.concat(pool.map(map_func,df_split))
        pool.close()
        pool.join()

        before_normalize = sumstats.loc[variants_to_check,[ea,nea]]
        changed_num = len(normalized_pd.loc[(before_normalize[ea]!=normalized_pd[ea]) | (before_normalize[nea]!=normalized_pd[nea]),:])
        if changed_num>0:
            if snpid in sumstats.columns:
                before_normalize_id = sumstats.loc[variants_to_check,snpid]
            elif rsid in sumstats.columns:
                before_normalize_id = sumstats.loc[variants_to_check,rsid]
            else:
                before_normalize_id = pd.DataFrame(sumstats.index[variants_to_check],index=sumstats.index[variants_to_check])
            
            log.write(" -Not normalized allele IDs:",end="", verbose=verbose)
            for i in before_normalize_id.loc[(before_normalize[ea]!=normalized_pd[ea]) | (before_normalize[nea]!=normalized_pd[nea])].head().values:
                log.write(i,end=" ",show_time=False)
            log.write("... \n",end="",show_time=False, verbose=verbose) 
        
            log.write(" -Not normalized allele:",end="", verbose=verbose)
            for i in before_normalize.loc[(before_normalize[ea]!=normalized_pd[ea]) | (before_normalize[nea]!=normalized_pd[nea]),[ea,nea]].head().values:
                log.write(i,end="",show_time=False, verbose=verbose)
            log.write("... \n",end="",show_time=False, verbose=verbose)     
            log.write(" -Modified "+str(changed_num) +" variants according to parsimony and left alignment principal.", verbose=verbose)
        else:
            log.write(" -All variants are already normalized..", verbose=verbose)
        ###################################################################################################################
    
    categories = set(sumstats[ea])|set(sumstats[nea]) |set(normalized_pd.loc[:,ea]) |set(normalized_pd.loc[:,nea])
    sumstats[ea]  = pd.Categorical(sumstats[ea],categories = categories) 
    sumstats[nea] = pd.Categorical(sumstats[nea],categories = categories ) 
    sumstats.loc[variants_to_check,[pos,nea,ea,status]] = normalized_pd.values
    
    try:
        sumstats[pos] = sumstats[pos].astype('Int64')
    except:
        sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
  
    finished(log,verbose,_end_line)
    return sumstats

def normalizeallele(sumstats,pos="POS" ,nea="NEA",ea="EA",status="STATUS"):
    #single df
    #normalized = sumstats.apply(lambda x: normalizevariant(x[0],x[1],x[2],x[3]),axis=1)
    normalized = sumstats.apply(lambda x: normalizevariant(x[pos],x[nea],x[ea],x[status]),axis=1)
    sumstats = pd.DataFrame(normalized.to_list(), columns=[pos,nea,ea,status],index=sumstats.index)
    return sumstats

def fastnormalizeallele(insumstats,pos="POS" ,nea="NEA",ea="EA",status="STATUS",chunk=3000000,log=Log(),verbose=False):
    log.write(" -Number of variants to check:{}".format(len(insumstats)), verbose=verbose)
    log.write(" -Chunk size:{}".format(chunk), verbose=verbose)
    log.write(" -Processing in chunks:",end="", verbose=verbose)
    changed_index = np.array([])
    for part_n in range(len(insumstats)//chunk+1):
        log.write(part_n, end=" ",show_time=False, verbose=verbose)
        insumstats["NEA"] = insumstats["NEA"].astype("string")
        insumstats["EA"] = insumstats["EA"].astype("string")
        insumstats.iloc[part_n*chunk:(part_n+1)*chunk,:],changed_index_single  = normalizae_chunk(insumstats.iloc[part_n*chunk:(part_n+1)*chunk,:].copy())
        changed_index = np.concatenate([changed_index,changed_index_single])
        gc.collect()
    log.write("\n",end="",show_time=False, verbose=verbose)   
    return insumstats, changed_index

def normalizae_chunk(sumstats,pos="POS" ,nea="NEA",ea="EA",status="STATUS"):
    # already normalized

    is_same = sumstats["NEA"] == sumstats["EA"]
    is_normalized = ((sumstats["NEA"].str.len()==1) | (sumstats["EA"].str.len()==1) ) & (~is_same)
    
    # a series to keep tracking of variants that are modified
    changed = sumstats["NEA"] != sumstats["NEA"]
    
    # right side
    ea_len = sumstats["NEA"].str.len()
    nea_len = sumstats["EA"].str.len()
    max_length=max(ea_len.max(), nea_len.max())

    for i in range(1, max_length):
        is_pop = (sumstats["NEA"].str[-1] == sumstats["EA"].str[-1]) & (~is_normalized)
        if sum(is_pop)==0:
            break
        if i ==1:
            changed = changed | is_pop
        nea_len[is_pop] = nea_len[is_pop] -1 
        ea_len[is_pop] = ea_len[is_pop] -1
        sumstats.loc[is_pop, "NEA"] = sumstats.loc[is_pop,"NEA"].str[:-1]
        sumstats.loc[is_pop, "EA"] = sumstats.loc[is_pop,"EA"].str[:-1]
        is_normalized = ((sumstats["NEA"].str.len()==1) | (sumstats["EA"].str.len()==1) ) & (~is_same)
        gc.collect()
    
    # left side 
    max_length=max(sumstats["NEA"].str.len().max(), sumstats["EA"].str.len().max())
    for i in range(1, max_length):
        is_pop = (sumstats["NEA"].str[0] == sumstats["EA"].str[0]) & (~is_normalized)
        if sum(is_pop)==0:
            break
        if i ==1:
            changed = changed | is_pop
        sumstats.loc[is_pop, "NEA"] = sumstats.loc[is_pop,"NEA"].str[1:]
        sumstats.loc[is_pop, "EA"] = sumstats.loc[is_pop,"EA"].str[1:]
        sumstats.loc[is_pop, "POS"] = sumstats.loc[is_pop,"POS"] + 1
        is_normalized = ((sumstats["NEA"].str.len()==1) | (sumstats["EA"].str.len()==1) ) & (~is_same)
        gc.collect()
    
    sumstats.loc[is_normalized,status]     = vchange_status(sumstats.loc[is_normalized, status],  5,"4","0")
    sumstats.loc[is_same,status]     = vchange_status(sumstats.loc[is_same, status],  5,"4","3") 
    changed_index = sumstats[changed].index
    return sumstats, changed_index.values

def normalizevariant(pos,a,b,status):
    # single record
    # https://genome.sph.umich.edu/wiki/Variant_Normalization
    # a - ref - nea    starting -> pos
    # b - alt - ea
    # status xx1xx /xx2xx /xx9xx
    
    status_pre=status[:4]
    status_end=status[5:]
    
    if len(a)==1 or len(b)==1:
        #return pos,a,b,change_status(status,5,0) #status_pre+"0"+status_end
        return pos,a,b,status_pre+"0"+status_end
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
            if len(a)==1 or len(b)==1:
                return pos,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"0"+status_end #change_status(status,5,0) #
            return pos,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"3"+status_end #change_status(status,5,3) #

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
    if len(a)==1 or len(b)==1:
        return pos+pos_change,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"0"+status_end #change_status(status,5,0) #
    return pos+pos_change,a[pointer_a_l:pointer_a_r+1],b[pointer_b_l:pointer_b_r+1],status_pre+"3"+status_end # change_status(status,5,3) #
""


###############################################################################################################
# 20220426
def add_tolerence(stats, float_tolerence, mode):
    if "l" in mode:
        stats = (stats[0] - float_tolerence if stats[0]!=float("Inf") else float("Inf"), stats[1])
    if "r" in mode:
        stats = (stats[0] , stats[1] + float_tolerence if stats[0]!=float("Inf") else float("Inf"))
    return stats


def check_range(sumstats, var_range, header, coltocheck, cols_to_check, log, verbose, dtype="Int64"):
    pre_number=len(sumstats)
    if header in coltocheck and header in sumstats.columns:
        cols_to_check.append(header)
        if header=="STATUS": 
            log.write(" -Checking STATUS and converting STATUS to categories....", verbose=verbose) 
            categories = {str(j+i) for j in [1300000,1800000,1900000,3800000,9700000,9800000,9900000] for i in range(0,100000)}
            sumstats[header] = pd.Categorical(sumstats[header],categories=categories)
            return sumstats
        
        if dtype in ["Int64","Int32","int","int32","in64"]:
            log.write(" -Checking if {} <= {} <= {} ...".format( var_range[0] ,header, var_range[1]), verbose=verbose) 
            sumstats[header] = np.floor(pd.to_numeric(sumstats[header], errors='coerce')).astype(dtype)
        
        elif dtype in ["Float64","Float32","float","float64","float32"]:
            log.write(" -Checking if {} < {} < {} ...".format( var_range[0] ,header, var_range[1]),verbose=verbose) 
            sumstats[header] = pd.to_numeric(sumstats[header], errors='coerce').astype(dtype)
        
        is_valid = (sumstats[header]>=var_range[0]) & (sumstats[header]<=var_range[1])
        is_valid = is_valid.fillna(False)

        if header=="P":
            is_low_p =  sumstats["P"] == 0 
            if sum(is_low_p) >0:
                log.warning("Extremely low P detected (P=0 or P < minimum positive value of float64) : {}".format(sum(is_low_p)))
                log.warning("Please consider using MLOG10P instead.")
        
        if header=="INFO":
            is_high_info =  sumstats["INFO"]>1 
            if sum(is_high_info) >0:
                log.warning("High INFO detected (INFO>1) : {}".format(sum(is_high_info)))
                log.warning("max(INFO): {}".format(sumstats["INFO"].max()))
                log.warning("Please check if this is as expected.")

        if sum(~is_valid)>0:
            try:
                if "SNPID" in sumstats.columns:
                    id_to_use = "SNPID"
                elif "rsID" in sumstats.columns:
                    id_to_use = "rsID"
                invalid_ids = sumstats.loc[~is_valid, id_to_use].head().astype("string")
                invalid_values = sumstats.loc[~is_valid, header].head().astype("string").fillna("NA")
                log.write("  -Examples of invalid variants({}): {} ...".format(id_to_use, ",".join(invalid_ids.to_list()) ), verbose=verbose) 
                log.write("  -Examples of invalid values ({}): {} ...".format(header, ",".join(invalid_values.to_list()) ), verbose=verbose) 
            except:
                pass

        sumstats = sumstats.loc[is_valid,:]
        after_number=len(sumstats)
        log.write(" -Removed {} variants with bad/na {}.".format(pre_number - after_number, header), verbose=verbose) 
    return sumstats

def sanitycheckstats(sumstats,
                     coltocheck=None,
                     n=(0,2**31-1),
                     ncase=(0,2**31-1),
                     ncontrol=(0,2**31-1),
                     eaf=(0,1),
                     mac=(0,2**31-1),
                     maf=(0,0.5),
                     chisq=(0,float("Inf")),
                     z=(-9999,9999),
                     t=(-99999,99999),
                     f=(0,float("Inf")),
                     p=(0,1),
                     mlog10p=(0,99999),
                     beta=(-100,100),
                     se=(0,float("Inf")),
                     OR=(-100,100),
                     OR_95L=(0,float("Inf")),
                     OR_95U=(0,float("Inf")),
                     HR=(-100,100),
                     HR_95L=(0,float("Inf")),
                     HR_95U=(0,float("Inf")),
                     info=(0,2),
                     float_tolerence = 1e-7,
                     verbose=True,
                     log=Log()):
    '''
    Sanity check (default; v3.4.33):
        N:      Int32    , N>0 , 
        EAF:    float32  , 0<= EAF <=1, 
        P:      float64  , 0< P < 1, 
        BETA:   float64  , abs(BETA) <10
        SE:     float64  , SE > 0
        OR:     float64  , -10 <log(OR) <10
        OR_95L: float64  , OR_95L >0
        OR_95U: float64  , OR_95L >0
        HR:     float64  , -10 <log(HR) <10
        HR_95L: float64  , HR_95L >0
        HR_95U: float64  , HR_95L >0
        INFO:   float32  , 1>=INFO>0
        Z       float64  , -9999 < Z < 9999
        T       float64  , -99999 < T < 99999
        F       float64  , F > 0 
    '''
    ##start function with col checking##########################################################
    _start_line = "perform sanity check for statistics"
    _end_line = "sanity check for statistics"
    _start_cols =[]
    _start_function = ".check_sanity()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    log.write(" -Comparison tolerance for floats: {}".format(float_tolerence), verbose=verbose) 
    eaf = add_tolerence(eaf, float_tolerence, "lr")
    maf = add_tolerence(maf, float_tolerence, "lr")
    beta = add_tolerence(beta, float_tolerence, "lr")
    se = add_tolerence(se, float_tolerence, "lr")
    mlog10p = add_tolerence(mlog10p, float_tolerence, "lr")
    OR = add_tolerence(OR, float_tolerence, "lr")
    OR_95L = add_tolerence(OR_95L, float_tolerence, "lr")
    OR_95U = add_tolerence(OR_95U, float_tolerence, "lr")
    HR = add_tolerence(HR, float_tolerence, "lr")
    HR_95L = add_tolerence(HR_95L, float_tolerence, "lr")
    HR_95U = add_tolerence(HR_95U, float_tolerence, "lr")
    info = add_tolerence(info, float_tolerence, "lr")
    z = add_tolerence(z, float_tolerence, "lr")
    p = add_tolerence(p, float_tolerence, "lr")
    f = add_tolerence(f, float_tolerence, "lr")
    chisq = add_tolerence(chisq, float_tolerence, "lr")
    ############################################################################################
    ## add direction
    if coltocheck is None:
        coltocheck = ["P","MLOG10P","INFO","Z","BETA","SE","EAF","CHISQ","F","N","N_CASE","N_CONTROL","OR","OR_95L","OR_95U","HR","HR_95L","HR_95U","STATUS"]
    
    cols_to_check=[]
    oringinal_number=len(sumstats)
    sumstats = sumstats.copy()
    
    ###Int64 ################################################################################################################################################
    sumstats = check_range(sumstats, var_range=n, header="N", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")
    sumstats = check_range(sumstats, var_range=ncase, header="N_CASE", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")
    sumstats = check_range(sumstats, var_range=ncontrol, header="N_CONTROL", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="Int64")

    ###float32 ################################################################################################################################################
    sumstats = check_range(sumstats, var_range=eaf, header="EAF", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float32")
    sumstats = check_range(sumstats, var_range=maf, header="MAF", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float32")
    sumstats = check_range(sumstats, var_range=info, header="INFO", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float32")

    ###float64 ################################################################################################################################################
    sumstats = check_range(sumstats, var_range=chisq, header="CHISQ", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=z, header="Z", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=t, header="T", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=f, header="F", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=p, header="P", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=mlog10p, header="MLOG10P", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=beta, header="BETA", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=se, header="SE", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=OR, header="OR", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=OR_95L, header="OR_95L", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=OR_95U, header="OR_95U", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=HR, header="HR", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=HR_95L, header="HR_95L", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    sumstats = check_range(sumstats, var_range=HR_95U, header="HR_95U", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="float64")
    ###STATUS ###############################################################################################################################################
    sumstats = check_range(sumstats, var_range=None, header="STATUS", coltocheck=coltocheck, cols_to_check=cols_to_check, log=log, verbose=verbose, dtype="category")

    after_number=len(sumstats)
    log.write(" -Removed "+str(oringinal_number - after_number)+" variants with bad statistics in total.",verbose=verbose) 
    log.write(" -Data types for each column:",verbose=verbose)
    check_datatype(sumstats,verbose=verbose, log=log)
    finished(log,verbose,_end_line)
    return sumstats

### check consistency #############################################################################################################################################

def _check_data_consistency(sumstats, beta="BETA", se="SE", p="P",mlog10p="MLOG10P",rtol=1e-3, atol=1e-3, equal_nan=True, verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "check data consistency across columns"
    _end_line = "checking data consistency across columns"
    _start_cols =[]
    _start_function = ".check_data_consistency()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    log.write(" -Tolerance: {} (Relative) and {} (Absolute)".format(rtol, atol),verbose=verbose)
    check_status = 0
    
    if "SNPID" in sumstats.columns:
        id_to_use = "SNPID"
    elif "rsID" in sumstats.columns:
        id_to_use = "rsID"
    else:
        log.write(" -SNPID/rsID not available...SKipping",verbose=verbose)
        log.write("Finished checking data consistency across columns.",verbose=verbose) 
        return 0
    
    
    if "BETA" in sumstats.columns and "SE" in sumstats.columns:
        if "MLOG10P" in sumstats.columns:
            log.write(" -Checking if BETA/SE-derived-MLOG10P is consistent with MLOG10P...",verbose=verbose)
            betase_derived_mlog10p =  _convert_betase_to_mlog10p(sumstats["BETA"], sumstats["SE"])
            is_close = np.isclose(betase_derived_mlog10p, sumstats["MLOG10P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
            diff = betase_derived_mlog10p - sumstats["MLOG10P"]
            if sum(~is_close)>0:
                log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
                log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
            else:
                log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
            check_status=1
        
        if "P" in sumstats.columns:
            log.write(" -Checking if BETA/SE-derived-P is consistent with P...",verbose=verbose)
            betase_derived_p =  _convert_betase_to_p(sumstats["BETA"], sumstats["SE"])
            is_close = np.isclose(betase_derived_p, sumstats["P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
            diff = betase_derived_p - sumstats["P"]
            if sum(~is_close)>0:
                log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
                log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
            else:
                log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
            check_status=1
    
    if "MLOG10P" in sumstats.columns and "P" in sumstats.columns:
        log.write(" -Checking if MLOG10P-derived-P is consistent with P...",verbose=verbose)
        mlog10p_derived_p = _convert_mlog10p_to_p(sumstats["MLOG10P"])
        is_close = np.isclose(mlog10p_derived_p, sumstats["P"], rtol=rtol, atol=atol, equal_nan=equal_nan)
        diff = mlog10p_derived_p - sumstats["P"]
        if sum(~is_close)>0:
            log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
            log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
        else:
            log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
        check_status=1

    if "N" in sumstats.columns and "N_CONTROL" in sumstats.columns and "N_CASE" in sumstats.columns:
        log.write(" -Checking if N is consistent with N_CASE + N_CONTROL ...", verbose=verbose) 
        is_close = sumstats["N"] == sumstats["N_CASE"] + sumstats["N_CONTROL"] 
        #is_close = np.isclose(sumstats["N"], sumstats["N_CASE"] + sumstats["N_CONTROL"] , rtol=rtol, atol=atol, equal_nan=equal_nan)
        diff = abs(sumstats["N"] - (sumstats["N_CASE"] + sumstats["N_CONTROL"] ))
        if sum(~is_close)>0:
            log.write("  -Not consistent: {} variant(s)".format(sum(~is_close)),verbose=verbose)
            log.write("  -Variant {} with max difference: {} with {}".format(id_to_use, sumstats.loc[diff.idxmax(),id_to_use], diff.max()),verbose=verbose)
        else:
            log.write("  -Variants with inconsistent values were not detected." ,verbose=verbose)
        check_status=1
        
    if check_status==1:
        log.write(" -Note: if the max difference is greater than expected, please check your original sumstats.",verbose=verbose)
    else:
        log.write(" -No availalbe columns for data consistency checking...Skipping...",verbose=verbose)
    finished(log,verbose,_end_line)

###############################################################################################################
# 20220426
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

def flip_by_swap(sumstats, matched_index, log, verbose):
    if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
        log.write(" -Swapping column: NEA <=> EA...", verbose=verbose) 
        sumstats.loc[matched_index,['NEA','EA']] = sumstats.loc[matched_index,['EA','NEA']].values
    return sumstats

def flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1):
    if "OR" in sumstats.columns:
        log.write(" -Flipping column: OR = 1 / OR...", verbose=verbose) 
        sumstats.loc[matched_index,"OR"] =   factor / sumstats.loc[matched_index,"OR"].values
    if "OR_95L" in sumstats.columns:
        log.write(" -Flipping column: OR_95U = 1 / OR_95L...", verbose=verbose) 
        sumstats.loc[matched_index,"OR_95U"] =   factor / sumstats.loc[matched_index,"OR_95L"].values
    if "OR_95U" in sumstats.columns:
        log.write(" -Flipping column: OR_95L = 1 / OR_95U...", verbose=verbose) 
        sumstats.loc[matched_index,"OR_95L"] =   factor / sumstats.loc[matched_index,"OR_95U"].values
    if "HR" in sumstats.columns:
        log.write(" -Flipping column: HR = 1 / HR...", verbose=verbose) 
        sumstats.loc[matched_index,"HR"] =   factor / sumstats.loc[matched_index,"HR"].values
    if "HR_95L" in sumstats.columns:
        log.write(" -Flipping column: HR_95U = 1 / HR_95L...", verbose=verbose) 
        sumstats.loc[matched_index,"HR_95U"] =   factor / sumstats.loc[matched_index,"HR_95L"].values
    if "HR_95U" in sumstats.columns:
        log.write(" -Flipping column: HR_95L = 1 / HR_95U...", verbose=verbose) 
        sumstats.loc[matched_index,"HR_95L"] =   factor / sumstats.loc[matched_index,"HR_95U"].values
    return sumstats

def flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1):
    if "EAF" in sumstats.columns:
        log.write(" -Flipping column: EAF = 1 - EAF...", verbose=verbose) 
        sumstats.loc[matched_index,"EAF"] =   factor - sumstats.loc[matched_index,"EAF"].values
    return sumstats

def flip_by_sign(sumstats, matched_index, log, verbose, cols=None):
    if "BETA" in sumstats.columns:
        log.write(" -Flipping column: BETA = - BETA...", verbose=verbose) 
        sumstats.loc[matched_index,"BETA"] =     - sumstats.loc[matched_index,"BETA"].values
    if "BETA_95L" in sumstats.columns:
        log.write(" -Flipping column: BETA_95U = - BETA_95L...", verbose=verbose) 
        sumstats.loc[matched_index,"BETA_95U"] =     - sumstats.loc[matched_index,"BETA_95L"].values
    if "BETA_95U" in sumstats.columns:
        log.write(" -Flipping column: BETA_95L = - BETA_95U...", verbose=verbose) 
        sumstats.loc[matched_index,"BETA_95L"] =     - sumstats.loc[matched_index,"BETA_95U"].values
    if "Z" in sumstats.columns:
        log.write(" -Flipping column: Z = - Z...", verbose=verbose) 
        sumstats.loc[matched_index,"Z"] =     - sumstats.loc[matched_index,"Z"].values
    if "T" in sumstats.columns:
        log.write(" -Flipping column: T = - T...", verbose=verbose) 
        sumstats.loc[matched_index,"Z"] =     - sumstats.loc[matched_index,"T"].values
    if "DIRECTION" in sumstats.columns:
        log.write(" -Flipping column: DIRECTION +-?0 <=> -+?0 ...", verbose=verbose) 
        sumstats.loc[matched_index,"DIRECTION"] =   sumstats.loc[matched_index,"DIRECTION"].apply(flip_direction)
    return sumstats

def flipallelestats(sumstats,status="STATUS",verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "adjust statistics based on STATUS code"
    _end_line = "adjusting statistics based on STATUS code"
    _start_cols =[]
    _start_function = ".flip_allele_stats()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    if_stats_flipped = False
    ###################get reverse complementary####################
    pattern = r"\w\w\w\w\w[45]\w"  
    #matched_index = status_match(sumstats[status],6,[4,5]) #
    matched_index = sumstats[status].str[5].str.match(r"4|5")
    if sum(matched_index)>0:
        log.write("Start to convert alleles to reverse complement for SNPs with status xxxxx[45]x...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(sum(matched_index)) +" variants...", verbose=verbose) 
        if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns) :
            log.write(" -Converting to reverse complement : EA and NEA...", verbose=verbose) 
            reverse_complement_nea = sumstats.loc[matched_index,'NEA'].apply(lambda x :get_reverse_complementary_allele(x)) 
            reverse_complement_ea = sumstats.loc[matched_index,'EA'].apply(lambda x :get_reverse_complementary_allele(x)) 
            categories = set(sumstats['EA'])|set(sumstats['NEA']) |set(reverse_complement_ea) |set(reverse_complement_nea)
            sumstats['EA']=pd.Categorical(sumstats['EA'],categories = categories) 
            sumstats['NEA']=pd.Categorical(sumstats['NEA'],categories = categories ) 
            sumstats.loc[matched_index,['NEA']] = reverse_complement_nea
            sumstats.loc[matched_index,['EA']] = reverse_complement_ea
            sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 6, "4","2")
            log.write(" -Changed the status for flipped variants : xxxxx4x -> xxxxx2x", verbose=verbose)
        if_stats_flipped = True
    ###################flip ref####################
    pattern = r"\w\w\w\w\w[35]\w"  
    #matched_index = status_match(sumstats[status],6,[3,5]) #sumstats[status].str.match(pattern)
    matched_index = sumstats[status].str[5].str.match(r"3|5")
    if sum(matched_index)>0:
        log.write("Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(sum(matched_index)) +" variants...", verbose=verbose) 
        
        flip_by_swap(sumstats, matched_index, log, verbose)
        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants : xxxxx[35]x -> xxxxx[12]x", verbose=verbose) 
        sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 6, "35","12")
        if_stats_flipped = True
        
    ###################flip ref for undistingushable indels####################
    pattern = r"\w\w\w\w[123][67]6"  
    #matched_index = status_match(sumstats[status],6,[1,2,3])|status_match(sumstats[status],6,[6,7])|status_match(sumstats[status],7,6) #sumstats[status].str.match(pattern)
    matched_index = sumstats[status].str[4:].str.match(r"[123][67]6")
    if sum(matched_index)>0:
        log.write("Start to flip allele-specific stats for standardized indels with status xxxx[123][67][6]: ALT->EA , REF->NEA...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(sum(matched_index)) +" variants...", verbose=verbose) 
        
        flip_by_swap(sumstats, matched_index, log, verbose)
        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4", verbose=verbose) 
        sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 7, "6","4")
        if_stats_flipped = True
         # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    pattern = r"\w\w\w\w\w[012]5"  
    #matched_index = status_match(sumstats[status],6,[0,1,2]) | status_match(sumstats[status],7,[5])#sumstats[status].str.match(pattern)
    matched_index = sumstats[status].str[5:].str.match(r"05|15|25")
    if sum(matched_index)>0:
        log.write("Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...{}".format(_get_version()), verbose=verbose) 
        log.write(" -Flipping "+ str(sum(matched_index)) +" variants...", verbose=verbose) 

        flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
        flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
        flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)
        
        #change status    
        log.write(" -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2", verbose=verbose) 
        sumstats.loc[matched_index,status] = vchange_status(sumstats.loc[matched_index,status], 7, "5","2")
        if_stats_flipped = True

    if if_stats_flipped != True:
        log.write(" -No statistics have been changed.")

    finished(log, verbose, _end_line)
    return sumstats


###############################################################################################################
# 20220426
def liftover_snv(row,chrom,converter,to_build):
    status_pre=""
    status_end=row.iloc[1][2]+"9"+row.iloc[1][4]+"99"  
    pos_0_based = int(row.iloc[0]) - 1
    results = converter[chrom][pos_0_based]
    if converter[chrom][pos_0_based]:
        # return chrom, pos_1_based
        if results[0][0].strip("chr")!=chrom:
            return pd.NA,pd.NA,"97"+status_end
        else:
            return results[0][0].strip("chr"),results[0][1]+1,to_build+status_end
    else:
        return pd.NA,pd.NA,"97"+status_end

def liftover_variant(sumstats, 
             chrom="CHR", 
             pos="POS", 
             status="STATUS",
             from_build="19", 
             to_build="38",
             chain=None):
    
    try:
        if chain is None:
            converter = get_lifter("hg{}".format(from_build),"hg{}".format(to_build),one_based=True)
        else:
            converter = ChainFile(chain,target="",query="", one_based=True)
    except:
        if chain is None:
            converter = get_lifter("hg{}".format(from_build),"hg{}".format(to_build))
        else:
            converter = ChainFile(chain, target="",query="")

    dic= get_number_to_chr(in_chr=False,xymt=["X","Y","M"])
    dic2= get_chr_to_number(out_chr=False)
    for i in sumstats[chrom].unique():
        chrom_to_convert = dic[i]
        variants_on_chrom_to_convert = sumstats[chrom]==i
        lifted = sumstats.loc[variants_on_chrom_to_convert,[pos,status]].apply(lambda x: liftover_snv(x[[pos,status]],chrom_to_convert,converter,to_build),axis=1)
        sumstats.loc[variants_on_chrom_to_convert,pos]     =   lifted.str[1]
        sumstats.loc[variants_on_chrom_to_convert,status]   =  lifted.str[2]
        sumstats.loc[variants_on_chrom_to_convert,chrom]    =  lifted.str[0].map(dic2).astype("Int64")
    return sumstats

def parallelizeliftovervariant(sumstats,n_cores=1,chrom="CHR", pos="POS", from_build="19", to_build="38",status="STATUS",remove=True,chain=None, verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "perform liftover"
    _end_line = "liftover"
    _start_cols =[chrom,pos,status]
    _start_function = ".liftover()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            n_cores=n_cores,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################
    
    lifter_from_build = _process_build(from_build,log=log,verbose=False)
    lifter_to_build = _process_build(to_build,log=log,verbose=False)

    if chain is not None:
        log.write(" -Creating converter using ChainFile: {}".format(chain), verbose=verbose)
    else:
        try:
            chain = get_chain(from_build=from_build, to_build=to_build)
            if chain is None or chain==False:
                raise ValueError("")
            log.write(" -Creating converter using provided ChainFile: {}".format(chain), verbose=verbose)
        except:
            chain = None
            lifter_from_build=from_build
            lifter_to_build=to_build
            log.write(" -Try creating converter using liftover package", verbose=verbose)

    log.write(" -Creating converter : {} -> {}".format(lifter_from_build, lifter_to_build), verbose=verbose)
    # valid chr and pos
    pattern = r"\w\w\w0\w\w\w"  
    to_lift = sumstats[status].str.match(pattern)
    sumstats = sumstats.loc[to_lift,:].copy()
    log.write(" -Converting variants with status code xxx0xxx :"+str(len(sumstats))+"...", verbose=verbose)
    ###########################################################################
    if sum(to_lift)>0:
        if sum(to_lift)<10000:
            n_cores=1
    
        #df_split = np.array_split(sumstats[[chrom,pos,status]], n_cores)
        df_split = _df_split(sumstats[[chrom,pos,status]], n_cores)
        pool = Pool(n_cores)
        #df = pd.concat(pool.starmap(func, df_split))
        func=liftover_variant
        sumstats[[chrom,pos,status]] = pd.concat(pool.map(partial(func,chrom=chrom,pos=pos,from_build=from_build,to_build=to_build,status=status,chain=chain),df_split))
        pool.close()
        pool.join()
    ############################################################################
    unmap_num = len(sumstats.loc[sumstats[pos].isna(),:])    
    
    if remove is True:
        log.write(" -Removed unmapped variants: "+str(unmap_num), verbose=verbose)
        sumstats = sumstats.loc[~sumstats[pos].isna(),:]
    
    # after liftover check chr and pos
    sumstats = fixchr(sumstats,chrom=chrom,add_prefix="",remove=remove, verbose=True)
    sumstats = fixpos(sumstats,pos=pos,remove=remove, verbose=True)
    
    finished(log,verbose,_end_line)
    return sumstats

###############################################################################################################
# 20220426
def sortcoordinate(sumstats,chrom="CHR",pos="POS",reindex=True,verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "sort the genome coordinates"
    _end_line = "sorting coordinates"
    _start_cols =[chrom,pos]
    _start_function = ".sort_coordinate()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################
    
    try:
        if sumstats[pos].dtype == "Int64":
            pass
        else:
            log.write(" -Force converting POS to Int64...", verbose=verbose)
            sumstats[pos]  = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    except:
        pass
    sumstats = sumstats.sort_values(by=[chrom,pos],ascending=True,ignore_index=True)
    
    finished(log,verbose,_end_line)
    return sumstats
###############################################################################################################
# 20230430 added HR HR_95 BETA_95 N_CASE N_CONTROL
def sortcolumn(sumstats,verbose=True,log=Log(),order = None):
    ##start function with col checking##########################################################
    _start_line = "reorder the columns"
    _end_line = "reordering the columns"
    _start_cols =[]
    _start_function = ".sort_column()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    if order is None:
        order = [
        "SNPID","rsID", "CHR", "POS", "EA", "NEA", "EAF", "MAF", "BETA", "SE","BETA_95L","BETA_95U", "Z","T","F",
        "CHISQ", "P", "MLOG10P", "OR", "OR_95L", "OR_95U","HR", "HR_95L", "HR_95U","INFO", "N","N_CASE","N_CONTROL","DIRECTION","I2","P_HET","DOF","SNPR2","STATUS"]
    output_columns = []
    for i in order:
        if i in sumstats.columns: output_columns.append(i)
    for i in sumstats.columns:
        if i not in order: output_columns.append(i)
    log.write(" -Reordering columns to    :", ",".join(output_columns), verbose=verbose)
    sumstats = sumstats[ output_columns]

    finished(log,verbose,_end_line)
    return sumstats


###############################################################################################################
def start_to(sumstats, 
             log, 
             verbose, 
             start_line,
             end_line,
             start_cols,
             start_function,
             ref_vcf=None,
             ref_fasta=None,
             n_cores=None,
             ref_tsv=None,
             **kwargs
             ):
    
    log.write("Start to {}...{}".format(start_line,_get_version()), verbose=verbose)
    
    if sumstats is not None:
        check_dataframe_shape(sumstats=sumstats, 
                            log=log, 
                            verbose=verbose)  
        is_enough_col = check_col(sumstats.columns, 
                                verbose=verbose, 
                                log=log, 
                                cols=start_cols, 
                                function=start_function)
        if is_enough_col==True:
            if n_cores is not None:
                log.write(" -Number of threads/cores to use: {}".format(n_cores))
            if ref_vcf is not None:
                log.write(" -Reference VCF: {}".format(ref_vcf))
            if ref_fasta is not None:
                log.write(" -Reference FASTA: {}".format(ref_fasta))
            if ref_tsv is not None:
                log.write(" -Reference TSV: {}".format(ref_tsv))
            
            is_args_valid = True
            for key, value in kwargs.items():
                is_args_valid = is_args_valid & check_arg(log, verbose, key, value, start_function)
            is_enough_col = is_args_valid & is_enough_col

        if  is_enough_col == False:
            skipped(log, verbose, end_line)
        return is_enough_col
    else:
        return None

def finished(log, verbose, end_line):
    log.write("Finished {}.".format(end_line), verbose=verbose)
    gc.collect()

def skipped(log, verbose, end_line):
    log.write("Skipped {}.".format(end_line), verbose=verbose)
    gc.collect()

def check_arg(log, verbose, key, value, function):
    if value is None:
        log.warning("Necessary argument {} for {} is not provided!".format(key, function))
        return False
    return True

def check_col(df_col_names, verbose=True, log=Log(), cols=None, function=None):
    not_in_df=[]
    for i in cols:
        if type(i) is str:
            # single check
            if i in df_col_names:
                continue
            else:
                not_in_df.append(i)
        else:
            # paried check
            count=0
            for j in i:
                if j not in df_col_names:
                    not_in_df.append(j)
                    count+=1

    if len(not_in_df)>0:
        if function is None:
            to_show_title=" "
        else:
            to_show_title = " for {} ".format(function)
        log.warning("Necessary columns{}were not detected:{}".format(to_show_title, ",".join(not_in_df)))
        skipped(log, verbose, end_line=function)
        return False
    
    return True

###############################################################################################################
def _df_split(dataframe, n):
    k, m = divmod(len(dataframe), n)
    return [dataframe.iloc[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]