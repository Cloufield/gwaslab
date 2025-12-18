import re
import gc
import pandas as pd
import numpy as np
from multiprocessing import  Pool
from liftover import get_lifter
from liftover import ChainFile
from functools import partial

from gwaslab.g_vchange_status import vchange_status
from gwaslab.g_Log import Log
from gwaslab.g_version import _get_version
from gwaslab.g_vchange_status import STATUS_CATEGORIES
from gwaslab.g_vchange_status import match_status

from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chain
from gwaslab.bd.bd_common_data import NA_STRINGS

from gwaslab.qc.qc_check_datatype import check_datatype
from gwaslab.qc.qc_pattern import RSID_PATTERN, SNPID_PATTERN_EXTRACT, CHR_PATTERN_EXTRACT, FLAGS, SNPID_SEP_PATTERN, CHR_PREFIX_PATTERN, SNPID_PATTERN_STRIP
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.util.util_in_fill_data import fill_extreme_mlog10p

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
#flipallelestats
#sortcoordinate
#sortcolumn

@with_logging(
        start_to_msg= "check SNPID/rsID",
        finished_msg= "checking SNPID/rsID",
        start_function= ".fix_id()",
        start_cols=[]
)
def fixID(sumstats,
       snpid="SNPID",rsid="rsID",chrom="CHR",pos="POS",nea="NEA",ea="EA",status="STATUS",fixprefix=False,
       fixchrpos=False,fixid=False,fixeanea=False,fixeanea_flip=False,fixsep=False, reversea=False,
       overwrite=False,verbose=True,forcefixid=False,log=Log()):  
    
    '''
    Fix various aspects of genomic data including SNPID, rsID, chromosome positions, and allele information.
    
    This function performs multiple data quality checks and fixes: (1) validates and fixes SNPID format
    (CHR:POS:NEA:EA pattern), (2) extracts and fixes chromosome and position from SNPID or rsID, (3) validates
    rsID format (rsxxxxxx pattern), (4) extracts and fixes EA and NEA from SNPID, (5) standardizes separators
    and removes prefixes in SNPID, and (6) generates new SNPID from available CHR, POS, EA, NEA data.
    
    Parameters
    ----------
    fixprefix : bool
        Whether to remove 'chr' prefix in SNPID.
    fixchrpos : bool
        Whether to fix chromosome and position from SNPID.
    fixid : bool
        Whether to generate new SNPID from available data.
    fixeanea : bool
        Whether to fix EA and NEA from SNPID.
    fixeanea_flip : bool
        Whether to flip EA and NEA during fixing.
    fixsep : bool
        Whether to standardize separators in SNPID.
    reversea : bool
        Whether to reverse alleles in SNPID.
    overwrite : bool
        Whether to overwrite existing values.
    verbose : bool, optional
        Whether to print progress.
    forcefixid : bool
        Whether to force fix even without status check.

    Returns
    -------
    pd.DataFrame
        Modified sumstats with fixed data.
    '''

    ############################  checking datatype ###################################################  
    if rsid in sumstats.columns:
        log.write(" -Checking rsID data type...", verbose=verbose)
        if sumstats[rsid].dtype != "string":
            log.write(" -Converting rsID to pd.string data type...", verbose=verbose)
            sumstats[rsid] = sumstats[rsid].astype("string")

    if snpid in sumstats.columns:
        log.write(" -Checking SNPID data type...", verbose=verbose)
        if sumstats[snpid].dtype != "string":
            log.write(" -Converting SNPID to pd.string data type...", verbose=verbose)
            sumstats[snpid] = sumstats[snpid].astype("string")
            
    ############################  checking string NA ###################################################
    log.write(" -Checking NA strings :{}".format(",".join(NA_STRINGS)),verbose=verbose)  
    if snpid in sumstats.columns:  
        log.write(" -Checking if SNPID contains NA strings...",verbose=verbose)
        is_snpid_string_na = sumstats[snpid].isin(NA_STRINGS)
        if is_snpid_string_na.sum() > 0:
            log.write("  -Converting {} NA strings in SNPID to pd.NA...".format(is_snpid_string_na.sum()), verbose=verbose)
            sumstats.loc[is_snpid_string_na ,snpid] = pd.NA

    if rsid in sumstats.columns: 
        log.write(" -Checking if rsID contains NA strings...",verbose=verbose)
        is_rsid_string_na = sumstats[rsid].isin(NA_STRINGS)
        if is_rsid_string_na.sum() > 0:
            log.write("  -Converting {} NA strings in rsID to pd.NA...".format(is_rsid_string_na.sum()), verbose=verbose)
            sumstats.loc[is_rsid_string_na ,rsid] = pd.NA
    ############################  checking ###################################################  
    if snpid in sumstats.columns:  
        log.write(" -Checking if SNPID is CHR:POS:NEA:EA...(separator: - ,: , _)",verbose=verbose)
        # check if SNPID is CHR:POS:EA:NEA
        is_chrposrefalt = sumstats[snpid].str.match(SNPID_PATTERN_EXTRACT, flags=FLAGS, na=False)
        # check if SNPID is NA
        is_snpid_na = sumstats[snpid].isna()

        # change STATUS code
        sumstats.loc[ is_chrposrefalt,status] = vchange_status(sumstats.loc[ is_chrposrefalt,status],3 ,"975" ,"630")
        sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status] = vchange_status(sumstats.loc[(~is_chrposrefalt)&(~is_snpid_na),status],3 ,"975" ,"842")
        
    if rsid in sumstats.columns: 
        log.write(" -Checking if rsID is rsxxxxxx...", verbose=verbose)
        is_rsid = sumstats[rsid].str.match(RSID_PATTERN, flags=FLAGS, na=False)
        
        sumstats.loc[ is_rsid,status] = vchange_status(sumstats.loc[ is_rsid,status], 3, "986","520")
        sumstats.loc[~is_rsid,status] = vchange_status(sumstats.loc[~is_rsid,status], 3, "986","743")
        
        log.write(" -Checking if CHR:POS:NEA:EA is mixed in rsID column ...", verbose=verbose)
        is_rs_chrpos = sumstats[rsid].str.match(SNPID_PATTERN_EXTRACT, flags=FLAGS, na=False)
        
        log.write(" -Number of CHR:POS:NEA:EA mixed in rsID column :", is_rs_chrpos.sum(), verbose=verbose)
        log.write(" -Number of Unrecognized rsID :", len(sumstats) - is_rs_chrpos.sum() - is_rsid.sum(), verbose=verbose)
        log.write(" -A look at the unrecognized rsID :",set(sumstats.loc[(~is_rsid)&(~is_rs_chrpos),rsid].head()),"...", verbose=verbose) 
      
    ############################  fixing chr pos###################################################  
    if reversea == True:
        if snpid in sumstats.columns: 
            log.write(" -Reversing Alleles in SNPID...", verbose=verbose)
            to_fix = is_chrposrefalt
            to_fix_num = to_fix.sum()
            if to_fix_num>0 and verbose: log.write(" -Number of variants could be reversed: "+str(to_fix_num)+" ...")
            extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
            sumstats.loc[to_fix, snpid] = extracted["CHR"] + ":" + extracted["POS"] + ":" + extracted["EA"] + ":" + extracted["NEA"]
            
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
                to_fix_num = to_fix.sum()
                if to_fix_num and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)
            
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                log.write(" -Initiating CHR columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna()
                to_fix_num = to_fix.sum()
                if to_fix_num>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)
            
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                log.write(" -Initiating CHR and POS column...", verbose=verbose)
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt & sumstats[chrom].isna() & sumstats[pos].isna() 
                to_fix_num = to_fix.sum()
                if to_fix_num>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)     
                
            else:
                log.write(" -Initiating CHR and POS columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_chrposrefalt
                to_fix_num = to_fix.sum()
                if to_fix_num>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix_num)+" ...")
                else:
                    log.write(" -No fixable variants. ...", verbose=verbose)   
                    
            if to_fix.sum()>0:
                log.write(" -Filling CHR and POS columns using valid SNPID's chr:pos...", verbose=verbose)
                # format and qc filled chr and pos
                extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
                sumstats.loc[to_fix, chrom] = extracted["CHR"]
                sumstats.loc[to_fix, pos] = extracted["POS"]
                
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
                if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)
            elif (chrom not in sumstats.columns) and (pos in sumstats.columns):
                log.write(" -Initiating CHR columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna()
                if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)
            elif (chrom in sumstats.columns) and (pos not in sumstats.columns):
                log.write(" -Initiating CHR and POS column...", verbose=verbose)
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos & sumstats[chrom].isna() & sumstats[pos].isna() 
                if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)
            else:
                log.write(" -Initiating CHR and POS columns...", verbose=verbose)
                sumstats[chrom]=pd.Series(dtype="string")   
                sumstats[pos]=pd.Series(dtype="Int64") 
                to_fix = is_rs_chrpos
                if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
                else:
                    log.write(" -No fixable variants ...", verbose=verbose)   
            
            if to_fix.sum()>0:

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
            if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
        elif (nea in sumstats.columns) and (ea not in sumstats.columns):
            log.write(" -Initiating EA columns...", verbose=verbose)
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
        elif (nea not in sumstats.columns) and (ea in sumstats.columns):
            log.write(" -Initiating NEA columns...", verbose=verbose)
            sumstats[nea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt&(sumstats[nea].isna()|sumstats[ea].isna())
            if to_fix.sum()>0 and verbose: log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...")
        else:
            log.write(" -Initiating EA and NEA columns...", verbose=verbose)
            sumstats[nea]=pd.Series(dtype="string")
            sumstats[ea]=pd.Series(dtype="string")
            to_fix = is_chrposrefalt
            if to_fix.sum()>0:
                log.write(" -Number of variants could be fixed: "+str(to_fix.sum())+" ...", verbose=verbose)
    #                
        if to_fix.sum()>0:
            log.write(" -Filling "+str(to_fix.sum())+" EA and NEA columns using SNPID's CHR:POS:NEA:EA...", verbose=verbose)
    #        
            if fixeanea_flip == True:
                log.write(" -Flipped : CHR:POS:NEA:EA -> CHR:POS:EA:NEA ", verbose=verbose)
                extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
                sumstats.loc[to_fix, ea] = extracted["NEA"]
                sumstats.loc[to_fix, nea] = extracted["EA"]
            else:
                log.write(" -Chr:pos:a1:a2...a1->EA , a2->NEA ", verbose=verbose)
                extracted = sumstats.loc[to_fix, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
                sumstats.loc[to_fix, ea] = extracted["EA"]
                sumstats.loc[to_fix, nea] = extracted["NEA"]
    
    #        #to_change_status = sumstats[status].str.match(r"\w\w\w[45]\w\w\w")
    #        #sumstats.loc[to_fix&to_change_status,status] = vchange_status(sumstats.loc[to_fix&to_change_status,status],4,"2")  
    #        #sumstats.loc[to_fix,snpid].apply(lambda x:re.split(':|_|-',x)[1]).astype("string")
    #        #sumstats.loc[to_fix,rsid].apply(lambda x:re.split(':|_|-',x)[1]).astype("Int64")
    
    ############################  fixing id ################################################### 
    if fixsep == True:
        if snpid in sumstats.columns: 
            log.write(' -Replacing separators in SNPID with ":" ...', verbose=verbose)
            sumstats[snpid] = sumstats[snpid].str.replace(SNPID_SEP_PATTERN, ":", regex=True, flags=FLAGS)
    
    if fixprefix == True:
        if snpid in sumstats.columns: 
            log.write(' -Removing /^chr/ in SNPID ...', verbose=verbose)
            sumstats[snpid] = sumstats[snpid].str.replace(CHR_PREFIX_PATTERN, '', regex=True, flags=FLAGS)

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
                to_fix = sumstats[snpid].isna() & match_status(sumstats[status], r"\w\w\w[0]\w\w\w") 
                #status_match(sumstats[status],4,0)
                #
            else:
                #fix all
                to_fix = match_status(sumstats[status], r"\w\w\w[0]\w\w\w") 
                # status_match(sumstats[status],4,0)
                #
            
            if (ea in sumstats.columns) and (nea in sumstats.columns):
            # when ea and nea is available  -> check status -> fix to chr:pos:nea:ea 
                
                pattern = r"\w\w\w[0]\w\w\w"  
                matched_index = match_status(sumstats[status], pattern)
                #matched_index = status_match(sumstats[status],4,0) #
                to_part_fix = matched_index & to_fix 
                
                #pattern = r"\w\w\w[0][01267][01234]\w"  
                pattern = r"\w\w\w\w[0123][01267][01234]"
                matched_index = match_status(sumstats[status], pattern)
                #status_match(sumstats[status],5,[0,1,2,3]) | status_match(sumstats[status],6,[0,1,2,6,7])| status_match(sumstats[status],7,[0,1,2,3,4])
                #
                if forcefixid is True:
                    matched_index = to_fix
                to_full_fix = matched_index & to_fix 
                
                
                if to_part_fix.sum() > 0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                if to_full_fix.sum() > 0:
                    sumstats.loc[to_full_fix,snpid] = sumstats.loc[to_full_fix,chrom].astype("string") + ":"+sumstats.loc[to_full_fix,pos].astype("string") +":"+ sumstats.loc[to_full_fix,nea].astype("string") +":"+ sumstats.loc[to_full_fix,ea].astype("string")
                log.write(" -Filling "+str(to_part_fix.sum() - to_full_fix.sum()) +" SNPID using CHR:POS...", verbose=verbose)
                log.write(" -Filling "+str(to_full_fix.sum()) +" SNPID using CHR:POS:NEA:EA...", verbose=verbose)
                sumstats.loc[(to_full_fix),status] = vchange_status(sumstats.loc[(to_full_fix),status],3,"975","630") 
                sumstats.loc[(to_part_fix),status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")  
                
            else:
            #when these is no ea or ena, just fix to chr:pos
                to_part_fix = to_fix & sumstats[chrom].notnull() & sumstats[pos].notnull()
                log.write(" -Filling "+str(to_part_fix.sum()) +" SNPID using CHR POS...", verbose=verbose)
                if to_part_fix.sum() > 0:
                    sumstats.loc[to_part_fix,snpid] = sumstats.loc[to_part_fix,chrom].astype("string") + ":"+sumstats.loc[to_part_fix,pos].astype("string")
                    sumstats.loc[to_part_fix,status] = vchange_status(sumstats.loc[(to_part_fix),status],3,"975","842")
                    
            after_number = sumstats[snpid].isna().sum()
            log.write(" -Fixed "+ str(pre_number - after_number) +" variants ID...", verbose=verbose)
        else:
            log.write(" -ID unfixable: no CHR and POS columns or no SNPID. ", verbose=verbose)

    return sumstats

@with_logging(
        start_to_msg= "strip SNPID",
        finished_msg= "stripping SNPID",
        start_function= ".strip_snpid()",
        start_cols=["SNPID"]
)
def stripSNPID(sumstats,snpid="SNPID",overwrite=False,verbose=True,log=Log()):  
    '''
    Strip non-standard characters from SNPID values to standardize format.
    
    Removes leading and trailing non-standard characters from SNPID values that match
    the pattern (xxx:)CHR:POS:ATCG_Allele:ATCG_Allele(:xxx), keeping only the core
    CHR:POS:NEA:EA format.

    Parameters
    ----------
    overwrite : bool
        Whether to overwrite existing values.
    verbose : bool, optional
        Whether to print progress.
    Returns
    -------
    pd.DataFrame
        Modified sumstats with stripped SNPIDs.
    '''
    log.write(" -Checking if SNPID is (xxx:)CHR:POS:ATCG_Allele:ATCG_Allele(:xxx)...(separator: - ,: , _)",verbose=verbose)
    is_chrposrefalt = sumstats[snpid].str.match(SNPID_PATTERN_STRIP, flags=FLAGS, na=False)
    # check if SNPID is NA
    log.write(" -Stripping {} non-NA fixable SNPIDs...".format(sum(is_chrposrefalt)),verbose=verbose)

    # flip 
    extracted = sumstats.loc[is_chrposrefalt, snpid].str.extract(SNPID_PATTERN_STRIP, flags=FLAGS)
    sumstats.loc[is_chrposrefalt, snpid] = (
        extracted["CHR"].astype("string") + ":" +
        extracted["POS"].astype("string") + ":" +
        extracted["NEA"].astype("string") + ":" +
        extracted["EA"].astype("string")
    )

    return sumstats

@with_logging(
        start_to_msg= "flip SNPID from CHR:POS:A1:A2 to CHR:POS:A2:A1",
        finished_msg= "flipping SNPID",
        start_function= ".flip_snpid()",
        start_cols=["SNPID"]
)
def flipSNPID(sumstats,snpid="SNPID",overwrite=False,verbose=True,log=Log()):  
    '''
    Flip alleles in SNPID values without changing status codes or statistics.
    
    Converts SNPID from CHR:POS:EA:NEA format to CHR:POS:NEA:EA format by swapping
    the last two allele components. This function only modifies the SNPID column and
    does not affect EA, NEA, STATUS, or any statistics.

    Parameters
    ----------
    snpid : str
        Column name for SNPID.
    overwrite : bool
        Whether to overwrite existing values.
    verbose : bool, optional
        Whether to print progress.

    Returns
    -------
    pd.DataFrame
        Modified sumstats with flipped alleles.
    '''
    ##start function with col checking##########################################################
    log.warning("This function only flips alleles in SNPID without changing EA, NEA, STATUS or any statistics.")
    log.write(" -Checking if SNPID is CHR:POS:ATCG_Allele:ATCG_Allele...(separator: - ,: , _)",verbose=verbose)
    is_chrposrefalt = sumstats[snpid].str.match(SNPID_PATTERN_EXTRACT, flags=FLAGS, na=False)
    
    log.write(" -Flipping {} non-NA fixable SNPIDs...".format(sum(is_chrposrefalt)),verbose=verbose)

    # flip 
    parts = sumstats.loc[is_chrposrefalt, snpid].str.extract(SNPID_PATTERN_EXTRACT, flags=FLAGS)
    parts = parts.astype({"CHR":"string", "POS":"string", "NEA":"string", "EA":"string"})
    sumstats.loc[is_chrposrefalt, snpid] = parts["CHR"].str.cat([parts["POS"], parts["EA"], parts["NEA"]], sep=":")

    return sumstats

###############################################################################################################
# 20230128 
@with_logging(
        start_to_msg= "remove duplicated/multiallelic variants",
        finished_msg= "removing duplicated/multiallelic variants",
        start_function= ".remove_dup()",
        start_cols=None
)
def removedup(sumstats,mode="dm",chrom="CHR",pos="POS",snpid="SNPID",ea="EA",nea="NEA",rsid="rsID",keep='first',keep_col="P",remove_na=False,keep_ascend=True,verbose=True,log=Log()):
    """
    Remove duplicate or multiallelic variants based on user-selected criteria.

    Supports multiple duplicate-identification strategies depending on variant identifiers
    (SNPID, rsID) or allele and coordinate combinations. Can also collapse multi-allelic
    sites by retaining a single representative variant. Variants are sorted by a specified
    column (e.g., P-value) before removal to ensure the best variant is kept.

    Parameters
    ----------
    mode : str
        String encoding the deduplication rules; may include one or more of:
        - 'ds' : Identify duplicates using SNPID.
        - 'dr' : Identify duplicates using rsID.
        - 'dc' : Identify duplicates using chromosome, position, effect allele, and non-effect allele.
        - 'm' : Identify multi-allelic variants (same chromosome + position).
    keep : {'first', 'last', False}, default 'first'
        Which record to retain when duplicates are detected.
    keep_col : str, optional
        Column to sort by prior to duplicate removal; used only when `keep` is
        not False.
    remove_na : bool, default False
        If True, remove rows containing missing values in deduplication-relevant columns.
    keep_ascend : bool, default True
        If True, sort in ascending order when determining which duplicate to keep.
    verbose : bool, default False
        If True, print progress and summary information.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with duplicates and multi-allelic variants removed
        according to the specified mode.

    """

    log.write(" -Removing mode:{}".format(mode), verbose=verbose)
    if keep_col == "P" and ("P" in sumstats.columns):
        zero_count = (sumstats["P"] == 0).sum()
        if zero_count > 1:
            log.write(" -Detected {} variants with P=0; converting to MLOG10P and sorting descending...".format(zero_count), verbose=verbose)
            status, _ = fill_extreme_mlog10p(sumstats, None, log, verbose=verbose, filled_count=0)
            if status == 1 and ("MLOG10P" in sumstats.columns):
                keep_col = "MLOG10P"
                keep_ascend = False
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
    
    if "n" in mode or remove_na==True:
        # if remove==True, remove NAs
        log.write(" -Removing NAs...", verbose=verbose)
        pre_number =len(sumstats) 
        specified_columns = []
        if "d" in mode:
            if rsid in sumstats.columns: specified_columns.append(rsid)
            if snpid in sumstats.columns: specified_columns.append(snpid)
            if chrom in sumstats.columns: specified_columns.append(chrom)
            if pos in sumstats.columns: specified_columns.append(pos)
            if ea in sumstats.columns: specified_columns.append(ea)
            if nea in sumstats.columns: specified_columns.append(nea)
        if "r" in mode:
            if rsid in sumstats.columns:specified_columns.append(rsid)
        if "s" in mode:
            if snpid in sumstats.columns:specified_columns.append(snpid)
        if "m" in mode:
            if chrom in sumstats.columns:specified_columns.append(chrom)
            if pos in sumstats.columns:specified_columns.append(pos)
        if "c" in mode:
            if chrom in sumstats.columns:specified_columns.append(chrom)
            if pos in sumstats.columns:specified_columns.append(pos)
            if ea in sumstats.columns:specified_columns.append(ea)
            if nea in sumstats.columns:specified_columns.append(nea)
        specified_columns = list(set(specified_columns))
        sumstats = sumstats.loc[~sumstats[specified_columns].isna().any(axis=1),:]
        after_number=len(sumstats) 
        log.write(" -Removed ",pre_number -after_number," variants with NA values in {} .".format(specified_columns), verbose=verbose)  

    return sumstats

###############################################################################################################
# 20230128
@with_logging(
        start_to_msg= "fix chromosome notation (CHR)",
        finished_msg= "fixing chromosome notation (CHR)",
        start_function= ".fix_chr()",
        start_cols=["CHR","STATUS"]
)
def fixchr(sumstats,chrom="CHR",status="STATUS",add_prefix="",x=("X",23),y=("Y",24),mt=("MT",25), remove=False, verbose=True, chrom_list = None, minchr=1,log=Log()):
    """
    Standardize chromosome notation and handle special chromosome cases (X, Y, MT).
    
    All chromosome notations are converted to string type first. After fix, all chromosome
    notations will be int. This function normalizes chromosome labels to a consistent format,
    extracts chromosome numbers from various formats (e.g., "chr1", "1", "chrX"), maps special
    chromosomes (X, Y, mitochondrial) to standardized numeric identifiers, and optionally
    removes invalid chromosome values.

    Parameters
    ----------
    add_prefix : str, optional, default=""
        Prefix to prepend to chromosome labels (e.g., "chr").
    x : list of [str, int], optional
        Mapping for the X chromosome, given as [label, numeric_value]. Default is ["X",23]
    y : list of [str, int], optional
        Mapping for the Y chromosome, given as [label, numeric_value]. Default is ["Y",24]
    mt : list of [str, int], optional
        Mapping for the mitochondrial chromosome, given as [label, numeric_value]. Default is ["MY",25]
    remove : bool, default False
        If True, remove records with invalid or unrecognized chromosome labels.
    verbose : bool, default False
        If True, print progress or diagnostic messages.
    minchr : int, default 1
        Minimum allowed chromosome number when interpreting numeric labels.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table with standardized chromosome identifiers.

    Less used parameters
    -------------------------
    chrom_list : list of str or int, optional
        List of valid chromosome labels or numeric values. Default is ["1", ..., "25", "X", "Y", "MT"]
    """

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
    ########################################################################################################################################
    # check if CHR is numeric
    is_chr_fixed = sumstats[chrom].str.isnumeric()
    # fill NAs with False
    is_chr_fixed = is_chr_fixed.fillna(False)
    log.write(" -Variants with standardized chromosome notation:", is_chr_fixed.sum(), verbose=verbose)  
    
    # if there are variants whose CHR need to be fixed
    if is_chr_fixed.sum() < len(sumstats):
        
        #extract the CHR number or X Y M MT
        chr_extracted = sumstats.loc[~is_chr_fixed, chrom].str.extract(CHR_PATTERN_EXTRACT,  flags=FLAGS)[1]

        is_chr_fixable = ~chr_extracted.isna()
        log.write(" -Variants with fixable chromosome notations:", is_chr_fixable.sum(), verbose=verbose)  

        # For not fixed variants, check if na
        is_chr_na  = sumstats.loc[~is_chr_fixed, chrom].isna()
        if is_chr_na.sum()>0 and verbose: 
            log.write(" -Variants with NA chromosome notations:", is_chr_na.sum())  
        
        # Check variants with CHR being not NA and not fixable
        is_chr_invalid = (~is_chr_fixable)&(~is_chr_na)
        if is_chr_invalid.sum()>0 and verbose: 
            log.write(" -Variants with invalid chromosome notations:", is_chr_invalid.sum(), verbose=verbose) 
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
        if sex_chr.sum()>0:
            log.write(" -Identifying non-autosomal chromosomes : {}, {}, and {} ...".format(x[0],y[0],mt[0]), verbose=verbose)
            log.write(" -Identified ", str(sex_chr.sum()), " variants on sex chromosomes...", verbose=verbose)
            
            # convert "X, Y, MT" to numbers
            convert_num_to_xymt={}
            if sumstats[chrom].isin([x[0].lower(), x[0].upper()]).any():
                convert_num_to_xymt[x[0].lower()] = str(x[1])
                convert_num_to_xymt[x[0].upper()] = str(x[1])
                log.write(" -Standardizing sex chromosome notations: {} to {}...".format(x[0], x[1]), verbose=verbose)
            if sumstats[chrom].isin([y[0].lower(), y[0].upper()]).any():
                convert_num_to_xymt[y[0].lower()] = str(y[1])
                convert_num_to_xymt[y[0].upper()] = str(y[1])
                log.write(" -Standardizing sex chromosome notations: {} to {}...".format(y[0], y[1]), verbose=verbose)
            if sumstats[chrom].isin([mt[0].lower(), mt[0].upper()]).any():
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
        unrecognized_num = (~sumstats[chrom].isin(chrom_list)).sum()
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

    ########################################################################################################################################
    # Convert string to int
    try:
        sumstats[chrom] = sumstats[chrom].astype('Int64')
    except:
    #    # force convert
        sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    ########################################################################################################################################
    # filter out variants with CHR <=0
    out_of_range_chr = sumstats[chrom] < minchr
    out_of_range_chr = out_of_range_chr.fillna(False)
    if out_of_range_chr.sum()>0:
        log.write(" -Sanity check for CHR...", verbose=verbose) 
        log.write(" -Removed {} variants with CHR < {}...".format(out_of_range_chr.sum(), minchr), verbose=verbose)
        sumstats = sumstats.loc[~out_of_range_chr,:]

    return sumstats

###############################################################################################################    
# 20230128
@with_logging(
        start_to_msg= "fix basepair positions (POS)",
        finished_msg= "fixing basepair positions (POS)",
        start_function= ".fix_pos()",
        start_cols=["POS","STATUS"]
)
def fixpos(sumstats,pos="POS",status="STATUS",remove=False, verbose=True, lower_limit=0 , upper_limit=None , limit=250000000, log=Log()):
    '''
    Standardize and validate genomic base-pair positions.
    
    This function checks that reported genomic positions fall within valid chromosomal bounds
    and optionally removes invalid entries. It handles string-formatted positions with thousands
    separators, converts positions to Int64 type, and filters out positions outside the specified
    range. If explicit limits are not provided, a default maximum bound is applied.

    Parameters
    ----------
    remove : bool, default False
        If True, remove records with invalid or out-of-range positions.
    verbose : bool, default False
        If True, print progress or diagnostic messages.
    lower_limit : int, optional
        Minimum acceptable genomic position. Deafult is 0.
    upper_limit : int, optional
        Maximum acceptable genomic position.
    limit : int, default 3_000_000_000
        Default upper limit applied when `upper_limit` is not provided.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with standardized and validated base-pair positions.
    '''

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
 
    return sumstats

###############################################################################################################    
# 20220514
@with_logging(
        start_to_msg= "fix alleles (EA and NEA)",
        finished_msg= "fixing alleles (EA and NEA)",
        start_function= ".fix_allele()",
        start_cols=["EA","NEA","STATUS"]
)
def fixallele(sumstats,ea="EA", nea="NEA",status="STATUS",remove=False,verbose=True,log=Log()):
    """
    Validate and standardize allele representations.
    
    This function checks allele fields for valid nucleotide characters (A, T, C, G), converts
    all alleles to uppercase, standardizes their format using categorical data types, and
    classifies variants as SNPs, indels, normalized, or not normalized. Optionally, rows with
    invalid allele values can be removed.

    Parameters
    ----------
    remove : bool, default False
        If True, remove variants with invalid allele representations.
    verbose : bool, default False
        If True, print progress or warning messages.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table with validated and standardized allele values.
    """

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

    log.write(" -Variants with bad EA  : {}".format(bad_ea.sum()), verbose=verbose)
    log.write(" -Variants with bad NEA : {}".format(bad_nea.sum()), verbose=verbose)
    
    ## check NA
    is_eanea_na = sumstats[ea].isna() |  sumstats[nea].isna()
    log.write(" -Variants with NA for EA or NEA: {}".format(is_eanea_na.sum()), verbose=verbose)
    
    ## check same alleles
    not_variant = sumstats[nea] == sumstats[ea]
    log.write(" -Variants with same EA and NEA: {}".format(not_variant.sum()), verbose=verbose)

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
    
    ea_len = sumstats[ea].str.len()
    nea_len = sumstats[nea].str.len()
    is_eanea_fixed = good_ea | good_nea
    is_snp = (ea_len==1) & (nea_len==1)
    is_indel = (ea_len!=nea_len)
    is_not_normalized = (ea_len>1) & (nea_len>1)
    is_normalized = (ea_len==1) ^ (nea_len==1)
    
    #is_snp = (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()==1)
    #is_indel = (sumstats[ea].str.len()!=sumstats[nea].str.len())
    #is_not_normalized = (sumstats[ea].str.len()>1) &(sumstats[nea].str.len()>1)
    #is_normalized = is_indel &( (sumstats[ea].str.len()==1) &(sumstats[nea].str.len()>1) | (sumstats[ea].str.len()>1) &(sumstats[nea].str.len()==1) )

    if is_invalid.sum()>0:
        sumstats.loc[is_invalid, status]                      = vchange_status(sumstats.loc[is_invalid,status],                5,"9","6") 
    if is_eanea_na.sum()>0:
        sumstats.loc[is_eanea_na,status]                      = vchange_status(sumstats.loc[is_eanea_na, status],              5,"9","7")
    if (is_eanea_fixed & is_not_normalized).sum()>0:
        sumstats.loc[is_eanea_fixed&is_not_normalized,status] = vchange_status(sumstats.loc[is_eanea_fixed&is_not_normalized,status], 5,"9","5")
    if (is_eanea_fixed & is_snp).sum()>0:
        sumstats.loc[is_eanea_fixed&is_snp, status]           = vchange_status(sumstats.loc[is_eanea_fixed&is_snp,status],        5,"9","0")
    if (is_eanea_fixed & is_indel).sum()>0:
        sumstats.loc[is_eanea_fixed&is_indel,status]          = vchange_status(sumstats.loc[is_eanea_fixed&is_indel, status],      5,"9","4")
    if (is_eanea_fixed & is_normalized).sum()>0:
        sumstats.loc[is_eanea_fixed&is_normalized,status]     = vchange_status(sumstats.loc[is_eanea_fixed&is_normalized, status],  5,"4","3")

    return sumstats

###############################################################################################################   
# 20220721

@with_logging(
        start_to_msg= "normalize indels",
        finished_msg= "normalizing indels",
        start_function= ".normalize()",
        start_cols=["EA","NEA","STATUS"]
)
def parallelnormalizeallele(sumstats,mode="s",snpid="SNPID",rsid="rsID",pos="POS",nea="NEA",ea="EA" ,status="STATUS",chunk=3000000,n_cores=1,verbose=True,log=Log()):
    '''
    Normalize indels in parallel using left-alignment and parsimony principles.
    
    This function standardizes allele representations for insertion/deletion variants by
    left-aligning and trimming shared sequence context. It removes common suffixes and
    prefixes from both alleles and adjusts positions accordingly, following the VCF
    normalization standard.

    Parameters
    ----------
    chunk : int, default 10000
        Size of chunks for parallel processing.
    n_cores : int, default 1
        Number of CPU cores used for parallel processing.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with normalized indel allele representations.
    '''
    ############################################################################################

    #variants_to_check = status_match(sumstats[status],5,[4,5]) #
    #r'\w\w\w\w[45]\w\w'
    variants_to_check = match_status(sumstats[status], r"\w\w\w\w[45]\w\w")
    if variants_to_check.sum() == 0:
        log.write(" -No available variants to normalize..", verbose=verbose)
        return sumstats
    ###############################################################################################################
    if mode=="v":
        if variants_to_check.sum() < 100000:
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
        if variants_to_check.sum() < 10000:
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
        if is_pop.sum() == 0:
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
        if is_pop.sum() == 0:
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

@with_logging(
        start_to_msg= "adjust statistics based on STATUS code",
        finished_msg= "adjusting statistics based on STATUS code",
        start_function= ".flip_allele_stats()",
        start_cols=None
)
def flipallelestats(sumstats,status="STATUS",verbose=True,log=Log()):
    '''
    Adjust statistics when allele direction has changed based on STATUS codes.
    
    This function adjusts effect sizes and allele-specific statistics when variants have been
    flipped or converted to reverse complement. It handles multiple scenarios: reverse
    complement conversion for SNPs, allele swapping for REF/ALT mismatches, flipping for
    standardized indels, and strand flipping for palindromic variants. Run after checking
    with reference sequence.

    Parameters
    ----------
    verbose : bool, default False
        If True, print progress messages during processing.

    Returns
    -------
    pandas.DataFrame
        Summary statistics with effect sizes and alleles flipped where required.
    '''

    if_stats_flipped = False
    ###################get reverse complementary####################
    pattern = r"\w\w\w\w\w[45]\w"  
    #matched_index = status_match(sumstats[status],6,[4,5]) #
    matched_index = sumstats[status].str[5].str.match(r"4|5")
    if matched_index.sum() > 0:
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
    matched_index = match_status(sumstats[status], r"\w\w\w\w\w[35]\w")
    if matched_index.sum() > 0:
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
    matched_index = match_status(sumstats[status], r"\w\w\w\w[123][67]6")
    if matched_index.sum() > 0:
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
    matched_index = match_status(sumstats[status], r"\w\w\w\w\w[012]5")
    if matched_index.sum() > 0:
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

@with_logging(
        start_to_msg= "perform liftover",
        finished_msg= "liftover",
        start_function= ".liftover()",
        start_cols=["CHR","POS","STATUS"]
)
def parallelizeliftovervariant(sumstats,n_cores=1,chrom="CHR", pos="POS", from_build="19", to_build="38",status="STATUS",remove=True,chain=None, verbose=True,log=Log()):
    '''
    Perform parallelized liftover of variants to a new genome build.
    
    Converts genomic coordinates from one genome build (e.g., hg19/GRCh37) to another
    (e.g., hg38/GRCh38) using UCSC chain files. Only processes variants with valid
    coordinates (status code xxx0xxx). After liftover, validates and fixes chromosome
    and position columns.

    Parameters
    ----------
    n_cores : int
        Number of CPU cores to use for parallel processing.
    from_build : str
        Name of the original genome build (e.g., "19").
    to_build : str
        Name of the target genome build (e.g., "38").
    remove : bool, default False
        If True, remove variants that fail to map.
    verbose : bool, default False
        If True, print progress messages during processing.

    Returns
    -------
    pandas.DataFrame
        Summary statistics table with updated genomic coordinates.
    '''
    
    lifter_from_build = _process_build(from_build,log=log,verbose=False)
    lifter_to_build = _process_build(to_build,log=log,verbose=False)

    if chain is not None:
        log.write(" -Creating converter using ChainFile: {}".format(chain), verbose=verbose)
    else:
        try:
            chain = get_chain(from_build=from_build, to_build=to_build)
            if chain is None or chain==False:
                raise ValueError("No available chain file for {} -> {}".format(lifter_from_build, lifter_to_build))
            log.write(" -Creating converter using provided ChainFile: {}".format(chain), verbose=verbose)
        except:
            chain = None
            log.write(" -Try creating converter using liftover package", verbose=verbose)

    log.write(" -Creating converter : {} -> {}".format(lifter_from_build, lifter_to_build), verbose=verbose)
    # valid chr and pos
    pattern = r"\w\w\w0\w\w\w"  
    to_lift = match_status(sumstats[status], pattern)
    sumstats = sumstats.loc[to_lift,:].copy()
    log.write(" -Converting variants with status code xxx0xxx :"+str(len(sumstats))+"...", verbose=verbose)
    ###########################################################################
    if to_lift.sum() > 0:
        if to_lift.sum() < 10000:
            n_cores=1
    
        #df_split = np.array_split(sumstats[[chrom,pos,status]], n_cores)
        df_split = _df_split(sumstats[[chrom,pos,status]], n_cores)
        pool = Pool(n_cores)
        #df = pd.concat(pool.starmap(func, df_split))
        func=liftover_variant
        sumstats[[chrom,pos,status]] = pd.concat(pool.map(partial(func,chrom=chrom,pos=pos,from_build=lifter_from_build,to_build=lifter_to_build,status=status,chain=chain),df_split))
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

    return sumstats

###############################################################################################################
# 20220426
@with_logging(
        start_to_msg= "sort the genome coordinates",
        finished_msg= "sorting coordinates",
        start_function= ".sort_coordinate()",
        start_cols=["CHR","POS"],
        show_shape=False
)
def sortcoordinate(sumstats,chrom="CHR",pos="POS",reindex=True,verbose=True,log=Log()):
    '''
    Sort variants by genomic coordinates (chromosome, then position).
    
    Sorts the dataframe first by chromosome number, then by position in ascending order.
    The index is reset to sequential integers after sorting.

    Parameters
    ----------
    verbose : bool, default False
        If True, print progress messages.

    Returns
    -------
    pandas.DataFrame
        DataFrame with sorted genomic coordinates.
    '''
    
    try:
        if sumstats[pos].dtype == "Int64":
            pass
        else:
            log.write(" -Force converting POS to Int64...", verbose=verbose)
            sumstats[pos]  = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    except:
        pass
    sumstats = sumstats.sort_values(by=[chrom,pos],ascending=True,ignore_index=True)
    return sumstats
###############################################################################################################
# 20230430 added HR HR_95 BETA_95 N_CASE N_CONTROL
@with_logging(
        start_to_msg= "reorder the columns",
        finished_msg= "reordering the columns",
        start_function= ".sort_column()",
        start_cols=None,
        show_shape=False
)
def sortcolumn(sumstats,verbose=True,log=Log(),order = None):
    '''
    Reorder columns according to a specified order.
    
    Reorders the dataframe columns to match a predefined standard order, placing standard
    GWAS columns first (SNPID, rsID, CHR, POS, EA, NEA, statistics, etc.) followed by
    any additional columns not in the standard list.

    Parameters
    ----------
    verbose : bool, optional
        Whether to print progress. Default is True.

    Returns
    -------
    pd.DataFrame
        Modified sumstats with reordered columns.
    '''

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

    return sumstats


###############################################################################################################
def _df_split(dataframe, n):
    k, m = divmod(len(dataframe), n)
    return [dataframe.iloc[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]
