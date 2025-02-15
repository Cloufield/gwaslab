import pandas as pd
import numpy as np
from pysam import VariantFile
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool
from functools import partial
import re
import os
import gc
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import fixchr
from gwaslab.qc_fix_sumstats import fixpos
from gwaslab.qc_fix_sumstats import sortcolumn
from gwaslab.qc_fix_sumstats import _df_split
from gwaslab.qc_fix_sumstats import check_col
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.qc_fix_sumstats import skipped
from gwaslab.qc_check_datatype import check_dataframe_shape
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_chr_list
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_NC
from gwaslab.bd_common_data import _maketrans
from gwaslab.g_vchange_status import vchange_status
from gwaslab.g_version import _get_version
from gwaslab.cache_manager import CacheManager, PALINDROMIC_INDEL, NON_PALINDROMIC

#rsidtochrpos
#checkref
#parallelizeassignrsid
#inferstrand
#parallelecheckaf

### CONSTANTS AND MAPPINGS ###

PADDING_VALUE = 100

# chr(0) should not be used in the mapping dict because it's a reserved value.
# Instead of starting from chr(1), we start from chr(2) because this could be useful in the future
# to compute the complementary allele with a simple XOR operation (e.g. 2 ^ 1 = 3, 3 ^ 1 = 2, 4 ^ 1 = 5, 5 ^ 1 = 4, ...)
MAPPING = {
    "A": chr(2),
    "T": chr(3),
    "C": chr(4),
    "G": chr(5),
    "N": chr(6),
}
assert all(value != chr(0) for value in MAPPING.values()), "Mapping in the dictionary should not be equal to chr(0). This is a reserved value"

_COMPLEMENTARY_MAPPING = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
}
COMPLEMENTARY_MAPPING = {k: MAPPING[v] for k,v in _COMPLEMENTARY_MAPPING.items()}

TRANSLATE_TABLE = _maketrans(MAPPING)
TRANSLATE_TABLE_COMPL = _maketrans(COMPLEMENTARY_MAPPING)

#20220808
#################################################################################################################

###~!!!!
def rsidtochrpos(sumstats,
         path=None, ref_rsid_to_chrpos_tsv=None, snpid="SNPID",
         rsid="rsID", chrom="CHR",pos="POS",ref_rsid="rsID",ref_chr="CHR",ref_pos="POS", build="19",
              overwrite=False,remove=False,chunksize=5000000,verbose=True,log=Log()):
    '''
    assign chr:pos based on rsID
    '''
    ##start function with col checking##########################################################
    _start_line = "assign CHR and POS using rsIDs"
    _end_line = "assigning CHR and POS using rsIDs"
    _start_cols = [rsid]
    _start_function = ".rsid_to_chrpos()"
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

    log.write(" -rsID dictionary file: "+ path,verbose=verbose)  
    
    if ref_rsid_to_chrpos_tsv is not None:
        path = ref_rsid_to_chrpos_tsv

    if snpid in sumstats.columns and sum(sumstats[rsid].isna())>0:
        log.write(" -Filling na in rsID columns with SNPID...",verbose=verbose)  
        sumstats.loc[sumstats[rsid].isna(),rsid] = sumstats.loc[sumstats[rsid].isna(),snpid]
    
    if sum(sumstats[rsid].isna())>0:
        log.write(" -Filling na in rsID columns with NA_xxx for {} variants...".format(sum(sumstats[rsid].isna())),verbose=verbose)  
        sumstats.loc[sumstats[rsid].isna(),rsid] = ["NA_" + str(x+1) for x in range(len(sumstats.loc[sumstats[rsid].isna(),rsid]))]

    dic_chuncks = pd.read_csv(path,sep="\t",usecols=[ref_rsid,ref_chr,ref_pos],
                      chunksize=chunksize,index_col = ref_rsid,
                      dtype={ref_rsid:"string",ref_chr:"Int64",ref_pos:"Int64"})
    
    sumstats = sumstats.set_index(rsid)
    
    #if chr or pos columns not in sumstats
    if chrom not in sumstats.columns:
        sumstats[chrom] =pd.Series(dtype="Int64")
    if pos not in sumstats.columns:    
        sumstats[pos] =pd.Series(dtype="Int64")
    
    log.write(" -Setting block size: ",chunksize,verbose=verbose)
    log.write(" -Loading block: ",end="",verbose=verbose)     
    for i,dic in enumerate(dic_chuncks): 
        dic_to_update = dic[dic.index.notnull()]
        log.write(i," ",end=" ",show_time=False)  
        dic_to_update = dic_to_update.rename(index={ref_rsid:rsid})
        dic_to_update = dic_to_update.rename(columns={ref_chr:chrom,ref_pos:pos})  
        dic_to_update = dic_to_update[~dic_to_update.index.duplicated(keep='first')]
        sumstats.update(dic_to_update,overwrite="True")
        gc.collect()
    
    log.write("\n",end="",show_time=False,verbose=verbose) 
    sumstats = sumstats.reset_index()
    sumstats = sumstats.rename(columns = {'index':rsid})
    log.write(" -Updating CHR and POS finished.Start to re-fixing CHR and POS... ",verbose=verbose)
    sumstats = fixchr(sumstats,verbose=verbose)
    sumstats = fixpos(sumstats,verbose=verbose)
    sumstats = sortcolumn(sumstats,verbose=verbose)

    finished(log,verbose,_end_line)
    return sumstats
    #################################################################################################### 


####################################################################################################################
    
def merge_chrpos(sumstats_part,all_groups_max,path,build,status):
    group=str(sumstats_part["group"].mode(dropna=True)[0])
    if group in [str(i) for i in range(all_groups_max+1)]:
        try:
            to_merge=pd.read_hdf(path, key="group_"+str(group)).drop_duplicates(subset="rsn")
            to_merge = to_merge.set_index("rsn")
            is_chrpos_fixable = sumstats_part.index.isin(to_merge.index)
            sumstats_part.loc[is_chrpos_fixable,status] = vchange_status(sumstats_part.loc[is_chrpos_fixable, status],  1,"139",3*build[0])
            sumstats_part.loc[is_chrpos_fixable,status] = vchange_status(sumstats_part.loc[is_chrpos_fixable, status],  2,"987",3*build[1])
            sumstats_part.update(to_merge)
        except:
            pass
    return sumstats_part


def parallelrsidtochrpos(sumstats, rsid="rsID", chrom="CHR",pos="POS", path=None, ref_rsid_to_chrpos_vcf = None, ref_rsid_to_chrpos_hdf5 = None, build="99",status="STATUS",
                         n_cores=4,block_size=20000000,verbose=True,log=Log()):

    ##start function with col checking##########################################################
    _start_line = "assign CHR and POS using rsIDs"
    _end_line = "assigning CHR and POS using rsIDs"
    _start_cols = [rsid]
    _start_function = ".rsid_to_chrpos2()"
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

    if ref_rsid_to_chrpos_hdf5 is not None:
        path = ref_rsid_to_chrpos_hdf5
    elif ref_rsid_to_chrpos_vcf is not None:
        vcf_file_name = os.path.basename(ref_rsid_to_chrpos_vcf)
        vcf_dir_path = os.path.dirname(ref_rsid_to_chrpos_vcf)
        path = "{}/{}.rsID_CHR_POS_groups_{}.h5".format(vcf_dir_path,vcf_file_name,int(block_size))
    
    if path is None:
        raise ValueError("Please provide path to hdf5 file.")
    
    sumstats["rsn"] = pd.to_numeric(sumstats[rsid].str.strip("rs"),errors="coerce").astype("Int64")
    
    log.write(" -Source hdf5 file: ",path,verbose=verbose)
    log.write(" -Cores to use : ",n_cores,verbose=verbose)
    log.write(" -Blocksize (make sure it is the same as hdf5 file ): ",block_size,verbose=verbose)
    
    input_columns= sumstats.columns
    sumstats_nonrs = sumstats.loc[sumstats["rsn"].isna()|sumstats["rsn"].duplicated(keep='first') ,:].copy()
    sumstats_rs  = sumstats.loc[sumstats["rsn"].notnull(),:].copy()
    
    log.write(" -Non-Valid rsIDs: ",sum(sumstats["rsn"].isna()),verbose=verbose)
    log.write(" -Duplicated rsIDs except for the first occurrence: ",sum(sumstats.loc[~sumstats["rsn"].isna(), "rsn"].duplicated(keep='first')),verbose=verbose)
    log.write(" -Valid rsIDs: ", len(sumstats_rs),verbose=verbose)
    
    del sumstats
    gc.collect()
    
    # assign group number
    sumstats_rs.loc[:,"group"]= sumstats_rs.loc[:,"rsn"]//block_size
    
    # all groups
    
    
    # set index
    sumstats_rs = sumstats_rs.set_index("rsn")
    
    #
    pool = Pool(n_cores)
    if chrom not in input_columns:
        log.write(" -Initiating CHR ... ",verbose=verbose)
        sumstats_rs[chrom]=pd.Series(dtype="Int64") 
        
    if pos not in input_columns:
        log.write(" -Initiating POS ... ",verbose=verbose)
        sumstats_rs[pos]=pd.Series(dtype="Int64") 
    
    df_split=[y for x, y in sumstats_rs.groupby('group', as_index=False)]
    log.write(" -Divided into groups: ",len(df_split),verbose=verbose)
    log.write("  -",set(sumstats_rs.loc[:,"group"].unique()),verbose=verbose)
    
    # check keys
    store = pd.HDFStore(path, 'r')
    all_groups = store.keys()
    all_groups_len = len(all_groups)
    store.close()
    all_groups_max = max(map(lambda x: int(x.split("_")[1]), all_groups))
    log.write(" -Number of groups in HDF5: ",all_groups_len,verbose=verbose)
    log.write(" -Max index of groups in HDF5: ",all_groups_max,verbose=verbose)

    # update CHR and POS using rsID with multiple threads
    sumstats_rs = pd.concat(pool.map(partial(merge_chrpos,all_groups_max=all_groups_max,path=path,build=build,status=status),df_split),ignore_index=True)
    sumstats_rs[["CHR","POS"]] = sumstats_rs[["CHR","POS"]].astype("Int64")
    del df_split
    gc.collect()
    log.write(" -Merging group data... ",verbose=verbose)
    # drop group and rsn
    sumstats_rs = sumstats_rs.drop(columns=["group"])
    sumstats_nonrs = sumstats_nonrs.drop(columns=["rsn"])
    
    # merge back
    log.write(" -Append data... ",verbose=verbose)
    sumstats = pd.concat([sumstats_rs,sumstats_nonrs],ignore_index=True)
    
    del sumstats_rs
    del sumstats_nonrs
    gc.collect()
    
    # check
    sumstats = fixchr(sumstats,verbose=verbose)
    sumstats = fixpos(sumstats,verbose=verbose)
    sumstats = sortcolumn(sumstats,verbose=verbose)

    pool.close()
    pool.join()

    finished(log, verbose, _end_line)
    return sumstats
####################################################################################################################
# old version
def _old_check_status(row,record):
    #pos,ea,nea
    # status 
    #0 /  ----->  match
    #1 /  ----->  Flipped Fixed
    #2 /  ----->  Reverse_complementary Fixed
    #3 /  ----->  flipped
    #4 /  ----->  reverse_complementary 
    #5 / ------>  reverse_complementary + flipped
    #6 /  ----->  both allele on genome + unable to distinguish
    #7 /  ----> reverse_complementary + both allele on genome + unable to distinguish
    #8 / -----> not on ref genome
    #9 / ------> unchecked
    
    status_pre=row.iloc[3][:5]
    status_end=row.iloc[3][6:]
    
    ## nea == ref
    if row.iloc[2] == record[row.iloc[0]-1: row.iloc[0]+len(row.iloc[2])-1].seq.upper():
        ## ea == ref
        if row.iloc[1] == record[row.iloc[0]-1: row.iloc[0]+len(row.iloc[1])-1].seq.upper():
            ## len(nea) >len(ea):
            if len(row.iloc[2])!=len(row.iloc[1]):
                # indels both on ref, unable to identify
                return status_pre+"6"+status_end 
        else:
            #nea == ref & ea != ref
            return status_pre+"0"+status_end 
    ## nea!=ref
    else:
        # ea == ref_seq -> need to flip
        if row.iloc[1] == record[row.iloc[0]-1: row.iloc[0]+len(row.iloc[1])-1].seq.upper():
            return status_pre+"3"+status_end 
        # ea !=ref
        else:
            #_reverse_complementary
            row.iloc[1] = get_reverse_complementary_allele(row.iloc[1])
            row.iloc[2] = get_reverse_complementary_allele(row.iloc[2])
            ## nea == ref
            if row.iloc[2] == record[row.iloc[0]-1: row.iloc[0]+len(row.iloc[2])-1].seq.upper():
                ## ea == ref
                if row.iloc[1] == record[row.iloc[0]-1: row.iloc[0]+len(row.iloc[1])-1].seq.upper():
                    ## len(nea) >len(ea):
                    if len(row.iloc[2])!=len(row.iloc[1]):
                        return status_pre+"8"+status_end  # indel reverse complementary
                else:
                    return status_pre+"4"+status_end 
            else:
                # ea == ref_seq -> need to flip
                if row.iloc[1] == record[row.iloc[0]-1: row.iloc[0]+len(row.iloc[1])-1].seq.upper():
                    return status_pre+"5"+status_end 
            # ea !=ref
            return status_pre+"8"+status_end

def oldcheckref(sumstats,ref_seq,chrom="CHR",pos="POS",ea="EA",nea="NEA",status="STATUS",chr_dict=get_chr_to_number(),remove=False,verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "check if NEA is aligned with reference sequence"
    _end_line = "checking if NEA is aligned with reference sequence"
    _start_cols = [chrom,pos,ea,nea,status]
    _start_function = ".check_ref()"
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
    log.write(" -Reference genome FASTA file: "+ ref_seq,verbose=verbose)  
    log.write(" -Checking records: ", end="",verbose=verbose)  
    chromlist = get_chr_list(add_number=True)
    records = SeqIO.parse(ref_seq, "fasta")
    for record in records:
        #record = next(records)
        if record is not None:
            record_chr = str(record.id).strip("chrCHR").upper()
            if record_chr in chr_dict.keys():
                i = chr_dict[record_chr]
            else:
                i = record_chr
            if i in chromlist:
                log.write(record_chr," ", end="",show_time=False,verbose=verbose) 
                to_check_ref = (sumstats[chrom]==i) & (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna())
                sumstats.loc[to_check_ref,status] = sumstats.loc[to_check_ref,[pos,ea,nea,status]].apply(lambda x:_old_check_status(x,record),axis=1)
    
    log.write("\n",end="",show_time=False,verbose=verbose) 
        
    CATEGORIES = {str(j+i) for j in [1300000,1800000,1900000,3800000,9700000,9800000,9900000] for i in range(0,100000)}
    sumstats[status] = pd.Categorical(sumstats[status],categories=CATEGORIES)
    #sumstats[status] = sumstats[status].astype("string")


    available_to_check =sum( (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna()))
    status_0=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[0]\w", case=False, flags=0, na=False))
    status_3=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[3]\w", case=False, flags=0, na=False))
    status_4=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[4]\w", case=False, flags=0, na=False))
    status_5=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[5]\w", case=False, flags=0, na=False))
    status_6=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[6]\w", case=False, flags=0, na=False))
    #status_7=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[7]\w", case=False, flags=0, na=False))
    status_8=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w", case=False, flags=0, na=False))
    
    log.write(" -Variants allele on given reference sequence : ",status_0,verbose=verbose)
    log.write(" -Variants flipped : ",status_3,verbose=verbose)
    raw_matching_rate = (status_3+status_0)/available_to_check
    flip_rate = status_3/available_to_check
    log.write("  -Raw Matching rate : ","{:.2f}%".format(raw_matching_rate*100),verbose=verbose)
    if raw_matching_rate <0.8:
        log.warning("Matching rate is low, please check if the right reference genome is used.")
    if flip_rate > 0.85 :
        log.write("  -Flipping variants rate > 0.85, it is likely that the EA is aligned with REF in the original dataset.",verbose=verbose)
    
    log.write(" -Variants inferred reverse_complement : ",status_4,verbose=verbose)
    log.write(" -Variants inferred reverse_complement_flipped : ",status_5,verbose=verbose)
    log.write(" -Both allele on genome + unable to distinguish : ",status_6,verbose=verbose)
    #log.write(" -Reverse_complementary + both allele on genome + unable to distinguish: ",status_7)
    log.write(" -Variants not on given reference sequence : ",status_8,verbose=verbose)
    
    if remove is True:
        sumstats = sumstats.loc[~sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w"),:]
        log.write(" -Variants not on given reference sequence were removed.",verbose=verbose)
    
    finished(log, verbose, _end_line)
    return sumstats

#20240320 check if non-effect allele is aligned with reference genome         
def _fast_check_status(x: pd.DataFrame, record: np.array, starting_positions: np.array, records_len: np.array):
    # starting_positions and records_len must be 1D arrays containing data only for the chromosomes contained in x,
    # and these arrays must be ordered in the same way as the chromosomes in np.unique(x['CHR'].values).

    # status 
    #0 /  ----->  match
    #1 /  ----->  Flipped Fixed
    #2 /  ----->  Reverse_complementary Fixed
    #3 /  ----->  flipped
    #4 /  ----->  reverse_complementary 
    #5 / ------>  reverse_complementary + flipped
    #6 /  ----->  both allele on genome + unable to distinguish
    #7 /  ----> reverse_complementary + both allele on genome + unable to distinguish
    #8 / -----> not on ref genome
    #9 / ------> unchecked
    if x.empty:
        return np.array([])
    
    # x is expected to be a DataFrame with these columns in that order: ['CHR', 'POS', 'EA', 'NEA', 'STATUS']
    # In this way, we don't need to specify the columns names
    _chrom = x.iloc[:, 0]
    _pos = x.iloc[:, 1]
    _ea = x.iloc[:, 2]
    _nea = x.iloc[:, 3]
    _status = x.iloc[:, 4]

    # position of the status (i.e. x['STATUS']) that will be modified
    status_flip_idx = 5 

    pos = _pos.values.astype(np.int64) # convert to int64 because they could be of type 'object'

    # Rebase the chromosome numbers to 0-based indexing
    # e.g. ['1', '2', '4', '2'] -> [0, 1, 2, 1]
    # This is needed because record is a single 1D array containing all the records for all the selected chromosomes,
    # so for instance if record contains the records for chr1, chr2, chr4 ([...chr1...chr2...chr4...]), we need to 
    # rebase the chromosome numbers to 0-based indexing to index the correct record portion when we do starting_positions[chrom]
    # Note that in x there are only the rows for the same chromosomes for which we have the records in record
    # (i.e. we don't have rows for chr3 if we don't have the record for chr3). This filtering is done in the caller function
    _chrom = _chrom.values
    unique_values, _ = np.unique(_chrom, return_inverse=True) # Get the sorted unique values and their indices
    chrom = np.searchsorted(unique_values, _chrom) # Replace each value in '_chrom' with its corresponding index in the sorted unique values

    max_len_nea = _nea.str.len().max()
    max_len_ea = _ea.str.len().max()

    ########################################## mask for variants with out of range POS
    mask_outlier = pos > records_len[chrom]
    #########################################

    # Let's apply the same magic used for the fasta records (check build_fasta_records() for details) to convert the NEA and EA to
    # a numpy array of integers in a very fast way.
    # In that case we start from a pd.Series to we can apply some built-in methods.
    # Also, when doing nea.view('<u4'), each row will be automatically right-padded with zeros to reach the max_len_nea.
    # For this reason, we then replace the zeros with out padding value
    # (and that's why the mapping dict can't have chr(0) as a value, otherwise we would have zeros for both padding and a character)
    # Reshaping is needed because .view('<u4') will create a flattened array    
    nea = _nea.str.translate(TRANSLATE_TABLE).to_numpy().astype(f'<U{max_len_nea}')
    nea = nea.view('<u4').reshape(-1, max_len_nea).astype(np.uint8)
    nea[nea == 0] = PADDING_VALUE # padding value
    ###########################################
    
    ###########################################
    # Create a mask holding True at the position of non-padding values
    mask_nea = nea != PADDING_VALUE

    # Create the reverse complement of NEA
    # In this case, we manually left-pad the translated string with the padding value, since the padding done by view('<u4') would be right-padded
    # and that will make hard the reverse operation (because we would have e.g. [2, 2, 4, 100, ..., 100] which will be hard to convert into [4, 2, 2, 100, ..., 100])
    rev_nea = _nea.str.translate(TRANSLATE_TABLE_COMPL).str.pad(max_len_nea, 'left', chr(PADDING_VALUE)).to_numpy().astype(f'<U{max_len_nea}')
    rev_nea = rev_nea.view('<u4').reshape(-1, max_len_nea).astype(np.uint8)
    rev_nea = rev_nea[:, ::-1]


    # Let's do everything again for EA
    ea = _ea.str.translate(TRANSLATE_TABLE).to_numpy().astype(f'<U{max_len_ea}')
    ea = ea.view('<u4').reshape(-1, max_len_ea).astype(np.uint8)
    ea[ea == 0] = PADDING_VALUE # padding value
    ###########################################
    
    ###########################################
    mask_ea = ea != PADDING_VALUE

    rev_ea = _ea.str.translate(TRANSLATE_TABLE_COMPL).str.pad(max_len_ea, 'left', chr(PADDING_VALUE)).to_numpy().astype(f'<U{max_len_ea}')
    rev_ea = rev_ea.view('<u4').reshape(-1, max_len_ea).astype(np.uint8)
    rev_ea = rev_ea[:, ::-1]


    # Convert the status (which are integers represented as strings) to a numpy array of integers.
    # Again, use the same concept as before to do this in a very fast way.
    # e.g. ["9999999", "9939999", "9929999"] -> [[9, 9, 9, 9, 9, 9, 9], [9, 9, 3, 9, 9, 9, 9], [9, 9, 2, 9, 9, 9, 9]]
    assert _status.str.len().value_counts().nunique() == 1 # all the status strings should have the same length, let's be sure of that.
    status_len = len(_status.iloc[0])
    mapping_status = {str(v): chr(v) for v in range(10)}
    table_stats = _maketrans(mapping_status)
    status = _status.str.translate(table_stats).to_numpy().astype(f'<U{status_len}')
    status = status.view('<u4').reshape(-1, status_len).astype(np.uint8)


    # Expand the position to a 2D array and subtract 1 to convert to 0-based indexing
    # e.g. [2, 21, 46] -> [[1], [20], [45]]
    pos = np.expand_dims(pos, axis=-1) - 1

    # Create a modified indices array specifying the starting position of each chromosome in the concatenated record array
    modified_indices = starting_positions[chrom]
    modified_indices = modified_indices[:, np.newaxis] # Add a new axis to modified_indices to align with the dimensions of pos

    # Create the range of indices: [0, ..., max_len_nea-1]
    indices_range = np.arange(max_len_nea)

    # Add the range of indices to the starting indices
    # e.g. pos = [[1], [20], [45]], indices_range = [0, 1, 2], indices = [[1, 2, 3], [20, 21, 22], [45, 46, 47]]
    indices = pos + indices_range

    # Modify indices to select the correct absolute position in the concatenated record array
    indices = indices + modified_indices

    # Let's pad the fasta records array because if there is a (pos, chrom) for which (pos+starting_position[chrom]+max_len_nea > len(record) we get out of bounds error.
    # This basically happens if there is a pos for the last chromosome for which pos+max_len_nea > len(record for that chrom).
    # This is very unlikely to happen but we should handle this case.
    record = np.pad(record, (0, max_len_nea), constant_values=PADDING_VALUE)
    
    # Index the record array using the computed indices.
    # Since we use np.take, indices must all have the same length, and this is why we added the padding to NEA
    # and we create the indices using max_len_nea (long story short, we can't obtain a scattered/ragged array)
    output_nea = np.take(record, indices, mode="clip")
    ##################################################################
    output_nea[mask_outlier] = PADDING_VALUE
    ##################################################################
    
    # Check if the NEA is equal to the reference sequence at the given position
    # In a non-matrix way, this is equivalent (for one single element) to:
    # nea == record[pos-1: pos+len(nea)-1]
    # where for example:
    #  a) nea = "AC", record = "ACTG", pos = 1 -> True
    #  b) nea = "T", record = "ACTG", pos = 3 -> True
    #  c) nea = "AG", record = "ACTG", pos = 1 -> False
    # Since we want to do everything in a vectorized way, we will compare the padded NEA with the output 
    # and then we use the mask to focus only on the non-padded elements
    # Pseudo example (X represents the padding value):
    #  nea = ['AC', 'T'], record = 'ACTGAAG', pos = [1, 3]
    #  -> nea = ['AC', 'TX'], indices = [[1, 2], [3, 4]], mask = [[True, True], [True, False]], output_nea = [['A', 'C'], ['T', 'G']]
    #  -> nea == output_nea: [[True, True], [True, False]], mask: [[True, True], [True, False]]
    #  -> nea == output_nea + ~mask: [[True, True], [True, True]]
    #  -> np.all(nea == output_nea + ~mask, 1): [True, True]
    nea_eq_ref = np.all((nea == output_nea) + ~mask_nea, 1)
    rev_nea_eq_ref = np.all((rev_nea == output_nea) + ~mask_nea, 1)

    # Let's do everything again for EA
    indices_range = np.arange(max_len_ea)
    indices = pos + indices_range
    indices = indices + modified_indices
    output_ea = np.take(record, indices, mode="clip")
    ##################################################################
    output_ea[mask_outlier] = PADDING_VALUE
    ##################################################################

    ea_eq_ref = np.all((ea == output_ea) + ~mask_ea, 1)
    rev_ea_eq_ref = np.all((rev_ea == output_ea) + ~mask_ea, 1)

    masks_max_len = max(mask_nea.shape[1], mask_ea.shape[1])

    len_nea_eq_len_ea = np.all(
        np.pad(mask_nea, ((0,0),(0, masks_max_len-mask_nea.shape[1])), constant_values=False) == 
        np.pad(mask_ea, ((0,0),(0, masks_max_len-mask_ea.shape[1])), constant_values=False)
        , axis=1) # pad masks with False to reach same shape
    len_rev_nea_eq_rev_len_ea = len_nea_eq_len_ea

    # The following conditions replicates the if-else statements of the original check_status function:
    # https://github.com/Cloufield/gwaslab/blob/f6b4c4e58a26e5d67d6587141cde27acf9ce2a11/src/gwaslab/hm_harmonize_sumstats.py#L238

    # nea == ref && ea == ref && len(nea) != len(ea)
    status[nea_eq_ref * ea_eq_ref * ~len_nea_eq_len_ea, status_flip_idx] = 6

    # nea == ref && ea != ref
    status[nea_eq_ref * ~ea_eq_ref, status_flip_idx] = 0

    # nea != ref && ea == ref
    status[~nea_eq_ref * ea_eq_ref, status_flip_idx] = 3

    # nea != ref && ea != ref && rev_nea == ref && rev_ea == ref && len(rev_nea) != len(rev_ea)
    status[~nea_eq_ref * ~ea_eq_ref * rev_nea_eq_ref * rev_ea_eq_ref * ~len_rev_nea_eq_rev_len_ea, status_flip_idx] = 8

    # nea != ref && ea != ref && rev_nea == ref && rev_ea != ref
    status[~nea_eq_ref * ~ea_eq_ref * rev_nea_eq_ref * ~rev_ea_eq_ref, status_flip_idx] = 4

    # nea != ref && ea != ref && rev_nea != ref && rev_ea == ref
    status[~nea_eq_ref * ~ea_eq_ref * ~rev_nea_eq_ref * rev_ea_eq_ref, status_flip_idx] = 5

    # nea != ref && ea != ref && rev_nea != ref && rev_ea != ref
    status[~nea_eq_ref * ~ea_eq_ref * ~rev_nea_eq_ref * ~rev_ea_eq_ref, status_flip_idx] = 8

    # Convert back the (now modified) 2D status array to a numpy array of strings in a very fast way.
    # Since 'status' is a 2D array of integers ranging from 0 to 9, we can build the integer representation
    # of each row using the efficent operation below (e.g. [1, 2, 3, 4, 5] -> [12345]).
    # Then we convert this integer to a string using the f'<U{status.shape[1]}' dtype (e.g. 12345 -> '12345')
    # The "naive" way would be:
    #   status_str = [''.join(map(str, l)) for l in status]
    #   status_arr = np.array(status_str)
    status_flat = np.sum(status * 10**np.arange(status.shape[1]-1, -1, -1), axis=1)
    status_arr = status_flat.astype(f'<U{status.shape[1]}')

    return status_arr


def check_status(sumstats: pd.DataFrame, fasta_records_dict, log=Log(), verbose=True):

    chrom,pos,ea,nea,status = sumstats.columns

    # First, convert the fasta records to a single numpy array of integers
    record, starting_positions_dict, records_len_dict = build_fasta_records(fasta_records_dict, pos_as_dict=True, log=log, verbose=verbose)

    # In _fast_check_status(), several 2D numpy arrays are created and they are padded to have shape[1] == max_len_nea or max_len_ea
    # Since most of the NEA and EA strings are short, we perform the check first on the records having short NEA and EA strings,
    # and then we perform the check on the records having long NEA and EA strings. In this way we can speed up the process (since the 
    # arrays are smaller) and save memory.
    max_len = 4 # this is a chosen value, we could compute it using some stats about the length and count of NEA and EA strings
    condition = (sumstats[nea].str.len() <= max_len) & (sumstats[ea].str.len() <= max_len)

    log.write(f"   -Checking records for ( len(NEA) <= {max_len} and len(EA) <= {max_len} )", verbose=verbose)
    sumstats_cond = sumstats[condition]
    unique_chrom_cond = sumstats_cond[chrom].unique()
    starting_pos_cond = np.array([starting_positions_dict[k] for k in unique_chrom_cond])
    records_len_cond = np.array([records_len_dict[k] for k in unique_chrom_cond])
    sumstats.loc[condition, status] = _fast_check_status(sumstats_cond, record=record, starting_positions=starting_pos_cond, records_len=records_len_cond)

    log.write(f"   -Checking records for ( len(NEA) > {max_len} or len(EA) > {max_len} )", verbose=verbose)
    sumstats_not_cond = sumstats[~condition]
    unique_chrom_not_cond = sumstats_not_cond[chrom].unique()
    starting_not_pos_cond = np.array([starting_positions_dict[k] for k in unique_chrom_not_cond])
    records_len_not_cond = np.array([records_len_dict[k] for k in unique_chrom_not_cond])
    sumstats.loc[~condition, status] = _fast_check_status(sumstats_not_cond, record=record, starting_positions=starting_not_pos_cond, records_len=records_len_not_cond)

    return sumstats[status].values
        

def checkref(sumstats,ref_seq,chrom="CHR",pos="POS",ea="EA",nea="NEA",status="STATUS",chr_dict=get_chr_to_number(),remove=False,verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "check if NEA is aligned with reference sequence"
    _end_line = "checking if NEA is aligned with reference sequence"
    _start_cols = [chrom,pos,ea,nea,status]
    _start_function = ".check_ref()"
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
    log.write(" -Reference genome FASTA file: "+ ref_seq,verbose=verbose)  
    log.write(" -Loading fasta records:",end="", verbose=verbose)
    chromlist = get_chr_list(add_number=True)
    records = SeqIO.parse(ref_seq, "fasta")

    all_records_dict = {}
    chroms_in_sumstats = sumstats[chrom].unique() # load records from Fasta file only for the chromosomes present in the sumstats
    for record in records:
        #record = next(records)
        if record is not None:
            record_chr = str(record.id).strip("chrCHR").upper()
            if record_chr in chr_dict.keys():
                i = chr_dict[record_chr]
            else:
                i = record_chr
            if (i in chromlist) and (i in chroms_in_sumstats):
                log.write(record_chr," ", end="",show_time=False,verbose=verbose)
                all_records_dict.update({i: record})
    log.write("",show_time=False,verbose=verbose)

    if len(all_records_dict) > 0:
        log.write(" -Checking records", verbose=verbose)
        all_records_dict = dict(sorted(all_records_dict.items())) # sort by key in case the fasta records are not already ordered by chromosome
        to_check_ref = (sumstats[chrom].isin(list(all_records_dict.keys()))) & (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna())
        sumstats_to_check = sumstats.loc[to_check_ref,[chrom,pos,ea,nea,status]]
        sumstats.loc[to_check_ref,status] = check_status(sumstats_to_check, all_records_dict, log=log, verbose=verbose)
        log.write(" -Finished checking records", verbose=verbose) 
    
    CATEGORIES = {str(j+i) for j in [1300000,1800000,1900000,3800000,9700000,9800000,9900000] for i in range(0,100000)}
    sumstats[status] = pd.Categorical(sumstats[status],categories=CATEGORIES)
    #sumstats[status] = sumstats[status].astype("string")

    available_to_check =sum( (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna()))
    status_0=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[0]\w", case=False, flags=0, na=False))
    status_3=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[3]\w", case=False, flags=0, na=False))
    status_4=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[4]\w", case=False, flags=0, na=False))
    status_5=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[5]\w", case=False, flags=0, na=False))
    status_6=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[6]\w", case=False, flags=0, na=False))
    #status_7=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[7]\w", case=False, flags=0, na=False))
    status_8=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w", case=False, flags=0, na=False))
    
    log.write(" -Variants allele on given reference sequence : ",status_0,verbose=verbose)
    log.write(" -Variants flipped : ",status_3,verbose=verbose)
    raw_matching_rate = (status_3+status_0)/available_to_check
    flip_rate = status_3/available_to_check
    log.write("  -Raw Matching rate : ","{:.2f}%".format(raw_matching_rate*100),verbose=verbose)
    if raw_matching_rate <0.8:
        log.warning("Matching rate is low, please check if the right reference genome is used.")
    if flip_rate > 0.85 :
        log.write("  -Flipping variants rate > 0.85, it is likely that the EA is aligned with REF in the original dataset.",verbose=verbose)
    
    log.write(" -Variants inferred reverse_complement : ",status_4,verbose=verbose)
    log.write(" -Variants inferred reverse_complement_flipped : ",status_5,verbose=verbose)
    log.write(" -Both allele on genome + unable to distinguish : ",status_6,verbose=verbose)
    #log.write(" -Reverse_complementary + both allele on genome + unable to distinguish: ",status_7)
    log.write(" -Variants not on given reference sequence : ",status_8,verbose=verbose)
    
    if remove is True:
        sumstats = sumstats.loc[~sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w"),:]
        log.write(" -Variants not on given reference sequence were removed.",verbose=verbose)

    
    finished(log, verbose, _end_line)
    return sumstats

def build_fasta_records(fasta_records_dict, pos_as_dict=True, log=Log(), verbose=True):
    log.write("   -Building numpy fasta records from dict", verbose=verbose)

    # Let's do some magic to convert the fasta record to a numpy array of integers in a very fast way.
    # fasta_record.seq._data is a byte-string, so we can use the bytes.maketrans to apply a translation.
    # Here we map the bytes to the unicode character representing the desired integer as defined in the mapping dict
    # (i.e. b'A' -> '\x02', b'T' -> '\x03', b'C' -> '\x04', b'G' -> '\x05', b'N' -> '\x06')
    # Then, using np.array(... dtype=<U..) we convert the string to a numpy array of unicode characters.
    # Then, we do a magic with view('<u4') to convert the unicode characters to 4-byte integers, so we obtain the actual integer representation of the characters
    # Lastly, we cast the array to np.uint8 to convert the 4-byte integers to 1-byte integers to save memory
    # Full example:
    # fasta_record.seq._data = b'ACTGN' -> b'\x02\x04\x03\x05\x06' -> np.array(['\x02\x04\x03\x05\x06'], dtype='<U5') -> np.array([2, 4, 3, 5, 6], dtype=uint32) -> np.array([2, 4, 3, 5, 6], dtype=uint8)
    all_r = []
    for r in fasta_records_dict.values():
        r = r.seq._data.translate(TRANSLATE_TABLE)
        r = np.array([r], dtype=f'<U{len(r)}').view('<u4').astype(np.uint8)
        all_r.append(r)

    # We've just created a list of numpy arrays, so we can concatenate them to obtain a single numpy array
    # Then we keep track of the starting position of each record in the concatenated array. This will be useful later
    # to index the record array depending on the position of the variant and the chromosome
    records_len = np.array([len(r) for r in all_r])
    starting_positions = np.cumsum(records_len) - records_len
    if pos_as_dict:
        starting_positions = {k: v for k, v in zip(fasta_records_dict.keys(), starting_positions)}
        records_len_dict =  {k: v for k, v in zip(fasta_records_dict.keys(), records_len)}
    record = np.concatenate(all_r)
    del all_r # free memory

    return record, starting_positions,records_len_dict

#######################################################################################################################################

#20220721
def chrposref_rsid(chr,end,ref,alt,vcf_reader,chr_dict=get_number_to_chr()):
    ## single record assignment
    start=end-1
    if chr_dict is not None: chr=chr_dict[chr]
    
    try:
        chr_seq = vcf_reader.fetch(chr,start,end)
    except:
        return pd.NA
    
    for record in chr_seq:
        if record.pos==end: 
            if record.alts is None:
                return pd.NA
            if record.ref==ref and (alt in record.alts):
                return record.id
            elif (ref in record.alts) and record.ref==alt:
                return record.id
    return pd.NA

def assign_rsid_single(sumstats,path,rsid="rsID",chr="CHR",pos="POS",ref="NEA",alt="EA",chr_dict=get_number_to_chr()):
    ## single df assignment
    vcf_reader = VariantFile(path)
    def rsid_helper(x,vcf_reader,chr_dict):
         return chrposref_rsid(x.iloc[0],x.iloc[1],x.iloc[2],x.iloc[3],vcf_reader,chr_dict)
    map_func=partial(rsid_helper,vcf_reader=vcf_reader,chr_dict=chr_dict)
    rsID = sumstats.apply(map_func,axis=1)
    return rsID

def parallelizeassignrsid(sumstats, path, ref_mode="vcf",snpid="SNPID",rsid="rsID",chr="CHR",pos="POS",ref="NEA",alt="EA",status="STATUS",
                          n_cores=1,chunksize=5000000,ref_snpid="SNPID",ref_rsid="rsID",
                          overwrite="empty",verbose=True,log=Log(),chr_dict=None):
    '''
    overwrite mode : 
    all ,    overwrite rsid for all availalbe rsid 
    invalid,  only assign rsid for variants with invalid rsid
    empty    only assign rsid for variants with na rsid
    '''

    if ref_mode=="vcf":
        ###################################################################################################################
        ##start function with col checking##########################################################
        _start_line = "assign rsID using reference VCF"
        _end_line = "assign rsID using reference file"
        _start_cols = [chr,pos,ref,alt,status]
        _start_function = ".assign_rsid()"
        _must_args ={}

        is_enough_info = start_to(sumstats=sumstats,
                                log=log,
                                verbose=verbose,
                                start_line=_start_line,
                                end_line=_end_line,
                                start_cols=_start_cols,
                                start_function=_start_function,
                                n_cores=n_cores,
                                ref_vcf=path,
                                **_must_args)
        if is_enough_info == False: return sumstats
        ############################################################################################
        chr_dict = auto_check_vcf_chr_dict(path, chr_dict, verbose, log)
        log.write(" -Assigning rsID based on CHR:POS and REF:ALT/ALT:REF...",verbose=verbose)
        ##############################################
        if rsid not in sumstats.columns:
            sumstats[rsid]=pd.Series(dtype="string")

        ###############################################
        total_number= len(sumstats)
        pre_number = sum(~sumstats[rsid].isna())

        ##################################################################################################################
        standardized_normalized = sumstats["STATUS"].str.match("\w\w\w[0][01234]\w\w", case=False, flags=0, na=False)
        if overwrite=="all":
            to_assign = standardized_normalized
        if overwrite=="invalid":
            to_assign = (~sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)) & standardized_normalized
        if overwrite=="empty":
            to_assign = sumstats[rsid].isna()& standardized_normalized
        ##################################################################################################################
        # multicore arrangement

        if sum(to_assign)>0:
            if sum(to_assign)<10000: n_cores=1
            #df_split = np.array_split(sumstats.loc[to_assign, [chr,pos,ref,alt]], n_cores)
            df_split = _df_split(sumstats.loc[to_assign, [chr,pos,ref,alt]], n_cores)
            pool = Pool(n_cores)
            map_func = partial(assign_rsid_single,path=path,chr=chr,pos=pos,ref=ref,alt=alt,chr_dict=chr_dict) 
            assigned_rsid = pd.concat(pool.map(map_func,df_split))
            sumstats.loc[to_assign,rsid] = assigned_rsid.values 
            pool.close()
            pool.join()
        gc.collect()
        ##################################################################################################################

        after_number = sum(~sumstats[rsid].isna())
        log.write(" -rsID Annotation for "+str(total_number - after_number) +" need to be fixed!",verbose=verbose)
        log.write(" -Annotated "+str(after_number - pre_number) +" rsID successfully!",verbose=verbose)
    
    ##################################################################################################################
    elif ref_mode=="tsv":
        '''
        assign rsID based on chr:pos
        '''
        ##start function with col checking##########################################################
        _start_line = "assign rsID by matching SNPID with CHR:POS:REF:ALT in the reference TSV"
        _end_line = "assign rsID using reference file"
        _start_cols = [snpid,status]
        _start_function = ".assign_rsid()"
        _must_args ={}

        is_enough_info = start_to(sumstats=sumstats,
                                log=log,
                                verbose=verbose,
                                start_line=_start_line,
                                end_line=_end_line,
                                start_cols=_start_cols,
                                start_function=_start_function,
                                n_cores=n_cores,
                                ref_tsv=path,
                                **_must_args)
        if is_enough_info == False: return sumstats
        ############################################################################################
        
        #standardized_normalized = sumstats["STATUS"].str.match("\w\w\w[0][01234]\w\w", case=False, flags=0, na=False)
        standardized_normalized = sumstats["STATUS"] == sumstats["STATUS"]

        if rsid not in sumstats.columns:
            sumstats[rsid]=pd.Series(dtype="string")
            
        if overwrite == "empty":
            to_assign = sumstats[rsid].isna() & standardized_normalized
        if overwrite=="all":
            to_assign = standardized_normalized
        if overwrite=="invalid":
            to_assign = (~sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)) & standardized_normalized
            
        total_number= len(sumstats)
        pre_number = sum(~sumstats[rsid].isna())
        log.write(" -"+str(sum(to_assign)) +" rsID could be possibly fixed...",verbose=verbose)
        if sum(to_assign)>0: 
            sumstats = sumstats.set_index(snpid)  
            dic_chuncks = pd.read_csv(path,sep="\t",usecols=[ref_snpid,ref_rsid],
                              chunksize=chunksize,index_col=ref_snpid,
                              dtype={ref_snpid:"string",ref_rsid:"string"})

            log.write(" -Setting block size: ",chunksize,verbose=verbose)
            log.write(" -Loading block: ",end="",verbose=verbose)     
            for i,dic in enumerate(dic_chuncks):
                gc.collect()
                log.write(i," ",end=" ",show_time=False)  
                dic = dic.rename(index={ref_snpid:snpid})
                dic = dic.rename(columns={ref_rsid:rsid})  
                dic = dic.loc[~dic.index.duplicated(keep=False),:]
                sumstats.update(dic,overwrite=True)

            log.write("\n",end="",show_time=False,verbose=verbose) 
            sumstats = sumstats.reset_index()
            sumstats = sumstats.rename(columns = {'index':snpid})

            after_number = sum(~sumstats[rsid].isna())
            log.write(" -rsID annotation for "+str(total_number - after_number) +" needed to be fixed!",verbose=verbose)
            log.write(" -Annotated "+str(after_number - pre_number) +" rsID successfully!",verbose=verbose)
        else:
            log.write(" -No rsID can be fixed...skipping...",verbose=verbose)
        ################################################################################################################
            
    finished(log,verbose,_end_line)
    return sumstats
#################################################################################################################################################
#single record assignment

def check_strand_status(chr,start,end,ref,alt,eaf,vcf_reader,alt_freq,status,chr_dict=get_number_to_chr()):
    ### 0 : not palindromic
    ### 1 : palindromic +strand 
    ### 2 : palindromic -strand -> need to flip -> flipped
    ### 5 : palindromic -strand -> need to flip
    ### 8 : no ref data
    if chr_dict is not None: chr=chr_dict[chr]
    status_pre=status[:6]
    status_end=""
    try:
        chr_seq = vcf_reader.fetch(chr,start,end)
    except:
        return status_pre+"8"+status_end
        
    
    for record in chr_seq:
        if record.pos==end and record.ref==ref and (alt in record.alts):
            
            if  (record.info[alt_freq][0]<0.5) and (eaf<0.5):
                return status_pre+"1"+status_end
            elif (record.info[alt_freq][0]>0.5) and (eaf>0.5):
                return status_pre+"1"+status_end
            else:
                return status_pre+"5"+status_end
    return status_pre+"8"+status_end

def check_strand_status_cache(data,cache,ref_infer=None,ref_alt_freq=None,chr_dict=get_number_to_chr(),trust_cache=True,log=Log(),verbose=True):
    if not trust_cache:
        assert ref_infer is not None, "If trust_cache is False, ref_infer must be provided"
        log.warning("You are not trusting the cache, this will slow down the process. Please consider building a complete cache.")

    if ref_infer is not None and not trust_cache:
        vcf_reader = VariantFile(ref_infer)

    if isinstance(data, pd.DataFrame):
        data = data.values
    
    in_cache = 0
    new_statuses = []
    
    for i in range(data.shape[0]):
        _chrom, pos, ref, alt, eaf, status = data[i]
        chrom = _chrom
        start = pos - 1
        end = pos
        
        if chr_dict is not None: chrom=chr_dict[chrom]
        
        status_pre=status[:6]
        status_end=""
        
        new_status = status_pre+"8"+status_end # default value
        
        cache_key = f"{chrom}:{pos}:{ref}:{alt}"
        if cache_key in cache:
            in_cache += 1
            record = cache[cache_key]
            if record is None:
                new_status = status_pre+"8"+status_end
            else:
                if (record<0.5) and (eaf<0.5):
                    new_status = status_pre+"1"+status_end
                elif (record>0.5) and (eaf>0.5):
                    new_status = status_pre+"1"+status_end
                else:
                    new_status = status_pre+"5"+status_end
        else:
            if not trust_cache:
                # If we don't trust the cache as a not complete cache, we should perform the check reading from the VCF file
                new_status = check_strand_status(_chrom, start, end, ref, alt, eaf, vcf_reader, ref_alt_freq, status, chr_dict)
        
        new_statuses.append(new_status)
        
    log.write(f"  -Elements in cache: {in_cache}", verbose=verbose)
    return new_statuses


def check_unkonwn_indel(chr,start,end,ref,alt,eaf,vcf_reader,alt_freq,status,chr_dict=get_number_to_chr(),daf_tolerance=0.2):
    ### input : unknown indel, both on genome (xx1[45]x)
    ### 3 no flip
    ### 4 unknown indel,fixed   (6->5)
    ### 6 flip
    
    if chr_dict is not None: chr=chr_dict[chr]
    status_pre=status[:6]
    status_end=""
    
    try:
        chr_seq = vcf_reader.fetch(chr,start,end)
    except:
        return status_pre+"8"+status_end

    for record in chr_seq:
        if record.pos==end and record.ref==ref and (alt in record.alts):
            if  abs(record.info[alt_freq][0] - eaf)<daf_tolerance:
                return status_pre+"3"+status_end
   
        elif record.pos==end and record.ref==alt and (ref in record.alts):
            if  abs(record.info[alt_freq][0] - (1 - eaf))<daf_tolerance:
                return status_pre+"6"+status_end

    return status_pre+"8"+status_end


def check_unkonwn_indel_cache(data,cache,ref_infer=None,ref_alt_freq=None,chr_dict=get_number_to_chr(),daf_tolerance=0.2,trust_cache=True,log=Log(),verbose=True):
    if not trust_cache:
        assert ref_infer is not None, "If trust_cache is False, ref_infer must be provided"
        log.warning("You are not trusting the cache, this will slow down the process. Please consider building a complete cache.")

    if ref_infer is not None:
        vcf_reader = VariantFile(ref_infer)

    if isinstance(data, pd.DataFrame):
        data = data.values
    
    in_cache = 0
    new_statuses = []
    
    for i in range(data.shape[0]):
        _chrom, pos, ref, alt, eaf, status = data[i]
        chrom = _chrom
        
        if chr_dict is not None: chrom=chr_dict[chrom]
        start = pos - 1
        end = pos
        
        status_pre=status[:6]
        status_end=""
        
        new_status = status_pre+"8"+status_end # default value
        
        cache_key_ref_alt = f"{chrom}:{pos}:{ref}:{alt}"
        cache_key_alt_ref = f"{chrom}:{pos}:{alt}:{ref}"

        if cache_key_ref_alt in cache:
            in_cache += 1
            record = cache[cache_key_ref_alt]
            if record is None:
                new_status = status_pre+"8"+status_end
            else:
                if  abs(record - eaf)<daf_tolerance:
                    new_status = status_pre+"3"+status_end

        elif cache_key_alt_ref in cache:
            in_cache += 1
            record = cache[cache_key_alt_ref]
            if record is None:
                new_status = status_pre+"8"+status_end
            else:
                if  abs(record - (1 - eaf))<daf_tolerance:
                    new_status = status_pre+"6"+status_end

        else:
            if not trust_cache:
                # If we don't trust the cache as a not complete cache, we should perform the check reading from the VCF file
                new_status = check_unkonwn_indel(_chrom, start, end, ref, alt, eaf, vcf_reader, ref_alt_freq, status, chr_dict, daf_tolerance)
                
        new_statuses.append(new_status)
        
    log.write(f"  -Elements in cache: {in_cache}", verbose=verbose)
    return new_statuses

                                               
def get_reverse_complementary_allele(a):
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a[::-1].translate(dic)
                                                 
def is_palindromic(sumstats,a1="EA",a2="NEA"):
    gc= (sumstats[a1]=="G") & (sumstats[a2]=="C")
    cg= (sumstats[a1]=="C") & (sumstats[a2]=="G")
    at= (sumstats[a1]=="A") & (sumstats[a2]=="T")
    ta= (sumstats[a1]=="T") & (sumstats[a2]=="A")
    palindromic = gc | cg | at | ta 
    return palindromic
##################################################################################################################################################
#single df assignment

def check_strand(sumstats,ref_infer,ref_alt_freq=None,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",chr_dict=get_number_to_chr(),status="STATUS"):
    vcf_reader = VariantFile(ref_infer)
    status_part = sumstats.apply(lambda x:check_strand_status(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],x.iloc[4],vcf_reader,ref_alt_freq,x.iloc[5],chr_dict),axis=1) 
    return status_part

def check_strand_cache(sumstats,cache,ref_infer,ref_alt_freq=None,chr_dict=get_number_to_chr(),trust_cache=True,log=Log(),verbose=True):
    assert cache is not None, "Cache must be provided"
    status_part = check_strand_status_cache(sumstats,cache,ref_infer,ref_alt_freq,chr_dict,trust_cache,log,verbose)
    return status_part

def check_indel(sumstats,ref_infer,ref_alt_freq=None,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",chr_dict=get_number_to_chr(),status="STATUS",daf_tolerance=0.2):
    vcf_reader = VariantFile(ref_infer)
    status_part = sumstats.apply(lambda x:check_unkonwn_indel(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],x.iloc[4],vcf_reader,ref_alt_freq,x.iloc[5],chr_dict,daf_tolerance),axis=1)
    return status_part

def check_indel_cache(sumstats,cache,ref_infer,ref_alt_freq=None,chr_dict=get_number_to_chr(),daf_tolerance=0.2,trust_cache=True,log=Log(),verbose=True):
    assert cache is not None, "Cache must be provided"
    status_part = check_unkonwn_indel_cache(sumstats,cache,ref_infer,ref_alt_freq,chr_dict,daf_tolerance,trust_cache,log,verbose)
    return status_part

##################################################################################################################################################

def parallelinferstrand(sumstats,ref_infer,ref_alt_freq=None,maf_threshold=0.40,daf_tolerance=0.20,remove_snp="",mode="pi",n_cores=1,remove_indel="",
                       chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",status="STATUS",
                       chr_dict=None,cache_options={},verbose=True,log=Log()):
    '''
    Args:
    cache_options : A dictionary with the following keys:
        - cache_manager: CacheManager object or None. If any between cache_loader and cache_process is not None, or use_cache is True, a CacheManager object will be created automatically.
        - trust_cache: bool (optional, default: True). Whether to completely trust the cache or not. Trusting the cache means that any key not found inside the cache will be considered as a missing value even in the VCF file.
        - cache_loader: Object with a get_cache() method or None.
        - cache_process: Object with an apply_fn() method or None.
        - use_cache: bool (optional, default: False). If any of the cache_manager, cache_loader or cache_process is not None, this will be set to True automatically.
                     If set to True and all between cache_manager, cache_loader and cache_process are None, the cache will be loaded (or built) on the spot.

        The usefulness of a cache_loader or cache_process object is to pass a custom object which already has the cache loaded. This can be useful if the cache is loaded in background in another thread/process while other operations are performed.
        The cache_manager is a CacheManager object is used to expose the API to interact with the cache.
    '''

    ##start function with col checking##########################################################
    _start_line = "infer strand for palindromic SNPs/align indistinguishable indels"
    _end_line = "inferring strand for palindromic SNPs/align indistinguishable indels"
    _start_cols = [chr,pos,ref,alt,eaf,status]
    _start_function = ".infer_strand()"
    _must_args ={"ref_alt_freq":ref_alt_freq}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            n_cores=n_cores,
                            ref_vcf=ref_infer,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    chr_dict = auto_check_vcf_chr_dict(ref_infer, chr_dict, verbose, log)
    
    # Setup cache variables
    cache_manager = cache_options.get("cache_manager", None)
    if cache_manager is not None:
        assert isinstance(cache_manager, CacheManager), "cache_manager must be a CacheManager object"
    trust_cache = cache_options.get("trust_cache", True)
    cache_loader = cache_options.get("cache_loader", None)
    cache_process = cache_options.get("cache_process", None)
    use_cache = any(c is not None for c in [cache_manager, cache_loader, cache_process]) or cache_options.get('use_cache', False)
    _n_cores = n_cores # backup n_cores
    
    log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)  

    if "p" in mode:
        ## checking \w\w\w\w[0]\w\w -> standardized and normalized snp
        good_chrpos =  sumstats[status].str.match(r'\w\w\w[0][0]\w\w', case=False, flags=0, na=False) 
        palindromic = good_chrpos & is_palindromic(sumstats[[ref,alt]],a1=ref,a2=alt)   
        not_palindromic_snp = good_chrpos & (~palindromic)

        ##not palindromic : change status
        sumstats.loc[not_palindromic_snp,status] = vchange_status(sumstats.loc[not_palindromic_snp,status], 7 ,"9","0")  
        log.write(" -Identified ", sum(palindromic)," palindromic SNPs...",verbose=verbose)
        
        #palindromic but can not infer
        maf_can_infer   = (sumstats[eaf] < maf_threshold) | (sumstats[eaf] > 1 - maf_threshold)
        
        sumstats.loc[palindromic&(~maf_can_infer),status] = vchange_status(sumstats.loc[palindromic&(~maf_can_infer),status],7,"9","7")
        
        #palindromic WITH UNKNWON OR UNCHECKED STATUS
        unknow_palindromic = sumstats[status].str.match(r'\w\w\w\w\w[012][89]', case=False, flags=0, na=False) 

        unknow_palindromic_to_check = palindromic & maf_can_infer & unknow_palindromic
        
        log.write(" -After filtering by MAF< {} , {} palindromic SNPs with unknown strand will be inferred...".format(maf_threshold, sum(unknow_palindromic_to_check)),verbose=verbose)

        ######################################################################################### 
        if sum(unknow_palindromic_to_check)>0:
            if sum(unknow_palindromic_to_check)<10000: 
                n_cores=1

            if use_cache and cache_manager is None:
                cache_manager = CacheManager(base_path=ref_infer, cache_loader=cache_loader, cache_process=cache_process,
                                             ref_alt_freq=ref_alt_freq, category=PALINDROMIC_INDEL,
                                             n_cores=_n_cores, log=log, verbose=verbose)

            log.write(" -Starting strand inference for palindromic SNPs...",verbose=verbose)
            df_to_check = sumstats.loc[unknow_palindromic_to_check,[chr,pos,ref,alt,eaf,status]]
           
            if use_cache and cache_manager.cache_len > 0:
                log.write("  -Using cache for strand inference",verbose=verbose)
                status_inferred = cache_manager.apply_fn(check_strand_cache, sumstats=df_to_check, ref_infer=ref_infer, ref_alt_freq=ref_alt_freq, chr_dict=chr_dict, trust_cache=trust_cache, log=log, verbose=verbose)
                sumstats.loc[unknow_palindromic_to_check,status] = status_inferred
            else:
                #df_split = np.array_split(df_to_check, n_cores)
                df_split = _df_split(df_to_check, n_cores)
                pool = Pool(n_cores)
                map_func = partial(check_strand,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,status=status,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,chr_dict=chr_dict) 
                status_inferred = pd.concat(pool.map(map_func,df_split))
                sumstats.loc[unknow_palindromic_to_check,status] = status_inferred.values
                pool.close()
                pool.join()
            log.write(" -Finished strand inference.",verbose=verbose)
        else:
            log.warning("No palindromic variants available for checking.")
        #########################################################################################
        #0 Not palindromic SNPs
        #1 Palindromic +strand  -> no need to flip
        #2 palindromic -strand  -> need to flip -> fixed
        #3 Indel no need flip
        #4 Unknown Indel -> fixed
        #5 Palindromic -strand -> need to flip
        #6 Indel need flip
        #7 indistinguishable
        #8 Not matching or No information
        #9 Unchecked

        status0 = sumstats[status].str.match(r'\w\w\w\w\w\w[0]', case=False, flags=0, na=False)  
        status1 = sumstats[status].str.match(r'\w\w\w\w\w\w[1]', case=False, flags=0, na=False)  
        status5 = sumstats[status].str.match(r'\w\w\w\w\w\w[5]', case=False, flags=0, na=False)  
        status7 = sumstats[status].str.match(r'\w\w\w\w\w\w[7]', case=False, flags=0, na=False)  
        status8 = sumstats[status].str.match(r'\w\w\w\w\w[123][8]', case=False, flags=0, na=False)  

        log.write("  -Non-palindromic : ",sum(status0),verbose=verbose)
        log.write("  -Palindromic SNPs on + strand: ",sum(status1),verbose=verbose)
        log.write("  -Palindromic SNPs on - strand and needed to be flipped:",sum(status5),verbose=verbose)   
        log.write("  -Palindromic SNPs with MAF not available to infer : ",sum(status7),verbose=verbose)  
        log.write("  -Palindromic SNPs with no macthes or no information : ",sum(status8),verbose=verbose)  

        if ("7" in remove_snp) and ("8" in remove_snp) :
            log.write("  -Palindromic SNPs with MAF not available to infer and with no macthes or no information will will be removed",verbose=verbose) 
            sumstats = sumstats.loc[~(status7 | status8),:].copy()
        elif "8" in remove_snp:
            log.write("  -Palindromic SNPs with no macthes or no information will be removed",verbose=verbose)
            sumstats = sumstats.loc[~status8,:].copy()
        elif "7" in remove_snp:
            log.write("  -Palindromic SNPs with MAF not available to infer will be removed",verbose=verbose) 
            sumstats = sumstats.loc[~status7,:].copy()

    ### unknow_indel
    if "i" in mode:
        unknow_indel = sumstats[status].str.match(r'\w\w\w\w\w[6][89]', case=False, flags=0, na=False)   
        log.write(" -Identified ", sum(unknow_indel)," indistinguishable Indels...",verbose=verbose)
        if sum(unknow_indel)>0:
            log.write(" -Indistinguishable indels will be inferred from reference vcf REF and ALT...",verbose=verbose)
            #########################################################################################  
            #with maf can not infer
            #maf_can_infer   = (sumstats[eaf] < maf_threshold) | (sumstats[eaf] > 1 - maf_threshold) 
            #sumstats.loc[unknow_indel&(~maf_can_infer),status] = vchange_status(sumstats.loc[unknow_indel&(~maf_can_infer),status],7,"9","8") 
            log.write(" -Difference in allele frequency (DAF) tolerance: {}".format(daf_tolerance),verbose=verbose)
                         
            if sum(unknow_indel)>0:
                if sum(unknow_indel)<10000: 
                    n_cores=1

                if use_cache and cache_manager is None:
                    cache_manager = CacheManager(base_path=ref_infer, cache_loader=cache_loader, cache_process=cache_process,
                                                ref_alt_freq=ref_alt_freq, category=PALINDROMIC_INDEL,
                                                n_cores=_n_cores, log=log, verbose=verbose)

                log.write(" -Starting indistinguishable indel inference...",verbose=verbose)
                df_to_check = sumstats.loc[unknow_indel,[chr,pos,ref,alt,eaf,status]]
            
                if use_cache and cache_manager.cache_len > 0:
                    log.write("  -Using cache for indel inference",verbose=verbose)
                    status_inferred = cache_manager.apply_fn(check_indel_cache, sumstats=df_to_check, ref_infer=ref_infer, ref_alt_freq=ref_alt_freq, chr_dict=chr_dict, daf_tolerance=daf_tolerance, trust_cache=trust_cache, log=log, verbose=verbose)
                    sumstats.loc[unknow_indel,status] = status_inferred
                else:
                    #df_split = np.array_split(sumstats.loc[unknow_indel, [chr,pos,ref,alt,eaf,status]], n_cores)
                    df_split = _df_split(sumstats.loc[unknow_indel, [chr,pos,ref,alt,eaf,status]], n_cores)
                    pool = Pool(n_cores)
                    map_func = partial(check_indel,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,status=status,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,chr_dict=chr_dict,daf_tolerance=daf_tolerance) 
                    status_inferred = pd.concat(pool.map(map_func,df_split))
                    sumstats.loc[unknow_indel,status] = status_inferred.values 
                    pool.close()
                    pool.join()
                log.write(" -Finished indistinguishable indel inference.",verbose=verbose)

            #########################################################################################

            status3 =  sumstats[status].str.match(r'\w\w\w\w\w\w[3]', case=False, flags=0, na=False)  
            status6 =  sumstats[status].str.match(r'\w\w\w\w\w\w[6]', case=False, flags=0, na=False)  
            status8 =  sumstats[status].str.match(r'\w\w\w\w\w[6][8]', case=False, flags=0, na=False)  

            log.write("  -Indels ea/nea match reference : ",sum(status3),verbose=verbose)
            log.write("  -Indels ea/nea need to be flipped : ",sum(status6),verbose=verbose)
            log.write("  -Indels with no macthes or no information : ",sum(status8),verbose=verbose)
            if "8" in remove_indel:
                log.write("  -Indels with no macthes or no information will be removed",verbose=verbose)
                sumstats = sumstats.loc[~status8,:].copy()    
        else:
            log.warning("No indistinguishable indels available for checking.") 
    
    finished(log,verbose,_end_line)
    return sumstats




















################################################################################################################
def parallelecheckaf(sumstats,ref_infer,ref_alt_freq=None,maf_threshold=0.4,column_name="DAF",suffix="",n_cores=1, chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",status="STATUS",chr_dict=None,force=False, verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "check the difference between EAF and reference VCF ALT frequency"
    _end_line = "checking the difference between EAF and reference VCF ALT frequency"
    _start_cols = [chr,pos,ref,alt,eaf,status]
    _start_function = ".check_daf()"
    _must_args ={"ref_alt_freq":ref_alt_freq}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            n_cores=n_cores,
                            ref_vcf=ref_infer,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################

    chr_dict = auto_check_vcf_chr_dict(ref_infer, chr_dict, verbose, log)

    column_name = column_name + suffix
    


    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)  
        if not force:
            good_chrpos =  sumstats[status].str.match(r'\w\w\w[0]\w\w\w', case=False, flags=0, na=False)  
        log.write(" -Checking variants:", sum(good_chrpos),verbose=verbose) 
        sumstats[column_name]=np.nan
    
    ########################  
        if sum(~sumstats[eaf].isna())<10000: 
            n_cores=1       
        #df_split = np.array_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]], n_cores)
        df_split = _df_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]], n_cores)
        pool = Pool(n_cores)
        if sum(~sumstats[eaf].isna())>0:
            map_func = partial(checkaf,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,column_name=column_name,chr_dict=chr_dict) 
            sumstats.loc[good_chrpos,[column_name]] = pd.concat(pool.map(map_func,df_split))
        pool.close()
        pool.join()
    ###########################
        #status_inferred = sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]].apply(lambda x:check_daf(x[0],x[1]-1,x[1],x[2],x[3],x[4],vcf_reader,ref_alt_freq,chr_dict),axis=1)
        
        #sumstats.loc[good_chrpos,"DAF"] = status_inferred.values
        #sumstats["DAF"]=sumstats["DAF"].astype("float")     
        log.write(" - {} max:".format(column_name), np.nanmax(sumstats[column_name]),verbose=verbose)
        log.write(" - {} min:".format(column_name), np.nanmin(sumstats[column_name]),verbose=verbose)
        log.write(" - {} sd:".format(column_name), np.nanstd(sumstats[column_name]),verbose=verbose)
        log.write(" - abs({}) min:".format(column_name), np.nanmin(np.abs(sumstats[column_name])),verbose=verbose) 
        log.write(" - abs({}) max:".format(column_name), np.nanmax(np.abs(sumstats[column_name])),verbose=verbose)
        log.write(" - abs({}) sd:".format(column_name), np.nanstd(np.abs(sumstats[column_name])),verbose=verbose) 
        log.write("Finished allele frequency checking!") 
    return sumstats

def checkaf(sumstats,ref_infer,ref_alt_freq=None,column_name="DAF",chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",chr_dict=None):
    #vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    vcf_reader = VariantFile(ref_infer)
    def afapply(x,vcf,alt_freq,chr_dict):
            return check_daf(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],x.iloc[4],vcf_reader,ref_alt_freq,chr_dict)
    map_func = partial(afapply,vcf=vcf_reader,alt_freq=ref_alt_freq,chr_dict=chr_dict)
    status_inferred = sumstats.apply(map_func,axis=1)
    sumstats[column_name] = status_inferred.values
    sumstats[column_name]=sumstats[column_name].astype("float") 
    return sumstats

def check_daf(chr,start,end,ref,alt,eaf,vcf_reader,alt_freq,chr_dict=None):
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    
    for record in chr_seq:
        if record.pos==end:
            if record.ref==ref and (alt in record.alts):
                return eaf - record.info[alt_freq][0]
    return np.nan
################################################################################################################

def paralleleinferaf(sumstats,ref_infer,ref_alt_freq=None,n_cores=1, chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",status="STATUS",chr_dict=None,force=False, verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "infer EAF using reference VCF ALT frequency"
    _end_line = "inferring EAF using reference VCF ALT frequency"
    _start_cols = [chr,pos,ref,alt,status]
    _start_function = ".infer_af()"
    _must_args ={"ref_alt_freq":ref_alt_freq}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            n_cores=n_cores,
                            ref_vcf=ref_infer,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################
    chr_dict = auto_check_vcf_chr_dict(ref_infer, chr_dict, verbose, log)
    
    if eaf not in sumstats.columns:
        sumstats[eaf]=np.nan
    
    prenumber = sum(sumstats[eaf].isna())

    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)   
        if not force:
            good_chrpos =  sumstats[status].str.match(r'\w\w\w[0]\w\w\w', case=False, flags=0, na=False)  
        log.write(" -Checking variants:", sum(good_chrpos),verbose=verbose) 
    
    ########################  
        if sum(sumstats[eaf].isna())<10000: 
            n_cores=1       
        #df_split = np.array_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], n_cores)
        df_split = _df_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], n_cores)
        pool = Pool(n_cores)
        map_func = partial(inferaf,chr=chr,pos=pos,ref=ref,alt=alt,eaf=eaf,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,chr_dict=chr_dict) 
        sumstats.loc[good_chrpos,[eaf]] = pd.concat(pool.map(map_func,df_split))
        pool.close()
        pool.join()
    ###########################
        
        afternumber = sum(sumstats[eaf].isna())
        log.write(" -Inferred EAF for {} variants.".format(prenumber - afternumber),verbose=verbose) 
        log.write(" -EAF is still missing for {} variants.".format(afternumber),verbose=verbose) 
    
    finished(log,verbose,_end_line)
    return sumstats

def inferaf(sumstats,ref_infer,ref_alt_freq=None,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",chr_dict=None):
    #vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    vcf_reader = VariantFile(ref_infer)
    def afapply(x,vcf,alt_freq,chr_dict):
            return infer_af(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],vcf_reader,ref_alt_freq,chr_dict)
    map_func = partial(afapply,vcf=vcf_reader,alt_freq=ref_alt_freq,chr_dict=chr_dict)
    status_inferred = sumstats.apply(map_func,axis=1)
    sumstats[eaf] = status_inferred.values
    sumstats[eaf]=sumstats[eaf].astype("float") 
    return sumstats

def infer_af(chr,start,end,ref,alt,vcf_reader,alt_freq,chr_dict=None):
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    
    for record in chr_seq:
        if record.pos==end:
            if record.ref==ref and (alt in record.alts):
                return record.info[alt_freq][0]
            elif record.ref==alt and (ref in record.alts):
                return 1 - record.info[alt_freq][0]
    return np.nan
##############################################################################################################################################################################################

################################################################################################################

def _paralleleinferafwithmaf(sumstats,ref_infer,ref_alt_freq=None,n_cores=1, chr="CHR",pos="POS",ref="NEA",alt="EA",
                            eaf="EAF",maf="MAF",ref_eaf="_REF_EAF",status="STATUS",chr_dict=None,force=False, verbose=True,log=Log()):
    ##start function with col checking##########################################################
    _start_line = "infer sumstats EAF from sumstats MAF using reference VCF ALT frequency"
    _end_line = "inferring sumstats EAF from sumstats MAF using reference VCF ALT frequency"
    _start_cols = [chr,pos,ref,alt,status]
    _start_function = ".infer_af()"
    _must_args ={"ref_alt_freq":ref_alt_freq}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            n_cores=n_cores,
                            ref_vcf=ref_infer,
                            **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################
    chr_dict = auto_check_vcf_chr_dict(ref_infer, chr_dict, verbose, log)
    
    if eaf not in sumstats.columns:
        sumstats[eaf]=np.nan
    if ref_eaf not in sumstats.columns:
        sumstats[ref_eaf]=np.nan

    prenumber = sum(sumstats[eaf].isna())

    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        log.write(" -Field for alternative allele frequency in VCF INFO: {}".format(ref_alt_freq), verbose=verbose)   
        if not force:
            good_chrpos =  sumstats[status].str.match(r'\w\w\w[0]\w\w\w', case=False, flags=0, na=False)  
        log.write(" -Checking variants:", sum(good_chrpos),verbose=verbose) 
    
        ########################  
        #extract ref af
        if sum(sumstats[eaf].isna())<10000: 
            n_cores=1       
        #df_split = np.array_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], n_cores)
        df_split = _df_split(sumstats.loc[good_chrpos,[chr,pos,ref,alt]], n_cores)
        pool = Pool(n_cores)
        map_func = partial(inferaf,chr=chr,pos=pos,ref=ref,alt=alt,eaf=ref_eaf,ref_infer=ref_infer,ref_alt_freq=ref_alt_freq,chr_dict=chr_dict) 
        sumstats.loc[good_chrpos,[ref_eaf]] = pd.concat(pool.map(map_func,df_split))
        pool.close()
        pool.join()
        
        ###########################
        # infer sumstats EAF 
        # based on sumstats MAF and reference EAF
        is_filpped = ((sumstats[ref_eaf]>=0.5)&(sumstats[maf]<=0.5)) |((sumstats[ref_eaf]<0.5)&(sumstats[maf]>0.5)) 
        sumstats[eaf] = sumstats[maf]
        log.write(" -Flipping MAF to obtain EAF for {} variants".format(sum(is_filpped)),verbose=verbose) 
        sumstats.loc[is_filpped,eaf] = 1 - sumstats.loc[is_filpped,maf]

        ###########################    
        afternumber = sum(sumstats[eaf].isna())
        log.write(" -Inferred EAF for {} variants.".format(prenumber - afternumber),verbose=verbose) 
        log.write(" -EAF is still missing for {} variants.".format(afternumber),verbose=verbose) 
        sumstats = sumstats.drop(columns=[ref_eaf])
    
    finished(log,verbose,_end_line)
    return sumstats

def inferaf(sumstats,ref_infer,ref_alt_freq=None,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",chr_dict=None):
    #vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    vcf_reader = VariantFile(ref_infer)
    def afapply(x,vcf,alt_freq,chr_dict):
            return infer_af(x.iloc[0],x.iloc[1]-1,x.iloc[1],x.iloc[2],x.iloc[3],vcf_reader,ref_alt_freq,chr_dict)
    map_func = partial(afapply,vcf=vcf_reader,alt_freq=ref_alt_freq,chr_dict=chr_dict)
    status_inferred = sumstats.apply(map_func,axis=1)
    sumstats[eaf] = status_inferred.values
    sumstats[eaf]=sumstats[eaf].astype("float") 
    return sumstats

def infer_af(chr,start,end,ref,alt,vcf_reader,alt_freq,chr_dict=None):
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    
    for record in chr_seq:
        if record.pos==end:
            if record.ref==ref and (alt in record.alts):
                return record.info[alt_freq][0]
            elif record.ref==alt and (ref in record.alts):
                return 1 - record.info[alt_freq][0]
    return np.nan

##############################################################################################################################################################################################
def auto_check_vcf_chr_dict(vcf_path, vcf_chr_dict, verbose, log):    
    if vcf_path is not None:
        if vcf_chr_dict is None:
            log.write(" -Checking chromosome notations in VCF/BCF files..." ,verbose=verbose)  
            vcf_chr_dict = check_vcf_chr_NC(vcf_path, log, verbose)
            if vcf_chr_dict is not None:
                return vcf_chr_dict
            log.write(" -Checking prefix for chromosomes in VCF/BCF files..." ,verbose=verbose)  
            prefix = check_vcf_chr_prefix(vcf_path, log,verbose) 
            if prefix is not None:
                log.write(" -Prefix for chromosomes: ",prefix) 
                vcf_chr_dict = get_number_to_chr(prefix=prefix)
            else:
                log.write(" -No prefix for chromosomes in the VCF/BCF files." ,verbose=verbose)  
                vcf_chr_dict = get_number_to_chr()
    return vcf_chr_dict

def check_vcf_chr_prefix(vcf_bcf_path,log,verbose):
    vcf_bcf = VariantFile(vcf_bcf_path)
    for i in list(vcf_bcf.header.contigs):
        m = re.search('(chr|Chr|CHR)([0-9xXyYmM]+)', i)
        if m is not None:
            return m.group(1)
    else:
        return None

def check_vcf_chr_NC(vcf_bcf_path,log,verbose):
    vcf_bcf = VariantFile(vcf_bcf_path)
    for i in list(vcf_bcf.header.contigs):
        if i in get_number_to_NC(build="19").values():
            log.write("  -RefSeq ID detected (hg19) in VCF/BCF...",verbose=verbose) 
            return get_number_to_NC(build="19")
        elif i in get_number_to_NC(build="38").values():
            log.write("  -RefSeq ID detected (hg38) in VCF/BCF...",verbose=verbose) 
            return get_number_to_NC(build="38")
    else:
        return None

