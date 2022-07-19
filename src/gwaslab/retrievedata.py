import gwaslab as gl
import pandas as pd
import numpy as np
import vcf
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool
from functools import partial
import re
import os
import psutil
#rsidtochrpos
#checkref
#parallelizeassignrsid
#inferstrand
#20220428 
#################################################################################################################

###~!!!!
def rsidtochrpos(sumstats,
         path="",
         rsid="rsID", chrom="CHR",pos="POS",ref_snp="SNP",ref_chr="CHR",ref_pos="POS", build="19",
              overwrite=False,remove=False,chunksize=5000000,verbose=True,log=gl.Log()):
    '''
    assign chr:pos based on rsID
    '''
    if verbose:  log.write("Start to update chromosome and position information based on rsID...")  
    if verbose:  log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose:  log.write(" -rsID dictionary file: "+ path)  
    dic_chuncks = pd.read_csv(path,"\s+",usecols=[ref_snp,ref_chr,ref_pos],
                      chunksize=chunksize,index_col=ref_snp,
                      dtype={ref_snp:"string",ref_chr:"string",ref_pos:"Int64"})
    sumstats = sumstats.set_index(rsid)
    
    #if chr or pos columns not in sumstats
    if chrom not in sumstats.columns:
        sumstats[chrom] =pd.Series(dtype="string")
    if pos not in sumstats.columns:    
        sumstats[pos] =pd.Series(dtype="string")
    
    if verbose:  log.write(" -Setting block size: ",chunksize)
    if verbose:  log.write(" -Loading block: ",end="")     
    for i,dic in enumerate(dic_chuncks):    
        log.write(i," ",end=" ",show_time=False)  
        dic = dic.rename(index={ref_snp:rsid})
        dic = dic.rename(columns={ref_chr:chrom,ref_pos:pos})  
        sumstats.update(dic,overwrite=overwrite)
    
    if verbose:  log.write("\n",end="",show_time=False) 
    sumstats = sumstats.reset_index()
    sumstats = sumstats.rename(columns = {'index':'rsID'})
    if verbose:  log.write(" -Updating CHR and POS finished.Start to re-fixing CHR and POS... ")
    sumstats = fixchr(sumstats,remove=remove,verbose=verbose)
    sumstats = fixpos(sumstats,remove=remove,verbose=verbose)
    return sumstats
####################################################################################################################
def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem
    
def merge_chrpos(sumstats_part,path,build,status):
    group=str(sumstats_part["group"].mode(dropna=True)[0])
    if group in [str(i) for i in range(75)]:
        to_merge=pd.read_hdf(path, key="part"+str(group))
        is_chrpos_fixable = sumstats_part.index.isin(to_merge.index)
        sumstats_part.loc[is_chrpos_fixable,status] = vchange_status(sumstats_part.loc[is_chrpos_fixable, status],  1,"139",3*build[0])
        sumstats_part.loc[is_chrpos_fixable,status] = vchange_status(sumstats_part.loc[is_chrpos_fixable, status],  2,"987",3*build[1])
        sumstats_part.update(to_merge)
    return sumstats_part



def parallelrsidtochrpos(sumstats, rsid="rsID", chrom="CHR",pos="POS", path=None,build="19",status="STATUS",
                         n_cores=4,block_size=20000000,verbose=True,log=gl.Log()):
    if verbose:  log.write(" -Start to assign CHR and POS using rsIDs... ")
    if path is None:
        raise ValueError("Please provide path to hdf5 file.")
    
    sumstats["rsn"] = pd.to_numeric(sumstats[rsid].str.strip("rs"),errors="coerce").astype("Int64")
    
    if verbose:  log.write(" -Source hdf5 file: ",path)
    if verbose:  log.write(" -Cores to use : ",n_cores)
    if verbose:  log.write(" -Blocksize (make sure it is the same as hdf5 file ): ",block_size)
    
    input_columns= sumstats.columns
    sumstats_nonrs = sumstats.loc[sumstats["rsn"].isna(),:].copy()
    sumstats_rs  = sumstats.loc[sumstats["rsn"].notnull(),:].copy()

    del sumstats
    
    if verbose:  log.write(" -Non-Valid rsIDs: ",len(sumstats_nonrs))
    if verbose:  log.write(" -Valid rsIDs: ",len(sumstats_rs))
    
    sumstats_rs.loc[:,"group"]= sumstats_rs.loc[:,"rsn"]//block_size
    if verbose:  log.write(" -Groups : ",set(sumstats_rs.loc[:,"group"].unique()))
    sumstats_rs = sumstats_rs.set_index("rsn")
 
    pool = Pool(n_cores)
    if chrom not in input_columns:
        if verbose:  log.write(" -Initiating CHR ... ")
        sumstats_rs[chrom]=pd.Series(dtype="string") 
        
    if pos not in input_columns:
        if verbose:  log.write(" -Initiating POS ... ")
        sumstats_rs[pos]=pd.Series(dtype="Int64") 
    
    df_split=[y for x, y in sumstats_rs.groupby('group', as_index=False)]
    if verbose:  log.write(" -Divided into groups: ",len(df_split))
    
    sumstats_rs = pd.concat(pool.map(partial(merge_chrpos,path=path,build=build,status=status),df_split))
    del df_split
   
    if verbose:  log.write(" -Merging group data... ")
    
    sumstats_rs = sumstats_rs.drop(columns=["group"])
    sumstats_nonrs = sumstats_nonrs.drop(columns=["rsn"])
  
    if verbose:  log.write(" -Append data... ")
    sumstats = pd.concat([sumstats_rs,sumstats_nonrs],ignore_index=True)
    del sumstats_rs
    del sumstats_nonrs
    sumstats = gl.fixchr(sumstats,verbose=False)
    sumstats = gl.fixpos(sumstats,verbose=False)
    pool.close()

    pool.join()
    return sumstats
####################################################################################################################
#20220426 check if mom-effect allele is aligned with reference genome 
def check_status(row,record):
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
    
    status_pre=row[3][:5]
    status_end=row[3][6:]
    
    ## nea == ref
    if row[2] == record[row[0]-1: row[0]+len(row[2])-1].seq:
        ## ea == ref
        if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
            ## len(nea) >len(ea):
            if len(row[2])!=len(row[1]):
                # indels both on ref, unable to identify
                return status_pre+"6"+status_end 
        else:
            #nea == ref & ea != ref
            return status_pre+"0"+status_end 
    ## nea!=ref
    else:
        # ea == ref_seq -> need to flip
        if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
            return status_pre+"3"+status_end 
        # ea !=ref
        else:
            #_reverse_complementary
            row[1] = get_reverse_complementary_allele(row[1])
            row[2] = get_reverse_complementary_allele(row[2])
            ## nea == ref
            if row[2] == record[row[0]-1: row[0]+len(row[2])-1].seq:
                ## ea == ref
                if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
                    ## len(nea) >len(ea):
                    if len(row[2])!=len(row[1]):
                        return status_pre+"8"+status_end  # indel reverse complementary
                else:
                    return status_pre+"4"+status_end 
            else:
                # ea == ref_seq -> need to flip
                if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
                    return status_pre+"5"+status_end 
            # ea !=ref
            return status_pre+"8"+status_end
        
    
def get_reverse_complementary_allele(a1,a2):
    dic = str.maketrans({
       "A":"T",
       "T":"A",
       "C":"G",
       "G":"C"})
    return a1[::-1].translate(dic),a2[::-1].translate(dic)

def checkref(sumstats,ref_path,chrom="CHR",pos="POS",ea="EA",nea="NEA",status="STATUS",remove=False,verbose=True,log=gl.Log()):
    if verbose: log.write("Start to check if NEA is aligned with reference sequence...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns)) 
    if verbose:  log.write(" -Reference genome fasta file: "+ ref_path)  
    if verbose:  log.write(" -Checking records: ", end="")  
        
    records = SeqIO.parse(ref_path, "fasta")
    for i in range(1,23):
        record = next(records)
        chromlist = gl.get_chr_list()
        record_chr = str(record.id).strip("chrCHR").upper()
        if record_chr in chromlist:
            if verbose:  log.write(record_chr," ", end="",show_time=False)  
            to_check_ref = (sumstats[chrom]==str(i)) & (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna())
            sumstats.loc[to_check_ref,status] = sumstats.loc[to_check_ref,[pos,ea,nea,status]].apply(lambda x:check_status(x,record),axis=1)
    
    if verbose:  log.write("\n",end="",show_time=False) 
        
    sumstats.loc[:,status] = sumstats.loc[:,status].astype("string")
    available_to_check =sum( (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna()))
    status_0=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[0]\w", case=False, flags=0, na=False))
    status_3=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[3]\w", case=False, flags=0, na=False))
    status_4=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[4]\w", case=False, flags=0, na=False))
    status_5=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[5]\w", case=False, flags=0, na=False))
    status_6=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[6]\w", case=False, flags=0, na=False))
    #status_7=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[7]\w", case=False, flags=0, na=False))
    status_8=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w", case=False, flags=0, na=False))
    
    if verbose: log.write(" -Variants allele on given reference sequence : ",status_0)
    if verbose: log.write(" -Variants flipped : ",status_3)
    raw_matching_rate = (status_3+status_0)/available_to_check
    flip_rate = status_3/available_to_check
    if verbose: log.write("  -Raw Matching rate : ","{:.2f}%".format(raw_matching_rate*100))
    if raw_matching_rate <0.8:
        if verbose: log.write("  -!!!Warning, matching rate is low, please check if the right reference genome is used.")
    if flip_rate > 0.85 :
        if verbose: log.write("  -Flipping variants rate > 0.85, it is likely that the EA is aligned with REF in the original dataset.")
    if verbose: log.write(" -Variants inferred reverse_complement : ",status_4)
    if verbose: log.write(" -Variants inferred reverse_complement_flipped : ",status_5)
    if verbose: log.write(" -Both allele on genome + unable to distinguish : ",status_6)
    #if verbose: log.write(" -Reverse_complementary + both allele on genome + unable to distinguish: ",status_7)
    if verbose: log.write(" -Variants not on given reference sequence : ",status_8)
    
    if remove is True:
        sumstats = sumstats.loc[~sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w"),:]
        if verbose: log.write(" -Variants not on given reference sequence were removed.")
    return sumstats

#######################################################################################################################################

#20220426
def chrposref_rsid(chr,start,end,ref,alt,vcf_reader,chr_dict=None):
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    for record in chr_seq:
        if record.POS==end and record.REF==ref and (alt in record.ALT):
            return record.ID
        elif record.POS==end and (ref in record.ALT) and record.REF==alt:
            return record.ID
    return pd.NA

def assign_rsid_single(sumstats,path,rsid="rsID",chr="CHR",pos="POS",ref="NEA",alt="EA",status="STATUS",overwrite="empty",chr_dict=None):
    vcf_reader = vcf.Reader(open(path, 'rb'))
    standardized_normalized = sumstats["STATUS"].str.match("\w\w\w[0][01234]\w\w", case=False, flags=0, na=False)
    if overwrite=="all":
        no_rsID = standardized_normalized
    if overwrite=="invalid":
        no_rsID = (~sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)) & standardized_normalized
    if overwrite=="empty":
        no_rsID = sumstats[rsid].isna()& standardized_normalized
    if sum(no_rsID)>0:
        rsID = sumstats.loc[no_rsID,[chr,pos,ref,alt]].apply(lambda x:chrposref_rsid(x[0],x[1]-1,x[1],x[2],x[3],vcf_reader,chr_dict),axis=1)
        sumstats.loc[no_rsID,rsid] = rsID.values
    return sumstats

def parallelizeassignrsid(sumstats, path ,rsid="rsID",chr="CHR",pos="POS",ref="NEA",alt="EA",status="STATUS",n_cores=1,overwrite="empty",verbose=True,log=gl.Log(),chr_dict=None):
    '''
    overwrite mode : 
    all ,    overwrite rsid for all availalbe rsid 
    invalid,  only assign rsid for variants with invalid rsid
    empty    only assign rsid for variants with na rsid
    '''  
    if verbose: log.write("Start to assign rsID ...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose: log.write(" -CPU Cores to use :",n_cores)
    if verbose: log.write(" -Reference VCF file:", path)
    if verbose: log.write(" -Assigning rsID based on chr:pos:ref:alt...")
    
    if rsid not in sumstats.columns:
        sumstats["rsID"]=pd.Series(dtype="string")
    
    total_number= len(sumstats)
    pre_number = sum(~sumstats["rsID"].isna())
    
    if overwrite == "empty":
        if sum(sumstats["rsID"].isna())<10000: n_cores=1
    
    ##########################        
    func=assign_rsid_single
    df_split = np.array_split(sumstats, n_cores)
    pool = Pool(n_cores)
    #df = pd.concat(pool.starmap(func, df_split))
    if sum(sumstats["rsID"].isna())>0 or overwrite!="empty":
        sumstats = pd.concat(pool.map(partial(func,rsid=rsid,chr=chr,pos=pos,ref=ref,overwrite=overwrite,status=status,alt=alt,path=path,chr_dict=chr_dict),df_split))
    pool.close()
    pool.join()
    ###########################
    
    after_number = sum(~sumstats["rsID"].isna())
    if verbose: log.write(" -rsID Annotation for "+str(total_number - after_number) +" need to be fixed!")
    if verbose: log.write(" -Annotated "+str(after_number - pre_number) +" rsID successfully!")
    return sumstats
#################################################################################################################################################


def check_strand_status(chr,start,end,ref,alt,eaf,vcf_reader,alt_freq,status,chr_dict=None):
    ### 0 : not palindromic
    ### 1 : palindromic +strand 
    ### 2 : palindromic -strand -> need to flip -> flipped
    ### 5 : palindromic -strand -> need to flip
    ### 8 : no ref data
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    status_pre=status[:6]
    status_end=""
    for record in chr_seq:
        if alt_freq not in record.INFO.keys():
            continue
        if record.POS==end and record.REF==ref and (alt in record.ALT):
            if (record.INFO[alt_freq][0]<0.5) and (eaf<0.5):
                return status_pre+"1"+status_end
            elif (record.INFO[alt_freq][0]>0.5) and (eaf>0.5):
                return status_pre+"1"+status_end
            else:
                return status_pre+"5"+status_end
    return status_pre+"8"+status_end


def check_unkonwn_indel(chr,start,end,ref,alt,vcf_reader,alt_freq,status,chr_dict=None):
    ### input : unknown indel, both on genome (xx1[45]x)
    ### 3 no flip
    ### 4 unknown indel,fixed   (6->5)
    ### 6 flip
    ### 9 noinfo or not matching
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    status_pre=status[:6]
    status_end=""
    for record in chr_seq:
        if record.POS==end and record.REF==ref and (alt in record.ALT):
            return status_pre+"3"+status_end
        elif record.POS==end and record.REF==alt and (ref in record.ALT):
            return status_pre+"6"+status_end
    return status_pre+"8"+status_end

                                               
#def get_reverse_complementary_allele(a):
#    dic = str.maketrans({
#       "A":"T",
#       "T":"A",
#       "C":"G",
#       "G":"C"})
#    return a[::-1].translate(dic)
                                                 
def is_palindromic(sumstats,a1="EA",a2="NEA"):
    gc= (sumstats[a1]=="G") & (sumstats[a2]=="C")
    cg= (sumstats[a1]=="C") & (sumstats[a2]=="G")
    at= (sumstats[a1]=="A") & (sumstats[a2]=="T")
    ta= (sumstats[a1]=="T") & (sumstats[a2]=="A")
    palindromic = gc | cg | at | ta 
    return palindromic

def inferstrand(sumstats,ref_infer,ref_alt_freq=None,maf_threshold=0.43,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",status="STATUS",chr_dict=None,verbose=True,log=gl.Log()):
    
    if verbose: log.write("Start to infer strand for palindromic SNPs...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose: log.write(" -Reference vcf file:", ref_infer)   
     
    vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    
    # check if the columns are complete
    if not ((chr in sumstats.columns) and (pos in sumstats.columns) and (ref in sumstats.columns) and (alt in sumstats.columns) and (status in sumstats.columns)):
        raise ValueError("Not enough information: CHR, POS, NEA , EA, ALT, STATUS...")
    
    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        if verbose: log.write(" -Alternative allele frequency in INFO:", ref_alt_freq)  
        
        ## checking \w\w\w\w[0]\w\w -> standardized and normalized snp
        palindromic = sumstats[status].str.match(r'\w\w\w\w[0]\w\w', case=False, flags=0, na=False) & is_palindromic(sumstats,a1="EA",a2="NEA")           
        
        not_palindromic_snp = sumstats[status].str.match(r'\w\w\w\w[0]\w\w', case=False, flags=0, na=False) & (~palindromic)
        ##
        good_chrpos =  sumstats[status].str.match(r'\w\w\w[0]\w\w\w', case=False, flags=0, na=False)  
        
        ##not palindromic
        sumstats.loc[not_palindromic_snp,status] = vchange_status(sumstats.loc[not_palindromic_snp,status], 7 ,"9","0")
        
        if verbose: log.write(" -Identified ", sum(palindromic)," palindromic SNPs...")
        
        palindromic = palindromic&good_chrpos
        
        #palindromic but can not infer
        maf_can_infer   = (sumstats.loc[:,"EAF"] < maf_threshold) | (sumstats.loc[:,"EAF"] > 1 - maf_threshold) 
        sumstats.loc[palindromic&(~maf_can_infer),status] = vchange_status(sumstats.loc[palindromic&(~maf_can_infer),status],7,"9","7")
        
        
        if verbose: log.write(" -After filtering by MAF< ", maf_threshold ," , the strand of ", sum(palindromic & maf_can_infer)," palindromic SNPs will be inferred...")
        if sum(palindromic & maf_can_infer)>0:
            status_inferred = sumstats.loc[(palindromic & maf_can_infer),
                                       [chr,pos,ref,alt,eaf,status]].apply(
                                        lambda x:check_strand_status(x[0],x[1]-1,x[1],x[2],x[3],x[4],vcf_reader,ref_alt_freq,x[5],chr_dict),axis=1)
            sumstats.loc[(palindromic & maf_can_infer),status] = status_inferred.values
        
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
        
        status0 =  sum( sumstats[status].str.match(r'\w\w\w\w\w\w[0]', case=False, flags=0, na=False)  )
        status1 =  sum( sumstats[status].str.match(r'\w\w\w\w\w\w[1]', case=False, flags=0, na=False)  )
        status5 =  sum( sumstats[status].str.match(r'\w\w\w\w\w\w[5]', case=False, flags=0, na=False)  )
        status7 =  sum( sumstats[status].str.match(r'\w\w\w\w\w\w[7]', case=False, flags=0, na=False)  )
        status8 =  sum( sumstats[status].str.match(r'\w\w\w\w\w[123][8]', case=False, flags=0, na=False)  )
        
        if verbose: log.write("  -Non-palindromic : ",status0)
        if verbose: log.write("  -Palindromic SNPs on + strand: ",status1)
        if verbose: log.write("  -Palindromic SNPs on - strand and need to be flipped:",status5)   
        if verbose: log.write("  -Palindromic SNPs with maf not availble to infer : ",status7)  
        if verbose: log.write("  -Palindromic SNPs with no macthes or no information : ",status8)   
    ### unknow_indel
    
    unknow_indel = sumstats[status].str.match(r'\w\w\w\w\w[6][89]', case=False, flags=0, na=False)   
    if sum(unknow_indel)>0:
        if verbose: log.write(" -Identified ", sum(unknow_indel)," indistinguishable Indels...")
        #if verbose: log.write(" -After filtering by MAF< ", maf_threshold ," , the strand of ", sum(unknow_indel & maf_can_infer ),"  SNPs not on reference genome will be inferred...")
        if verbose: log.write(" -Indistinguishable Indels will be inferred from reference vcf ref and alt...")
        status_inderred = sumstats.loc[unknow_indel, [chr,pos,ref,alt,status]].apply(lambda x:check_unkonwn_indel(x[0],x[1]-1,x[1],x[2],x[3],vcf_reader,ref_alt_freq,x[4],chr_dict),axis=1)
        sumstats.loc[unknow_indel,status] = status_inderred.values 
        
        status3 =  sum( sumstats[status].str.match(r'\w\w\w\w\w\w[3]', case=False, flags=0, na=False)  )
        status6 =  sum( sumstats[status].str.match(r'\w\w\w\w\w\w[6]', case=False, flags=0, na=False)  )
        status8 =  sum( sumstats[status].str.match(r'\w\w\w\w\w[6][8]', case=False, flags=0, na=False)  )
        
        if verbose: log.write("  -Indels ea/nea match reference : ",status3)
        if verbose: log.write("  -Indels ea/nea need to be flipped : ",status6)
        if verbose: log.write("  -Indels with no macthes or no information : ",status8)
                       
    return sumstats
################################################################################################################
def checkaf(sumstats,ref_infer,ref_alt_freq=None,maf_threshold=0.43,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",status="STATUS",chr_dict=None,verbose=True,log=gl.Log()):
        
    if verbose: log.write("Start to check the difference between EAF and refence vcf alt frequency ...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose: log.write(" -Reference vcf file:", ref_infer)   
    vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    
    # check if the columns are complete
    if not ((chr in sumstats.columns) and (pos in sumstats.columns) and (ref in sumstats.columns) and (alt in sumstats.columns) and (status in sumstats.columns)):
        raise ValueError("Not enough information: CHR, POS, NEA , EA, ALT, STATUS...")
    
    # ref_alt_freq INFO in vcf was provided
    if ref_alt_freq is not None:
        if verbose: log.write(" -Alternative allele frequency in INFO:", ref_alt_freq)  
        good_chrpos =  sumstats[status].str.match(r'\w\w\w[0]\w\w\w', case=False, flags=0, na=False)  
        if verbose: log.write(" -Checking variants:", sum(good_chrpos)) 
        sumstats["DAF"]=np.nan
        status_inferred = sumstats.loc[good_chrpos,[chr,pos,ref,alt,eaf]].apply(lambda x:check_daf(x[0],x[1]-1,x[1],x[2],x[3],x[4],vcf_reader,ref_alt_freq,chr_dict),axis=1)
        sumstats.loc[good_chrpos,"DAF"] = status_inferred.values
        sumstats.loc[:,"DAF"]=sumstats.loc[:,"DAF"].astype("float")
        if verbose: log.write(" - DAF min:", np.nanmax(sumstats.loc[:,"DAF"])) 
        if verbose: log.write(" - DAF max:", np.nanmin(sumstats.loc[:,"DAF"])) 
        if verbose: log.write(" - abs(DAF) min:", np.nanmax(np.abs(sumstats.loc[:,"DAF"]))) 
        if verbose: log.write(" - abs(DAF) max:", np.nanmin(np.abs(sumstats.loc[:,"DAF"])))
        if verbose: log.write(" - DAF sd:", np.nanstd(sumstats.loc[:,"DAF"])) 
        if verbose: log.write(" - abs(DAF) sd:", np.nanstd(np.abs(sumstats.loc[:,"DAF"]))) 
        
    return sumstats

def check_daf(chr,start,end,ref,alt,eaf,vcf_reader,alt_freq,chr_dict=None):
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    for record in chr_seq:
        if alt_freq not in record.INFO.keys():
            continue
        if record.POS==end and record.REF==ref and (alt in record.ALT):
            return eaf - record.INFO[alt_freq][0]
    return np.nan
################################################################################################################
def change_status(status_value,digit,to):
    if digit>1:
        status_pre=status_value[:digit-1]
    else:
        status_pre=""
        
    status_end=status_value[digit:]
    return status_pre+str(to)+status_end

################################################################################################################
def vchange_status(status,digit,before,after):
    for i in range(len(before)):
        if digit>1:
            pattern= (digit-1) * r'\w' + before[i] + (7 - digit)* r'\w'
        else:
            pattern=before[i]+r'\w\w\w\w\w\w'        
        
        to_change = status.str.match(pattern, case=False, flags=0, na=False)  
        if sum(to_change)>0:
            if digit>1:
                status_pre = status[to_change].str[:digit-1]
            else:
                status_pre = ""

            status_end=status[to_change].str[digit:]

            status[to_change] = status_pre+after[i]+status_end
    return status