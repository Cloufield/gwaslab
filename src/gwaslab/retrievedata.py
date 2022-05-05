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
def rsidtochrpos(sumstats,
         path="",
         rsid="rsID", chrom="CHR",pos="POS",ref_snp="SNP",ref_chr="CHR",ref_pos="POS",  
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
    
def merge_chrpos(sumstats_part,path):
  
    group=str(sumstats_part["group"].mode(dropna=True)[0])
    if group in [str(i) for i in range(75)]:
        to_merge=pd.read_hdf(path, key="part"+str(group))
        sumstats_part.update(to_merge)
    return sumstats_part



def parallelrsidtochrpos(sumstats, rsid="rsID", chrom="CHR",pos="POS", path=None,
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
    
    sumstats_rs = pd.concat(pool.map(partial(merge_chrpos,path=path),df_split))
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
    #xxx1x /  ----->  match
    #xxx2x /  ----->  8 + flipped => 2
    #xxx4x /   ----->  both allele on genome + unable to distinguish
    #xxxx5 /  ----->  reverse_complementary + both allele on genome + unable to distinguish
    #xxxx6 /  ----->  reverse_complementary 
    #xxxx7 / ------>  reverse_complementary + flipped
    #xxx8x /  ----->  flipped
    #xxx9x / ----->  not on ref genome
    status_pre=row[3][:3]
    status_end=row[3][4:]
    ## nea == ref
    if row[2] == record[row[0]-1: row[0]+len(row[2])-1].seq:
        ## ea == ref
        if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
            ## len(nea) >len(ea):
            if len(row[2])!=len(row[1]):
                # indels both on ref, unable to identify
                return status_pre+"4"+status_end 
        else:
            #nea == ref & ea != ref
            return status_pre+"1"+status_end 
    ## nea!=ref
    else:
        # ea == ref_seq -> need to flip, assign 8
        if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
            return status_pre+"8"+status_end 
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
                        return status_pre+"5"+status_end 
                else:
                    return status_pre+"6"+status_end 
            else:
                # ea == ref_seq -> need to flip, assign 8
                if row[1] == record[row[0]-1: row[0]+len(row[1])-1].seq:
                    return status_pre+"7"+status_end 
            # ea !=ref
            return status_pre+"9"+status_end
        
    
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
    records = SeqIO.parse(ref_path, "fasta")
    if verbose:  log.write(" -Checking records: ", end="")  
    for i in range(1,23):
        record = next(records)
        chromlist = gl.get_chr_list()
        record_chr = str(record.id).strip("chrCHR").upper()
        if record_chr in chromlist:
            if verbose:  log.write(record_chr," ", end="",show_time=False)  
            sumstats.loc[sumstats[chrom]==str(i),status] = sumstats.loc[sumstats[chrom] == str(i),[pos,ea,nea,status]].apply(lambda x:check_status(x,record),axis=1)
    if verbose:  log.write("\n",end="",show_time=False) 
    sumstats.loc[:,status] = sumstats.loc[:,status].astype("string")
    status_1=len(sumstats.loc[sumstats["STATUS"].str.match("[0-9][0-9][0-9][1][0-9]"),:])
    status_8=len(sumstats.loc[sumstats["STATUS"].str.match("[0-9][0-9][0-9][8][0-9]"),:])
    status_6=len(sumstats.loc[sumstats["STATUS"].str.match("[0-9][0-9][0-9][6][0-9]"),:])
    status_7=len(sumstats.loc[sumstats["STATUS"].str.match("[0-9][0-9][0-9][7][0-9]"),:])
    status_9=len(sumstats.loc[sumstats["STATUS"].str.match("[0-9][0-9][0-9][9][0-9]"),:])
    if verbose: log.write(" -Variants allele on given reference sequence : ",status_1)
    if verbose: log.write(" -Variants flipped : ",status_8)
    if verbose: log.write(" -Variants inferred reverse_complement : ",status_6)
    if verbose: log.write(" -Variants inferred reverse_complement_flipped : ",status_7)
    if verbose: log.write(" -Variants not on given reference sequence : ",status_9)
    if remove is True:
        sumstats = sumstats.loc[~sumstats["STATUS"].str.match("[0-9][0-9][0-9][9][0-9][0-9]"),:]
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

def assign_rsid_single(sumstats,path,rsid="rsID",chr="CHR",pos="POS",ref="NEA",alt="EA",overwrite=False,chr_dict=None):
    vcf_reader = vcf.Reader(open(path, 'rb'))
    no_rsID = sumstats[rsid].isna()| ~sumstats[rsid].str.match(r'rs([0-9]+)', case=False, flags=0, na=False)
    if overwrite is True: no_rsID = sumstats[rsid].isna() | sumstats[rsid].notnull()
    rsID = sumstats.loc[no_rsID,[chr,pos,ref,alt]].apply(lambda x:chrposref_rsid(x[0],x[1]-1,x[1],x[2],x[3],vcf_reader,chr_dict),axis=1)
    sumstats.loc[no_rsID,rsid] = rsID.values
    return sumstats

def parallelizeassignrsid(sumstats, path,rsid="rsID",chr="CHR",pos="POS",ref="NEA",alt="EA",n_cores=1,overwrite=False,verbose=True,log=gl.Log(),chr_dict=None):
    if verbose: log.write("Start to assign rsID ...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose: log.write(" -CPU Cores to use :",n_cores)
    if verbose: log.write(" -Reference VCF file:", path)
    if verbose: log.write(" -Assigning rsID based on chr:pos:ref:alt...")
    if rsid not in sumstats.columns:
        sumstats["rsID"]=pd.Series(dtype="string")
    total_number= len(sumstats)
    pre_number = sum(~sumstats["rsID"].isna())
    
    if overwrite is False:
        if sum(sumstats["rsID"].isna())<10000: n_cores=1
    func=assign_rsid_single
    df_split = np.array_split(sumstats, n_cores)
    pool = Pool(n_cores)
    #df = pd.concat(pool.starmap(func, df_split))
    if sum(sumstats["rsID"].isna())>0:
        sumstats = pd.concat(pool.map(partial(func,rsid=rsid,chr=chr,pos=pos,ref=ref,overwrite=overwrite,alt=alt,path=path,chr_dict=chr_dict),df_split))
    pool.close()
    pool.join()
    after_number = sum(~sumstats["rsID"].isna())
    if verbose: log.write(" -rsID Annotation for "+str(total_number - after_number) +" were not found!")
    if verbose: log.write(" -Annotated "+str(after_number - pre_number) +" rsID successfully!")
    return sumstats
#################################################################################################################################################


def check_strand_status(chr,start,end,ref,alt,eaf,vcf_reader,alt_freq,status,chr_dict=None):
    ### 1 : not palindromic
    ### 2 : palindromic +strand 
    ### 3 : palindromic -strand -> need to flip -> flipped
    ### 8 : palindromic -strand -> need to flip
    ### 9 : no ref data
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    status_pre=status[:4]
    status_end=""
    for record in chr_seq:
        if alt_freq not in record.INFO.keys():
            continue
        if record.POS==end and record.REF==ref and (alt in record.ALT):
            if (record.INFO[alt_freq][0]<0.5) and (eaf<0.5):
                return status_pre+"2"+status_end
            elif (record.INFO[alt_freq][0]>0.5) and (eaf>0.5):
                return status_pre+"2"+status_end
            else:
                return status_pre+"8"+status_end
    return status_pre+"9"+status_end


def check_unkonwn_indel(chr,start,end,ref,alt,vcf_reader,alt_freq,status,chr_dict=None):
    ### input : unknown indel, both on genome (xx1[45]x)
    ### 4 no flip
    ### 5 unknown indel,fixed   (6->5)
    ### 6 flip
    ### 9 noinfo or not matching
    if chr_dict is not None: chr=chr_dict[chr]
    chr_seq = vcf_reader.fetch(chr,start,end)
    status_pre=status[:4]
    status_end=""
    for record in chr_seq:
        if record.POS==end and record.REF==ref and (alt in record.ALT):
            return status_pre+"4"+status_end
        elif record.POS==end and record.REF==alt and (ref in record.ALT):
            return status_pre+"6"+status_end
    return status_pre+"9"+status_end

                                               
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

def inferstrand(sumstats,ref_infer,ref_alt_freq=None,maf_threshold=0.43,chr="CHR",pos="POS",ref="NEA",alt="EA",eaf="EAF",status="STATUS",chr_dict=None,verbose=True,log=gl.Log()):
    if verbose: log.write("Start to infer strand for palindromic SNPs...")
    if verbose: log.write(" -Current Dataframe shape :",len(sumstats)," x ", len(sumstats.columns))   
    if verbose: log.write(" -Reference vcf file:", ref_infer)   
     
    vcf_reader = vcf.Reader(open(ref_infer, 'rb'))
    
    if not ((chr in sumstats.columns) and (pos in sumstats.columns) and (ref in sumstats.columns) and (alt in sumstats.columns) and (status in sumstats.columns)):
        raise ValueError("Not enough information: CHR, POS, NEA , EA, ALT, STATUS...")
    
    if ref_alt_freq is not None:
        
        
        if verbose: log.write(" -Alternative allele frequency in INFO:", ref_alt_freq)  
        ## checking
        palindromic = is_palindromic(sumstats,a1="EA",a2="NEA")                   
        ##
        good_chrpos = ~(sumstats[chr].isna()|sumstats[pos].isna())
        ##not palindromic
        sumstats.loc[~palindromic,status] = sumstats.loc[~palindromic,status].apply(lambda x:change_status(x,5,"1")).astype("string") 

        if verbose: log.write(" -Identified ", sum(palindromic)," palindromic SNPs...")
        palindromic = palindromic&good_chrpos
        
        maf_can_infer   = (sumstats.loc[:,"EAF"] < maf_threshold) | (sumstats.loc[:,"EAF"] > 1 - maf_threshold) 
        sumstats.loc[palindromic&(~maf_can_infer),status] = sumstats.loc[palindromic&(~maf_can_infer),status].apply(lambda x:change_status(x,5,"7")).astype("string")
        if verbose: log.write(" -After filtering by MAF< ", maf_threshold ," , the strand of ", sum(palindromic & maf_can_infer)," palindromic SNPs will be inferred...")
        if sum(palindromic & maf_can_infer)>0:
            status_inderred = sumstats.loc[(palindromic & maf_can_infer),
                                       [chr,pos,ref,alt,eaf,status]].apply(
                                        lambda x:check_strand_status(x[0],x[1]-1,x[1],x[2],x[3],x[4],vcf_reader,ref_alt_freq,x[5],chr_dict),axis=1)
            sumstats.loc[(palindromic & maf_can_infer),status] = status_inderred.values
            
    
    ### unknow_indel
    unknow_indel = sumstats[status].str.match(r'[\w][\w][1][45][\w]', case=False, flags=0, na=False)   
    if sum(unknow_indel)>0:
        if verbose: log.write(" -Identified ", sum(unknow_indel)," indistinguishable Indels...")
        #if verbose: log.write(" -After filtering by MAF< ", maf_threshold ," , the strand of ", sum(unknow_indel & maf_can_infer ),"  SNPs not on reference genome will be inferred...")
        if verbose: log.write(" -Indistinguishable Indels will be inferred from reference vcf ref and alt...")
        status_inderred = sumstats.loc[unknow_indel,
                                       [chr,pos,ref,alt,status]].apply(
                                        lambda x:check_unkonwn_indel(x[0],x[1]-1,x[1],x[2],x[3],vcf_reader,ref_alt_freq,x[4],chr_dict),axis=1)
        sumstats.loc[unknow_indel,status] = status_inderred.values                    
                          
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