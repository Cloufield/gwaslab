import pandas as pd
import yaml
import hashlib
from gwaslab.Log import Log
from gwaslab.CommonData import get_format_dict
from gwaslab.CommonData import get_number_to_chr
from gwaslab.version import gwaslab_info
from pysam import tabix_compress 
from pysam import tabix_index
from gwaslab.CommonData import get_formats_list
from datetime import datetime
from datetime import date
from gwaslab.preformat_input import print_format_info
# to vcf
# to fmt
    ## vcf
    ## vep
    ## bed
    ## annovar
    ## general : ldsc, plink, plink2, saige, regenie 
###################################################################################################################################################
    
###################################################################################################################################################
def tofmt(sumstats,
          meta,
          path=None,
          suffix=None,
          fmt=None,
          cols=[],
          xymt_number=False,
          xymt=["X","Y","MT"],
          build="19",
          chr_prefix="",
          ssfmeta=False,
          md5sum=False,
          bgzip=False,
          tabix=False,
          verbose=True,
          log=Log(),
          to_csvargs=None):
    
    if to_csvargs is None:
        to_csvargs=dict()
    
    if fmt in ["ssf"]: 
        xymt_number=True
        if "SNPID" in sumstats.columns:
            if verbose: log.write(' - Replacing SNPID separator from ":" to "_"...')
            sumstats["SNPID"] = sumstats["SNPID"].str.replace(":","_")
    if verbose: log.write(" - Start outputting sumstats in "+fmt+" format...")
    
    if "CHR" in sumstats.columns:
        if xymt_number is False and sumstats["CHR"].dtype in ["Int64","int"]:
            sumstats["CHR"]= sumstats["CHR"].map(get_number_to_chr(xymt=xymt,prefix=chr_prefix))
        elif chr_prefix is not None:
            sumstats["CHR"]= chr_prefix + sumstats["CHR"].astype("string")
    
    ### calculate meta data
    if "EAF" in sumstats.columns:
        min_maf = sumstats["EAF"].min()
    else:
        min_maf = "Unknown"
    
    if "N" in sumstats.columns:
        n_median =  sumstats["N"].median()
        n_max = sumstats["N"].max()
        n_min = sumstats["N"].min()
    else:
        n_median = "Unknown"
        n_max = "Unknown"
        n_min = "Unknown"


    if fmt=="bed":
        # bed-like format, 0-based, 
        # first 3 columns : chromosome, start, end
        # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
        is_snp = (sumstats["EA"].str.len() == sumstats["NEA"].str.len())
        is_insert = (sumstats["EA"].str.len()>1) &(sumstats["NEA"].str.len()==1)
        is_delete = (sumstats["EA"].str.len()==1) &(sumstats["NEA"].str.len()>1)
        
        if verbose: log.write(" -Number of SNPs :",sum(is_snp))
        if verbose: log.write(" -Number of Insertions :",sum(is_insert)) 
        if verbose: log.write(" -Number of Deletions :",sum(is_delete)) 
        
        if verbose: log.write(" -formatting to 0-based bed-like file...")
        # for snp  
        # start = pos - 1 ; end = pos
        # A/G
        # AT/CG
        sumstats.loc[is_snp,"START"]  = sumstats.loc[is_snp,"POS"]-1 
        sumstats.loc[is_snp,"END"]    = sumstats.loc[is_snp,"POS"]-1 + sumstats.loc[is_snp,"NEA"].str.len()
        sumstats.loc[is_snp,"NEA/EA"] = sumstats.loc[is_snp,"NEA"].astype("string")+"/"+sumstats.loc[is_snp,"EA"].astype("string")
        
        # for insertion
        # start = pos : end = pos
        # A/ATC -> -/TC
        sumstats.loc[is_insert,"START"]  = sumstats.loc[is_insert,"POS"]
        sumstats.loc[is_insert,"END"]    = sumstats.loc[is_insert,"POS"]
        sumstats.loc[is_insert,"NEA/EA"] = "-/"+sumstats.loc[is_insert,"EA"].str.slice(start=1)
        
        # for deletion 
        # start = pos - 1 +1; end = pos -1 +1+ len(Ref)
        # ATC/A -> TC/-
        sumstats.loc[is_delete,"START"]  = sumstats.loc[is_delete,"POS"]
        sumstats.loc[is_delete,"END"]    = sumstats.loc[is_delete,"POS"] + sumstats.loc[is_delete,"NEA"].str.len() - 1
        sumstats.loc[is_delete,"NEA/EA"] = sumstats.loc[is_delete,"NEA"].str.slice(start=1)+"/-"
        
        sumstats["STRAND"]="+"

        sumstats["START"] = sumstats["START"].astype("Int64")
        sumstats["END"] = sumstats["END"].astype("Int64")
        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"] + cols

        sumstats = sumstats.loc[:,ouput_cols ]
        path = path + "."+suffix
        if verbose: log.write(" -Output columns:",sumstats.columns)
        if verbose: log.write(" -Output path:",path) 
        
        sumstats.to_csv(path,sep="\t",index=None,header=None,**to_csvargs)
        #tabix_compress
        #tabix_index
        if bgzip is True:
            if verbose: log.write(" -bgzip compressing ...") 
            tabix_compress(path, path+".gz",force=True)
        if tabix is True:
            if verbose: log.write(" -tabix indexing...") 
            tabix_index(path+".gz" ,preset="bed",force=True) 
    ####################################################################################################################   
    elif fmt=="vep":
        # bed-like format, 1-based
        # first 6 columns : chromosome, start, end, allele, strand, identifier
        # https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html
        
        is_snp = (sumstats["EA"].str.len() == sumstats["NEA"].str.len())
        is_insert = (sumstats["EA"].str.len()>1) &(sumstats["NEA"].str.len()==1)
        is_delete = (sumstats["EA"].str.len()==1) &(sumstats["NEA"].str.len()>1)
        
        if verbose: log.write(" -Number of SNPs :",sum(is_snp)) 
        if verbose: log.write(" -Number of Insertions :",sum(is_insert)) 
        if verbose: log.write(" -Number of Deletions :",sum(is_delete)) 
            
        if verbose: log.write(" -formatting to 1-based bed-like file (for vep)...")
        # for snp  
        # start = pos ; end = pos
        sumstats.loc[is_snp,"START"]  = sumstats.loc[is_snp,"POS"] + (sumstats.loc[is_snp,"NEA"].str.len() - 1 )
        sumstats.loc[is_snp,"END"]    = sumstats.loc[is_snp,"POS"] + (sumstats.loc[is_snp,"NEA"].str.len() - 1 )
        sumstats.loc[is_snp,"NEA/EA"] = sumstats.loc[is_snp,"NEA"].astype("string")+"/"+sumstats.loc[is_snp,"EA"].astype("string")
        
        # for insertion
        # start = pos+1 ; end = pos
        # A/ATC -> -/TC
        sumstats.loc[is_insert,"START"]  = sumstats.loc[is_insert,"POS"] + 1
        sumstats.loc[is_insert,"END"]    = sumstats.loc[is_insert,"POS"]
        sumstats.loc[is_insert,"NEA/EA"] = "-/" + sumstats.loc[is_insert,"EA"].str.slice(start=1)
        
        # for deletion 
        # start = pos ; end = pos + len(Ref) -1
        # ATC/A -> TC/-
        sumstats.loc[is_delete,"START"]  = sumstats.loc[is_delete,"POS"] + 1
        sumstats.loc[is_delete,"END"]    = sumstats.loc[is_delete,"POS"] + (sumstats.loc[is_delete,"NEA"].str.len() -1)
        sumstats.loc[is_delete,"NEA/EA"] = sumstats.loc[is_delete,"NEA"].str.slice(start=1)+"/-"
           
        sumstats["STRAND"]="+"
        
        sumstats["START"] = sumstats["START"].astype("Int64")
        sumstats["END"] = sumstats["END"].astype("Int64")
        
        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"]+ cols
        sumstats = sumstats.loc[:,ouput_cols]
        path = path + "."+suffix+".gz"
        if verbose: log.write(" -Output columns:",sumstats.columns)
        if verbose: log.write(" -Output path:",path) 

        sumstats.to_csv(path,sep="\t",index=None,header=None,**to_csvargs)
        if md5sum is True: md5sum_file(path,log,verbose)
    ####################################################################################################################
    elif fmt=="annovar":
        # bed-like format, 1-based, 
        # first 3 columns : Chromosome ("chr" prefix is optional), Start, End, Reference Allelel, Alternative Allele
        # https://annovar.openbioinformatics.org/en/latest/user-guide/input/
        is_snp = (sumstats["EA"].str.len() == sumstats["NEA"].str.len())
        is_insert = (sumstats["EA"].str.len()>1) &(sumstats["NEA"].str.len()==1)
        is_delete = (sumstats["EA"].str.len()==1) &(sumstats["NEA"].str.len()>1)
        
        if verbose: log.write(" -Number of SNPs :",sum(is_snp))
        if verbose: log.write(" -Number of Insertions :",sum(is_insert)) 
        if verbose: log.write(" -Number of Deletions :",sum(is_delete)) 
        
        if verbose: log.write(" -formatting to 1-based bed-like file...")
        # for snp  
        # start = pos ; end = pos
        # A/G
        # AT/CG
        sumstats.loc[is_snp,"START"]  = sumstats.loc[is_snp,"POS"]
        sumstats.loc[is_snp,"END"]    = sumstats.loc[is_snp,"POS"]-1 + sumstats.loc[is_snp,"NEA"].str.len()
        sumstats.loc[is_snp,"NEA_out"] = sumstats.loc[is_snp,"NEA"].astype("string")
        sumstats.loc[is_snp,"EA_out"] = sumstats.loc[is_snp,"EA"].astype("string")
        
        # for insertion
        # start = pos : end = pos
        # A/ATC -> -/TC
        sumstats.loc[is_insert,"START"]  = sumstats.loc[is_insert,"POS"]+1
        sumstats.loc[is_insert,"END"]   = sumstats.loc[is_insert,"POS"]+1
        sumstats.loc[is_insert,"NEA_out"] = "-"
        sumstats.loc[is_insert,"EA_out"] = sumstats.loc[is_insert,"EA"].str.slice(start=1)
        
        # for deletion 
        # start = pos - 1 +1; end = pos -1 +1+ len(Ref)
        # ATC/A -> TC/-
        sumstats.loc[is_delete,"START"] = sumstats.loc[is_delete,"POS"]
        sumstats.loc[is_delete,"END"]  = sumstats.loc[is_delete,"POS"]- 1 + sumstats.loc[is_delete,"NEA"].str.len() 
        sumstats.loc[is_delete,"NEA_out"] = sumstats.loc[is_delete,"NEA"].str.slice(start=1)
        sumstats.loc[is_delete,"EA_out"] = "-"
        
        sumstats["START"] = sumstats["START"].astype("Int64")
        sumstats["END"] = sumstats["END"].astype("Int64")

        ouput_cols=["CHR","START","END","NEA_out","EA_out","SNPID"]+ cols
        sumstats = sumstats.loc[:,ouput_cols]
        path = path + "."+suffix
        if verbose: log.write(" -Output columns:",sumstats.columns)
        if verbose: log.write(" -Output path:",path) 
        
        sumstats.to_csv(path,sep="\t",index=None,header=None,**to_csvargs)
        #tabix_compress
        #tabix_index
        if bgzip is True:
            if verbose: log.write(" -bgzip compressing ...") 
            tabix_compress(path, path+".gz",force=True)
            if md5sum is True: md5sum_file(path+".gz",log,verbose)
        if tabix is True:
            if verbose: log.write(" -tabix indexing...") 
            tabix_index(path+".gz" ,preset="bed",force=True) 
    ####################################################################################################################       
    elif fmt=="vcf":
        if verbose: log.write(" -"+fmt+" format will be loaded...")
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        #loading format data
        if verbose:             
            log.write(" -"+fmt+" format meta info:")   
            for key,value in meta_data.items():
                if key not in ["format_fixed_header","format_contig_19","format_contig_38"]:
                    log.write("  -",key," : ",value)
        
        # determine which ID to use
        if "rsID" in sumstats.columns:
            rename_dictionary["rsID"]="ID"
        else:
            rename_dictionary["SNPID"]="ID" 

        # logging
        if verbose: 
            log.write(" -gwaslab to "+fmt+" format dictionary:")  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            log.write("  - gwaslab keys:",','.join(keys)) 
            log.write("  - "+fmt+" values:",','.join(values)) 
        
        # get the columns to output
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        ouput_cols = ouput_cols +["STATUS"]+ cols
        sumstats = sumstats.loc[:,ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        
        # calculate meta data
        harmonised = sum(sumstats["STATUS"].str.match( r"\w\w\w[0][0123][012][01234]", case=False, flags=0, na=False ) )
        switchedalleles = sum(sumstats["STATUS"].str.match( r"\w\w\w[0][0123][12][24]", case=False, flags=0, na=False ) )
        sumstats["ID"] = sumstats["ID"].str.replace(":","_")
        
        # process Allele frequency data
        if "AF" in sumstats.columns:
            sumstats["INFO"] = "AF="+sumstats["AF"].astype("string")
        else:
            sumstats["INFO"] = "."

        # sumstats columns placed in vcf-FORMAT column     
        output_format=[]
        for i in sumstats.columns:
            if i in meta_data["format_format"]:
                output_format.append(i)  

        # Create vcf header
        vcf_header= meta_data["format_fixed_header"] +"\n"+ meta_data["format_contig_"+str(build)]+"\n"
        # Create sample header
        vcf_header+="##SAMPLE=<ID={},TotalVariants={},VariantsNotRead=0,HarmonisedVariants={},VariantsNotHarmonised={},SwitchedAlleles={},StudyType={}>\n".format(meta["gwaslab"]["study_name"],len(sumstats),harmonised,len(sumstats)-harmonised,switchedalleles,meta["gwaslab"]["study_type"])
        vcf_header+="##gwaslab_version="+gwaslab_info()["version"]+"\n"
        
        
         #StudyID=meta["Name"]
        #otalVariants = len(sumstats)
        #HarmonisedVariants = 
        #VariantsNotHarmonised = 
        #StudyType=
        ##SAMPLE=<ID=IEU-b-1,TotalVariants=9851866,VariantsNotRead=0,HarmonisedVariants=9851866,VariantsNotHarmonised=0,SwitchedAlleles=9851866,StudyType=Continuous>
        
        
        # determine path
        path = path + "."+suffix
        if verbose: log.write(" -Output path:",path) 
        if verbose: log.write(" -vcf header contig build:"+str(build)) 
        
        # output header
        with open(path,"w") as file:
            file.write(vcf_header)
        
        with open(path,"a") as file:
            if verbose: log.write(" -Output columns:"," ".join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]]))
            file.write("\t".join(meta_data["format_fixed"]+[meta["gwaslab"]["study_name"]])+"\n")
            if verbose: log.write(" -Outputing data...")
            counter=0
            QUAL="."
            FILTER="PASS"
            for index,row in sumstats.iterrows():
                CHROM=str(row["#CHROM"])
                POS=str(row["POS"])
                ID=str(row["ID"])
                REF=str(row["REF"])
                ALT=str(row["ALT"])
                INFO=str(row["INFO"])
                FORMAT=":".join(output_format)
                DATA=":".join(row[output_format].astype("string"))
                file.write(CHROM+"\t"+POS+"\t"+ID+"\t"+REF+"\t"+ALT+"\t"+QUAL+"\t"+FILTER+"\t"+INFO+"\t"+FORMAT+"\t"+DATA+"\n")

        if bgzip==True:
            if verbose: log.write(" -bgzip compressing ...") 
            tabix_compress(path, path+".gz",force=True)
            if md5sum is True: md5sum_file(path+".gz",log,verbose)
        if tabix==True:
            if verbose: log.write(" -tabix indexing...") 
            tabix_index(path+".gz" ,preset="vcf",force=True) 
    ####################################################################################################################
    elif fmt in get_formats_list():       
        if verbose: log.write(" -"+fmt+" format will be loaded...")
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        print_format_info(fmt=fmt, meta_data=meta_data,rename_dictionary=rename_dictionary,verbose=verbose, log=log, output=True)
        #if verbose:             
        #    log.write(" -"+fmt+" format meta info:")   
        #    for key,value in meta_data.items():
        #        if type(value) is list:
        #            log.write("  -",key," : ",','.join(value))
        #        else:
        #            log.write("  -",key," : ",value)
        #if verbose: 
        #    log.write(" -gwaslab to "+fmt+" format dictionary:",)  
        #    keys=[]
        #    values=[]
        #    for key,value in rename_dictionary.items():
        #        keys.append(key)
        #        values.append(value)
        #    log.write("  - gwaslab keys:",  ','.join(keys)) 
        #    log.write("  - "+fmt+" values:"  , ','.join(values))
        
        # grab format cols that exist in sumstats
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        # + additional cols
        ouput_cols = ouput_cols + cols
        
        sumstats = sumstats.loc[:,ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        

        ymal_path = path + "."+suffix+".tsv-meta.ymal"
        path = path + "."+suffix+".tsv.gz"
        
        if verbose: log.write(" -Output path:",path) 
        
        if path is not None: 
            if "format_separator" in meta_data.keys():
                to_csvargs["sep"] = meta_data["format_separator"]
            else:
                to_csvargs["sep"]="\t"
            if "format_na" in meta_data.keys():
                to_csvargs["na_rep"] = meta_data["format_na"]
            if "format_col_order" in meta_data.keys():
                fixed_col =[]
                other_col=[]
                for i in meta_data["format_col_order"]:
                    if i in sumstats.columns:
                        fixed_col.append(i)
                for i in sumstats.columns:
                    if i not in meta_data["format_col_order"]:
                        other_col.append(i)
                
                sumstats = sumstats.loc[:,fixed_col + other_col]
                if verbose: log.write(" -Reordering columns...") 

            if verbose: log.write(" -Output columns:",','.join(sumstats.columns))
            sumstats.to_csv(path, index=None,**to_csvargs)

        if md5sum is True: 
            md5_value = md5sum_file(path,log,verbose)
        else:
            md5_value = calculate_md5sum_file(path)
        
        ## update ssf-style meta data and export to yaml file
        if ssfmeta==True:
            meta_copy = meta.copy()
            if "format_cite_name" in meta_data.keys():
                meta_copy["file_type"] = meta_data["format_cite_name"]
            else:
                meta_copy["file_type"] = fmt
            meta_copy["minor_allele_freq_lower_limit"] = min_maf
            meta_copy["data_file_name"] = path
            meta_copy["data_file_md5sum"] = md5_value
            meta_copy["date_last_modified"] = get_format_date_and_time()
            meta_copy["samples"]["sample_size"] = n_max
            meta_copy["gwaslab"]["samples"]["sample_size_min"] = n_min
            meta_copy["gwaslab"]["samples"]["sample_size_median"] = n_median
            meta_copy["gwaslab"]["variants"]["variant_number"] = len(sumstats)
            if verbose: log.write(" -Exporting SSF-style meta data to {}".format(ymal_path)) 
            with open(ymal_path, 'w') as outfile:
                yaml.dump(meta_copy, outfile)
        return sumstats  
    
def md5sum_file(filename,log,verbose):
    if verbose: log.write(" -md5sum hashing for the file:",filename) 
    md5_hash = hashlib.md5()
    with open(filename,"rb") as f:
        # Read and update hash in chunks
        for byte_block in iter(lambda: f.read(4096*1000),b""):
            md5_hash.update(byte_block)
    with open(filename+".md5sum","w") as f:
        out = str(md5_hash.hexdigest())
        f.write(out+"\n")
        if verbose: log.write(" -md5sum path:",filename+".md5sum")    
        return out

def calculate_md5sum_file(filename):
    md5_hash = hashlib.md5()
    with open(filename,"rb") as f:
        # Read and update hash in chunks
        for byte_block in iter(lambda: f.read(4096*1000),b""):
            md5_hash.update(byte_block)
        out = str(md5_hash.hexdigest())
        return out    

def get_format_date_and_time():
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d-%H:%M:%S")
    return dt_string