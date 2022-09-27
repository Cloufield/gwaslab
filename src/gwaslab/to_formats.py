import pandas as pd
from gwaslab.Log import Log
from gwaslab.CommonData import get_format_dict
from gwaslab.CommonData import get_number_to_chr
from gwaslab.CommonData import gwaslab_info
from pysam import tabix_compress 
from pysam import tabix_index
from gwaslab.CommonData import get_formats_list

# to vcf
# to fmt
    ## general : ldsc, plink, plink2, saige, regenie ...
    ## vep
    ## bed
    ## annovar

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
          chr_prefix=None,
          bgzip=False,
          tabix=False,
          verbose=True,
          log=Log(),
          to_csvargs={}):
    
    if verbose: log.write(" - Start outputting sumstats in "+fmt+" format...")
    if xymt_number is False:
        sumstats["CHR"]= sumstats["CHR"].map(get_number_to_chr(xymt=xymt))
    if chr_prefix is not None:
        sumstats["CHR"]= chr_prefix + sumstats["CHR"].astype("string")
    
    ### calculate meta data

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

        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"]
        
        sumstats["START"] = sumstats["START"].astype("Int64")
        sumstats["END"] = sumstats["END"].astype("Int64")
        
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
        
        ouput_cols=["CHR","START","END","NEA/EA","STRAND","SNPID"]
        sumstats = sumstats.loc[:,ouput_cols]
        path = path + "."+suffix+".gz"
        if verbose: log.write(" -Output columns:",sumstats.columns)
        if verbose: log.write(" -Output path:",path) 

        sumstats.to_csv(path,sep="\t",index=None,header=None,**to_csvargs)
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
        
        ouput_cols=["CHR","START","END","NEA_out","EA_out","SNPID"]
        
        sumstats["START"] = sumstats["START"].astype("Int64")
        sumstats["END"] = sumstats["END"].astype("Int64")
        
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
        if tabix is True:
            if verbose: log.write(" -tabix indexing...") 
            tabix_index(path+".gz" ,preset="bed",force=True) 
    ####################################################################################################################       
    elif fmt=="vcf":
        if verbose: log.write(" -"+fmt+" format will be loaded...")
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        if verbose:             
            log.write(" -"+fmt+" format meta info:")   
            for key,value in meta_data.items():
                if key not in ["format_fixed_header","format_contig_19","format_contig_38"]:
                    log.write("  -",key," : ",value)
        if "rsID" in sumstats.columns:
            rename_dictionary["rsID"]="ID"
        else:
            rename_dictionary["SNPID"]="ID"
        
        if verbose: 
            log.write(" -gwaslab to "+fmt+" format dictionary:")  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            log.write("  - gwaslab keys:",keys) 
            log.write("  - "+fmt+" values:",values) 
            
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        
        ouput_cols = ouput_cols +["STATUS"]+ cols
        sumstats = sumstats.loc[:,ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        
        harmonised = sum(sumstats["STATUS"].str.match( r"\w\w\w[0][0123][012][01234]", case=False, flags=0, na=False ) )
        switchedalleles = sum(sumstats["STATUS"].str.match( r"\w\w\w[0][0123][12][24]", case=False, flags=0, na=False ) )
        sumstats["ID"] = sumstats["ID"].str.replace(":","_")
        
        if "AF" in sumstats.columns:
            sumstats["INFO"] = "AF="+sumstats["AF"].astype("string")
        else:
            sumstats["INFO"] = "."
            
        output_format=[]
        for i in sumstats.columns:
            if i in meta_data["format_format"]:
                output_format.append(i)  
                
        path = path + "."+suffix
        if verbose: log.write(" -vcf header contig build:"+str(build)) 
        # Create vcf header
        
        vcf_header= meta_data["format_fixed_header"] +"\n"+ meta_data["format_contig_"+str(build)]+"\n"
        # Create sample header
        vcf_header+="##SAMPLE=<ID={},TotalVariants={},VariantsNotRead=0,HarmonisedVariants={},VariantsNotHarmonised={},SwitchedAlleles={},StudyType={}>\n".format(meta["Name"],len(sumstats),harmonised,len(sumstats)-harmonised,switchedalleles,meta["StudyType"])
        vcf_header+="##gwaslab_version="+gwaslab_info()["version"]+"\n"
        
        
         #StudyID=meta["Name"]
        #otalVariants = len(sumstats)
        #HarmonisedVariants = 
        #VariantsNotHarmonised = 
        #StudyType=
        ##SAMPLE=<ID=IEU-b-1,TotalVariants=9851866,VariantsNotRead=0,HarmonisedVariants=9851866,VariantsNotHarmonised=0,SwitchedAlleles=9851866,StudyType=Continuous>
        
        # output header
        with open(path,"w") as file:
            file.write(vcf_header)
        
        with open(path,"a") as file:
            if verbose: log.write(" -Output columns:"," ".join(meta_data["format_fixed"]+[meta["Name"]]))
            file.write("\t".join(meta_data["format_fixed"]+[meta["Name"]])+"\n")
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

        if verbose: log.write(" -Output path:",path) 
#        if path is not None: sumstats.to_csv(path,sep="\t",index=None,mode='a',**to_csvargs)
        if bgzip is True:
            if verbose: log.write(" -bgzip compressing ...") 
            tabix_compress(path, path+".gz",force=True)
        if tabix is True:
            if verbose: log.write(" -tabix indexing...") 
            tabix_index(path+".gz" ,preset="vcf",force=True) 
    ####################################################################################################################
    elif fmt in get_formats_list():       
        if verbose: log.write(" -"+fmt+" format will be loaded...")
        meta_data,rename_dictionary = get_format_dict(fmt,inverse=True)
        if verbose:             
            log.write(" -"+fmt+" format meta info:")   
            for key,value in meta_data.items():
                log.write("  -",key," : ",value)
        if verbose: 
            log.write(" -gwaslab to "+fmt+" format dictionary:",)  
            keys=[]
            values=[]
            for key,value in rename_dictionary.items():
                keys.append(key)
                values.append(value)
            log.write("  - gwaslab keys:",keys) 
            log.write("  - "+fmt+" values:",values) 
            
        ouput_cols=[]
        for i in sumstats.columns:
            if i in rename_dictionary.keys():
                ouput_cols.append(i)  
        ouput_cols = ouput_cols + cols
        sumstats = sumstats.loc[:,ouput_cols]
        sumstats = sumstats.rename(columns=rename_dictionary) 
        path = path + "."+suffix+".tsv.gz"
        if verbose: log.write(" -Output columns:",sumstats.columns)
        if verbose: log.write(" -Output path:",path) 
        if path is not None: sumstats.to_csv(path,sep="\t",index=None,**to_csvargs)
        return sumstats  
    
    
    
    
    
    
    
    
    
    