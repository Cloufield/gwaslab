import pandas as pd
from gwaslab.Log import Log
from gwaslab.CommonData import get_format_dict
from gwaslab.CommonData import get_number_to_chr
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
def tovcf(sumstats,path=None,snpid="MARKERNAME", chrom="CHR", pos="POS", ea="EA", nea="NEA", beta="BETA", se="SE" ,p="P" ,n="N",verbose=True):
    #https://github.com/mrcieu/gwasvcf/blob/master/R/manipulate.r
    #df['QUAL'] = "NA"
    #df['FILTER'] = "PASS"
    #df['INFO'] = ""
    #df = df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']] 

    #header = """##fileformat=VCFv4.1
    ##fileDate=20090805
    ##source=myImputationProgramV3.1
    ##reference=file:///seq/references/
    #CHROM POS ID REF ALT QUAL FILTER INFO
    #"""
    #output_VCF = "myfile.vcf"
    #with open(output_VCF, 'w') as vcf:
    #    vcf.write(header)
    #df.to_csv(output_VCF, sep="\t", mode='a', index=False)
    ##fileformat=VCFv4.2
    ##fileDate=20090805
    ##source=myImputationProgramV3.1
    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##phasing=partial
    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
    ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
    ##FILTER=<ID=q10,Description="Quality below 10">
    ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
    #20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
    #20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
    #20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
    #20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
    #20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G
    pass
    
###################################################################################################################################################
def tofmt(sumstats,
          path=None,
          suffix=None,
          fmt=None,
          cols=[],
          xymt_number=False,
          xymt=["X","Y","MT"],
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
    if fmt in get_formats_list():       
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
    
    elif fmt=="bed":
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
    
    
    
    
    
    
    
    
    
    
    