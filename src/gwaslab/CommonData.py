import pandas as pd
from os import path
import json
from pyensembl import Genome

#hard-coded data
def get_chr_NC_dict(build,inverse=False):
    #https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
    if build =="19":
        dic={
        "1":"NC_000001.10",
        "2":"NC_000002.11",
        "3":"NC_000003.11",
        "4":"NC_000004.11",
        "5":"NC_000005.9",
        "6":"NC_000006.11",
        "7":"NC_000007.13",
        "8":"NC_000008.10",
        "9":"NC_000009.11",
        "10":"NC_000010.10",
        "11":"NC_000011.9",
        "12":"NC_000012.11",
        "13":"NC_000013.10",
        "14":"NC_000014.8",
        "15":"NC_000015.9",
        "16":"NC_000016.9",
        "17":"NC_000017.10",
        "18":"NC_000018.9",
        "19":"NC_000019.9",
        "20":"NC_000020.10",
        "21":"NC_000021.8",
        "22":"NC_000022.10",
        "X":"NC_000023.10",
        "Y":"NC_000024.9"}
    elif build=="38":
        dic={
        "1":"NC_000001.11",
        "2":"NC_000002.12",
        "3":"NC_000003.12",
        "4":"NC_000004.12",
        "5":"NC_000005.10",
        "6":"NC_000006.12",
        "7":"NC_000007.14",
        "8":"NC_000008.11",
        "9":"NC_000009.12",
        "10":"NC_000010.11",
        "11":"NC_000011.10",
        "12":"NC_000012.12",
        "13":"NC_000013.11",
        "14":"NC_000014.9",
        "15":"NC_000015.10",
        "16":"NC_000016.10",
        "17":"NC_000017.11",
        "18":"NC_000018.10",
        "19":"NC_000019.10",
        "20":"NC_000020.11",
        "21":"NC_000021.9",
        "22":"NC_000022.11",
        "X":"NC_000023.11",
        "Y":"NC_000024.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic

def get_chr_list():
    chrom_list=["1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18",
            "19","20","21","22","23","24","25","X","Y","M","MT"]
    return chrom_list

def get_chr_to_number(out_chr=False,xymt=["X","Y","MT"]):
    if out_chr is True:
        dic= {str(i):str(i) for i in range(1,26)}
        dic[xymt[0]]="23"
        dic[xymt[1]]="24"
        dic[xymt[2]]="25"
    
    else:
        dic= {str(i):i for i in range(1,26)}
        dic[xymt[0]]=23
        dic[xymt[1]]=24
        dic[xymt[2]]=25
    return dic

def get_number_to_chr(in_chr=False,xymt=["X","Y","MT"]):
    if in_chr is True:
        dic= {str(i):str(i) for i in range(1,23)}
        dic["23"]=xymt[0]
        dic["24"]=xymt[1]
        dic["25"]=xymt[2]
    else:
        dic= {i:str(i) for i in range(1,23)}
        dic[23]=xymt[0]
        dic[24]=xymt[1]
        dic[25]=xymt[2]
    return dic


# reading from files    
###################################################################################################################    
def get_high_ld(build="19"):
    if build=="19":
        data_path =  path.dirname(__file__) + '/data/high_ld/high_ld_hla_hg19.bed.gz'
    elif build=="38":
        data_path =  path.dirname(__file__) + '/data/high_ld/high_ld_hla_hg38.bed.gz'
    return data_path

def get_format_dict(fmt,inverse=False):
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    dicts = json.load(open(data_path))
    dic_meta = dicts[fmt]["meta_data"]
    dic_dict = dicts[fmt]["format_dict"]
    if inverse is True:
        inv_dic = {v: k for k, v in dic_dict.items()}
        return dic_meta,inv_dic
    return dic_meta, dic_dict

def get_formats_list():
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    dicts = json.load(open(data_path))
    format_list = list(dicts.keys())
    return format_list

def get_recombination_rate(chrom, build="19"):
    #http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110106_recombination_hotspots/
    if build=="19":
        data_path =  path.dirname(__file__) + '/data/recombination/hg19/genetic_map_GRCh37_chr'+str(chrom)+'.txt.gz'
        recombination_rate = pd.read_csv(data_path,sep="\t")
    return recombination_rate


####################################################################################################################
def get_gtf(chrom, build="19",source="ensembl"):
    if source=="ensembl":
        if build=="19":
            data_path =  path.dirname(__file__) + '/data/Ensembl/release75/Homo_sapiens.GRCh37.75.gtf.gz'
            gtf = pd.read_csv(data_path,sep="\t",header=None,comment="#",dtype={0:"string"})
            gtf = gtf.loc[gtf[0]==chrom,:]
        if build=="38":
            data_path =  path.dirname(__file__) + '/data/Ensembl/release107/Homo_sapiens.GRCh38.107.chr.gtf.gz'
            gtf = pd.read_csv(data_path,sep="\t",header=None,comment="#",dtype={0:"string"})
            gtf = gtf.loc[gtf[0]==chrom,:]
    if source=="refseq":
        chrom_NC = get_chr_NC_dict(build=build)[str(chrom)]
        if build=="19":
            data_path =  path.dirname(__file__) + '/data/RefSeq/GRCh37/GRCh37_latest_genomic.gtf.gz'
            gtf = pd.read_csv(data_path,sep="\t",header=None,comment="#")
            gtf = gtf.loc[gtf[0]==chrom_NC,:]
        if build=="38":
            data_path =  path.dirname(__file__) + '/data/RefSeq/GRCh38/GRCh38_latest_genomic.gtf.gz'
            gtf = pd.read_csv(data_path,sep="\t",header=None,comment="#")
            gtf = gtf.loc[gtf[0]==chrom_NC,:]
        gtf[0] = str(chrom)
    return gtf

#def get_gtf_chr(chrom, build="19"):
#    if build=="19":
#        data_path =  path.dirname(__file__) + '/data/Ensembl/release75/Homo_sapiens.GRCh37.75.chr'+str(chrom)+'.gtf.gz'
#        gtf = pd.read_csv(data_path,sep="\t",header=None,comment="#")
#    if build=="38":
#        data_path =  path.dirname(__file__) + '/data/Ensembl/release107/Homo_sapiens.GRCh38.107.chr'+str(chrom)+'.gtf.gz'
#        gtf = pd.read_csv(data_path,sep="\t",header=None,comment="#")
#    return gtf

def gtf_index():
    if build=="19":
        data_path =  path.dirname(__file__) + '/data/Ensembl/release75/Homo_sapiens.GRCh37.75.gtf.gz'
        data = Genome(  reference_name='GRCh37',
                        annotation_name='hg19_genes',
                        gtf_path_or_url=data_path)
        data.index()
    if build=="38":
        data_path =  path.dirname(__file__) + '/data/Ensembl/release107/Homo_sapiens.GRCh38.107.chr.gtf.gz'
        data = Genome(  reference_name='GRCh38',
                        annotation_name='hg38_genes',
                        gtf_path_or_url=data_path)
        data.index()
    
####################################################################################################################
    
def get_vcf_format(build="19",inverse=True):
    #Lyon, M.S., Andrews, S.J., Elsworth, B. et al. The variant call format provides efficient and robust storage of GWAS summary statistics. Genome Biol 22, 32 (2021). https://doi.org/10.1186/s13059-020-02248-0
    meta={}
    meta["format_name"]="vcf"
    meta["format_source"]="https://github.com/MRCIEU/gwas-vcf-specification/tree/1.0.0"
    meta["fixed_header"] = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">
##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">
##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">
##FORMAT=<ID=SS,Number=A,Type=Float,Description="Sample size used to estimate genetic effect">
##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">
##FORMAT=<ID=SI,Number=A,Type=Float,Description="Accuracy score of summary data imputation">
##FORMAT=<ID=NC,Number=A,Type=Float,Description="Number of cases used to estimate genetic effect">
##FORMAT=<ID=ID,Number=1,Type=String,Description="Study variant identifier">
##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">
##META=<ID=VariantsNotRead,Number=1,Type=Integer,Description="Number of variants that could not be read">
##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">
##META=<ID=VariantsNotHarmonised,Number=1,Type=Integer,Description="Total number of variants that could not be harmonised">
##META=<ID=SwitchedAlleles,Number=1,Type=Integer,Description="Total number of variants strand switched">
##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">
##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">
##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">
"""        
    format_dict={"ID":"SNPID",
        "ID":"rsID",
        "#CHROM":"CHR",
        "POS":"POS",
        "REF":"NEA",
        "ALT":"EA",
        "SS":"N",
        "AF":"EAF",
        "ES":"BETA",
        "SE":"SE",
        "LP":"MLOG10P",
        "SI":"INFO",
        "EZ":"Z"}
    if inverse is True:
        format_dict={"SNPID":"ID",
        "rsID":"ID",
        "CHR":"#CHROM",
        "POS":"POS",
        "NEA":"REF",
        "EA":"ALT",
        "N":"SS",
        "EAF":"AF",
        "BETA":"ES",
        "SE":"SE",
        "MLOG10P":"LP",
        "INFO":"SI",
        "Z":"EZ"}

    meta["FIXED"]=["ID","#CHROM","POS","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
    meta["FORMAT"]=["ID","SS","ES","SE","LP","SI","EZ"]
    if build=="19":
        meta["contig"]="""##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>
/n##contig=<ID=2,length=243199373,assembly=HG19/GRCh37>
/n##contig=<ID=3,length=198022430,assembly=HG19/GRCh37>
/n##contig=<ID=4,length=191154276,assembly=HG19/GRCh37>
/n##contig=<ID=5,length=180915260,assembly=HG19/GRCh37>
/n##contig=<ID=6,length=171115067,assembly=HG19/GRCh37>
/n##contig=<ID=7,length=159138663,assembly=HG19/GRCh37>
/n##contig=<ID=8,length=146364022,assembly=HG19/GRCh37>
/n##contig=<ID=9,length=141213431,assembly=HG19/GRCh37>
/n##contig=<ID=10,length=135534747,assembly=HG19/GRCh37>
/n##contig=<ID=11,length=135006516,assembly=HG19/GRCh37>
/n##contig=<ID=12,length=133851895,assembly=HG19/GRCh37>
/n##contig=<ID=13,length=115169878,assembly=HG19/GRCh37>
/n##contig=<ID=14,length=107349540,assembly=HG19/GRCh37>
/n##contig=<ID=15,length=102531392,assembly=HG19/GRCh37>
/n##contig=<ID=16,length=90354753,assembly=HG19/GRCh37>
/n##contig=<ID=17,length=81195210,assembly=HG19/GRCh37>
/n##contig=<ID=18,length=78077248,assembly=HG19/GRCh37>
/n##contig=<ID=19,length=59128983,assembly=HG19/GRCh37>
/n##contig=<ID=20,length=63025520,assembly=HG19/GRCh37>
/n##contig=<ID=21,length=48129895,assembly=HG19/GRCh37>
/n##contig=<ID=22,length=51304566,assembly=HG19/GRCh37>
/n##contig=<ID=X,length=155270560,assembly=HG19/GRCh37>
/n##contig=<ID=Y,length=59373566,assembly=HG19/GRCh37>
/n##contig=<ID=MT,length=16569,assembly=HG19/GRCh37>
"""
    if build=="38":
        meta["contig"]="""##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>
##contig=<ID=2,length=243199373,assembly=HG19/GRCh37>
##contig=<ID=3,length=198022430,assembly=HG19/GRCh37>
##contig=<ID=4,length=191154276,assembly=HG19/GRCh37>
##contig=<ID=5,length=180915260,assembly=HG19/GRCh37>
##contig=<ID=6,length=171115067,assembly=HG19/GRCh37>
##contig=<ID=7,length=159138663,assembly=HG19/GRCh37>
##contig=<ID=8,length=146364022,assembly=HG19/GRCh37>
##contig=<ID=9,length=141213431,assembly=HG19/GRCh37>
##contig=<ID=10,length=135534747,assembly=HG19/GRCh37>
##contig=<ID=11,length=135006516,assembly=HG19/GRCh37>
##contig=<ID=12,length=133851895,assembly=HG19/GRCh37>
##contig=<ID=13,length=115169878,assembly=HG19/GRCh37>
##contig=<ID=14,length=107349540,assembly=HG19/GRCh37>
##contig=<ID=15,length=102531392,assembly=HG19/GRCh37>
##contig=<ID=16,length=90354753,assembly=HG19/GRCh37>
##contig=<ID=17,length=81195210,assembly=HG19/GRCh37>
##contig=<ID=18,length=78077248,assembly=HG19/GRCh37>
##contig=<ID=19,length=59128983,assembly=HG19/GRCh37>
##contig=<ID=20,length=63025520,assembly=HG19/GRCh37>
##contig=<ID=21,length=48129895,assembly=HG19/GRCh37>
##contig=<ID=22,length=51304566,assembly=HG19/GRCh37>
##contig=<ID=X,length=155270560,assembly=HG19/GRCh37>
##contig=<ID=Y,length=59373566,assembly=HG19/GRCh37>
##contig=<ID=MT,length=16569,assembly=HG19/GRCh37>
"""
    return meta,format_dict
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    