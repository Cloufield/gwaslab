import pandas as pd
from os import path
import json
from pyensembl import Genome
import requests
from gwaslab.Log import Log

#hard-coded data
def get_chr_to_NC(build,inverse=False):
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
        "Y":"NC_000024.9",
        "MT":"NC_012920.1"}
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
        "Y":"NC_000024.1",
        "MT":"NC_012920.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic

def get_NC_to_chr(build):
    return get_chr_to_NC(build=build,inverse=True)


def get_number_to_NC(build,inverse=False):
    #https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
    if build =="19":
        dic={
        1:"NC_000001.10",
        2:"NC_000002.11",
        3:"NC_000003.11",
        4:"NC_000004.11",
        5:"NC_000005.9",
        6:"NC_000006.11",
        7:"NC_000007.13",
        8:"NC_000008.10",
        9:"NC_000009.11",
        10:"NC_000010.10",
        11:"NC_000011.9",
        12:"NC_000012.11",
        13:"NC_000013.10",
        14:"NC_000014.8",
        15:"NC_000015.9",
        16:"NC_000016.9",
        17:"NC_000017.10",
        18:"NC_000018.9",
        19:"NC_000019.9",
        20:"NC_000020.10",
        21:"NC_000021.8",
        22:"NC_000022.10",
        23:"NC_000023.10",
        24:"NC_000024.9",
        25:"NC_012920.1"}
    elif build=="38":
        dic={
        1:"NC_000001.11",
        2:"NC_000002.12",
        3:"NC_000003.12",
        4:"NC_000004.12",
        5:"NC_000005.10",
        6:"NC_000006.12",
        7:"NC_000007.14",
        8:"NC_000008.11",
        9:"NC_000009.12",
        10:"NC_000010.11",
        11:"NC_000011.10",
        12:"NC_000012.12",
        13:"NC_000013.11",
        14:"NC_000014.9",
        15:"NC_000015.10",
        16:"NC_000016.10",
        17:"NC_000017.11",
        18:"NC_000018.10",
        19:"NC_000019.10",
        20:"NC_000020.11",
        21:"NC_000021.9",
        22:"NC_000022.11",
        23:"NC_000023.11",
        24:"NC_000024.1",
        25:"NC_012920.1"
        }
    if inverse is True:
        inv_dic = {v: k for k, v in dic.items()}
        return inv_dic
    return dic


def get_NC_to_number(build):
    return get_number_to_NC(build=build,inverse=True)



def get_chr_list(add_number=False):
    chrom_list=["1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18",
            "19","20","21","22","23","24","25","X","Y","M","MT"]
    if add_number is True:
        chrom_list = ["1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18",
            "19","20","21","22","23","24","25","X","Y","M","MT"] +[i for i in range(1,26)]
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
        dic["M"]=25
    return dic

def get_number_to_chr(in_chr=False,xymt=["X","Y","MT"],prefix=""):
    if in_chr is True:
        dic= {str(i):prefix+str(i) for i in range(1,23)}
        dic["23"]=prefix+xymt[0]
        dic["24"]=prefix+xymt[1]
        dic["25"]=prefix+xymt[2]
    else:
        dic= {i:prefix+str(i) for i in range(1,23)}
        dic[23]=prefix+xymt[0]
        dic[24]=prefix+xymt[1]
        dic[25]=prefix+xymt[2]
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
        chrom_NC = get_chr_to_NC(build=build)[str(chrom)]
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

      
def update_formatbook(log=Log()):
    url = 'https://raw.github.com/Cloufield/formatbook/main/formatbook.json'
    log.write("Updating formatbook from:",url)
    r = requests.get(url, allow_redirects=True)
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    with open(data_path, 'wb') as file:
        file.write(r.content)
    book=json.load(open("formatbook.json"))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))
    log.write("Formatbook has been updated!")
        
         
def list_formats(log=Log()):
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    book=json.load(open("formatbook.json"))
    available_formats = list(book.keys())
    available_formats.sort()
    log.write("Available formats:",",".join(available_formats))    

def check_format(fmt,log=Log()):
    data_path =  path.dirname(__file__) + '/data/formatbook.json'
    book=json.load(open("formatbook.json"))
    log.write("Available formats:",end="")
    for i in book[fmt].keys():
        log.write(i,end="")
    log.write("") 
    for i in book[fmt].values():
        log.write(i,end="")
        
####################################################################################################################       
        
def gwaslab_info():
    dic={
   "version":"3.3.3",
   "release_date":"20221029"
    }
    return dic       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

    