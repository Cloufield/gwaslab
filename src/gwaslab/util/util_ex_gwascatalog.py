from typing import Optional, Union
import requests
import json
import pandas as pd
import gwaslab as gl
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
from datetime import datetime
import os

def find_efo_cache(efo: str, path: str) -> Union[str, bool]:
    for root, dirs, files in os.walk(path):
        for file in files:
            if efo in file:
                return os.path.join(root, file)
    return False

@with_logging(
        start_to_msg="retrieve data from GWASCatalog",
        finished_msg="retrieving data from GWASCatalog"
)
def gwascatalog_trait(efo: str,
                      source: str = "NCBI",
                      sig_level: float = 5e-8,
                      use_cache: bool = True,
                      cache_dir: str = "./",
                      verbose: bool = True,
                      log: Log = Log()) -> pd.DataFrame:
    
    #https://www.ebi.ac.uk/gwas/rest/docs/api
    
    base_url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"+efo
    
    
    if use_cache==True:
        log.write("searching cache in : {}".format(cache_dir))
        cache = find_efo_cache(efo, cache_dir)
        if cache==False:
            log.write(" -Cache not found for {}... Downloading from GWASCatalog...".format(cache), verbose=verbose)
    else:
        cache = False

    if cache==False:
        #log.write(" -Please make sure your sumstats is based on GRCh38...", verbose=verbose)
        log.write(" -Requesting (GET) trait information through the GWASCatalog API...", verbose=verbose)
        log.write(" -EFO trait api: "+ base_url, verbose=verbose)
        text = requests.get(base_url)

        log.write(" -Status code: {}".format(text.status_code), verbose=verbose) 
        if text.status_code!=200:
            log.write(" -Status code is not 200. Access failed. Please check your internet or the GWAS Catalog sever status.", verbose=verbose) 
            log.write(" -Message:{}".format(text.text), verbose=verbose) 
            return 0

        api_response = json.loads(text.text)
        log.write(" -Trait Name:",api_response["trait"], verbose=verbose)
        log.write(" -Trait URL:",api_response["uri"], verbose=verbose) 
            
        base_url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"+efo+"/associations?projection=associationByEfoTrait"    
        log.write(" -Requesting (GET) GWAS associations through the GWASCatalog API...", verbose=verbose)
        log.write(" -associationsByTraitSummary API: "+ base_url, verbose=verbose)   
        log.write(" -Note: this step might take a while...", verbose=verbose)   
        
        # get request and check status code of response
        raw_data = requests.get(base_url)
        
        # whether to proceed based on status code
        is_proceed = check_request_status_code(raw_data.status_code,verbose=verbose,log=log)
        if is_proceed is False: return False
        

        log.write(" -Loading json ...", verbose=verbose)
        # Transform API response from JSON into Python dictionary
        api_response = json.loads(raw_data.text)

        now = datetime.now() # current date and time
        datestring = now.strftime("%Y%m%d")
        json_path = cache_dir + "GWASCatalog_{}_associationsByTraitSummary_text_{}.json".format(efo, datestring)
        
        try:
            log.write(" -Saving json to: {} ...".format(json_path), verbose=verbose) 
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(api_response, f, ensure_ascii=False, indent=4)
        except:
            pass
    else:
        log.write(" -Loading cache for {}: {} ...".format(efo, cache), verbose=verbose) 
        with open(cache) as f:
            api_response = json.load(f)
        
    log.write(" -Parsing json ...", verbose=verbose)        
    # An 
    records=list()
    log.write(" -Number of reported associations for "+ efo +" in GWASCatalog:",len( api_response["_embedded"]["associations"]), verbose=verbose)
   
    for association in api_response["_embedded"]["associations"]:
        #association statistics:       
        p=float(association["pvalue"])    
        # filter association by p value
        if p < sig_level:
            # obtain statistics
            try:
                function_class=association["functionalClass"] 
            except:
                function_class=None
            try:
                eaf=association["riskFrequency"]
            except:
                eaf= None
            try:
                study=association["study"]['publicationInfo']["title"]
                pubmedid=association["study"]['publicationInfo']["pubmedId"]
                author=association["study"]['publicationInfo']["author"]["fullname"]
            except:
                study= None
                pubmedid = None
                author=None
            try:
                gene = association["loci"][0]["authorReportedGenes"][0]["geneName"]
            except:
                gene = None
            try:
                OR=association["orPerCopyNum"]
            except:
                OR=None
            try:
                beta=association["betaNum"]
                se=association["standardError"]
            except:
                beta=None
                se=None    
            #########################################################
            #obtain snp information
            for snp in association["snps"]:
                if len(snp["locations"])>0:
                    for record_num in range(len(snp["locations"])):
                        if snp["locations"][record_num]["chromosomeName"] in [str(i) for i in range(1,26)]+["x","X","y","Y","mt","MT"]:
                            if len(snp["genomicContexts"])>0:
                                closegenes=set()
                                distances=set()
                                ingene=set()
                                for gene_num in range(len(snp["genomicContexts"])):
                                    if snp["genomicContexts"][gene_num]["source"]==source:
                                        distance= str(snp["genomicContexts"][gene_num]["distance"])
                                        ingene_name = snp["genomicContexts"][gene_num]["gene"]["geneName"]
                                        if distance==0:
                                            ingene.add(ingene_name)
                                            continue
                                        if snp["genomicContexts"][gene_num]["isClosestGene"] is True:
                                            closegene_name = snp["genomicContexts"][gene_num]["gene"]["geneName"]
                                            closegenes.add(closegene_name)
                                            distances.add(distance)

                                if len(ingene)>0:
                                    autogenes =  ",".ingene
                                else:
                                    autogenes = ",".join(closegenes)

                                row=[ snp["rsId"],
                                      snp["locations"][record_num]["chromosomeName"],
                                      snp["locations"][record_num]["chromosomePosition"],
                                      gene,
                                      autogenes,
                                      function_class,
                                      OR,
                                      beta,
                                      se,
                                      p,
                                      association["study"]["diseaseTrait"]["trait"],
                                      study,
                                      pubmedid,
                                      author
                                    ]
                                records.append(row)
            #rsid locations
    gwascatalog_lead_snps = pd.DataFrame(records,columns=["SNPID","CHR","POS","REPORT_GENENAME","CLOSEST_GENENAMES","FUNCTION_CLASS","OR","BETA","SE","P","TRAIT","STUDY","PUBMEDID","AUTHOR"])
    log.write(" -Loading retrieved data into gwaslab Sumstats object ...", verbose=verbose)  
    sigs = gl.Sumstats(gwascatalog_lead_snps.copy(),fmt="gwaslab",other=['REPORT_GENENAME', 'CLOSEST_GENENAMES','TRAIT', 'STUDY', 'PUBMEDID','AUTHOR'],verbose=False)
    sigs.fix_pos(verbose=False)
    sigs.fix_chr(verbose=False)
    sigs.sort_coordinate(verbose=False)
    log.write("Finished retrieving data from GWASCatalog...", verbose=verbose)
    #return gwaslab Sumstats object
    return sigs


###### helper ##################################################################################################
def check_request_status_code(
    request_code: int,
    verbose: bool = True,
    log: Log = Log()
) -> bool:
    
    is_proceed=False
    
    if request_code == 200:
        log.write(" -Status code 200 OK: Retrieved data from GWASCatalog successffully ...", verbose=verbose)
        is_proceed=True
    elif request_code == 404:
        log.write(" -Status code 404 Not Found: The requested resource did not exist ...", verbose=verbose)
    elif request_code == 301:
        log.write(" -Status code 301 Moved Permanently: The requested resource did not exist ...", verbose=verbose)
    elif request_code == 400:
        log.write(" -Status code 400 Bad Request: The requested resource did not exist ...", verbose=verbose)
    
    return is_proceed

