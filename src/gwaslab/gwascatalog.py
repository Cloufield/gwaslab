import requests
import json
import pandas as pd
import gwaslab as gl
from gwaslab.Log import Log

def gwascatalog_trait(efo,source="NCBI",sig_level=5e-8,verbose=True,log=Log()):
    
    #https://www.ebi.ac.uk/gwas/rest/docs/api
    
    base_url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"+efo
    if verbose: log.write("Start to retrieve data from GWASCatalog...")
    if verbose: log.write(" -Please make sure your sumstats is based on GRCh38...")
    if verbose: log.write(" -Requesting (GET) trait information through the GWASCatalog API...")
    if verbose: log.write(" -EFO trait api: "+ base_url)
    api_response = json.loads(requests.get(base_url).text)
    if verbose: log.write(" -Trait Name:",api_response["trait"])
    if verbose: log.write(" -Trait URL:",api_response["uri"]) 
        
    base_url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/"+efo+"/associations?projection=associationByEfoTrait"    
    if verbose: log.write(" -Requesting (GET) GWAS associations through the GWASCatalog API...")
    if verbose: log.write(" -associationsByTraitSummary API: "+ base_url)   
    if verbose: log.write(" -Note: this step might take a while...")   
    
    # get request and check status code of response
    raw_data = requests.get(base_url)
    
    # whether to proceed based on status code
    is_proceed = check_request_status_code(raw_data.status_code,verbose=verbose,log=log)
    if is_proceed is False: return False
    
    if verbose: log.write(" -Loading json ...")
    # Transform API response from JSON into Python dictionary
    api_response = json.loads(raw_data.text)
    if verbose: log.write(" -Parsing json ...")        
    # An 
    records=list()
    if verbose: log.write(" -Number of reported associations for "+ efo +" in GWASCatalog:",len( api_response["_embedded"]["associations"]))
   
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
    if verbose: log.write(" -Loading retrieved data into gwaslab Sumstats object ...")  
    sigs = gl.Sumstats(gwascatalog_lead_snps,fmt="gwaslab",other=['REPORT_GENENAME', 'CLOSEST_GENENAMES','TRAIT', 'STUDY', 'PUBMEDID','AUTHOR'],verbose=False)
    sigs.fix_pos(verbose=False)
    sigs.fix_chr(verbose=False)
    sigs.sort_coordinate(verbose=False)
    if verbose: log.write("Finished retrieving data from GWASCatalog...")
    #return gwaslab Sumstats object
    return sigs


###### helper ##################################################################################################
def check_request_status_code(request_code,verbose=True,log=Log()):
    
    is_proceed=False
    
    if request_code == 200:
        if verbose: log.write(" -Status code 200 OK: Retrieved data from GWASCatalog successffully ...")
        is_proceed=True
    elif request_code == 404:
        if verbose: log.write(" -Status code 404 Not Found: The requested resource did not exist ...")
    elif request_code == 301:
        if verbose: log.write(" -Status code 301 Moved Permanently: The requested resource did not exist ...")
    elif request_code == 400:
        if verbose: log.write(" -Status code 400 Bad Request: The requested resource did not exist ...")
    
    return is_proceed

