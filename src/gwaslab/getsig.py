import pandas as pd
import numpy as np
import scipy as sp
from gwaslab.Log import Log
from gwaslab.CommonData import get_chr_to_number
from gwaslab.CommonData import get_number_to_chr
from gwaslab.CommonData import get_chr_to_NC
from gwaslab.download import check_and_download
from gwaslab.gwascatalog import gwascatalog_trait
from pyensembl import EnsemblRelease
from pyensembl import Genome
from os import path
from gwaslab.fill import fill_p
import gc

# getsig
# closest_gene
# annogene
# getnovel

def getsig(insumstats,
           id,
           chrom,
           pos,
           p,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           build="19",
           source="ensembl",
           verbose=True):
    
    if verbose: log.write("Start to extract lead variants...")
    if verbose: log.write(" -Processing "+str(len(insumstats))+" variants...")
    if verbose: log.write(" -Significance threshold :", sig_level)
    if verbose: log.write(" -Sliding window size:", str(windowsizekb) ," kb")
    #load data
    
    sumstats=insumstats.loc[~insumstats[id].isna(),:].copy()
    
    #convert chrom to int
    if sumstats[chrom].dtype in ["object",str,pd.StringDtype]:
        chr_to_num = get_chr_to_number(out_chr=True,xymt=["X","Y","MT"])
        sumstats[chrom]=sumstats[chrom].map(chr_to_num)
    
    sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    if p not in sumstats.columns:
        if "MLOG10P" in sumstats.columns: 
            fill_p(sumstats,log)
    sumstats[p] = pd.to_numeric(sumstats[p], errors='coerce')
    
    #extract all significant variants
    sumstats_sig = sumstats.loc[sumstats[p]<sig_level,:]
    if verbose:log.write(" -Found "+str(len(sumstats_sig))+" significant variants in total...")
    
    #sort the coordinates
    sumstats_sig = sumstats_sig.sort_values([chrom,pos])
    
    if sumstats_sig is None:
        if verbose:log.write(" -No lead snps at given significance threshold!")
        return None
    
    sig_index_list=[]
    current_sig_index = False
    current_sig_p = 1
    current_sig_pos = 0
    current_sig_chr = 0
    
    #iterate through all significant snps
    for line_number,(index, row) in enumerate(sumstats_sig.iterrows()):
        #when finished one chr 
        if row[chrom]!=current_sig_chr:
            #add the current lead variants id to lead variant list
            if current_sig_index is not False:sig_index_list.append(current_sig_index)
            
            #update lead vairant info to the new variant
            current_sig_chr=row[chrom]
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
            
            # only one significant variant on a new chromsome and this is the last sig variant
            if  line_number == len(sumstats_sig)-1:
                sig_index_list.append(current_sig_index)
            continue
        
        # next loci : gap > windowsizekb*1000
        if row[pos]>current_sig_pos + windowsizekb*1000:
            sig_index_list.append(current_sig_index)
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
            continue
        
        # update current pos and p
        if row[p]<current_sig_p:
            current_sig_pos=row[pos]
            current_sig_p=row[p]
            current_sig_index=row[id]
        else:
            current_sig_pos=row[pos]
            
        #when last line in sig_index_list
        if  line_number == len(sumstats_sig)-1:
            sig_index_list.append(current_sig_index)
            continue
    
    if verbose:log.write(" -Identified "+str(len(sig_index_list))+" lead variants!")
    
    #num_to_chr = get_number_to_chr(in_chr=True,xymt=xymt)
    #sumstats_sig.loc[:,chrom] = sumstats_sig[chrom].astype("string")
    #sumstats_sig.loc[:,chrom] = sumstats_sig.loc[:,chrom].map(num_to_chr)
    output = sumstats_sig.loc[sumstats_sig[id].isin(sig_index_list),:].copy()
    if anno is True and len(output)>0:
        output = annogene(
               output,
               id=id,
               chrom=chrom,
               pos=pos,
               log=log,
               xymt=xymt,
               build=build,
               source=source,
               verbose=verbose)
    if verbose: log.write("Finished extracting lead variants successfully!")
    gc.collect()
    return output.copy()


def closest_gene(x,data,chrom="CHR",pos="POS",maxiter=20000,step=50,source="ensembl",build="19"):
        #
        # data
        # check snp position
        #convert 23,24,25 back to X,Y,MT for EnsemblRelease query
        if source=="ensembl":
            contig = get_number_to_chr()[x[chrom]]
        elif source=="refseq":
            contig = get_number_to_chr()[x[chrom]]
            contig = get_chr_to_NC(build=build)[contig]
        
        position = int(x[pos])
        # query
        gene = data.gene_names_at_locus(contig=contig, position=position)
        
        if len(gene)==0:
            # if not in any gene
            i=0
            while i<=maxiter:
                # using distance to check upstram and downstream region
                distance = i*step
                # upstream
                gene_u = data.gene_names_at_locus(contig=contig, position=position-distance)
                
                # downstream
                gene_d = data.gene_names_at_locus(contig=contig, position=position+distance)
                
                if len(gene_u)>0 and len(gene_d)>0:
                    # if found gene uptream and downstream at the same time 
                    # go back to last step
                    distance = (i-1)*step
                    for j in range(0,step,1):
                        # use small step to finemap                        
                        gene_u = data.gene_names_at_locus(contig=contig, position=position-distance-j)
                        gene_d = data.gene_names_at_locus(contig=contig, position=position+distance+j)
                        if len(gene_u)>0:
                            return -distance-j,",".join(gene_u)
                        elif len(gene_d)>0:
                            return distance+j,",".join(gene_d)
                elif len(gene_u)>0:                    
                    # if found gene uptream
                    distance = (i-1)*step
                    for j in range(0,step,1):
                        gene_u2 = data.gene_names_at_locus(contig=contig, position=position-distance-j)
                        if len(gene_u2)>0:
                            return -distance-j,",".join(gene_u)
                elif len(gene_d)>0:
                    # if found gene downstream
                    distance = (i-1)*step
                    for j in range(0,step,1):
                        gene_d2 = data.gene_names_at_locus(contig=contig, position=position+distance+j)
                        if len(gene_d2)>0:
                            return distance+j,",".join(gene_d)
                i+=1
                # increase i by 1
            return distance,"intergenic"
        else:
            return 0,",".join(gene)


def annogene(
           insumstats,
           id,
           chrom,
           pos,
           log=Log(),
           xymt=["X","Y","MT"],
           build="19",
           source="ensembl",
           verbose=True):
    
    if verbose: log.write("Start to annotate variants with nearest gene name(s)...")
    output = insumstats.copy()
    
    if source == "ensembl":
        if build=="19":
            #data = EnsemblRelease(75)
            if verbose:log.write(" -Assigning Gene name using Ensembl Release hg19")
            #zcat Homo_sapiens.GRCh37.75.gtf.gz| 
            #grep -E 'processed_transcript|protein_coding|_gene' 
            #| gzip >Homo_sapiens.GRCh37.75.processed.chr.gtf.gz     
            
            gtf_path = check_and_download("ensembl_hg19_gtf_protein_coding")
            gtf_db_path = gtf_path[:-2]+"db"
            
            data = Genome(
                reference_name='GRCh37',
                annotation_name='hg19_gene',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,source=source), axis=1)), 
                index=output.index).values
        elif build=="38":
            if verbose:log.write(" -Assigning Gene name using built-in Ensembl Release hg38")
            gtf_path = check_and_download("ensembl_hg38_gtf_protein_coding")
            gtf_db_path = gtf_path[:-2]+"db"
            data = Genome(
                reference_name='GRCh38',
                annotation_name='hg38_gene',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,source=source), axis=1)), 
                index=output.index).values
    
    if source == "refseq":
        if build=="19":
            if verbose:log.write(" -Assigning Gene name using NCBI refseq latest GRCh37")
            gtf_path = check_and_download("refseq_hg19_gtf_protein_coding")
            gtf_db_path = gtf_path[:-2]+"db"
            data = Genome(
                reference_name='GRCh37',
                annotation_name='hg19_gene',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,source=source,build=build), axis=1)), 
                index=output.index).values
        elif build=="38":
            if verbose:log.write(" -Assigning Gene name using built-in NCBI refseq latest GRCh38")
            gtf_path = check_and_download("refseq_hg38_gtf_protein_coding")
            gtf_db_path = gtf_path[:-2]+"db"
            data = Genome(
                reference_name='GRCh38',
                annotation_name='hg38_gene',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,source=source,build=build), axis=1)), 
                index=output.index).values
    if verbose: log.write("Finished annotating variants with nearest gene name(s) successfully!")
    return output

def getnovel(insumstats,
           id,
           chrom,
           pos,
           p,
           known=False,
           efo=False,
           only_novel=False,
           windowsizekb_for_novel=1000,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           build="19",
           source="ensembl",
           gwascatalog_source="NCBI",
           output_known=False,
           verbose=True):
    if verbose: log.write("Start to check if lead variants are known...")
    allsig = getsig(insumstats=insumstats,
           id=id,chrom=chrom,pos=pos,p=p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
           xymt=xymt,anno=anno,build=build, source=source,verbose=verbose)
    
    big_number = 1000000000
    for i in range(7):
        if insumstats["POS"].max()*10 >  big_number:
            big_number = int(big_number * 10)
        else:
            break

    # create helper column TCHR+POS for allsig
    allsig["TCHR+POS"]=allsig[chrom]*big_number + allsig[pos]
    
    knownsig = pd.DataFrame()
    if efo != False:
        known_Sumstats = gwascatalog_trait(efo,source=gwascatalog_source,sig_level=sig_level,verbose=verbose,log=log)
        knownsig = known_Sumstats.data.copy()
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
        if verbose: log.write(" -Retrieved {} associations from GWAS catalog.".format(len(knownsig)))
    if type(known) is pd.DataFrame:
        knownsig_2 = known.copy()
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
        if "SNPID" not in knownsig.columns:
            knownsig["SNPID"] =knownsig["CHR"].astype("string") + ":" + knownsig["POS"].astype("string")
    elif type(known) is str:
        knownsig_2 = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
        if "SNPID" not in knownsig.columns:
            knownsig["SNPID"] =knownsig["CHR"].astype("string") + ":" + knownsig["POS"].astype("string")
    if len(knownsig)<1:
        raise ValueError("Please input a dataframe of known loci or valid efo code")

    # create helper column TCHR+POS for knownsig
    knownsig["TCHR+POS"]=knownsig[chrom]*big_number + knownsig[pos]
    
    if verbose: log.write(" -Lead variants in known loci:",len(knownsig))
    if verbose: log.write(" -Checking the minimum distance between identified lead variants and provided known variants...")
    
    #sorting
    allsig = allsig.sort_values(by="TCHR+POS",ignore_index=True)
    knownsig = knownsig.sort_values(by="TCHR+POS",ignore_index=True)
    
    if "SNPID" in knownsig.columns:
        knownids=knownsig["SNPID"].values
    if "PUBMEDID" in knownsig.columns:
        knownpubmedids=knownsig["PUBMEDID"].values
    if "AUTHOR" in knownsig.columns:
        knownauthor=knownsig["AUTHOR"].values
    
    # get distance
    lambda x:np.min(np.abs(knownsig["TCHR+POS"]-x))
    allsig["DISTANCE_TO_KNOWN"] = allsig["TCHR+POS"].apply(lambda x:min(knownsig["TCHR+POS"]-x, key=abs))
    
    # get other info 
    if "SNPID" in knownsig.columns:
        allsig["KNOWN_ID"] = allsig["TCHR+POS"].apply(lambda x:knownids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])    
    if "PUBMEDID" in knownsig.columns:
        allsig["KNOWN_PUBMED_ID"] = allsig["TCHR+POS"].apply(lambda x:knownpubmedids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    if "AUTHOR" in knownsig.columns:
        allsig["KNOWN_AUTHOR"] = allsig["TCHR+POS"].apply(lambda x:knownauthor[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    
    # determine if novel
    allsig["NOVEL"] = allsig["DISTANCE_TO_KNOWN"].abs() > windowsizekb_for_novel*1000
    
    # determine location
    allsig["LOCATION_OF_KNOWN"]="Unknown"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"]== 0,"LOCATION_OF_KNOWN"] = "Same"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"] > 0 ,"LOCATION_OF_KNOWN"] = "Upstream"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"] < 0 ,"LOCATION_OF_KNOWN"] = "Downstream"

    # if not on same chromosome, distance set to pd.NA
    if sum(allsig["DISTANCE_TO_KNOWN"].abs() > insumstats["POS"].max())>0:
        not_on_same_chromosome = allsig["DISTANCE_TO_KNOWN"].abs() > insumstats["POS"].max()
        allsig.loc[ not_on_same_chromosome ,"DISTANCE_TO_KNOWN"] = pd.NA
        allsig.loc[ not_on_same_chromosome ,"LOCATION_OF_KNOWN"] = "NoneOnThisChr"
        if "SNPID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_ID"] = pd.NA
        if "PUBMEDID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_PUBMED_ID"] = pd.NA
        if "AUTHOR" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_AUTHOR"] = pd.NA

    # drop helper column TCHR+POS
    allsig = allsig.drop(["TCHR+POS"], axis=1)

    if verbose: log.write(" -Identified ",len(allsig)-sum(allsig["NOVEL"])," known vairants in current sumstats...")
    if verbose: log.write(" -Identified ",sum(allsig["NOVEL"])," novel vairants in current sumstats...")
    if verbose: log.write("Finished checking known or novel successfully!")
    gc.collect()
    
    # how to return
    if only_novel is True:
        if output_known is True:
            return allsig.loc[allsig["NOVEL"],:], knownsig
        else:
            return allsig.loc[allsig["NOVEL"],:]
    else:
        if output_known is True:
            return allsig, knownsig
        else:
            return allsig
