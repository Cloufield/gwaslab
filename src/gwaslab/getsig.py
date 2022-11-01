import pandas as pd
import numpy as np
import scipy as sp
from gwaslab.Log import Log
from gwaslab.CommonData import get_chr_to_number
from gwaslab.CommonData import get_number_to_chr
from gwaslab.CommonData import get_chr_to_NC
from gwaslab.CommonData import gtf_index
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
    if anno is True:
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
            if verbose:log.write(" -Assigning Gene name using built-in Ensembl Release",75 , " (hg19)")
            #zcat Homo_sapiens.GRCh37.75.gtf.gz| 
            #grep -E 'processed_transcript|protein_coding|_gene' 
            #| gzip >Homo_sapiens.GRCh37.75.processed.chr.gtf.gz     
            gtf_path = path.dirname(__file__) + '/data/Ensembl/release75/Homo_sapiens.GRCh37.75.protein_coding.gtf.gz'
            gtf_db_path = path.dirname(__file__) + '/data/Ensembl/release75/Homo_sapiens.GRCh37.75.protein_coding.gtf.db'
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
            if verbose:log.write(" -Assigning Gene name using built-in Ensembl Release",107 , " (hg38)")
            gtf_path = path.dirname(__file__) + '/data/Ensembl/release107/Homo_sapiens.GRCh38.107.protein_coding.chr.gtf.gz'
            gtf_db_path = path.dirname(__file__) + '/data/Ensembl/release107/Homo_sapiens.GRCh38.107.protein_coding.chr.gtf.db'
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
            if verbose:log.write(" -Assigning Gene name using built-in NCBI refseq latest GRCh37")
            gtf_path = path.dirname(__file__) + '/data/RefSeq/GRCh37/GRCh37_latest_genomic.protein_coding.gtf.gz'
            gtf_db_path = path.dirname(__file__) + '/data/RefSeq/GRCh37/GRCh37_latest_genomic.protein_coding.gtf.db'
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
            gtf_path = path.dirname(__file__) + '/data/RefSeq/GRCh38/GRCh38_latest_genomic.protein_coding.gtf.gz'
            gtf_db_path = path.dirname(__file__) + '/data/RefSeq/GRCh38/GRCh38_latest_genomic.protein_coding.gtf.db'
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
           known=None,
           efo=None,
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
           verbose=True):
    
    allsig = getsig(insumstats=insumstats,
           id=id,chrom=chrom,pos=pos,p=p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
           xymt=xymt,anno=anno,build=build, source=source,verbose=verbose)
    
    allsig["TCHR+POS"]=allsig[chrom]*1000000000 + allsig[pos]
    
    if known is None and efo is not None:
        known = gwascatalog_trait(efo,source=gwascatalog_source,sig_level=sig_level,verbose=verbose,log=log)
        knownsig = known.data.copy()
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
    elif type(known) is pd.DataFrame:
        knownsig = known.copy()
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
    elif type(known) is str:
        knownsig = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
    else:
        raise ValueError("Please input a dataframe of known loci or efo code")
    
    if "SNPID" in knownsig.columns:
        knownids=knownsig["SNPID"].values
    if "PUBMEDID" in knownsig.columns:
        knownpubmedids=knownsig["PUBMEDID"].values
    if "AUTHOR" in knownsig.columns:
        knownauthor=knownsig["AUTHOR"].values
        
    if verbose: log.write("Start to check if lead variants are known...")
    knownsig["TCHR+POS"]=knownsig[chrom]*1000000000 + knownsig[pos]
    
    if verbose: log.write(" -Lead variants in known loci:",len(knownsig))
    if verbose: log.write(" -Checking the minimum distance between identified lead variants and provided known variants...")
    
    #sorting
    allsig = allsig.sort_values(by="TCHR+POS",ignore_index=True)
    knownsig = knownsig.sort_values(by="TCHR+POS",ignore_index=True)

    allsig["DISTANCE_TO_KNOWN"] = allsig["TCHR+POS"].apply(lambda x:np.min(np.abs(knownsig["TCHR+POS"]-x)))
    
    if "SNPID" in knownsig.columns:
        allsig["KNOWN_ID"] = allsig["TCHR+POS"].apply(lambda x:knownids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])    
    if "PUBMEDID" in knownsig.columns:
        allsig["KNOWN_PUBMED_ID"] = allsig["TCHR+POS"].apply(lambda x:knownpubmedids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    if "AUTHOR" in knownsig.columns:
        allsig["KNOWN_AUTHOR"] = allsig["TCHR+POS"].apply(lambda x:knownauthor[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    allsig["NOVEL"] = allsig["DISTANCE_TO_KNOWN"] > windowsizekb_for_novel*1000
    
    if verbose: log.write(" -Identified ",len(allsig)-sum(allsig["NOVEL"])," known vairants in current sumstats...")
    if verbose: log.write(" -Identified ",sum(allsig["NOVEL"])," novel vairants in current sumstats...")
    if verbose: log.write("Finished checking known or novel successfully!")
    gc.collect()
    if only_novel is True:
        return allsig.loc[allsig["NOVEL"],:]
    else:
        return allsig
