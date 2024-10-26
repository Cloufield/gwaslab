import pandas as pd
import numpy as np
import scipy as sp
import gc
from pyensembl import EnsemblRelease
from pyensembl import Genome
from os import path
from gwaslab.util_in_fill_data import fill_p
from gwaslab.g_Log import Log
from gwaslab.bd_common_data import get_chr_to_number
from gwaslab.bd_common_data import get_number_to_chr
from gwaslab.bd_common_data import get_chr_to_NC
from gwaslab.bd_common_data import gtf_to_protein_coding
from gwaslab.bd_common_data import gtf_to_all_gene
from gwaslab.bd_download import check_and_download
from gwaslab.util_ex_gwascatalog import gwascatalog_trait
from gwaslab.qc_fix_sumstats import check_dataframe_shape
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.util_in_correct_winnerscurse import wc_correct
# getsig
# closest_gene
# annogene
# getnovel

def getsig(insumstats,
           id,
           chrom="CHR",
           pos="POS",
           p="P",
           mlog10p="MLOG10P",
           scaled=False,
           use_p=False,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           wc_correction=False,
           build="19",
           source="ensembl",
           gtf_path=None,
           verbose=True):
    """
    Extract the lead variants using a sliding window. P or MLOG10P will be used and converted to SCALEDP for sorting. 
    """
    ##start function with col checking##########################################################
    _start_line = "extract lead variants"
    _end_line = "extracting lead variants"
    _start_cols = [chrom,pos]
    _start_function = ".get_lead()"
    _must_args ={}

    is_enough_info = start_to(sumstats=insumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None
    ############################################################################################

    log.write(" -Processing "+str(len(insumstats))+" variants...", verbose=verbose)
    log.write(" -Significance threshold :", sig_level, verbose=verbose)
    log.write(" -Sliding window size:", str(windowsizekb) ," kb", verbose=verbose)
    
    #load data
    sumstats=insumstats.loc[~insumstats[id].isna(),:].copy()
    
    #convert chrom to int
    if sumstats[chrom].dtype in ["object",str,pd.StringDtype]:
        chr_to_num = get_chr_to_number(out_chr=True,xymt=["X","Y","MT"])
        sumstats[chrom]=sumstats[chrom].map(chr_to_num)
    
    # make sure the dtype is integer
    sumstats[chrom] = np.floor(pd.to_numeric(sumstats[chrom], errors='coerce')).astype('Int64')
    sumstats[pos] = np.floor(pd.to_numeric(sumstats[pos], errors='coerce')).astype('Int64')
    
    #create internal uniqid
    sumstats["__ID"] = range(len(sumstats))
    id = "__ID"

    #extract all significant variants
    ## use mlog10p first
    if use_p==False and (mlog10p in sumstats.columns):
        log.write(" -Using {} for extracting lead variants...".format(mlog10p),verbose=verbose)
        sumstats_sig = sumstats.loc[sumstats[mlog10p]> -np.log10(sig_level),:].copy()
        sumstats_sig.loc[:,"__SCALEDP"] = -pd.to_numeric(sumstats_sig[mlog10p], errors='coerce')
    else:
        #use P
        log.write(" -Using {} for extracting lead variants...".format(p),verbose=verbose)         
        sumstats[p] = pd.to_numeric(sumstats[p], errors='coerce')
        sumstats_sig = sumstats.loc[sumstats[p]<sig_level,:].copy()
        sumstats_sig.loc[:,"__SCALEDP"] = pd.to_numeric(sumstats_sig[p], errors='coerce')
    log.write(" -Found "+str(len(sumstats_sig))+" significant variants in total...", verbose=verbose)

    #sort the coordinates
    sumstats_sig = sumstats_sig.sort_values([chrom,pos])
    if sumstats_sig is None:
        log.write(" -No lead snps at given significance threshold!", verbose=verbose)
        return None
    
    #init 
    sig_index_list=[]
    current_sig_index = False
    current_sig_p = 1
    current_sig_pos = 0
    current_sig_chr = 0
    
    p="__SCALEDP"
    
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
            if  line_number == len(sumstats_sig)-1:
                sig_index_list.append(current_sig_index)
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
    
    log.write(" -Identified "+str(len(sig_index_list))+" lead variants!", verbose=verbose)
    
    # drop internal __SCALEDP
    sumstats_sig = sumstats_sig.drop("__SCALEDP",axis=1)
    
    # extract the lead variants
    output = sumstats_sig.loc[sumstats_sig[id].isin(sig_index_list),:].copy()

    # annotate GENENAME
    if anno is True and len(output)>0:
        log.write(" -Annotating variants using references:{}".format(source), verbose=verbose)
        log.write(" -Annotating variants using references based on genome build:{}".format(build), verbose=verbose)
        
        output = annogene(
               output,
               id=id,
               chrom=chrom,
               pos=pos,
               log=log,
               xymt=xymt,
               build=build,
               source=source,
               gtf_path=gtf_path,
               verbose=verbose)
        
    # drop internal id
    output = output.drop("__ID",axis=1)
    if wc_correction==True:
        log.write(" -Conducting Winner's Curse correction for BETA...", verbose=verbose)
        log.write(" -Referece: Zhong, H., & Prentice, R. L. (2008). Bias-reduced estimators and confidence intervals for odds ratios in genome-wide association studies. Biostatistics, 9(4), 621-634.")
        output["BETA_WC"] = output[["BETA","SE"]].apply(lambda x: wc_correct(x["BETA"],x["SE"],sig_level),axis=1)
    finished(log,verbose,_end_line)
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
            # for refseq , gene names are stored as gene_id, using gene_ids_at_locus instead
            data.gene_names_at_locus = data.gene_ids_at_locus
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
                            return -distance-j,",".join(gene_u).strip(",")
                        elif len(gene_d)>0:
                            return distance+j,",".join(gene_d).strip(",")
                elif len(gene_u)>0:                    
                    # if found gene uptream
                    distance = (i-1)*step
                    for j in range(0,step,1):
                        gene_u2 = data.gene_names_at_locus(contig=contig, position=position-distance-j)
                        if len(gene_u2)>0:
                            return -distance-j,",".join(gene_u).strip(",")
                elif len(gene_d)>0:
                    # if found gene downstream
                    distance = (i-1)*step
                    for j in range(0,step,1):
                        gene_d2 = data.gene_names_at_locus(contig=contig, position=position+distance+j)
                        if len(gene_d2)>0:
                            return distance+j,",".join(gene_d).strip(",")
                i+=1
                # increase i by 1
            return distance,"intergenic"
        else:
            return 0,",".join(gene).strip(",")


def annogene(
           insumstats,
           id,
           chrom="CHR",
           pos="POS",
           log=Log(),
           xymt=["X","Y","MT"],
           build="19",
           source="ensembl",
           gtf_path=None,
           verbose=True):
    
    log.write("Start to annotate variants with nearest gene name(s)...", verbose=verbose)
    output = insumstats.copy()
    
    if source == "ensembl":
        if build=="19":
            #data = EnsemblRelease(75)
            log.write(" -Assigning Gene name using ensembl_hg19_gtf for protein coding genes", verbose=verbose)
            #zcat Homo_sapiens.GRCh37.75.gtf.gz| 
            #grep -E 'processed_transcript|protein_coding|_gene' 
            #| gzip >Homo_sapiens.GRCh37.75.processed.chr.gtf.gz     
            
            #gtf_path = check_and_download("ensembl_hg19_gtf_protein_coding")
            if gtf_path is None:
                gtf_path = check_and_download("ensembl_hg19_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)

            gtf_db_path = gtf_path[:-2]+"db"
            
            data = Genome(
                reference_name='GRCh37',
                annotation_name='Ensembl',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,chrom=chrom,pos=pos,source=source), axis=1)), 
                index=output.index).values
        elif build=="38":
            log.write(" -Assigning Gene name using ensembl_hg38_gtf for protein coding genes", verbose=verbose)
            #gtf_path = check_and_download("ensembl_hg38_gtf_protein_coding")
            if gtf_path is None:
                gtf_path = check_and_download("ensembl_hg38_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            
            gtf_db_path = gtf_path[:-2]+"db"
            data = Genome(
                reference_name='GRCh38',
                annotation_name='Ensembl',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,chrom=chrom,pos=pos,source=source), axis=1)), 
                index=output.index).values
    
    if source == "refseq":
        if build=="19":
            log.write(" -Assigning Gene name using NCBI refseq latest GRCh37 for protein coding genes", verbose=verbose)
            #gtf_path = check_and_download("refseq_hg19_gtf_protein_coding")
            if gtf_path is None:
                gtf_path = check_and_download("refseq_hg19_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            
            gtf_db_path = gtf_path[:-2]+"db"
            data = Genome(
                reference_name='GRCh37',
                annotation_name='Refseq',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,chrom=chrom,pos=pos,source=source,build=build), axis=1)), 
                index=output.index).values
        elif build=="38":
            log.write(" -Assigning Gene name using NCBI refseq latest GRCh38 for protein coding genes", verbose=verbose)
            #gtf_path = check_and_download("refseq_hg38_gtf_protein_coding")
            if gtf_path is None:
                gtf_path = check_and_download("refseq_hg38_gtf")
                gtf_path = gtf_to_protein_coding(gtf_path,log=log,verbose=verbose)
            else:
                log.write(" -Using user-provided gtf:{}".format(gtf_path))
                gtf_path = gtf_to_all_gene(gtf_path,log=log,verbose=verbose)
            
            gtf_db_path = gtf_path[:-2]+"db"
            data = Genome(
                reference_name='GRCh38',
                annotation_name='Refseq',
                gtf_path_or_url=gtf_path)
            if path.isfile(gtf_db_path) is False:
                data.index()
            output.loc[:,["LOCATION","GENE"]] = pd.DataFrame(
                list(output.apply(lambda x:closest_gene(x,data=data,chrom=chrom,pos=pos,source=source,build=build), axis=1)), 
                index=output.index).values
    log.write("Finished annotating variants with nearest gene name(s) successfully!", verbose=verbose)
    return output

def getnovel(insumstats,
           id,
           chrom,
           pos,
           p,
           use_p=False,
           known=False,
           efo=False,
           only_novel=False,
           group_key=None,
           if_get_lead = True,
           windowsizekb_for_novel=1000,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           wc_correction=False,
           build="19",
           source="ensembl",
           gwascatalog_source="NCBI",
           output_known=False,
           verbose=True):
    ##start function with col checking##########################################################
    _start_line = "check if lead variants are known"
    _end_line = "checking if lead variants are known"
    _start_cols = [chrom,pos]
    _start_function = ".get_novel()"
    _must_args ={}

    is_enough_info = start_to(sumstats=insumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None
    ############################################################################################
    
    if if_get_lead == True:
        allsig = getsig(insumstats=insumstats,
            id=id,chrom=chrom,pos=pos,p=p,use_p=use_p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
            xymt=xymt,anno=anno,build=build,wc_correction=wc_correction, source=source,verbose=verbose)
    else:
        allsig = insumstats.copy()

    ############################################################################################
    knownsig = pd.DataFrame()
    if efo != False:
        if type(efo) is not list:
            log.write("Start to retrieve data using EFO: {}...".format(efo), verbose=verbose)
            known_Sumstats = gwascatalog_trait(efo,source=gwascatalog_source,sig_level=sig_level,verbose=verbose,log=log)
            knownsig = known_Sumstats.data.copy()
        else:
            knownsig=pd.DataFrame()
            log.write("Start to retrieve data using {} EFOs: {}...".format(len(efo),efo), verbose=verbose)
            for single_efo in efo:
                known_Sumstats = gwascatalog_trait(single_efo,source=gwascatalog_source,sig_level=sig_level,verbose=verbose,log=log)
                known_Sumstats.data["EFOID"] = single_efo
                knownsig = pd.concat([known_Sumstats.data, knownsig],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["POS"] = knownsig["POS"].astype("Int64")
        log.write(" -Retrieved {} associations from GWAS catalog.".format(len(knownsig)), verbose=verbose)
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
    
    if group_key is not None:
        if (group_key not in allsig.columns) or (group_key not in knownsig.columns):
            raise ValueError("Please check if group_key is in both sumstats and list of known associations.")

    # create helper column TCHR+POS for knownsig and all sig
     ############################################################################################
    maxpos = insumstats["POS"].max()
    big_number = determine_big_number(maxpos)
    knownsig = add_tchr_pos(knownsig, chrom, pos, big_number)
    allsig = add_tchr_pos(allsig, chrom, pos, big_number)
    ############################################################################################   
    #sorting
    allsig = allsig.sort_values(by="TCHR+POS",ignore_index=True)
    knownsig = knownsig.sort_values(by="TCHR+POS",ignore_index=True)
    ############################################################################################
    if group_key is not None:
        number_of_groups_allsig = allsig[group_key].nunique()
        number_of_groups_known = knownsig[group_key].nunique()
        log.write(" -Number of groups in sumstats:{}".format(number_of_groups_allsig), verbose=verbose)
        log.write(" -Number of groups in reference:{}".format(number_of_groups_known), verbose=verbose)

    log.write(" -Lead variants in known loci:",len(knownsig), verbose=verbose)
    log.write(" -Checking the minimum distance between identified lead variants and provided known variants...", verbose=verbose)
    
    ############################################################################################
    if group_key is None:
        # get distance
        allsig = determine_distance(allsig, knownsig)
        # get other info 
        allsig = fill_meta_info_for_known(allsig, knownsig)
        ############################################################################################
        # determine if novel
        allsig = determine_novel(allsig, windowsizekb_for_novel)
        # determine location
        allsig = determine_location(allsig)
        # if not on same chromosome, distance set to pd.NA
        allsig = determine_if_same_chromosome(allsig, knownsig, maxpos)
        ############################################################################################
    else:
        #groups1 = set(allsig[group_key].unique())
        #groups2 = set(knownsig[group_key].unique())
        #common_group = groups1.intersection(groups2)
        
        #allsig_no_group = allsig.loc[~allsig[group_key].isin(common_group),:].copy()
        allsig_group = pd.DataFrame()

        for key in allsig[group_key].unique():
            allsig_single_group = allsig.loc[allsig[group_key]==key,:].copy()
            knownsig_single_group = knownsig.loc[knownsig[group_key]==key,:].copy()

            #if len(allsig_single_group) >0 and len(knownsig_single_group) >0:
            allsig_single_group = determine_distance(allsig_single_group, knownsig_single_group)
            # get other info 
            allsig_single_group = fill_meta_info_for_known(allsig_single_group, knownsig_single_group)
            
            # determine if novel
            allsig_single_group = determine_novel(allsig_single_group, windowsizekb_for_novel)
            
            # determine location
            allsig_single_group = determine_location(allsig_single_group)
            
            # if not on same chromosome, distance set to pd.NA
            allsig_single_group = determine_if_same_chromosome(allsig_single_group, knownsig_single_group, maxpos) 
            
            allsig_group = pd.concat([allsig_group, allsig_single_group], ignore_index=True)
        
        allsig = allsig_group
        #pd.concat([allsig_no_group, allsig_group], ignore_index=True)
        
    # drop helper column TCHR+POS
    allsig = allsig.drop(["TCHR+POS"], axis=1)
    
    try:
        allsig = allsig.where(~pd.isna(allsig), pd.NA)
    except:
        pass

    log.write(" -Identified ",len(allsig)-sum(allsig["NOVEL"])," known vairants in current sumstats...", verbose=verbose)
    log.write(" -Identified ",sum(allsig["NOVEL"])," novel vairants in current sumstats...", verbose=verbose)
    
    finished(log,verbose,_end_line)
    
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
##################################################################################################################################################################################################


def _check_cis(insumstats,
           id,
           chrom,
           pos,
           p,
           use_p=False,
           known=False,
           group_key=None,
           if_get_lead = False,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           build="19",
           source="ensembl",
           verbose=True):
    ##start function with col checking##########################################################
    _start_line = "check if variants are in cis or trans regions"
    _end_line = "checking if variants are in cis or trans regions"
    _start_cols = [chrom,pos, group_key]
    _start_function = ".check_cis()"
    _must_args ={}

    is_enough_info = start_to(sumstats=insumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None
    ############################################################################################
    
    if if_get_lead == True:
        allsig = getsig(insumstats=insumstats,
            id=id,chrom=chrom,pos=pos,p=p,use_p=use_p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
            xymt=xymt,anno=anno,build=build, source=source,verbose=verbose)
    else:
        allsig = insumstats.copy()

    ############################################################################################
    knownsig = pd.DataFrame()
    if type(known) is pd.DataFrame:
        knownsig_2 = known.copy()
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["START"] = knownsig["START"].astype("Int64")
        knownsig["END"] = knownsig["END"].astype("Int64")
    elif type(known) is str:
        knownsig_2 = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig["CHR"] = knownsig["CHR"].astype("Int64")
        knownsig["START"] = knownsig["START"].astype("Int64")
        knownsig["END"] = knownsig["END"].astype("Int64")
    
    if len(knownsig)<1:
        raise ValueError("Please input a dataframe of gene list with GENE, CHR, START, END.")
    
    if group_key is not None:
        if group_key not in knownsig.columns:
            raise ValueError("Please check if group_key is in both sumstats and list of known associations.")

    ############################################################################################
    if group_key is not None:
        number_of_groups_allsig = allsig[group_key].nunique()
        number_of_groups_known = knownsig[group_key].nunique()
        log.write(" -Number of groups in sumstats:{}".format(number_of_groups_allsig), verbose=verbose)
        log.write(" -Number of groups in reference:{}".format(number_of_groups_known), verbose=verbose)

    log.write(" -Checking if variants in cis/trans regions grouped by {}...".format(group_key), verbose=verbose)
    log.write(" -Window size in kb adding to start and end: {}...".format(windowsizekb), verbose=verbose)
    ############################################################################################
    #convert to  a dict
    reference_dict = {}
    for index,row in knownsig.iterrows():
        reference_dict[row[group_key]] = (row["CHR"], row["START"], row["END"] )
    ############################################################################################
    try:
        no_reference_avaialble = allsig.loc[~allsig[group_key].isin(reference_dict.keys()),group_key]
        if len(no_reference_avaialble)>0:
            log.write(" -Groups not in reference: {}".format( ",".join(no_reference_avaialble.unique())), verbose=verbose)
    except:
        pass

    #allsig["CIS/TRANS"] = allsig.apply(lambda x: determine_if_cis(x, group_key,windowsizekb, reference_dict), axis=1)
    cis_tuples = allsig.apply(lambda x: determine_if_cis2(x, group_key,windowsizekb, reference_dict), axis=1)
    allsig[["CIS/TRANS","REF_CHR","REF_START","REF_END"]] = pd.DataFrame(cis_tuples.tolist(), index=allsig.index)

    try:
        allsig = allsig.where(~pd.isna(allsig), pd.NA)
    except:
        pass
    
    try:
        number_of_cis = sum(allsig["CIS/TRANS"] == "Cis")
        number_of_trans = sum(allsig["CIS/TRANS"] == "Trans")
        number_of_noreference = sum(allsig["CIS/TRANS"] == "NoReference")
        log.write (" -Number of Cis variants: {}".format(number_of_cis),verbose=verbose)
        log.write (" -Number of Trans variants: {}".format(number_of_trans),verbose=verbose)
        log.write (" -Number of NoReference variants: {}".format(number_of_noreference),verbose=verbose)
    except:
        pass

    finished(log,verbose,_end_line)
    
    return allsig

###################################################################################################################################################################################################


def determine_big_number(maxpos, big_number = 1000000000):
    for i in range(7):
        if maxpos*10 >  big_number:
            big_number = int(big_number * 10)
        else:
            break
    return big_number

def add_tchr_pos(df, chrom, pos, big_number):
    df["TCHR+POS"]=df[chrom]*big_number + df[pos]
    return df

def fill_meta_info_for_known(allsig, knownsig):
    if len(allsig)==0 or len(knownsig)==0: return allsig
    if "SNPID" in knownsig.columns:
        knownids=knownsig["SNPID"].values
    if "PUBMEDID" in knownsig.columns:
        knownpubmedids=knownsig["PUBMEDID"].values
    if "AUTHOR" in knownsig.columns:
        knownauthor=knownsig["AUTHOR"].values
    if "EFOID" in knownsig.columns:
        knownefo=knownsig["EFOID"].values

    if "SNPID" in knownsig.columns:
        allsig["KNOWN_ID"] = allsig["TCHR+POS"].apply(lambda x:knownids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])    
    if "PUBMEDID" in knownsig.columns:
        allsig["KNOWN_PUBMED_ID"] = allsig["TCHR+POS"].apply(lambda x:knownpubmedids[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    if "AUTHOR" in knownsig.columns:
        allsig["KNOWN_AUTHOR"] = allsig["TCHR+POS"].apply(lambda x:knownauthor[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])
    if "EFOID" in knownsig.columns:
        allsig["KNOWN_EFOID"] = allsig["TCHR+POS"].apply(lambda x:knownefo[np.argmin(np.abs(knownsig["TCHR+POS"]-x))])    
    return allsig

def determine_if_cis(x, group_key,windowsizekb, reference_dict):
    if x[group_key] in reference_dict.keys():
        is_same_chr = str(reference_dict[x[group_key]][0]) == str(x["CHR"])
        is_large_than_start = int(reference_dict[x[group_key]][1]) - windowsizekb*1000 <= x["POS"]
        is_smaller_than_end = int(reference_dict[x[group_key]][2]) + windowsizekb*1000 >= x["POS"]               
        
        if  is_same_chr and is_large_than_start  and is_smaller_than_end:
            return "Cis"
        else: 
            return "Trans"
    else:
        return "NoReference"

def determine_if_cis2(x, group_key,windowsizekb, reference_dict):
    if x[group_key] in reference_dict.keys():
        is_same_chr = str(reference_dict[x[group_key]][0]) == str(x["CHR"])
        is_large_than_start = int(reference_dict[x[group_key]][1]) - windowsizekb*1000 <= x["POS"]
        is_smaller_than_end = int(reference_dict[x[group_key]][2]) + windowsizekb*1000 >= x["POS"]               
        
        if  is_same_chr and is_large_than_start  and is_smaller_than_end:
            return "Cis", int(reference_dict[x[group_key]][0]), int(reference_dict[x[group_key]][1]), int(reference_dict[x[group_key]][2])
        else: 
            return "Trans", int(reference_dict[x[group_key]][0]), int(reference_dict[x[group_key]][1]), int(reference_dict[x[group_key]][2])
    else:
        return "NoReference", pd.NA, pd.NA, pd.NA


def determine_distance(allsig, knownsig):
    if len(allsig)==0: 
        return allsig
    if len(knownsig)==0:
        allsig["DISTANCE_TO_KNOWN"] = pd.NA
        return allsig
    allsig["DISTANCE_TO_KNOWN"] = allsig["TCHR+POS"].apply(lambda x:min(knownsig["TCHR+POS"]-x, key=abs))
    return allsig

def determine_novel(allsig, windowsizekb_for_novel):
    if len(allsig)==0 or "DISTANCE_TO_KNOWN" not in allsig.columns:
        return allsig
    allsig["NOVEL"] = allsig["DISTANCE_TO_KNOWN"].abs() > windowsizekb_for_novel*1000
    allsig.loc[allsig["DISTANCE_TO_KNOWN"].isna(), "NOVEL"] = True
    return allsig

def determine_location(allsig):
    allsig["LOCATION_OF_KNOWN"]="NoReference"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"]== 0,"LOCATION_OF_KNOWN"] = "Same"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"] > 0 ,"LOCATION_OF_KNOWN"] = "Upstream"
    allsig.loc[ allsig["DISTANCE_TO_KNOWN"] < 0 ,"LOCATION_OF_KNOWN"] = "Downstream"
    return allsig

def determine_if_same_chromosome(allsig, knownsig, maxpos):
    if sum(allsig["DISTANCE_TO_KNOWN"].abs() > maxpos)>0:
        not_on_same_chromosome = allsig["DISTANCE_TO_KNOWN"].abs() > maxpos
        allsig.loc[ not_on_same_chromosome ,"DISTANCE_TO_KNOWN"] = pd.NA
        allsig.loc[ not_on_same_chromosome ,"LOCATION_OF_KNOWN"] = "NoneOnThisChr"
        if "SNPID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_ID"] = pd.NA
        if "PUBMEDID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_PUBMED_ID"] = pd.NA
        if "AUTHOR" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_AUTHOR"] = pd.NA
        if "EFOID" in knownsig.columns:
            allsig.loc[ not_on_same_chromosome ,"KNOWN_EFOID"] = pd.NA
    return allsig

def _check_novel_set(insumstats,
           id,
           chrom,
           pos,
           p,
           use_p=False,
           known=False,
           group_key=None,
           snpset="SNPSET",
           snpid="SNPID",
           if_get_lead = False,
           windowsizekb=500,
           sig_level=5e-8,
           log=Log(),
           xymt=["X","Y","MT"],
           anno=False,
           build="19",
           source="ensembl",
           verbose=True):
    
    ##start function with col checking##########################################################
    _start_line = "check if variant sets are overlapping with those in reference file"
    _end_line = "checking if variant sets are overlapping with those in reference file"
    _start_cols = [chrom,pos, group_key]
    _start_function = ".check_cis()"
    _must_args ={}

    is_enough_info = start_to(sumstats=insumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None
    ############################################################################################
    
    if if_get_lead == True:
        allsig = getsig(insumstats=insumstats,
            id=id,chrom=chrom,pos=pos,p=p,use_p=use_p,windowsizekb=windowsizekb,sig_level=sig_level,log=log,
            xymt=xymt,anno=anno,build=build, source=source,verbose=verbose)
    else:
        allsig = insumstats.copy()

    ############################################################################################
    knownsig = pd.DataFrame()
    if type(known) is pd.DataFrame:
        knownsig_2 = known.copy()
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig[snpid] = knownsig[snpid].astype("string")
        knownsig[snpset] = knownsig[snpset].astype("string")
        knownsig[group_key] = knownsig[group_key].astype("string")
    elif type(known) is str:
        knownsig_2 = pd.read_csv(known,sep="\s+",dtype={"CHR":"Int64","POS":"Int64"})
        knownsig = pd.concat([knownsig, knownsig_2],ignore_index=True)
        knownsig[snpid] = knownsig[snpid].astype("string")
        knownsig[snpset] = knownsig[snpset].astype("string")
        knownsig[group_key] = knownsig[group_key].astype("string")
    
    if len(knownsig)<1:
        raise ValueError("Please input a dataframe of gene list with GENE, CHR, START, END.")
    
    if group_key is not None:
        if group_key not in knownsig.columns:
            raise ValueError("Please check if group_key is in both sumstats and list of known associations.")

    ############################################################################################
    if group_key is not None:
        number_of_groups_allsig = allsig[group_key].nunique()
        number_of_groups_known = knownsig[group_key].nunique()
        log.write(" -Number of groups in sumstats:{}".format(number_of_groups_allsig), verbose=verbose)
        log.write(" -Number of groups in reference:{}".format(number_of_groups_known), verbose=verbose)

    log.write(" -Checking if variants in cis/trans regions grouped by {}...".format(group_key), verbose=verbose)
    
    ############################################################################################
    #convert to  a dict
    reference_dict = {}

    for index,row in knownsig.iterrows():
        if row[group_key] in reference_dict.keys():
            if row[snpset] in reference_dict[row[group_key]].keys():
                reference_dict[row[group_key]][row[snpset]].add(row[snpid])
            else:
                reference_dict[row[group_key]][row[snpset]] = set([row[snpid]])
        else:
            reference_dict[row[group_key]] = {row[snpset]:set([row[snpid]])}
    ############################################################################################
    
    try:
        no_reference_avaialble = allsig.loc[~allsig[group_key].isin(reference_dict.keys()),group_key]
        if len(no_reference_avaialble)>0:
            log.write(" -Groups not in reference: {}".format( ",".join(no_reference_avaialble)), verbose=verbose)
    except:
        pass

    log.write(" -Checking if variants are in reference variant sets...", verbose=verbose)    
    known_list = allsig.apply(lambda x: check_overlap(x,snpid, group_key,reference_dict), axis=1)
    
    allsig["KNOWN_SET"] = known_list.str[0]
    allsig["KNOWN_VARIANT"] = known_list.str[1]

    back_dict={}
    for i in allsig[group_key].unique():
        back_dict[i] ={}
        for j in allsig.loc[allsig[group_key]==i,snpset].unique():
            back_dict[i][j] =set()
            for index, row in allsig.loc[(allsig[group_key]==i) & (allsig[snpset]==j) & (~allsig["KNOWN_SET"].isna()),:].iterrows():
                back_dict[i][j].add("{}-{}-{}".format(row[group_key], row["KNOWN_SET"],row["KNOWN_VARIANT"]))

    allsig["KNOWN_SET_VARIANT"] = allsig.apply(lambda x: assign_set_variant(x,group_key,snpset,back_dict), axis=1)
    
    finished(log,verbose,_end_line)
    
    return allsig

def check_overlap(x,snpid, group_key,reference_dict):
    if x[group_key] in reference_dict.keys():
        for key, value in reference_dict[x[group_key]].items():
            if x[snpid] in value:
                return key, x[snpid]
    return pd.NA, pd.NA, 

def assign_set_variant(x,group_key,snpset,back_dict):
    if x[group_key] in back_dict.keys():
        if x[snpset] in back_dict[x[group_key]].keys():
            if len(back_dict[x[group_key]][x[snpset]]) >0:
                return back_dict[x[group_key]][x[snpset]]
    return pd.NA