import pandas as pd
import numpy as np
import subprocess
from gwaslab.g_Log import Log
import os
from gwaslab.g_version import _checking_plink_version
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished

def _process_plink_input_files(chrlist, 
                               bfile=None, 
                               pfile=None, 
                               vcf=None, 
                               bgen=None,  
                               sample=None, 
                               n_cores=1, 
                               plink_log="", 
                               log=Log(), 
                               overwrite=False, 
                               bgen_mode="ref-first", 
                               convert="bfile", 
                               memory=None, 
                               load_bim=False,
                               plink="plink",
                               plink2="plink2"):
    """
    Process input files (bfile,pfile,vcf,bgen) to either PLINK1 bed/bim/fam or PLINK2 pgen/psam/pvar. 
    
    Returns:
    
    ref_file_prefix : prefix for either bfile or pfile.
    plink_log : if plink was used, return the log. Otherwise, return an empty string.
    ref_bims : if load_bim is True, return bim files as a list of pd.DataFrame. Otherwise, empty list.
    filetype : either bfile or pfile.
    """

    ref_bims = []
    # No change: bfile, pfile
    # Need to convert to bfile or pfile
    # file prefix : single or with @
    #
    filetype, ref_file_prefix, is_wild_card = _check_file_type(bfile,  pfile,  vcf,  bgen, sample)

    if filetype == "bfile": 
        ref_file_prefix, ref_bims = _process_bfile(chrlist=chrlist,
                                                   ref_file_prefix=ref_file_prefix, 
                                                   ref_bims=ref_bims, 
                                                   is_wild_card=is_wild_card, 
                                                   log=log, 
                                                   load_bim=load_bim)

    elif filetype == "pfile": 
        ref_file_prefix, ref_bims = _process_pfile(chrlist=chrlist,
                                                   ref_file_prefix=ref_file_prefix, 
                                                   ref_bims=ref_bims, 
                                                   is_wild_card=is_wild_card, 
                                                   log=log, 
                                                   load_bim=load_bim)

    elif filetype == "vcf": 
        ref_file_prefix, plink_log, ref_bims = _process_vcf(ref_file_prefix=ref_file_prefix, 
                                                            chrlist=chrlist, 
                                                            ref_bims=ref_bims, 
                                                            is_wild_card=is_wild_card, 
                                                            log=log, 
                                                            plink_log=plink_log, 
                                                            n_cores=n_cores, 
                                                            convert=convert, 
                                                            memory=memory, 
                                                            overwrite=overwrite, 
                                                            load_bim=load_bim,
                                                            plink=plink,
                                                            plink2=plink2)
        filetype = convert
    elif filetype == "bgen": 
        ref_file_prefix, plink_log, ref_bims = _process_bgen(ref_file_prefix=ref_file_prefix, 
                                                            chrlist=chrlist, 
                                                            bgen_mode=bgen_mode,
                                                            sample=sample,
                                                            ref_bims=ref_bims, 
                                                            is_wild_card=is_wild_card, 
                                                            log=log, 
                                                            plink_log=plink_log, 
                                                            n_cores=n_cores, 
                                                            convert=convert, 
                                                            memory=memory, 
                                                            overwrite=overwrite, 
                                                            load_bim=load_bim,
                                                            plink=plink,
                                                            plink2=plink2)
        filetype = convert
    return ref_file_prefix, plink_log, ref_bims, filetype

def _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log):
    bim_path =bpfile_prefix+".bim"
    single_bim = pd.read_csv(bim_path,
                             sep="\s+",
                             usecols=[0,1,3,4,5],
                             header=None,
                             dtype={1:"string",0:"category", 3:"int", 4:"string", 5:"string"}).rename(columns={1:"SNPID",0:"CHR_bim",3:"POS_bim",4:"EA_bim",5:"NEA_bim"})
    log.write("   -Variants in ref file: {}".format(len(single_bim))) 
    ref_bims.append(single_bim)
    return ref_bims

def _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log):
    if os.path.exists(bpfile_prefix+".pvar"):
        bim_path =bpfile_prefix+".pvar"
    elif os.path.exists(bpfile_prefix+".pvar.zst"):
        bim_path =bpfile_prefix+".pvar.zst"
    single_bim = pd.read_csv(bim_path,
                             sep="\s+",
                             usecols=[0,1,2,3,4],
                             header=None,
                             comment="#",
                             dtype={2:"string",0:"category", 1:"int", 3:"string", 4:"string"}).rename(columns={2:"SNPID",0:"CHR_bim",1:"POS_bim",3:"EA_bim",4:"NEA_bim"})
    log.write("   -Variants in ref file: {}".format(len(single_bim))) 
    ref_bims.append(single_bim)
    return ref_bims

def _check_file_type(bfile=None,  
                     pfile=None,  
                     vcf=None,  
                     bgen=None, 
                     sample=None):
    
    is_wild_card = False
    if bfile is None:
        if pfile is None:
            if vcf is None:
                if (bgen is None) or (sample is None):
                    raise ValueError("You need to provide one from bfile, pfile, bgen, vcf; for bgen, sample file is required.")
                else:
                    if "@" in bgen:
                        is_wild_card = True
                    return "bgen",  bgen, is_wild_card
            else:
                if "@" in vcf:
                    is_wild_card = True
                return "vcf",  vcf, is_wild_card
        else:
            if "@" in pfile:
                is_wild_card = True
            return "pfile",  pfile, is_wild_card
    else:
        if "@" in bfile:
            is_wild_card = True
        return "bfile",  bfile, is_wild_card
    
def _process_bfile(chrlist ,ref_file_prefix, ref_bims, is_wild_card, log, load_bim=False):
    if is_wild_card==False:
        is_bim_exist = os.path.exists(ref_file_prefix+".bim")
        is_bed_exist = os.path.exists(ref_file_prefix+".bed")
        is_fam_exist = os.path.exists(ref_file_prefix+".fam")
        if not (is_bim_exist and is_bed_exist  and is_fam_exist):
            raise ValueError("PLINK bfiles are missing : {} ".format(ref_file_prefix))
        if load_bim==True:  
            ref_bims = _load_single_bim_to_ref_bims(ref_file_prefix, ref_bims, log)    
        log.write(" -Single PLINK bfile as LD reference panel: {}".format(ref_file_prefix))   

    else:
        for i in chrlist:
            single_chr_ref_file_prefix = ref_file_prefix.replace("@",str(i))
            is_bim_exist = os.path.exists(single_chr_ref_file_prefix+".bim")
            is_bed_exist = os.path.exists(single_chr_ref_file_prefix+".bed")
            is_fam_exist = os.path.exists(single_chr_ref_file_prefix+".fam")
            if not (is_bim_exist and is_bed_exist  and  is_fam_exist):
                raise ValueError("PLINK bfiles for CHR {} are missing...".format(i))
            if load_bim==True:  
                ref_bims = _load_single_bim_to_ref_bims(single_chr_ref_file_prefix, ref_bims, log)    
        log.write(" -Split PLINK bfiles for each CHR as LD reference panel: {}".format(ref_file_prefix))  

    return ref_file_prefix, ref_bims

def _process_pfile(chrlist, ref_file_prefix, ref_bims, is_wild_card, log, load_bim=False):
    if is_wild_card==False:
        is_bim_exist = os.path.exists(ref_file_prefix+".pvar")
        is_bed_exist = os.path.exists(ref_file_prefix+".pgen")
        is_fam_exist = os.path.exists(ref_file_prefix+".psam")
        if not (is_bim_exist and is_bed_exist  and is_fam_exist):
            raise ValueError("PLINK pfiles are missing : {} ".format(ref_file_prefix))
        if load_bim==True:  
            ref_bims = _load_single_pvar_to_ref_bims(ref_file_prefix, ref_bims, log) 
        log.write(" -Single PLINK bfile as LD reference panel: {}".format(ref_file_prefix))   
    else:
        for i in chrlist:
            single_chr_ref_file_prefix = ref_file_prefix.replace("@",str(i))
            is_bim_exist = os.path.exists(single_chr_ref_file_prefix+".pvar")
            is_bed_exist = os.path.exists(single_chr_ref_file_prefix+".pgen")
            is_fam_exist = os.path.exists(single_chr_ref_file_prefix+".psam")
            if not (is_bim_exist and is_bed_exist  and  is_fam_exist):
                raise ValueError("PLINK pfiles for CHR {} are missing...".format(i))
            if load_bim==True:  
                ref_bims = _load_single_pvar_to_ref_bims(ref_file_prefix, ref_bims, log)  
        log.write(" -Split PLINK pfiles for each CHR as LD reference panel: {}".format(ref_file_prefix))  
  
    return ref_file_prefix, ref_bims

def _process_vcf(ref_file_prefix, 
                 chrlist, 
                 ref_bims, 
                 is_wild_card, 
                 log, 
                 plink_log, 
                 n_cores=1, 
                 convert="bfile", 
                 memory=None, 
                 overwrite=False, 
                 load_bim=False,
                 plink="plink",
                 plink2="plink2"):
    log.write(" -Processing VCF : {}...".format(ref_file_prefix))
    
    #check plink version
    log = _checking_plink_version(plink2=plink2,log=log)
    
    # file path prefix to return
    if is_wild_card==True:
        ref_file_prefix_converted = ref_file_prefix.replace(".vcf.gz","")
    else:
        ref_file_prefix_converted = ref_file_prefix.replace(".vcf.gz","") + ".@"

    for i in chrlist:
    # for each chr    
        log.write("  -Processing VCF for CHR {}...".format(i))
        
        # if multiple bfiles
        if is_wild_card==True:
            vcf_to_load = ref_file_prefix.replace("@",str(i))
             # output bpfile prefix for each chr : chr to chr
            bpfile_prefix = vcf_to_load.replace(".vcf.gz","")
        else:
            vcf_to_load = ref_file_prefix
             # output bpfile prefix for each chr : all to chr
            bpfile_prefix = vcf_to_load.replace(".vcf.gz","") + ".{}".format(i)

        
        #check if the file exists
        if convert=="bfile":
            is_file_exist = os.path.exists(bpfile_prefix+".bed")
            make_flag = "--make-bed"
        else:
            is_file_exist = os.path.exists(bpfile_prefix+".pgen")
            make_flag = "--make-pgen vzs"
        
        # figure memory
        if memory is not None:
            memory_flag = "--memory {}".format(memory)  
        else:
            memory_flag = ""
        
        #if not existing or overwrite is True
        if (not is_file_exist) or overwrite:
            script_vcf_to_bfile = """
            {} \
                --vcf {} \
                --chr {} \
                {} \
                --rm-dup force-first \
                --threads {}{}\
                --out {}
            """.format(plink2,
                        vcf_to_load, 
                       i, 
                       make_flag,
                       n_cores, memory_flag,
                       bpfile_prefix)
        
            # execute 
            try:
                if convert=="bfile":
                    log.write("  -Converting VCF to bfile: {}.bim/bed/fam...".format(bpfile_prefix))
                else:
                    log.write("  -Converting VCF to pfile: {}.pgen/pvar/psam...".format(bpfile_prefix))
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
            except subprocess.CalledProcessError as e:
                log.write(e.output)    
        else:
            log.write("  -Plink {} for CHR {} exists: {}. Skipping...".format(convert ,i, bpfile_prefix))
        
        if load_bim == True:
            if convert == "bfile":
                ref_bims = _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log)
            else:
                ref_bims = _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log)
    return ref_file_prefix_converted, plink_log, ref_bims

def _process_bgen(ref_file_prefix, 
                  chrlist, 
                  ref_bims, 
                  is_wild_card, 
                  log=Log(), 
                  plink_log="",
                  sample=None, 
                  bgen_mode="ref-first", 
                  n_cores=1, 
                  convert="bfile", 
                  memory=None, 
                  overwrite=False, 
                  load_bim=False,
                  plink="plink",
                 plink2="plink2"):
    log.write(" -Processing BGEN files : {}...".format(ref_file_prefix))
    
    #check plink version
    log = _checking_plink_version(log=log,plink2=plink2)
    
    # file path prefix to return
    if is_wild_card==True:
        ref_file_prefix_converted = ref_file_prefix.replace(".bgen","")
    else:
        ref_file_prefix_converted = ref_file_prefix.replace(".bgen","") + ".@"

    for i in chrlist:
    # for each chr    
        log.write("  -Processing BGEN for CHR {}...".format(i))
        
        # if multiple bfiles
        if is_wild_card==True:
            bgen_to_load = ref_file_prefix.replace("@",str(i))
            bpfile_prefix = bgen_to_load.replace(".bgen","")
        else:
            bgen_to_load = ref_file_prefix
            bpfile_prefix = bgen_to_load.replace(".bgen","") + ".{}".format(i)
        
        #check if the file exists
        if convert=="bfile":
            is_file_exist = os.path.exists(bpfile_prefix+".bed")
            make_flag = "--make-bed"
        else:
            is_file_exist = os.path.exists(bpfile_prefix+".pgen")
            make_flag = "--make-pgen vzs"
        
        #figure out sample file
        if sample is not None:
            if sample=="auto":
                sample_flag = "--sample {}.sample".format(bpfile_prefix)
            else:
                sample_flag = "--sample {}".format(sample)
        else:
            sample_flag=""
        
        # figure memory
        if memory is not None:
            memory_flag = "--memory {}".format(memory)  
        else:
            memory_flag = ""

        #if not existing or overwrite is True
        if (not is_file_exist) or overwrite:
            script_vcf_to_bfile = """
            {} \
                --bgen {} {} {}\
                --chr {} \
                {} \
                --rm-dup force-first \
                --threads {}{}\
                --out {}
            """.format(plink2,bgen_to_load, bgen_mode, sample_flag,
                       i, 
                       make_flag,
                       n_cores, memory_flag,
                       bpfile_prefix)
            # execute 
            try:
                if convert=="bfile":
                    log.write("  -Converting VCF to bfile: {}.bim/bed/fam...".format(bpfile_prefix))
                else:
                    log.write("  -Converting VCF to pfile: {}.pgen/pvar/psam...".format(bpfile_prefix))
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
            except subprocess.CalledProcessError as e:
                log.write(e.output)    
        else:
            log.write("  -PLINK {} for CHR {} exists. Skipping...".format(convert ,i))
        
        if load_bim == True:
            if convert == "bfile":
                ref_bims = _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log)
            else:
                ref_bims = _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log)
    return ref_file_prefix_converted, plink_log, ref_bims