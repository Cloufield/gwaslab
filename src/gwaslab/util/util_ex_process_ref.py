from typing import TYPE_CHECKING, Optional, List, Tuple
import pandas as pd
import numpy as np
import subprocess
from gwaslab.info.g_Log import Log
import os
from gwaslab.extension import _checking_plink_version

if TYPE_CHECKING:
    pass

def _process_plink_input_files(chrlist: List[int], 
                               bfile: Optional[str] = None, 
                               pfile: Optional[str] = None, 
                               vcf: Optional[str] = None, 
                               bgen: Optional[str] = None,  
                               sample: Optional[str] = None, 
                               threads: int = 1, 
                               plink_log: str = "", 
                               log: Log = Log(), 
                               overwrite: bool = False, 
                               bgen_mode: str = "ref-first", 
                               convert: str = "bfile", 
                               memory: Optional[int] = None, 
                               load_bim: bool = False,
                               plink: str = "plink",
                               plink2: str = "plink2") -> Tuple[str, str, List[pd.DataFrame], str]:
    """
    Process input files (bfile,pfile,vcf,bgen) to either PLINK1 bed/bim/fam or PLINK2 pgen/psam/pvar. 
    
    Parameters:
    -----------
    load_bim : bool, default=False
        If True, load BIM/PVAR variant information into ref_bims list.
        If False (default), ref_bims will be an empty list.
        **Important**: Pass load_bim=True to populate ref_bims when using VCF or BGEN inputs.
    
    Returns:
    --------
    ref_file_prefix : prefix for either bfile or pfile.
    plink_log : if plink was used, return the log. Otherwise, return an empty string.
    ref_bims : if load_bim is True, return bim files as a list of pd.DataFrame. Otherwise, empty list.
    filetype : either bfile or pfile.
    """

    # Step 1: Initialize list to store BIM/PVAR dataframes if load_bim is True
    ref_bims = []
    
    # Step 2: Determine the input file type and check if it uses wildcard pattern (@)
    # Returns: filetype (bfile/pfile/vcf/bgen), ref_file_prefix, and is_wild_card flag
    # Note: bfile and pfile require no conversion, while vcf and bgen need to be converted
    # File prefix can be single file or wildcard pattern with @ for chromosome-specific files
    filetype, ref_file_prefix, is_wild_card = _check_file_type(bfile,  pfile,  vcf,  bgen, sample)

    # Step 3: Process bfile format (PLINK1 bed/bim/fam)
    # No conversion needed, just validate files exist and optionally load BIM data
    if filetype == "bfile": 
        ref_file_prefix, ref_bims = _process_bfile(chrlist=chrlist,
                                                   ref_file_prefix=ref_file_prefix, 
                                                   ref_bims=ref_bims, 
                                                   is_wild_card=is_wild_card, 
                                                   log=log, 
                                                   load_bim=load_bim)

    # Step 4: Process pfile format (PLINK2 pgen/pvar/psam)
    # No conversion needed, just validate files exist and optionally load PVAR data
    elif filetype == "pfile": 
        ref_file_prefix, ref_bims = _process_pfile(chrlist=chrlist,
                                                   ref_file_prefix=ref_file_prefix, 
                                                   ref_bims=ref_bims, 
                                                   is_wild_card=is_wild_card, 
                                                   log=log, 
                                                   load_bim=load_bim)

    # Step 5: Process VCF format - convert to bfile or pfile using PLINK2
    # VCF files need conversion, so filetype is updated to the target format (convert parameter)
    elif filetype == "vcf": 
        ref_file_prefix, plink_log, ref_bims = _process_vcf(ref_file_prefix=ref_file_prefix, 
                                                            chrlist=chrlist, 
                                                            ref_bims=ref_bims, 
                                                            is_wild_card=is_wild_card, 
                                                            log=log, 
                                                            plink_log=plink_log, 
                                                            threads=threads, 
                                                            convert=convert, 
                                                            memory=memory, 
                                                            overwrite=overwrite, 
                                                            load_bim=load_bim,
                                                            plink=plink,
                                                            plink2=plink2)
        # Update filetype to the converted format (bfile or pfile)
        filetype = convert
    
    # Step 6: Process BGEN format - convert to bfile or pfile using PLINK2
    # BGEN files need conversion, so filetype is updated to the target format (convert parameter)
    elif filetype == "bgen": 
        ref_file_prefix, plink_log, ref_bims = _process_bgen(ref_file_prefix=ref_file_prefix, 
                                                            chrlist=chrlist, 
                                                            bgen_mode=bgen_mode,
                                                            sample=sample,
                                                            ref_bims=ref_bims, 
                                                            is_wild_card=is_wild_card, 
                                                            log=log, 
                                                            plink_log=plink_log, 
                                                            threads=threads, 
                                                            convert=convert, 
                                                            memory=memory, 
                                                            overwrite=overwrite, 
                                                            load_bim=load_bim,
                                                            plink=plink,
                                                            plink2=plink2)
        # Update filetype to the converted format (bfile or pfile)
        filetype = convert
    
    # Step 7: Return processed file prefix, PLINK log (if conversion occurred), 
    #         list of BIM/PVAR dataframes (if load_bim=True), and final filetype
    return ref_file_prefix, plink_log, ref_bims, filetype

def _load_single_bim_to_ref_bims(bpfile_prefix: str, ref_bims: List[pd.DataFrame], log: Log) -> List[pd.DataFrame]:
    bim_path =bpfile_prefix+".bim"
    single_bim = pd.read_csv(bim_path,
                             sep=r"\s+",
                             usecols=[0,1,3,4,5],
                             header=None,
                             dtype={1:"string",0:"category", 3:"int", 4:"string", 5:"string"}).rename(columns={1:"SNPID",0:"CHR_bim",3:"POS_bim",4:"EA_bim",5:"NEA_bim"})
    log.write("   -Variants in ref file: {}".format(len(single_bim))) 
    ref_bims.append(single_bim)
    return ref_bims

def _load_single_pvar_to_ref_bims(bpfile_prefix: str, ref_bims: List[pd.DataFrame], log: Log) -> List[pd.DataFrame]:
    if os.path.exists(bpfile_prefix+".pvar"):
        bim_path =bpfile_prefix+".pvar"
    elif os.path.exists(bpfile_prefix+".pvar.zst"):
        bim_path =bpfile_prefix+".pvar.zst"
    single_bim = pd.read_csv(bim_path,
                             sep=r"\s+",
                             usecols=[0,1,2,3,4],
                             header=None,
                             comment="#",
                             dtype={2:"string",0:"category", 1:"int", 3:"string", 4:"string"}).rename(columns={2:"SNPID",0:"CHR_bim",1:"POS_bim",3:"EA_bim",4:"NEA_bim"})
    log.write("   -Variants in ref file: {}".format(len(single_bim))) 
    ref_bims.append(single_bim)
    return ref_bims

def _check_file_type(bfile: Optional[str] = None,  
                     pfile: Optional[str] = None,  
                     vcf: Optional[str] = None,  
                     bgen: Optional[str] = None, 
                     sample: Optional[str] = None) -> Tuple[str, str, bool]:
    
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
    
def _process_bfile(chrlist: List[int], ref_file_prefix: str, ref_bims: List[pd.DataFrame], is_wild_card: bool, log: Log, load_bim: bool = False) -> Tuple[str, List[pd.DataFrame]]:
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

def _process_pfile(chrlist: List[int], ref_file_prefix: str, ref_bims: List[pd.DataFrame], is_wild_card: bool, log: Log, load_bim: bool = False) -> Tuple[str, List[pd.DataFrame]]:
    if is_wild_card==False:
        is_bim_exist = os.path.exists(ref_file_prefix+".pvar") or os.path.exists(ref_file_prefix+".pvar.zst")
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
            is_bim_exist = os.path.exists(single_chr_ref_file_prefix+".pvar") or os.path.exists(single_chr_ref_file_prefix+".pvar.zst")
            is_bed_exist = os.path.exists(single_chr_ref_file_prefix+".pgen")
            is_fam_exist = os.path.exists(single_chr_ref_file_prefix+".psam")
            if not (is_bim_exist and is_bed_exist  and  is_fam_exist):
                raise ValueError("PLINK pfiles for CHR {} are missing...".format(i))
            if load_bim==True:  
                ref_bims = _load_single_pvar_to_ref_bims(ref_file_prefix, ref_bims, log)  
        log.write(" -Split PLINK pfiles for each CHR as LD reference panel: {}".format(ref_file_prefix))  
  
    return ref_file_prefix, ref_bims

def _process_vcf(ref_file_prefix: str, 
                 chrlist: List[int], 
                 ref_bims: List[pd.DataFrame], 
                 is_wild_card: bool, 
                 log: Log, 
                 plink_log: str, 
                 threads: int = 1, 
                 convert: str = "bfile", 
                 memory: Optional[int] = None, 
                 overwrite: bool = False, 
                 load_bim: bool = False,
                 plink: str = "plink",
                 plink2: str = "plink2") -> Tuple[str, str, List[pd.DataFrame]]:
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
        conversion_success = True
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
                       threads, memory_flag,
                       bpfile_prefix)
        
            # execute conversion
            try:
                if convert=="bfile":
                    log.write("  -Converting VCF to bfile: {}.bim/bed/fam...".format(bpfile_prefix))
                else:
                    log.write("  -Converting VCF to pfile: {}.pgen/pvar/psam...".format(bpfile_prefix))
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
                # Verify conversion succeeded by checking if output file exists
                if convert=="bfile":
                    conversion_success = os.path.exists(bpfile_prefix+".bed")
                else:
                    conversion_success = os.path.exists(bpfile_prefix+".pgen")
                # Read PLINK log file if it exists
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file, 'r') as f:
                        plink_log += f.read() + "\n"
            except subprocess.CalledProcessError as e:
                log.write(e.output)
                conversion_success = False
                # Try to read PLINK log file even on error
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file, 'r') as f:
                        plink_log += f.read() + "\n"
        else:
            log.write("  -Plink {} for CHR {} exists: {}. Skipping...".format(convert ,i, bpfile_prefix))
        
        # Load BIM/PVAR data if requested and conversion was successful
        # Note: load_bim must be True to populate ref_bims, otherwise it remains empty
        if load_bim == True:
            # Only load if file exists (either from successful conversion or pre-existing)
            if convert == "bfile":
                if os.path.exists(bpfile_prefix+".bim"):
                    ref_bims = _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log)
                else:
                    log.write("  -Warning: BIM file not found for CHR {}: {}.bim. Skipping load.".format(i, bpfile_prefix))
            else:
                if os.path.exists(bpfile_prefix+".pvar") or os.path.exists(bpfile_prefix+".pvar.zst"):
                    ref_bims = _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log)
                else:
                    log.write("  -Warning: PVAR file not found for CHR {}: {}.pvar. Skipping load.".format(i, bpfile_prefix))
    return ref_file_prefix_converted, plink_log, ref_bims

def _process_bgen(ref_file_prefix: str, 
                  chrlist: List[int], 
                  ref_bims: List[pd.DataFrame], 
                  is_wild_card: bool, 
                  log: Log = Log(), 
                  plink_log: str = "",
                  sample: Optional[str] = None, 
                  bgen_mode: str = "ref-first", 
                  threads: int = 1, 
                  convert: str = "bfile", 
                  memory: Optional[int] = None, 
                  overwrite: bool = False, 
                  load_bim: bool = False,
                  plink: str = "plink",
                  plink2: str = "plink2") -> Tuple[str, str, List[pd.DataFrame]]:
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
        conversion_success = True
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
                       threads, memory_flag,
                       bpfile_prefix)
            # execute conversion
            try:
                if convert=="bfile":
                    log.write("  -Converting BGEN to bfile: {}.bim/bed/fam...".format(bpfile_prefix))
                else:
                    log.write("  -Converting BGEN to pfile: {}.pgen/pvar/psam...".format(bpfile_prefix))
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
                # Verify conversion succeeded by checking if output file exists
                if convert=="bfile":
                    conversion_success = os.path.exists(bpfile_prefix+".bed")
                else:
                    conversion_success = os.path.exists(bpfile_prefix+".pgen")
                # Read PLINK log file if it exists
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file, 'r') as f:
                        plink_log += f.read() + "\n"
            except subprocess.CalledProcessError as e:
                log.write(e.output)
                conversion_success = False
                # Try to read PLINK log file even on error
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file, 'r') as f:
                        plink_log += f.read() + "\n"
        else:
            log.write("  -PLINK {} for CHR {} exists. Skipping...".format(convert ,i))
        
        # Load BIM/PVAR data if requested and conversion was successful
        # Note: load_bim must be True to populate ref_bims, otherwise it remains empty
        if load_bim == True:
            # Only load if file exists (either from successful conversion or pre-existing)
            if convert == "bfile":
                if os.path.exists(bpfile_prefix+".bim"):
                    ref_bims = _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log)
                else:
                    log.write("  -Warning: BIM file not found for CHR {}: {}.bim. Skipping load.".format(i, bpfile_prefix))
            else:
                if os.path.exists(bpfile_prefix+".pvar") or os.path.exists(bpfile_prefix+".pvar.zst"):
                    ref_bims = _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log)
                else:
                    log.write("  -Warning: PVAR file not found for CHR {}: {}.pvar. Skipping load.".format(i, bpfile_prefix))
    return ref_file_prefix_converted, plink_log, ref_bims
