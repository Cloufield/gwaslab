from typing import TYPE_CHECKING, Optional, List, Tuple, Union
import pandas as pd
import numpy as np
import subprocess
from gwaslab.info.g_Log import Log
import os
from gwaslab.extension import _checking_plink_version
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper

if TYPE_CHECKING:
    pass

# PLINK file suffix definitions
PLINK_SUFFIXES = ('.bim', '.map', '.pvar', '.pvar.zst')

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

def _process_vcf(
    ref_file_prefix: str,
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
    plink2: str = "plink2"
) -> Tuple[str, str, List[pd.DataFrame]]:
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

def _process_bgen(
    ref_file_prefix: str,
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
    plink2: str = "plink2"
) -> Tuple[str, str, List[pd.DataFrame]]:
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

def _load_bim_single(chrom: int, bfile_prefix: str, log: Log) -> pd.DataFrame:
    if "@" in bfile_prefix:
        bim_path = bfile_prefix.replace("@",str(chrom)) +".bim"
    else:
        bim_path = bfile_prefix+".bim"
    single_bim = pd.read_csv(bim_path,
                             sep=r"\s+",
                             usecols=[0, 1,3,4,5],
                             header=None,
                             dtype={1:"string",0:"category", 3:"int", 4:"string", 5:"string"}).rename(columns={1:"SNPID_bim",0:"CHR_bim",3:"POS_bim",4:"NEA_bim",5:"EA_bim"})
    single_bim = single_bim.loc[single_bim["CHR_bim"]==str(chrom),:]
    log.write("   -#variants in ref file: {}".format(len(single_bim))) 
    return single_bim

def _load_pvar_single(chrom: int, bpfile_prefix: str, log: Log) -> pd.DataFrame:
    if "@" in bpfile_prefix:
        bim_path = bpfile_prefix.replace("@",str(chrom)) +".pvar"
    else:
        bim_path = bpfile_prefix+".pvar"
    
    if os.path.exists(bpfile_prefix):
        bim_path =bpfile_prefix
    elif os.path.exists(bpfile_prefix+".zst"):
        bim_path =bpfile_prefix+".zst"
    
    single_bim = pd.read_csv(bim_path,
                             sep=r"\s+",
                             usecols=[0,1,2,3,4],
                             header=None,
                             comment="#",
                             dtype={2:"string",0:"category", 1:"int", 3:"string", 4:"string"}).rename(columns={2:"SNPID_bim",0:"CHR_bim",1:"POS_bim",3:"NEA_bim",4:"EA_bim"})
    single_bim = single_bim.loc[single_bim["CHR_bim"]==str(chrom),:]
    log.write("   -#variants in ref file: {}".format(len(single_bim))) 
    return single_bim

def _plink_chr_to_sumstats_chr(plink_chr: Union[str, int, float], 
                                mapper: ChromosomeMapper) -> Union[str, int, None]:
    """
    Convert PLINK chromosome (from BIM/PVAR file) to sumstats chromosome format using ChromosomeMapper.
    
    PLINK uses numeric codes for chromosomes:
    - 1-22: Autosomes
    - 23: X chromosome
    - 24: Y chromosome
    - 25: XY (pseudo-autosomal region)
    - 26: MT (mitochondrial)
    
    This notation is consistent across both PLINK 1.9 and PLINK 2.0. BIM files (PLINK 1.9)
    and PVAR files (PLINK 2.0) typically store chromosomes as strings like "1", "23", "24", "25", "26".
    
    By default, sumstats only have 1-25, and XY (PLINK 25) is treated as X (23).
    MT (PLINK 26) maps to 25 in sumstats.
    
    Note: The pseudo-autosomal region (PAR, PLINK 25) is mapped to chrX/chr23 in sumstats.
    
    Parameters
    ----------
    plink_chr : str, int, or float
        PLINK chromosome from BIM/PVAR file (typically string like "1", "23", "24", "25", "26")
    mapper : ChromosomeMapper
        ChromosomeMapper instance with detected sumstats format
    
    Returns
    -------
    str, int, or None
        Chromosome in sumstats format, or None if unconvertible
    
    References
    ----------
    PLINK 1.9 chromosome notation: https://www.cog-genomics.org/plink/1.9/filter#chr
    PLINK 2.0 chromosome notation: https://www.cog-genomics.org/plink/2.0/filter#chr
    """
    # Handle None, NaN, or empty values
    if pd.isna(plink_chr) or plink_chr == "":
        return None
    
    # Convert PLINK chromosome to numeric (PLINK notation)
    # PLINK files store chromosomes as numeric strings, so parse directly
    try:
        plink_chr_num = int(str(plink_chr).strip())
    except (ValueError, TypeError):
        return None
    
    # Map PLINK numeric codes to sumstats numeric codes
    # PLINK 25 (XY) -> 23 (X), PLINK 26 (MT) -> 25 (MT)
    if plink_chr_num == 25:
        # PAR (XY) maps to X (23) in sumstats
        sumstats_chr_num = 23
    elif plink_chr_num == 26:
        # MT (26) maps to 25 in sumstats
        sumstats_chr_num = 25
    elif 1 <= plink_chr_num <= 24:
        # Autosomes (1-22), X (23), Y (24) map directly
        sumstats_chr_num = plink_chr_num
    else:
        # Out of range
        return None
    
    # Use mapper to convert numeric code to sumstats format
    return mapper.number_to_sumstats(sumstats_chr_num)


def _match_sumstats_with_ref_bim(sumstats: pd.DataFrame, 
                                 ref_bim_all: pd.DataFrame, 
                                 has_allele_info: bool, 
                                 log: Log, 
                                 verbose: bool = True) -> int:
    """
    Match sumstats variants with reference BIM using CHR, POS, and optionally EA, NEA.
    
    Handles PLINK chromosome notation by converting PLINK chromosomes (from BIM files) 
    to match sumstats chromosome notation using ChromosomeMapper. PLINK uses numeric codes:
    - 1-22: Autosomes
    - 23: X chromosome
    - 24: Y chromosome
    - 25: XY (pseudo-autosomal region)
    - 26: MT (mitochondrial)
    
    By default, sumstats only have 1-25, and XY (PLINK 25) is treated as X (23).
    MT (PLINK 26) maps to 25 in sumstats.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe with CHR, POS, and optionally EA, NEA columns
    ref_bim_all : pd.DataFrame
        Reference BIM dataframe with CHR_bim, POS_bim, SNPID, EA_bim, NEA_bim columns.
        CHR_bim is typically stored as string/category in BIM files (e.g., "1", "23", "24", "25").
    has_allele_info : bool
        Whether to use EA/NEA for matching (True) or just CHR/POS (False)
    log : Log
        Logger instance
    verbose : bool
        Whether to log messages
        
    Returns
    -------
    int
        Number of matched variants
    """
    # Initialize SNPID_bim column to store the matched BIM ID (default to original SNPID)
    sumstats["SNPID_bim"] = sumstats["SNPID"].copy()
    
    if len(ref_bim_all) == 0:
        log.write(" -Warning: ref_bim is empty. Using original SNPIDs.",verbose=verbose)
        return 0
    
    # Initialize ChromosomeMapper and detect sumstats format
    mapper = ChromosomeMapper(species="homo sapiens", log=log, verbose=verbose)
    sumstats_format = mapper.detect_sumstats_format(sumstats["CHR"])
    log.write(" -Detected sumstats chromosome format: {}".format(sumstats_format),verbose=verbose)
    
    # Convert PLINK chromosomes (from BIM/PVAR) to sumstats format using mapper
    # PLINK 25 (XY) is treated as X (23) by default, PLINK 26 (MT) maps to 25
    ref_bim_all["_CHR_bim_sumstats"] = ref_bim_all["CHR_bim"].apply(
        lambda x: _plink_chr_to_sumstats_chr(x, mapper)
    )
    
    # Filter out unconvertible chromosomes
    ref_bim_valid = ref_bim_all["_CHR_bim_sumstats"].notna()
    n_invalid_bim = (~ref_bim_valid).sum()
    
    if n_invalid_bim > 0:
        log.write(" -Warning: {} BIM variants have unconvertible chromosome notation.".format(n_invalid_bim),verbose=verbose)
    
    # Filter to only valid chromosomes
    ref_bim_filtered = ref_bim_all[ref_bim_valid].copy()
    
    if len(ref_bim_filtered) == 0:
        log.write(" -Error: No valid chromosomes found in BIM after conversion.",verbose=verbose)
        return 0
    
    # Convert POS to int for consistent comparison (only if not already int)
    if not pd.api.types.is_integer_dtype(ref_bim_filtered["POS_bim"]):
        ref_bim_filtered["POS_bim"] = ref_bim_filtered["POS_bim"].astype(int)
    if not pd.api.types.is_integer_dtype(sumstats["POS"]):
        sumstats["POS"] = sumstats["POS"].astype(int)
    
    # Detect if type conversion is needed for CHR matching
    # Check if sumstats CHR and converted BIM CHR have compatible types
    sumstats_chr_sample = sumstats["CHR"].dropna().head(10)
    bim_chr_sample = ref_bim_filtered["_CHR_bim_sumstats"].dropna().head(10)
    
    if len(sumstats_chr_sample) > 0 and len(bim_chr_sample) > 0:
        sumstats_chr_type = type(sumstats_chr_sample.iloc[0])
        bim_chr_type = type(bim_chr_sample.iloc[0])
        needs_conversion = sumstats_chr_type != bim_chr_type
    else:
        needs_conversion = False
    
    # Convert BIM CHR to match sumstats CHR type if needed
    if needs_conversion:
        if pd.api.types.is_integer_dtype(sumstats["CHR"]):
            # Sumstats uses int, convert BIM CHR to int
            ref_bim_filtered["_CHR_bim_sumstats"] = pd.to_numeric(
                ref_bim_filtered["_CHR_bim_sumstats"], errors='coerce'
            )
        else:
            # Sumstats uses string, convert BIM CHR to string
            ref_bim_filtered["_CHR_bim_sumstats"] = ref_bim_filtered["_CHR_bim_sumstats"].astype(str)
    
    # Create set of (CHR, POS) tuples from sumstats for fast lookup
    # Sumstats CHR is already standardized, so use it directly
    sumstats_chr_pos = set(zip(sumstats["CHR"], sumstats["POS"]))
    
    # Create a boolean mask for filtering - more efficient than apply
    ref_bim_filtered["_chr_pos_tuple"] = list(zip(ref_bim_filtered["_CHR_bim_sumstats"], ref_bim_filtered["POS_bim"]))
    mask = ref_bim_filtered["_chr_pos_tuple"].isin(sumstats_chr_pos)
    ref_bim_matched = ref_bim_filtered[mask].drop(columns=["_chr_pos_tuple"]).copy()
    
    log.write(" -Filtered reference BIM from {} to {} variants matching sumstats CHR/POS...".format(
        len(ref_bim_all), len(ref_bim_matched)),verbose=verbose)
    
    if len(ref_bim_matched) == 0:
        log.write(" -No matching CHR/POS found in reference BIM.",verbose=verbose)
        return 0
    
    # Create mapping dictionary for efficient lookup
    # Map (CHR, POS) -> list of variant dicts with SNPID, EA, NEA
    pos_to_variants = {}
    for _, row in ref_bim_matched.iterrows():
        key = (row["_CHR_bim_sumstats"], int(row["POS_bim"]))
        if key not in pos_to_variants:
            pos_to_variants[key] = []
        pos_to_variants[key].append({
            "SNPID": row["SNPID"],
            "EA": str(row["EA_bim"]).upper() if pd.notna(row["EA_bim"]) else None,
            "NEA": str(row["NEA_bim"]).upper() if pd.notna(row["NEA_bim"]) else None
        })
    
    # Match variants using dictionary lookup
    n_matched = 0
    
    if has_allele_info:
        log.write(" -Matching sumstats with reference BIM using CHR, POS, EA, NEA...",verbose=verbose)
        
        for idx in sumstats.index:
            chr_val = sumstats.loc[idx, "CHR"]
            if pd.isna(chr_val):
                continue
            chr_pos = (chr_val, int(sumstats.loc[idx, "POS"]))
            
            if chr_pos not in pos_to_variants:
                continue
            
            # Get sumstats alleles
            ea_sum = str(sumstats.loc[idx, "EA"]).upper() if pd.notna(sumstats.loc[idx, "EA"]) else None
            nea_sum = str(sumstats.loc[idx, "NEA"]).upper() if pd.notna(sumstats.loc[idx, "NEA"]) else None
            
            if ea_sum == "nan" or nea_sum == "nan" or ea_sum is None or nea_sum is None:
                continue
            
            # Try to find match with same or swapped alleles
            for variant in pos_to_variants[chr_pos]:
                ea_bim = variant["EA"]
                nea_bim = variant["NEA"]
                
                if ea_bim is None or nea_bim is None:
                    continue
                
                # Check same allele order
                if ea_sum == ea_bim and nea_sum == nea_bim:
                    sumstats.loc[idx, "SNPID_bim"] = variant["SNPID"]
                    n_matched += 1
                    break
                # Check swapped allele order
                elif ea_sum == nea_bim and nea_sum == ea_bim:
                    sumstats.loc[idx, "SNPID_bim"] = variant["SNPID"]
                    n_matched += 1
                    break
        
        log.write(" -Matched {} variants using CHR, POS, EA, NEA...".format(n_matched),verbose=verbose)
    else:
        log.write(" -Warning: EA/NEA columns not found in sumstats. Using CHR, POS only for matching...",verbose=verbose)
        # CHR, POS only matching - use first match if multiple variants at same position
        for idx in sumstats.index:
            chr_val = sumstats.loc[idx, "CHR"]
            if pd.isna(chr_val):
                continue
            chr_pos = (chr_val, int(sumstats.loc[idx, "POS"]))
            
            if chr_pos in pos_to_variants:
                # Use first variant at this position
                sumstats.loc[idx, "SNPID_bim"] = pos_to_variants[chr_pos][0]["SNPID"]
                n_matched += 1
        
        log.write(" -Matched {} variants using CHR, POS only...".format(n_matched),verbose=verbose)
    
    return n_matched

