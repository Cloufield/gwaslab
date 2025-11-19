from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_NC
from gwaslab.bd.bd_common_data import _maketrans
from pysam import VariantFile
from Bio import SeqIO
import re

def auto_check_vcf_chr_dict(
    vcf_path: str | None, 
    vcf_chr_dict: dict | None, 
    verbose: bool, 
    log: "Log"
) -> dict:
    """
    Automatically determine chromosome naming convention used in VCF/BCF files.

    This function checks the chromosome naming convention used in VCF/BCF files and
    returns an appropriate chromosome dictionary mapping. It first checks if the
    chromosome IDs match RefSeq IDs (for hg19 or hg38 builds), then checks for
    common prefixes like "chr", and finally defaults to standard numeric chromosomes.

    Parameters
    ----------
    vcf_path : str or None
        Path to the VCF/BCF file to check. If None, the function returns vcf_chr_dict.
    vcf_chr_dict : dict or None
        Optional pre-defined chromosome dictionary. If None, the function will
        attempt to determine the appropriate dictionary.
    verbose : bool
        If True, print detailed progress messages.
    log : gwaslab.g_Log.Log
        Logging object for recording process information.

    Returns
    -------
    dict
        A chromosome dictionary mapping that matches the chromosome naming
        convention used in the VCF/BCF file. This dictionary maps standard
        chromosome numbers to the format used in the VCF/BCF file.

    Notes
    -----
    The function checks for chromosome naming conventions in the following order:
    1. RefSeq IDs (for hg19 or hg38 builds)
    2. Chromosome prefixes (e.g., "chr1", "Chr1", "CHR1")
    3. Standard numeric chromosomes (e.g., 1, 2, 3, ...)

    If no specific convention is detected, it defaults to standard numeric
    chromosomes without any prefix.
    """
    if vcf_path is not None and vcf_chr_dict is None:
        log.write(" -Checking chromosome notations in VCF/BCF files...", verbose=verbose)
        vcf_chr_dict = check_vcf_chr_NC(vcf_path, log, verbose)
        if vcf_chr_dict is not None:
            return vcf_chr_dict
        
        log.write(" -Checking prefix for chromosomes in VCF/BCF files...", verbose=verbose)
        prefix = check_vcf_chr_prefix(vcf_path, log, verbose)
        if prefix is not None:
            log.write(f" -Prefix for chromosomes: {prefix}")
            vcf_chr_dict = get_number_to_chr(prefix=prefix)
        else:
            log.write(" -No prefix for chromosomes in the VCF/BCF files.", verbose=verbose)
            vcf_chr_dict = get_number_to_chr()
    
    # Filter to only include contigs present in the VCF file
    if vcf_path is not None:
        vcf_bcf = VariantFile(vcf_path)
        valid_contigs = set(vcf_bcf.header.contigs)
        vcf_chr_dict_filtered = {k: v for k, v in vcf_chr_dict.items() if v in valid_contigs}
        if vcf_chr_dict_filtered:
            return vcf_chr_dict_filtered
    
    return vcf_chr_dict or get_number_to_chr()

def check_vcf_chr_prefix(
    vcf_bcf_path: str, 
    log: "Log", 
    verbose: bool
) -> str | None:
    """
    Check for chromosome prefix in VCF/BCF file headers.

    Parameters
    ----------
    vcf_bcf_path : str
        Path to VCF/BCF file.
    log : gwaslab.g_Log.Log
        Logging object.
    verbose : bool
        Whether to log detailed messages.

    Returns
    -------
    str or None
        Detected chromosome prefix (e.g., "chr", "Chr", "CHR") if found, otherwise None.
    """
    vcf_bcf = VariantFile(vcf_bcf_path)
    for contig in vcf_bcf.header.contigs:
        match = re.search(r'(chr|Chr|CHR)([0-9xXyYmM]+)', contig)
        if match:
            return match.group(1)
    return None

def is_vcf_file(path: str) -> bool:
    """Check if the given path is a VCF or BCF file by examining headers."""
    try:
        with VariantFile(path) as f:
            return True
    except Exception:
        return False

def check_vcf_chr_NC(
    vcf_bcf_path: str, 
    log: "Log", 
    verbose: bool
) -> dict | None:
    """
    Check for RefSeq chromosome IDs in VCF/BCF file headers.

    Parameters
    ----------
    vcf_bcf_path : str
        Path to VCF/BCF file.
    log : gwaslab.g_Log.Log
        Logging object.
    verbose : bool
        Whether to log detailed messages.

    Returns
    -------
    dict or None
        Chromosome mapping dictionary for detected build (hg19/hg38) if found, otherwise None.
    """
    vcf_bcf = VariantFile(vcf_bcf_path)
    for contig in vcf_bcf.header.contigs:
        if contig in get_number_to_NC(build="19").values():
            log.write("  -RefSeq ID detected (hg19) in VCF/BCF...", verbose=verbose)
            return get_number_to_NC(build="19")
        if contig in get_number_to_NC(build="38").values():
            log.write("  -RefSeq ID detected (hg38) in VCF/BCF...", verbose=verbose)
            return get_number_to_NC(build="38")
    return None
