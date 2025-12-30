from typing import TYPE_CHECKING, Union, Optional, Tuple, List
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_NC
from gwaslab.bd.bd_common_data import _maketrans
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from pysam import VariantFile
import re
import numpy as np
import pandas as pd
import copy
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
from gwaslab.info.g_Log import Log
from shutil import which

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

# VCF/BCF file suffix definitions
VCF_BCF_SUFFIXES = ('.vcf.gz', '.bcf', '.vcf')

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
            log.write(f" -Prefix for chromosomes: {prefix}", verbose=verbose)
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
    contigs_list = list(vcf_bcf.header.contigs)
    
    # Use ChromosomeMapper to detect format
    if len(contigs_list) > 0:
        contigs_series = pd.Series(contigs_list[:min(10, len(contigs_list))])
        mapper = ChromosomeMapper(verbose=False)
        detected_format = mapper.detect_format(contigs_series)
        
        if detected_format == "chr":
            # Extract prefix from first matching contig
            for contig in contigs_list:
                match = re.search(r'(chr|Chr|CHR)([0-9xXyYmM]+)', contig, re.IGNORECASE)
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
    contigs_list = list(vcf_bcf.header.contigs)
    
    # Use ChromosomeMapper to detect format
    if len(contigs_list) > 0:
        contigs_series = pd.Series(contigs_list[:min(10, len(contigs_list))])
        mapper = ChromosomeMapper(verbose=False)
        detected_format = mapper.detect_format(contigs_series)
        
        if detected_format == "nc":
            # Check which build matches
            for build in ["19", "38"]:
                nc_mapping = get_number_to_NC(build=build)
                if nc_mapping:
                    # Check if any contig matches this build's NC IDs
                    nc_values = set(nc_mapping.values())
                    if any(contig in nc_values for contig in contigs_list):
                        log.write(f"  -RefSeq ID detected (hg{build}) in VCF/BCF...", verbose=verbose)
                        return nc_mapping
    return None


def _get_ld_matrix_from_vcf(sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], 
                vcf_path: Optional[str] = None, 
                region: Optional[Tuple[Union[int, str], int, int]] = None,
                log: Log = Log(), 
                verbose: bool = True, 
                pos: str = "POS",
                nea: str = "NEA",
                ea: str = "EA", 
                mapper: Optional[ChromosomeMapper] = None,
                tabix: Optional[bool] = None,
                export_path: Optional[str] = None) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Calculate full LD matrix from VCF file and return both the LD matrix and corresponding sumstats.
    
    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    vcf_path : str
        Path to the VCF file
    region : tuple
        Region specification (chromosome, start, end)
    log : Log
        Logging object
    verbose : bool
        Verbose flag
    pos : str
        Position column name
    nea : str
        Non-effect allele column name
    ea : str
        Effect allele column name
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance to use for chromosome name conversion.
        If not provided, creates a default mapper with automatic format detection.
    tabix : bool
        Whether to use tabix indexing
    export_path : str, optional
        Path to export the results. If provided, sumstats and LD matrix will be saved.
        Files will be named based on the region (e.g., chr1_1000_2000_sumstats.tsv.gz).
    
    Returns
    -------
    tuple
        (matched_sumstats, ld_matrix, valid_indices) where:
        - matched_sumstats: Subset of sumstats with valid variants, ordered to match ld_matrix
        - ld_matrix: Full LD matrix as numpy array (r^2 values)
        - valid_indices: Indices of variants used in LD calculation
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data.copy()
    
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)
    log.write(" -region : chr{}:{}-{}".format(region[0], region[1], region[2]), verbose=verbose)
    from gwaslab.util.util_in_filter_value import _filter_region
    sumstats = _filter_region(sumstats.copy(), region, log=log, verbose=verbose)

    if tabix is None:
        tabix = which("tabix")
        log.write(" -tabix will be used: {}".format(tabix),verbose=verbose)
    
    # Get or create mapper
    if mapper is None:
        mapper = ChromosomeMapper(log=log, verbose=verbose)
        # Auto-detect reference format from VCF file
        mapper.detect_reference_format(vcf_path)
    
    # Convert region chromosome to reference format
    # as_string=True for region string formatting
    region_chr_ref = mapper.sumstats_to_reference(region[0], reference_file=vcf_path, as_string=True)
    # load genotype data of the targeted region
    region_str = f"{region_chr_ref}:{region[1]}-{region[2]}"
    log.write(" -loading VCF region: {}".format(region_str), verbose=verbose)
    ref_genotype = read_vcf(vcf_path,region=region_str,tabix=tabix)
    if ref_genotype is None:
        log.warning("No data was retrieved. Skipping ...")
        ref_genotype=dict()
        ref_genotype["variants/POS"]=np.array([],dtype="int64")
    log.write(" -Retrieving index...", verbose=verbose)
    log.write(" -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
    # match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
     #######################################################################################
    def match_variant(x):
        # x: "POS,NEA,EA"
        if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
            # position match
            matches = np.where(ref_genotype["variants/POS"] == x.iloc[0])[0]
            if len(matches) > 1:  
                # multiple position matches
                for j in matches:
                    # for each possible match, compare ref and alt
                    if x.iloc[1] == ref_genotype["variants/REF"][j]:
                        if x.iloc[2] in ref_genotype["variants/ALT"][j]:
                            return j
                    elif x.iloc[1] in ref_genotype["variants/ALT"][j]:
                        if x.iloc[2] == ref_genotype["variants/REF"][j]:
                            return j    
                return None
            else: 
                # single match
                #log.write("  -Single match found for position {}".format(x.iloc[0]), verbose=verbose)
                return matches[0]
        else:
            # no position match
            return None
    log.write(" -Matching variants using POS, NEA, EA ...", verbose=verbose)
    log.write(" -Total variants in sumstats: {}".format(len(sumstats)), verbose=verbose)
    #############################################################################################
    sumstats["REFINDEX"] = sumstats[[pos,nea,ea]].apply(lambda x: match_variant(x),axis=1)
    #############################################################################################
    
    # Report matching results
    matched_count = sumstats["REFINDEX"].notna().sum()
    log.write(" -Matched variants: {}".format(matched_count), verbose=verbose)
    log.write(" -Unmatched variants: {}".format(len(sumstats) - matched_count), verbose=verbose)

    # Get non-na SNP indices for LD calculation
    valid_indices = sumstats["REFINDEX"].dropna().astype("int").values
    log.write(" -Valid indices for LD calculation: {}".format(len(valid_indices)), verbose=verbose)
    
    # Calculate full LD matrix using allel's standard method
    if len(valid_indices) > 0:
        log.write(" -Calculating LD matrix for {} variants...".format(len(valid_indices)), verbose=verbose)
        # Get genotypes for all valid SNPs
        all_snp_genotypes = GenotypeArray(ref_genotype["calldata/GT"][valid_indices]).to_n_alt()
        log.write(" -Genotype array shape: {}".format(all_snp_genotypes.shape), verbose=verbose)
        
        # Calculate pairwise LD matrix using allel's rogers_huff_r function
        log.write(" -Computing pairwise LD values...", verbose=verbose)
        r_values = rogers_huff_r_between(all_snp_genotypes, all_snp_genotypes)
        ld_matrix = np.power(r_values, 2)  # Convert r to r^2
        log.write(" -LD matrix calculated with shape: {}".format(ld_matrix.shape), verbose=verbose)
        
        # Create a subset of sumstats that matches the order of the LD matrix
        # First get the indices of non-na entries
        non_na_indices = sumstats.index[sumstats["REFINDEX"].notna()]
        matched_sumstats = sumstats.loc[non_na_indices].copy()
        log.write(" -Matched sumstats before cleaning: {}".format(len(matched_sumstats)), verbose=verbose)
        
        # Remove rows with missing pos, nea, or ea values while preserving order
        # We need to do this carefully to maintain the correspondence with ld_matrix
        na_mask = matched_sumstats[[pos, nea, ea]].isna().any(axis=1)
        clean_indices = matched_sumstats.index[~na_mask]
        matched_sumstats = matched_sumstats.loc[clean_indices]
        log.write(" -Matched sumstats after cleaning: {}".format(len(matched_sumstats)), verbose=verbose)
        
        # Also need to filter ld_matrix and valid_indices to match
        clean_original_indices = np.where(~na_mask)[0]
        log.write(" -Cleaning LD matrix and indices...", verbose=verbose)
        if len(clean_original_indices) > 0:
            ld_matrix = ld_matrix[np.ix_(clean_original_indices, clean_original_indices)]
            # Update valid_indices to match the cleaned data
            valid_indices = valid_indices[clean_original_indices]
            log.write(" -Final LD matrix shape: {}".format(ld_matrix.shape), verbose=verbose)
        else:
            ld_matrix = np.array([]).reshape(0, 0)
            log.write(" -No valid variants after cleaning, LD matrix is empty", verbose=verbose)
    else:
        log.write(" -No valid indices found for LD calculation, returning empty results", verbose=verbose)
        ld_matrix = np.array([]).reshape(0, 0)
        matched_sumstats = sumstats.iloc[[]].copy()  # Empty DataFrame with same structure
    
    ####################################################################################################
    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    
    # Export results if export_path is provided
    if export_path is not None:
        import os
        # Create export directory if it doesn't exist
        os.makedirs(export_path, exist_ok=True)
        
        # Automatically construct prefix based on region
        if region is not None:
            prefix = "chr{}_{}_{}".format(region[0], region[1], region[2])
        else:
            prefix = "ld"
        
        # Export matched sumstats
        sumstats_export_path = os.path.join(export_path, "{}_sumstats.tsv.gz".format(prefix))
        log.write(" -Exporting matched sumstats to: {}".format(sumstats_export_path), verbose=verbose)
        matched_sumstats.to_csv(sumstats_export_path, sep="\t", index=False, compression="gzip")
        
        # Export LD matrix
        ld_matrix_export_path = os.path.join(export_path, "{}_ldmatrix.npy".format(prefix))
        log.write(" -Exporting LD matrix to: {}".format(ld_matrix_export_path), verbose=verbose)
        np.save(ld_matrix_export_path, ld_matrix)
        
        log.write(" -Export completed successfully!", verbose=verbose)
    
    # Return both the LD matrix and the matched sumstats
    return matched_sumstats, ld_matrix