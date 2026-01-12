import scipy.sparse as sparse
import numpy as np
import pandas as pd

import subprocess
import os
import re
import gc
from typing import TYPE_CHECKING, Optional, List, Dict, Any, Tuple, Union

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_plink_version

from gwaslab.hm.hm_casting import _merge_mold_with_sumstats_by_chrpos

from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.io.io_plink import _process_plink_input_files
from gwaslab.util.util_in_filter_value import _exclude_hla
from gwaslab.util.util_ex_calculate_ldmatrix import _extract_variants_in_locus, _align_sumstats_with_bim, _calculate_ld_r
from gwaslab.qc.qc_decorator import with_logging
from shutil import which
from gwaslab.io.io_load_ld import _load_ld_matrix, _load_ld_map

####################################################################################################
# Helper functions for unified LD source handling
####################################################################################################

def _normalize_parameters_to_lists(
    n_studies: int,
    bfile: Optional[Union[str, List[str]]],
    pfile: Optional[Union[str, List[str]]],
    vcf: Optional[Union[str, List[str]]],
    ld_paths: Optional[List[str]],
    ld_maps: Optional[List[str]],
    ld_types: Optional[List[str]],
    ld_map_dics: Optional[List[Dict[str, str]]],
    tabix: Optional[Union[bool, str, List[Union[bool, str]]]]
) -> Tuple[List[Optional[str]], List[Optional[str]], List[Optional[str]], List[Optional[str]], List[Optional[str]], List[str], List[Optional[Dict[str, str]]], List[Optional[Union[bool, str]]]]:
    """
    Normalize LD source parameters to lists (one per study).
    
    Converts single values to lists, or keeps lists as-is.
    Ensures all parameters have the same length as n_studies.
    
    Returns:
        Tuple of normalized parameter lists
    """
    def _normalize_to_list(value, n_studies, default=None):
        if value is None:
            return [default] * n_studies
        elif isinstance(value, str):
            return [value] * n_studies
        elif isinstance(value, list):
            if len(value) != n_studies:
                raise ValueError("List parameter has {} elements but {} studies expected".format(len(value), n_studies))
            return value
        else:
            return [value] * n_studies
    
    bfiles = _normalize_to_list(bfile, n_studies)
    pfiles = _normalize_to_list(pfile, n_studies)
    vcfs = _normalize_to_list(vcf, n_studies)
    ld_paths_list = _normalize_to_list(ld_paths, n_studies)
    ld_maps_list = _normalize_to_list(ld_maps, n_studies)
    ld_types_list = _normalize_to_list(ld_types, n_studies, default="npz")
    ld_map_dics_list = _normalize_to_list(ld_map_dics, n_studies)
    tabix_list = _normalize_to_list(tabix, n_studies)
    
    return bfiles, pfiles, vcfs, ld_paths_list, ld_maps_list, ld_types_list, ld_map_dics_list, tabix_list

def _detect_ld_source_types(
    n_studies: int,
    ld_paths_list: List[Optional[str]],
    ld_maps_list: List[Optional[str]],
    bfiles: List[Optional[str]],
    pfiles: List[Optional[str]],
    vcfs: List[Optional[str]],
    log: Log,
    verbose: bool
) -> List[str]:
    """
    Detect LD source type for each study.
    
    Priority: ld_paths+ld_maps > bfile > pfile > vcf
    
    Returns:
        List of source types, one per study
    """
    ld_source_types = []
    for i in range(n_studies):
        source_type = _detect_ld_source_type(
            ld_path=ld_paths_list[i],
            ld_map_path=ld_maps_list[i],
            bfile=bfiles[i],
            pfile=pfiles[i],
            vcf=vcfs[i]
        )
        ld_source_types.append(source_type)
        log.write(" -Study {} will use LD source type: {}".format(i+1, source_type), verbose=verbose)
    return ld_source_types

def _filter_and_prepare_sumstats(
    sumstats: pd.DataFrame,
    required_cols: List[str],
    log: Log,
    verbose: bool
) -> pd.DataFrame:
    """
    Filter NA values and remove duplicates from sumstats.
    
    Returns:
        Cleaned sumstats DataFrame
    """
    log.write(" -Filtering NA values from sumstats...", verbose=verbose)
    missing_cols = [col for col in required_cols if col not in sumstats.columns]
    if missing_cols:
        raise ValueError("Missing required columns in sumstats: {}".format(missing_cols))
    
    initial_count = len(sumstats)
    sumstats = sumstats.dropna(subset=required_cols).copy()
    filtered_count = len(sumstats)
    n_filtered = initial_count - filtered_count
    if n_filtered > 0:
        log.write(" -Filtered {} variants with NA values ({} -> {})".format(
            n_filtered, initial_count, filtered_count), verbose=verbose)
    
    if len(sumstats) == 0:
        raise ValueError("No variants remaining after filtering NA values")
    
    # Drop duplicate SNPIDs (keep first occurrence)
    log.write(" -Dropping duplicated SNPIDs...", verbose=verbose)
    sumstats = sumstats.drop_duplicates(subset=["SNPID"]).copy()
    
    return sumstats

def _extract_lead_variants(
    sumstats: pd.DataFrame,
    loci: Optional[List[str]],
    loci_chrpos: Optional[List[Tuple[int, int]]],
    locus: Optional[Tuple[int, int]],
    suffixes: List[str],
    getlead_kwargs: Dict[str, Any],
    exclude_hla: bool,
    log: Log,
    verbose: bool
) -> pd.DataFrame:
    """
    Extract lead variants (significant loci) from sumstats.
    
    Returns:
        DataFrame with lead variants
    """
    if loci is None and loci_chrpos is None and locus is None:
        log.write(" -Loci were not provided. All significant loci will be automatically extracted...", verbose=verbose)
        sig_df = _get_sig(sumstats, variant_id="SNPID", chrom="CHR", pos="POS", p="P"+suffixes[0], **getlead_kwargs)
    elif loci is not None:
        sig_df = sumstats.loc[sumstats["SNPID"].isin(loci),:]
    elif loci_chrpos is not None:
        # Extract by CHR:POS tuples
        sig_df = pd.DataFrame()
        for chr_pos in loci_chrpos:
            matches = sumstats[(sumstats["CHR"] == chr_pos[0]) & (sumstats["POS"] == chr_pos[1])]
            sig_df = pd.concat([sig_df, matches], ignore_index=True)
    else:  # locus is not None
        matches = sumstats[(sumstats["CHR"] == locus[0]) & (sumstats["POS"] == locus[1])]
        sig_df = matches
    
    if exclude_hla:
        sig_df = _exclude_hla(sig_df, log=log, verbose=verbose)
    
    sig_df = sig_df.reset_index(drop=True)
    return sig_df

def _get_all_ld_maps(
    n_studies: int,
    studies: List[str],
    ld_source_types: List[str],
    row: pd.Series,
    locus_sumstats: pd.DataFrame,
    ld_maps_list: List[Optional[str]],
    ld_map_dics_list: List[Optional[Dict[str, str]]],
    ld_map_kwargs: Dict[str, Any],
    bfiles: List[Optional[str]],
    pfiles: List[Optional[str]],
    vcfs: List[Optional[str]],
    windowsizekb: int,
    plink: str,
    plink2: str,
    threads: int,
    memory: Optional[int],
    overwrite: bool,
    tabix_list: List[Optional[Union[bool, str]]],
    log: Log,
    verbose: bool,
    **kwargs: Any
) -> List[pd.DataFrame]:
    """
    Get LD maps (variant information) for all studies.
    
    This is STEP 2: Extract variant information from each LD source
    without calculating LD matrices yet.
    
    Returns:
        List of LD maps, one per study
    """
    all_ld_maps = []
    log.write(" -Getting LD maps for all studies to determine common SNPs...", verbose=verbose)
    
    for i in range(n_studies):
        gc.collect()
        study_name = studies[i] if studies else "study_{}".format(i+1)
        
        log.write(" -Study {} ({}): Getting LD map from {}...".format(i+1, study_name, ld_source_types[i]), verbose=verbose)
        
        # Get LD map from source (variant information only, not LD matrix yet)
        ld_map = _get_ld_map_from_source(
            source_type=ld_source_types[i],
            study_index=i,
            row=row,
            locus_sumstats=locus_sumstats,
            ld_map_path=ld_maps_list[i],
            ld_map_rename_dic=ld_map_dics_list[i],
            ld_map_kwargs=ld_map_kwargs,
            bfile=bfiles[i],
            pfile=pfiles[i],
            vcf=vcfs[i],
            windowsizekb=windowsizekb,
            plink=plink,
            plink2=plink2,
            threads=threads,
            memory=memory,
            overwrite=overwrite,
            tabix=tabix_list[i],
            log=log,
            verbose=verbose,
            **kwargs
        )
        
        all_ld_maps.append(ld_map)
    
    return all_ld_maps

def _find_common_variants_across_studies(
    n_studies: int,
    studies: List[str],
    locus_sumstats: pd.DataFrame,
    all_ld_maps: List[pd.DataFrame],
    suffixes: List[str],
    row: pd.Series,
    log: Log,
    verbose: bool
) -> pd.DataFrame:
    """
    Find common variants across all studies and sumstats.
    
    This is STEP 3: Determine intersection of variants based on:
    - CHR and POS match
    - EA-NEA match (perfect) OR EA-NEA flipped
    
    Returns:
        DataFrame with common variants, including _INDEX_BIM_{i} and _FLIPPED_{i} columns
    """
    log.write(" -Determining common SNPs across all studies and sumstats...", verbose=verbose)
    
    # Start with sumstats variants in locus
    common_variants = locus_sumstats.copy()
    
    # For each study, find matching variants
    for i in range(n_studies):
        gc.collect()
        study_name = studies[i] if studies else "study_{}".format(i+1)
        ld_map = all_ld_maps[i]
        
        log.write(" -Matching with study {} ({} variants in LD map)...".format(i+1, len(ld_map)), verbose=verbose)
        
        # Match on CHR and POS first
        # CRITICAL: Ensure CHR and POS types are compatible for merging
        # Keep CHR as int for consistent matching
        if "CHR" in common_variants.columns and "CHR" in ld_map.columns:
            # Convert both to int for consistent matching
            if pd.api.types.is_string_dtype(common_variants["CHR"]):
                common_variants["CHR"] = pd.to_numeric(common_variants["CHR"], errors='coerce').astype('Int64')
            elif not pd.api.types.is_integer_dtype(common_variants["CHR"]):
                common_variants["CHR"] = common_variants["CHR"].astype('Int64')
            
            if pd.api.types.is_string_dtype(ld_map["CHR"]):
                ld_map["CHR"] = pd.to_numeric(ld_map["CHR"], errors='coerce').astype('Int64')
            elif not pd.api.types.is_integer_dtype(ld_map["CHR"]):
                ld_map["CHR"] = ld_map["CHR"].astype('Int64')
        
        # Ensure POS is integer for both
        if "POS" in common_variants.columns and "POS" in ld_map.columns:
            if not pd.api.types.is_integer_dtype(common_variants["POS"]):
                common_variants["POS"] = common_variants["POS"].astype(int)
            if not pd.api.types.is_integer_dtype(ld_map["POS"]):
                ld_map["POS"] = ld_map["POS"].astype(int)
        
        # Convert alleles to string for reliable comparison
        common_variants["EA"] = common_variants["EA"].astype("string")
        common_variants["NEA"] = common_variants["NEA"].astype("string")
        ld_map["EA_bim"] = ld_map["EA_bim"].astype("string")
        ld_map["NEA_bim"] = ld_map["NEA_bim"].astype("string")
        
        # Merge on CHR and POS (inner join - only keeps matching positions)
        combined = pd.merge(common_variants, ld_map, on=["CHR", "POS"], how="inner", suffixes=("", "_ref"))
        
        if len(combined) == 0:
            log.write(" -No variants match on CHR/POS for study {}. No common variants.".format(i+1), verbose=verbose)
            common_variants = pd.DataFrame()
            break
        
        # Check allele matching: perfect match or flipped
        index_suffix = "_{}".format(i+1)
        combined["_FLIPPED" + index_suffix] = False
        
        perfect_match = ((combined["EA"] == combined["EA_bim"]) & (combined["NEA"] == combined["NEA_bim"]))
        flipped_match = ((combined["EA"] == combined["NEA_bim"]) & (combined["NEA"] == combined["EA_bim"]))
        allele_match = perfect_match | flipped_match
        
        combined.loc[flipped_match, "_FLIPPED" + index_suffix] = True
        
        log.write("   -Study {}: {} perfect matches, {} flipped matches, {} total matches".format(
            i+1, perfect_match.sum(), flipped_match.sum(), allele_match.sum()), verbose=verbose)
        
        # Keep only variants with valid allele matches
        common_variants = combined[allele_match].copy()
        
        # Store reference SNPID and index for later use
        if "SNPID_bim" in common_variants.columns:
            common_variants["SNPID_bim" + index_suffix] = common_variants["SNPID_bim"]
        if "_INDEX_BIM" in common_variants.columns:
            common_variants["_INDEX_BIM" + index_suffix] = common_variants["_INDEX_BIM"]
        
        # Keep only essential columns for next iteration
        keep_cols = ["SNPID", "CHR", "POS", "EA", "NEA"] + [col for col in common_variants.columns if "_FLIPPED_" in col or "_INDEX_BIM_" in col or "SNPID_bim_" in col]
        # Also keep study-specific sumstats columns
        for suffix in suffixes:
            for col in ["BETA", "SE", "Z", "EAF", "N", "P"]:
                if col + suffix in common_variants.columns:
                    keep_cols.append(col + suffix)
        
        common_variants = common_variants[[col for col in keep_cols if col in common_variants.columns]].copy()
        
        if len(common_variants) == 0:
            log.write(" -No common variants after matching with study {}.".format(i+1), verbose=verbose)
            break
    
    return common_variants

def _process_ld_matrices_for_all_studies(
    n_studies: int,
    studies: List[str],
    ld_source_types: List[str],
    common_variants: pd.DataFrame,
    row: pd.Series,
    ld_paths_list: List[Optional[str]],
    ld_maps_list: List[Optional[str]],
    ld_types_list: List[str],
    ld_if_square: bool,
    ld_if_add_T: bool,
    bfiles: List[Optional[str]],
    pfiles: List[Optional[str]],
    vcfs: List[Optional[str]],
    windowsizekb: int,
    out: str,
    group: str,
    plink: str,
    plink2: str,
    threads: int,
    mode: str,
    memory: Optional[int],
    overwrite: bool,
    extra_plink_option: str,
    tabix_list: List[Optional[Union[bool, str]]],
    suffixes: List[str],
    log: Log,
    verbose: bool,
    **kwargs: Any
) -> Tuple[List[np.ndarray], str, pd.DataFrame]:
    """
    Process LD matrices for all studies.
    
    This is STEP 4: For each study, get LD matrix for common variants,
    extract sub-matrix, apply allele flip corrections, and export.
    
    Returns:
        Tuple of (list of LD matrices, plink_log, output_file_list)
    """
    all_ld_matrices = []
    plink_log = ""
    output_file_list = pd.DataFrame(columns=["SNPID","SNPID_LIST","LD_R_MATRIX","LOCUS_SUMSTATS"])
    
    log.write(" -Processing LD matrices for each study...", verbose=verbose)
    
    for i in range(n_studies):
        gc.collect()
        study_name = studies[i] if studies else "study_{}".format(i+1)
        
        log.write(" -Study {} ({}): Computing LD matrix from {}...".format(i+1, study_name, ld_source_types[i]), verbose=verbose)
        
        # Get LD matrix for common variants only
        r_matrix, study_plink_log, matched_sumstats_from_vcf = _get_ld_matrix_for_common_variants(
            source_type=ld_source_types[i],
            study_index=i,
            row=row,
            common_variants=common_variants,
            ld_path=ld_paths_list[i],
            ld_map_path=ld_maps_list[i],
            ld_fmt=ld_types_list[i],
            ld_if_square=ld_if_square,
            ld_if_add_T=ld_if_add_T,
            bfile=bfiles[i],
            pfile=pfiles[i],
            vcf=vcfs[i],
            windowsizekb=windowsizekb,
            out=out,
            group=group,
            study=study_name,
            plink=plink,
            plink2=plink2,
            threads=threads,
            mode=mode,
            memory=memory,
            overwrite=overwrite,
            extra_plink_option=extra_plink_option,
            tabix=tabix_list[i],
            log=log,
            verbose=verbose,
            **kwargs
        )
        
        plink_log += study_plink_log
        
        # CRITICAL: For VCF sources, update indices to sequential
        index_col = "_INDEX_BIM_{}".format(i + 1)
        if ld_source_types[i] == "vcf" and r_matrix.shape[0] == len(common_variants):
            if index_col in common_variants.columns:
                common_variants[index_col] = range(len(common_variants))
                log.write(" -Updated _INDEX_BIM for VCF source to sequential indices (0-{})".format(
                    len(common_variants) - 1), verbose=verbose)
            elif ld_source_types[i] in ["bfile", "pfile"] and r_matrix.shape[0] == len(common_variants):
                if index_col in common_variants.columns:
                    max_index = common_variants[index_col].max() if len(common_variants) > 0 else -1
                    if max_index >= r_matrix.shape[0]:
                        common_variants[index_col] = range(len(common_variants))
        
        all_ld_matrices.append(r_matrix)
        
        # Extract and export LD matrix for this study
        matched_ld_matrix_path = _extract_variants_from_ld_matrix_m(
            merged_sumstats=common_variants, 
            r_matrix=r_matrix, 
            out=out, 
            group=group, 
            row=row, 
            windowsizekb=windowsizekb, 
            index=i,
            log=log, 
            verbose=verbose
        )
        
        # Export sumstats for this study
        matched_sumstats_path = _export_locus_sumstats_for_study(
            common_variants=common_variants,
            study_index=i,
            out=out,
            group=group,
            row=row,
            windowsizekb=windowsizekb,
            suffixes=suffixes,
            log=log
        )
        
        # Add to file list
        row_dict = {
            "SNPID": row["SNPID"],
            "SNPID_LIST": "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), group, row["SNPID"], windowsizekb),
            "LD_R_MATRIX": matched_ld_matrix_path,
            "LOCUS_SUMSTATS": matched_sumstats_path,
            "LOCUS": row["SNPID"],
            "SUBSTUDY": i+1,
            "STUDY": study_name
        }
        file_row = pd.Series(row_dict).to_frame().T
        output_file_list = pd.concat([output_file_list, file_row], ignore_index=True)
    
    return all_ld_matrices, plink_log, output_file_list

def _detect_ld_source_type(
    ld_path: Optional[str] = None,
    ld_map_path: Optional[str] = None,
    bfile: Optional[str] = None,
    pfile: Optional[str] = None,
    vcf: Optional[str] = None
) -> str:
    """
    Detect the type of LD source based on provided parameters.
    
    Parameters
    ----------
    ld_path : Optional[str]
        Path to pre-computed LD matrix file
    ld_map_path : Optional[str]
        Path to LD map file (variant annotation file)
    bfile : Optional[str]
        Path to PLINK bfile prefix
    pfile : Optional[str]
        Path to PLINK pfile prefix
    vcf : Optional[str]
        Path to VCF file
    
    Returns
    -------
    str
        LD source type: "precomputed", "bfile", "pfile", or "vcf"
    """
    if ld_path is not None and ld_map_path is not None:
        return "precomputed"
    elif bfile is not None:
        return "bfile"
    elif pfile is not None:
        return "pfile"
    elif vcf is not None:
        return "vcf"
    else:
        raise ValueError("No LD source provided. Please provide one of: ld_path+ld_map_path, bfile, pfile, or vcf")

def _get_ld_map_from_source(
    source_type: str,
    study_index: int,
    row: pd.Series,
    locus_sumstats: pd.DataFrame,
    ld_map_path: Optional[str] = None,
    ld_map_rename_dic: Optional[Union[Dict[str, str], List[str]]] = None,
    ld_map_kwargs: Optional[Dict[str, Any]] = None,
    bfile: Optional[str] = None,
    pfile: Optional[str] = None,
    vcf: Optional[str] = None,
    windowsizekb: int = 1000,
    plink: str = "plink",
    plink2: str = "plink2",
    threads: int = 1,
    memory: Optional[int] = None,
    overwrite: bool = False,
    tabix: Optional[Union[bool, str]] = None,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs: Any
) -> pd.DataFrame:
    """
    Get LD map (variant information) from any source type without calculating LD matrix.
    
    This is used in the first pass to determine common variants across all studies.
    Only variant information (CHR, POS, EA, NEA) is extracted, not LD values.
    
    Parameters
    ----------
    source_type : str
        Type of LD source: "precomputed", "bfile", "pfile", or "vcf"
    study_index : int
        Index of the study (0-based)
    row : pd.Series
        Lead variant information (CHR, POS, SNPID)
    locus_sumstats : pd.DataFrame
        Summary statistics for the locus
    ld_map_path : Optional[str]
        Path to LD map file (for precomputed source)
    ld_map_rename_dic : Optional[Union[Dict[str, str], List[str]]]
        Dictionary mapping column names in LD map files
    ld_map_kwargs : Optional[Dict[str, Any]]
        Additional keyword arguments for loading LD map files
    bfile : Optional[str]
        Path to PLINK bfile prefix (for bfile source)
    pfile : Optional[str]
        Path to PLINK pfile prefix (for pfile source)
    vcf : Optional[str]
        Path to VCF file (for vcf source)
    windowsizekb : int
        Window size in kilobases
    plink : str
        Path to plink executable
    plink2 : str
        Path to plink2 executable
    threads : int
        Number of threads
    memory : Optional[int]
        Memory limit in MB
    overwrite : bool
        Whether to overwrite existing files
    tabix : Optional[Union[bool, str]]
        Tabix executable path or True to auto-detect. If None, auto-detects.
        Used for efficient VCF region queries when reading VCF files.
    log : Log
        Log object
    verbose : bool
        Verbose output
    **kwargs : Any
        Additional keyword arguments
    
    Returns
    -------
    pd.DataFrame
        LD map with columns: SNPID_bim, CHR, POS, EA_bim, NEA_bim, _INDEX_BIM
    """
    if source_type == "precomputed":
        if ld_map_path is None:
            raise ValueError("ld_map_path is required for precomputed LD source")
        log.write(" -Loading LD map for study {}...".format(study_index + 1), verbose=verbose)
        ld_map = _load_ld_map(ld_map_path, ld_map_rename_dic=ld_map_rename_dic, **ld_map_kwargs if ld_map_kwargs else {})
        if "_INDEX_BIM" not in ld_map.columns:
            ld_map["_INDEX_BIM"] = range(len(ld_map))
        return ld_map
    
    elif source_type in ["bfile", "pfile"]:
        # Logging already done at caller level
        
        # Process PLINK input files and get BIM/PVAR
        bfile_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=[row["CHR"]],
            bfile=bfile if source_type == "bfile" else None,
            pfile=pfile if source_type == "pfile" else None,
            vcf=None,
            plink_log="",
            threads=threads,
            log=log,
            load_bim=True,
            overwrite=overwrite,
            plink=plink,
            plink2=plink2,
            convert="bfile",
            **kwargs
        )
        
        if len(ref_bims) == 0:
            raise ValueError("No BIM/PVAR data loaded for study {}".format(study_index + 1))
        
        ref_bim = ref_bims[0]
        
        # Filter to locus window
        locus_ref_bim = ref_bim[
            (ref_bim["CHR_bim"] == str(row["CHR"])) &
            (ref_bim["POS_bim"] >= row["POS"] - windowsizekb * 1000) &
            (ref_bim["POS_bim"] <= row["POS"] + windowsizekb * 1000)
        ].copy()
        
        # Create LD map from BIM/PVAR
        ld_map = pd.DataFrame()
        ld_map["SNPID_bim"] = locus_ref_bim["SNPID"]
        # Keep CHR as int (convert from string if needed)
        if pd.api.types.is_string_dtype(locus_ref_bim["CHR_bim"]):
            ld_map["CHR"] = pd.to_numeric(locus_ref_bim["CHR_bim"], errors='coerce').astype('Int64')
        else:
            ld_map["CHR"] = locus_ref_bim["CHR_bim"].astype('Int64')
        ld_map["POS"] = locus_ref_bim["POS_bim"]
        ld_map["EA_bim"] = locus_ref_bim["EA_bim"]
        ld_map["NEA_bim"] = locus_ref_bim["NEA_bim"]
        ld_map["_INDEX_BIM"] = range(len(ld_map))
        
        return ld_map
    
    elif source_type == "vcf":
        from gwaslab.io.io_vcf import read_vcf
        from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
        
        # Define region for VCF extraction
        # Ensure start position is at least 1 (chromosomes start at position 1)
        region_start = max(1, row["POS"] - windowsizekb * 1000)
        region_end = row["POS"] + windowsizekb * 1000
        region = (row["CHR"], region_start, region_end)
        
        # Get mapper and convert chromosome
        mapper = ChromosomeMapper(log=log, verbose=False)  # Reduce verbosity
        mapper.detect_reference_format(vcf)
        region_chr_ref = mapper.sumstats_to_reference(region[0], reference_file=vcf, as_string=True)
        region_str = f"{region_chr_ref}:{region[1]}-{region[2]}"
        
        # Auto-detect tabix if not provided (only log once if not already detected)
        if tabix is None:
            tabix = which("tabix")
            if not tabix:
                log.write(" -tabix not found, VCF will be read without indexing (slower)", verbose=verbose)
        elif tabix is True:
            tabix = which("tabix")
        
        # Read VCF (tabix enables efficient region queries for indexed VCF files)
        ref_genotype = read_vcf(vcf, region=region_str, tabix=tabix)
        if ref_genotype is None or len(ref_genotype.get("variants/POS", [])) == 0:
            raise ValueError("No variants found in VCF for study {}".format(study_index + 1))
        
        # Create LD map from VCF variants
        # Keep CHR as int (region[0] should already be int)
        ld_map = pd.DataFrame()
        ld_map["CHR"] = pd.Series([int(region[0])] * len(ref_genotype["variants/POS"]), dtype='Int64')
        ld_map["POS"] = ref_genotype["variants/POS"]
        ld_map["EA_bim"] = [alt[0] if len(alt) > 0 else "N" for alt in ref_genotype["variants/ALT"]]
        ld_map["NEA_bim"] = ref_genotype["variants/REF"]
        ld_map["SNPID_bim"] = ref_genotype.get("variants/ID", ["{}:{}".format(region[0], pos) for pos in ld_map["POS"]])
        ld_map["_INDEX_BIM"] = range(len(ld_map))
        
        return ld_map
    
    else:
        raise ValueError("Unknown LD source type: {}".format(source_type))

def _get_ld_matrix_and_map_from_source(
    source_type: str,
    study_index: int,
    row: pd.Series,
    locus_sumstats: pd.DataFrame,
    ld_path: Optional[str] = None,
    ld_map_path: Optional[str] = None,
    ld_fmt: str = "npz",
    ld_if_square: bool = False,
    ld_if_add_T: bool = False,
    ld_map_rename_dic: Optional[Union[Dict[str, str], List[str]]] = None,
    ld_map_kwargs: Optional[Dict[str, Any]] = None,
    bfile: Optional[str] = None,
    pfile: Optional[str] = None,
    vcf: Optional[str] = None,
    windowsizekb: int = 1000,
    out: str = "./",
    group: Optional[str] = None,
    study: Optional[str] = None,
    plink: str = "plink",
    plink2: str = "plink2",
    threads: int = 1,
    mode: str = "r",
    memory: Optional[int] = None,
    overwrite: bool = False,
    extra_plink_option: str = "",
    tabix: Optional[Union[bool, str]] = None,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs: Any
) -> Tuple[np.ndarray, pd.DataFrame, str]:
    """
    Get LD matrix and map from any source type (pre-computed, bfile, pfile, or vcf).
    
    This is a unified interface that handles different LD sources by:
    - Loading pre-computed LD matrices and maps
    - Calculating LD from bfile/pfile using PLINK
    - Calculating LD from VCF using allel
    
    Parameters
    ----------
    source_type : str
        Type of LD source: "precomputed", "bfile", "pfile", or "vcf"
    study_index : int
        Index of the study (0-based)
    row : pd.Series
        Lead variant information (CHR, POS, SNPID)
    locus_sumstats : pd.DataFrame
        Summary statistics for the locus
    ld_path : Optional[str]
        Path to pre-computed LD matrix (for precomputed source)
    ld_map_path : Optional[str]
        Path to LD map file (for precomputed source)
    ld_fmt : str
        Format of pre-computed LD matrix: "npz" or "txt"
    bfile : Optional[str]
        Path to PLINK bfile prefix (for bfile source)
    pfile : Optional[str]
        Path to PLINK pfile prefix (for pfile source)
    vcf : Optional[str]
        Path to VCF file (for vcf source)
    windowsizekb : int
        Window size in kilobases
    out : str
        Output directory
    group : Optional[str]
        Group identifier for output naming
    study : Optional[str]
        Study name for output naming
    plink : str
        Path to plink executable
    plink2 : str
        Path to plink2 executable
    threads : int
        Number of threads
    mode : str
        PLINK mode for LD calculation
    memory : Optional[int]
        Memory limit in MB
    overwrite : bool
        Whether to overwrite existing files
    extra_plink_option : str
        Additional PLINK options
    tabix : Optional[Union[bool, str]]
        Tabix executable path or True to auto-detect. If None, auto-detects.
        Used for efficient VCF region queries when reading VCF files.
    log : Log
        Log object
    verbose : bool
        Verbose output
    **kwargs : Any
        Additional keyword arguments
    
    Returns
    -------
    Tuple[np.ndarray, pd.DataFrame, str]
        (ld_matrix, ld_map, plink_log) where:
        - ld_matrix: LD correlation matrix as numpy array
        - ld_map: DataFrame with variant information matching the LD matrix
        - plink_log: PLINK log output (empty string if PLINK not used)
    """
    plink_log = ""
    
    if source_type == "precomputed":
        # Load pre-computed LD matrix and map
        if ld_path is None or ld_map_path is None:
            raise ValueError("ld_path and ld_map_path are required for precomputed LD source")
        
        log.write(" -Loading pre-computed LD matrix and map for study {}...".format(study_index + 1), verbose=verbose)
        r_matrix = _load_ld_matrix(ld_path, fmt=ld_fmt, if_square=ld_if_square, if_add_T=ld_if_add_T, log=log, verbose=verbose)
        ld_map = _load_ld_map(ld_map_path, ld_map_rename_dic=ld_map_rename_dic, **ld_map_kwargs if ld_map_kwargs else {})
        
        # Ensure _INDEX_BIM is set to match LD matrix row/column indices
        # For precomputed sources, the LD map should be in the same order as the LD matrix
        if "_INDEX_BIM" not in ld_map.columns:
            ld_map["_INDEX_BIM"] = range(len(ld_map))
        
        # Verify dimensions match
        if len(ld_map) != r_matrix.shape[0]:
            log.warning(" -Warning: LD map has {} variants but LD matrix has {} rows. They should match.".format(
                len(ld_map), r_matrix.shape[0]), verbose=verbose)
        
        return r_matrix, ld_map, plink_log
    
    elif source_type in ["bfile", "pfile"]:
        # Calculate LD from bfile/pfile using PLINK
        log.write(" -Calculating LD from {} for study {}...".format(source_type, study_index + 1), verbose=verbose)
        
        # Process PLINK input files and get BIM/PVAR
        bfile_prefix, plink_log, ref_bims, filetype = _process_plink_input_files(
            chrlist=[row["CHR"]],
            bfile=bfile if source_type == "bfile" else None,
            pfile=pfile if source_type == "pfile" else None,
            vcf=None,
            plink_log=plink_log,
            threads=threads,
            log=log,
            load_bim=True,
            overwrite=overwrite,
            plink=plink,
            plink2=plink2,
            convert="bfile",  # Always convert to bfile for PLINK1 LD calculation
            **kwargs
        )
        
        if len(ref_bims) == 0:
            raise ValueError("No BIM/PVAR data loaded for study {}".format(study_index + 1))
        
        ref_bim = ref_bims[0]
        
        # Align sumstats with BIM
        matched_sumstats = _align_sumstats_with_bim(
            row=row,
            locus_sumstats=locus_sumstats,
            ref_bim=ref_bim,
            log=log,
            suffixes=None
        )
        
        if len(matched_sumstats) == 0:
            raise ValueError("No matching variants found between sumstats and reference for study {}".format(study_index + 1))
        
        # Export SNP list and sumstats for PLINK
        study_name = study if study else "study_{}".format(study_index + 1)
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), study_name, row["SNPID"], windowsizekb)
        matched_sumstats_path = "{}/{}_{}_{}.sumstats".format(out.rstrip("/"), study_name, row["SNPID"], windowsizekb)
        
        matched_sumstats["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)
        matched_sumstats[["SNPID", "CHR", "POS", "EA", "NEA"]].to_csv(matched_sumstats_path, sep="\t", index=None)
        
        # Calculate LD matrix using PLINK
        matched_ld_matrix_path, plink_log = _calculate_ld_r(
            study=study_name,
            matched_sumstats_snpid=matched_sumstats["SNPID"],
            row=row,
            bfile_prefix=bfile_prefix,
            threads=threads,
            windowsizekb=windowsizekb,
            out=out,
            plink_log=plink_log,
            log=log,
            memory=memory,
            mode=mode,
            filetype=filetype,
            plink=plink,
            plink2=plink2,
            ref_allele_path=matched_sumstats_path,
            extra_plink_option=extra_plink_option,
            verbose=verbose
        )
        
        # Load the calculated LD matrix
        # The matrix is square, with rows/cols in the same order as the SNP list provided to PLINK
        r_matrix = np.loadtxt(matched_ld_matrix_path)
        
        # Create LD map from matched sumstats (matching the order of LD matrix)
        # CRITICAL: The LD matrix from PLINK matches the order of variants in the snplist
        # We create the LD map in the same order, and set _INDEX_BIM = [0, 1, 2, ...]
        # This ensures _INDEX_BIM points to the correct row/column in the LD matrix
        ld_map = matched_sumstats[["SNPID", "CHR", "POS", "EA", "NEA"]].copy()
        ld_map = ld_map.rename(columns={"EA": "EA_bim", "NEA": "NEA_bim"})
        ld_map["SNPID_bim"] = ld_map["SNPID"]
        # Set indices to match matrix rows/cols (0-indexed, matching numpy array indexing)
        ld_map["_INDEX_BIM"] = range(len(ld_map))
        
        return r_matrix, ld_map, plink_log
    
    elif source_type == "vcf":
        # Calculate LD from VCF using allel
        log.write(" -Calculating LD from VCF for study {}...".format(study_index + 1), verbose=verbose)
        
        from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
        
        # Auto-detect tabix if not provided (only log if not found)
        if tabix is None:
            tabix = which("tabix")
            if not tabix:
                log.write(" -tabix not found, VCF will be read without indexing (slower)", verbose=verbose)
        elif tabix is True:
            tabix = which("tabix")
        
        # Define region for VCF extraction
        # Ensure start position is at least 1 (chromosomes start at position 1)
        region_start = max(1, row["POS"] - windowsizekb * 1000)
        region_end = row["POS"] + windowsizekb * 1000
        region = (row["CHR"], region_start, region_end)
        
        # Get LD matrix from VCF (tabix enables efficient region queries for indexed VCF files)
        matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
            sumstats_or_dataframe=locus_sumstats,
            vcf_path=vcf,
            region=region,
            log=log,
            verbose=verbose,
            pos="POS",
            nea="NEA",
            ea="EA",
            tabix=tabix
        )
        
        if len(matched_sumstats) == 0:
            raise ValueError("No matching variants found in VCF for study {}".format(study_index + 1))
        
        # Create LD map from matched sumstats
        # CRITICAL: The LD matrix from _get_ld_matrix_from_vcf has rows/cols in the same order
        # as matched_sumstats. We create the LD map in the same order and set _INDEX_BIM accordingly.
        ld_map = matched_sumstats[["SNPID", "CHR", "POS", "EA", "NEA"]].copy()
        ld_map = ld_map.rename(columns={"EA": "EA_bim", "NEA": "NEA_bim"})
        ld_map["SNPID_bim"] = ld_map["SNPID"]
        # Set indices to match matrix rows/cols (0-indexed, matching numpy array indexing)
        ld_map["_INDEX_BIM"] = range(len(ld_map))
        
        # Note: _get_ld_matrix_from_vcf returns r^2, but PLINK and other sources return r (correlation)
        # Convert r^2 to |r| for consistency with PLINK output format
        # We lose sign information when converting from r^2, but this is acceptable for fine-mapping
        # as we're primarily interested in LD magnitude
        ld_matrix = np.sqrt(np.abs(ld_matrix))
        
        return ld_matrix, ld_map, plink_log
    
    else:
        raise ValueError("Unknown LD source type: {}".format(source_type))

def _get_ld_matrix_for_common_variants(
    source_type: str,
    study_index: int,
    row: pd.Series,
    common_variants: pd.DataFrame,
    ld_path: Optional[str] = None,
    ld_map_path: Optional[str] = None,
    ld_fmt: str = "npz",
    ld_if_square: bool = False,
    ld_if_add_T: bool = False,
    bfile: Optional[str] = None,
    pfile: Optional[str] = None,
    vcf: Optional[str] = None,
    windowsizekb: int = 1000,
    out: str = "./",
    group: Optional[str] = None,
    study: Optional[str] = None,
    plink: str = "plink",
    plink2: str = "plink2",
    threads: int = 1,
    mode: str = "r",
    memory: Optional[int] = None,
    overwrite: bool = False,
    extra_plink_option: str = "",
    tabix: Optional[Union[bool, str]] = None,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs: Any
) -> Tuple[np.ndarray, str]:
    """
    Get LD matrix for common variants only.
    
    For precomputed: Load full matrix and extract sub-matrix for common variants.
    For bfile/pfile/vcf: Calculate LD matrix using only common variants.
    
    Returns:
        (ld_matrix, plink_log) where ld_matrix has shape (len(common_variants), len(common_variants))
    """
    plink_log = ""
    index_suffix = "_{}".format(study_index + 1)
    index_col = "_INDEX_BIM" + index_suffix
    
    if source_type == "precomputed":
        # Load full LD matrix
        r_matrix = _load_ld_matrix(ld_path, fmt=ld_fmt, if_square=ld_if_square, if_add_T=ld_if_add_T, log=log, verbose=verbose)
        # Extract sub-matrix for common variants using original indices
        index_col = "_INDEX_BIM_{}".format(study_index + 1)
        indices = common_variants[index_col].values.astype(int)
        r_matrix = r_matrix[np.ix_(indices, indices)]
        # For precomputed, we use the original indices, so no need to return matched_sumstats
        return r_matrix, plink_log, None
    
    elif source_type in ["bfile", "pfile"]:
        # Calculate LD matrix using only common variants
        study_name = study if study else "study_{}".format(study_index + 1)
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), study_name, row["SNPID"], windowsizekb)
        matched_sumstats_path = "{}/{}_{}_{}.sumstats".format(out.rstrip("/"), study_name, row["SNPID"], windowsizekb)
        
        # Export SNP list and sumstats for PLINK
        common_variants["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)
        common_variants[["SNPID", "CHR", "POS", "EA", "NEA"]].to_csv(matched_sumstats_path, sep="\t", index=None)
        
        # Get bfile prefix
        bfile_prefix, _, _, filetype = _process_plink_input_files(
            chrlist=[row["CHR"]],
            bfile=bfile if source_type == "bfile" else None,
            pfile=pfile if source_type == "pfile" else None,
            vcf=None,
            plink_log="",
            threads=threads,
            log=log,
            load_bim=False,
            overwrite=overwrite,
            plink=plink,
            plink2=plink2,
            convert="bfile",
            **kwargs
        )
        
        # Calculate LD matrix
        matched_ld_matrix_path, plink_log = _calculate_ld_r(
            study=study_name,
            matched_sumstats_snpid=common_variants["SNPID"],
            row=row,
            bfile_prefix=bfile_prefix,
            threads=threads,
            windowsizekb=windowsizekb,
            out=out,
            plink_log=plink_log,
            log=log,
            memory=memory,
            mode=mode,
            filetype=filetype,
            plink=plink,
            plink2=plink2,
            ref_allele_path=matched_sumstats_path,
            extra_plink_option=extra_plink_option,
            verbose=verbose
        )
        
        # Load the calculated LD matrix
        # For bfile/pfile, the LD matrix is calculated for the matched variants only,
        # so it's already in the correct order and size matches common_variants
        r_matrix = np.loadtxt(matched_ld_matrix_path)
        return r_matrix, plink_log, None
    
    elif source_type == "vcf":
        # Calculate LD matrix from VCF using only common variants
        from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
        
        # Auto-detect tabix if not provided (only log if not found)
        if tabix is None:
            tabix = which("tabix")
            if not tabix:
                log.write(" -tabix not found, VCF will be read without indexing (slower)", verbose=verbose)
        elif tabix is True:
            tabix = which("tabix")
        
        # CRITICAL: _get_ld_matrix_from_vcf filters sumstats to region using _filter_region,
        # which uses strict boundaries (POS > start and POS < end, excluding boundaries).
        # However, _extract_variants_in_locus uses (POS >= start - window and POS < start + window).
        # To ensure common_variants pass the filter, we expand the region slightly to account
        # for boundary exclusion. We add 1bp to start (so POS > start-1 includes POS >= start)
        # and subtract 1bp from end (so POS < end+1 includes POS < end).
        # Also ensure start position is at least 1 (chromosomes start at position 1)
        region_start = max(1, row["POS"] - windowsizekb * 1000 - 1)  # Expand by 1bp to account for > vs >=
        region_end = row["POS"] + windowsizekb * 1000 + 1     # Expand by 1bp to account for < vs <=
        region = (row["CHR"], region_start, region_end)
        
        # Keep CHR as int - _filter_region will handle normalization internally
        # Get LD matrix from VCF (will match common variants)
        # tabix enables efficient region queries for indexed VCF files
        matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
            sumstats_or_dataframe=common_variants,
            vcf_path=vcf,
            region=region,
            log=log,
            verbose=verbose,
            pos="POS",
            nea="NEA",
            ea="EA",
            tabix=tabix
        )
        
        # Check if we got valid results
        if len(matched_sumstats) == 0 or ld_matrix.size == 0:
            raise ValueError("No variants matched in VCF for study {} after region filtering. "
                           "This may indicate a region format mismatch or that all variants were filtered out. "
                           "Common variants count: {}, Region: {}".format(
                               study_index + 1, len(common_variants), region))
        
        # CRITICAL: For VCF sources, the LD matrix returned by _get_ld_matrix_from_vcf
        # is already filtered to matched variants and is in the same order as matched_sumstats.
        # The matrix shape (n, n) matches len(matched_sumstats), so we need to use sequential
        # indices (0, 1, 2, ...) rather than the original _INDEX_BIM from the LD map extraction.
        # 
        # We return the matched_sumstats so the caller can update _INDEX_BIM accordingly.
        # For now, we'll handle this in the caller by checking if matrix shape matches variant count.
        
        # Convert r^2 to |r|
        ld_matrix = np.sqrt(np.abs(ld_matrix))
        
        # Return matched_sumstats so caller can verify order and update indices if needed
        # The matched_sumstats should be in the same order as the LD matrix rows/cols
        return ld_matrix, plink_log, matched_sumstats
    
    else:
        raise ValueError("Unknown LD source type: {}".format(source_type))

def _export_locus_sumstats_for_study(
    common_variants: pd.DataFrame,
    study_index: int,
    out: str,
    group: str,
    row: pd.Series,
    windowsizekb: int,
    suffixes: List[str],
    log: Log
) -> str:
    """
    Export locus sumstats for a single study.
    
    Returns:
        Path to exported sumstats file
    """
    study_suffix = suffixes[study_index] if study_index < len(suffixes) else "_{}".format(study_index + 1)
    matched_sumstats_path = "{}/{}_{}_{}_{}.sumstats".format(out.rstrip("/"), group, row["SNPID"], windowsizekb, study_index + 1)
    
    to_export_columns = ["CHR", "POS", "EA", "NEA"]
    
    if "Z" + study_suffix in common_variants.columns:
        to_export_columns.append("Z" + study_suffix)
    if ("BETA" + study_suffix in common_variants.columns) and ("SE" + study_suffix in common_variants.columns):
        to_export_columns.append("BETA" + study_suffix)
        to_export_columns.append("SE" + study_suffix)
    if "EAF" + study_suffix in common_variants.columns:
        to_export_columns.append("EAF" + study_suffix)
    if "N" + study_suffix in common_variants.columns:
        to_export_columns.append("N" + study_suffix)
    
    # Export path is already clear from file path, no need for verbose message
    rename_dic = {
        "BETA" + study_suffix: "Beta",
        "SE" + study_suffix: "Se",
        "SNPID": "SNP"
    }
    
    export_df = common_variants[["SNPID"] + to_export_columns].rename(columns=rename_dic)
    export_df.to_csv(matched_sumstats_path, sep="\t", index=None)
    export_df.to_csv(matched_sumstats_path + ".gz", sep="\t", index=None)
    
    return matched_sumstats_path + ".gz"

@with_logging(
        start_to_msg="calculate LD matrix",
        finished_msg="calculating LD matrix",
        start_cols=["SNPID","CHR","POS","EA","NEA"],
        start_function=".calculate_ld_matrix()"
)
def tofinemapping_m(
    sumstats: pd.DataFrame, 
    studies: Optional[List[str]] = None, 
    group: Optional[str] = None,
    ld_paths: Optional[List[str]] = None,
    ld_types: Optional[List[str]] = None, 
    ld_maps: Optional[List[str]] = None,
    ld_map_dics: Optional[List[Dict[str, str]]] = None,
    bfile: Optional[Union[str, List[str]]] = None, 
    pfile: Optional[Union[str, List[str]]] = None,
    vcf: Optional[Union[str, List[str]]] = None, 
    locus: Optional[Tuple[int, int]] = None,
    loci: Optional[List[str]] = None,
    loci_chrpos: Optional[List[Tuple[int, int]]] = None,
    out: str = "./",
    plink: str = "plink",
    plink2: str = "plink2",
    windowsizekb: int = 1000,
    threads: int = 1, 
    mode: str = "r",
    exclude_hla: bool = False, 
    getlead_kwargs: Optional[Dict[str, Any]] = None, 
    memory: Optional[int] = None, 
    overwrite: bool = False,
    log: Log = Log(),
    suffixes: Optional[List[str]] = None,
    ld_map_kwargs: Optional[Dict[str, Any]] = None,
    ld_if_square: bool = False,
    ld_if_add_T: bool = False,
    extra_plink_option: str = "",
    tabix: Optional[Union[bool, str, List[Union[bool, str]]]] = None,
    verbose: bool = True,
    **kwargs: Any
) -> Tuple[Optional[str], pd.DataFrame, str]:
    """
    Prepare data for multi-study fine-mapping by matching sumstats with LD reference files
    and extracting LD matrices for specified loci.
    
    This function is designed for fine-mapping analyses involving multiple studies (typically two).
    It processes significant loci, matches variants between sumstats and LD reference files,
    extracts LD matrices for each study, and exports the necessary files for downstream
    fine-mapping tools (e.g., mesusie, coloc).
    
    The function handles multiple LD source types:
    - Pre-computed LD matrices (npz or txt format) with LD map files
    - PLINK bfile format (bed/bim/fam)
    - PLINK pfile format (pgen/pvar/psam)
    - VCF files (BCF supported via tabix indexing for fast region queries)
    
    Each study can use a different LD source type. The function automatically detects
    the LD source type for each study based on provided parameters.
    
    The function handles:
    - Automatic extraction of significant loci if not provided
    - Matching variants between sumstats and LD reference files by CHR, POS, and alleles
    - Handling allele flips between reference and sumstats
    - Extracting/calculating LD matrices from various sources for each study
    - Exporting SNP lists, LD matrices, and locus-specific sumstats files
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe with columns: SNPID, CHR, POS, EA, NEA, and study-specific
        columns (e.g., BETA_1, SE_1, Z_1, BETA_2, SE_2, Z_2).
    studies : Optional[List[str]], optional
        List of study names (typically 2 studies). Defaults to None.
    group : Optional[str], optional
        Group identifier for output file naming. Defaults to None.
    ld_paths : Optional[List[str]], optional
        List of paths to pre-computed LD matrix files (one per study). 
        Required if using pre-computed LD matrices. Defaults to None.
    ld_types : Optional[List[str]], optional
        List of LD matrix file types ("npz" or "txt", one per study). 
        Defaults to None (assumes "npz").
    ld_maps : Optional[List[str]], optional
        List of paths to LD map files (variant annotation files matching LD matrices,
        one per study). Required if using pre-computed LD matrices. Defaults to None.
    ld_map_dics : Optional[List[Dict[str, str]]], optional
        List of dictionaries mapping column names in LD map files to standard names.
        Defaults to None.
    bfile : Optional[Union[str, List[str]]], optional
        Path(s) to PLINK binary files (bfile prefix). Can be a single string (applied to all studies)
        or a list (one per study). Defaults to None.
    pfile : Optional[Union[str, List[str]]], optional
        Path(s) to PLINK pfile format (pgen/pvar/psam). Can be a single string (applied to all studies)
        or a list (one per study). Defaults to None.
    vcf : Optional[Union[str, List[str]]], optional
        Path(s) to VCF file(s). Can be a single string (applied to all studies)
        or a list (one per study). For best performance, VCF files should be 
        tabix-indexed (use `tabix -p vcf file.vcf.gz`). Defaults to None.
    locus : Optional[Tuple[int, int]], optional
        Single locus as (CHR, POS) tuple. Defaults to None.
    loci : Optional[List[str]], optional
        List of lead variant SNPIDs to process. If None, significant loci are auto-extracted.
        Defaults to None.
    loci_chrpos : Optional[List[Tuple[int, int]]], optional
        List of loci as (CHR, POS) tuples. Defaults to None.
    out : str, optional
        Output directory path. Defaults to "./".
    plink : str, optional
        Path to plink executable. Defaults to "plink".
    plink2 : str, optional
        Path to plink2 executable. Defaults to "plink2".
    windowsizekb : int, optional
        Window size in kilobases around each lead variant for extracting variants.
        Defaults to 1000.
    threads : int, optional
        Number of threads for parallel processing. Defaults to 1.
    mode : str, optional
        Processing mode. Defaults to "r".
    exclude_hla : bool, optional
        Whether to exclude HLA region variants. Defaults to False.
    getlead_kwargs : Optional[Dict[str, Any]], optional
        Additional keyword arguments for lead variant extraction. Defaults to None.
    memory : Optional[int], optional
        Memory limit in MB. Defaults to None.
    overwrite : bool, optional
        Whether to overwrite existing output files. Defaults to False.
    log : Log, optional
        Log object for recording progress. Defaults to Log().
    suffixes : Optional[List[str]], optional
        List of suffixes for study-specific columns (e.g., ["_1", "_2"]).
        If None, defaults to [""]. Defaults to None.
    ld_map_kwargs : Optional[Dict[str, Any]], optional
        Additional keyword arguments for loading LD map files. Defaults to None.
    ld_if_square : bool, optional
        Whether to square the LD matrix elements (for pre-computed matrices). Defaults to False.
    ld_if_add_T : bool, optional
        Whether to add transpose to LD matrix (for pre-computed matrices). Defaults to False.
    extra_plink_option : str, optional
        Additional options to pass to PLINK commands. Defaults to "".
    tabix : Optional[Union[bool, str, List[Union[bool, str]]]], optional
        Tabix executable path(s) or True to auto-detect. Can be a single value (applied to all studies)
        or a list (one per study). If None, auto-detects for each VCF file.
        Used for efficient VCF region queries when reading VCF files.
        Defaults to None.
    verbose : bool, optional
        Whether to print verbose output. Defaults to True.
    **kwargs : Any
        Additional keyword arguments passed to PLINK processing functions.
    
    Returns
    -------
    Tuple[Optional[str], pd.DataFrame, str]
        A tuple containing:
        - output_file_list_path (Optional[str]): Path to the file list CSV, or None if no
          loci were processed successfully.
        - output_file_list (pd.DataFrame): DataFrame with columns: SNPID, SNPID_LIST,
          LD_R_MATRIX, LOCUS_SUMSTATS, LOCUS, SUBSTUDY, STUDY, GROUP. Each row represents
          one study-locus combination.
        - plink_log (str): PLINK log output (empty string if PLINK was not used).
    
    Notes
    -----
    - This function supports multiple studies (typically 2, but can handle any number).
    - LD source type is automatically detected for each study based on provided parameters.
    - Priority: ld_paths+ld_maps > bfile > pfile > vcf
    - Variants are matched by CHR, POS, and alleles. Allele flips are automatically detected
      and handled by negating the corresponding LD matrix entries.
    - Output files are named with pattern: {group}_{nstudy}_{lead_snpid}_{windowsizekb}kb.*
    
    Examples
    --------
    # Using pre-computed LD matrices
    >>> filelist_path, filelist_df, plink_log = tofinemapping_m(
    ...     sumstats=sumstats_df,
    ...     studies=["Study1", "Study2"],
    ...     group="my_analysis",
    ...     ld_paths=["/path/to/study1.ld.npz", "/path/to/study2.ld.npz"],
    ...     ld_maps=["/path/to/study1.map", "/path/to/study2.map"],
    ...     windowsizekb=500,
    ...     out="./finemapping_output"
    ... )
    
    # Using bfile for both studies
    >>> filelist_path, filelist_df, plink_log = tofinemapping_m(
    ...     sumstats=sumstats_df,
    ...     studies=["Study1", "Study2"],
    ...     group="my_analysis",
    ...     bfile="/path/to/reference",
    ...     windowsizekb=500,
    ...     out="./finemapping_output"
    ... )
    
    # Using different LD sources for each study
    >>> filelist_path, filelist_df, plink_log = tofinemapping_m(
    ...     sumstats=sumstats_df,
    ...     studies=["Study1", "Study2"],
    ...     group="my_analysis",
    ...     ld_paths=["/path/to/study1.ld.npz", None],
    ...     ld_maps=["/path/to/study1.map", None],
    ...     bfile=[None, "/path/to/study2_reference"],
    ...     windowsizekb=500,
    ...     out="./finemapping_output"
    ... )
    """
    
    ############################################################################################
    # STEP 0: Initialize and validate parameters
    ############################################################################################
    if suffixes is None:
        suffixes=[""]
    if getlead_kwargs is None:
        getlead_kwargs={"windowsizekb":1000}
    if ld_map_kwargs is None:
        ld_map_kwargs={}
    
    # Determine number of studies
    if studies is None:
        # Try to infer from other parameters
        n_studies = 2  # Default to 2
        if ld_paths is not None:
            n_studies = max(n_studies, len(ld_paths))
        if ld_maps is not None:
            n_studies = max(n_studies, len(ld_maps))
        if isinstance(bfile, list):
            n_studies = max(n_studies, len(bfile))
        if isinstance(pfile, list):
            n_studies = max(n_studies, len(pfile))
        if isinstance(vcf, list):
            n_studies = max(n_studies, len(vcf))
        studies = ["Study_{}".format(i+1) for i in range(n_studies)]
    else:
        n_studies = len(studies)
    
    if n_studies < 1:
        raise ValueError("At least one study must be specified")
    
    log.write(" -Number of studies: {}".format(n_studies), verbose=verbose)
    
    # Normalize LD source parameters to lists (one per study)
    bfiles, pfiles, vcfs, ld_paths_list, ld_maps_list, ld_types_list, ld_map_dics_list, tabix_list = _normalize_parameters_to_lists(
        n_studies=n_studies,
        bfile=bfile,
        pfile=pfile,
        vcf=vcf,
        ld_paths=ld_paths,
        ld_maps=ld_maps,
        ld_types=ld_types,
        ld_map_dics=ld_map_dics,
        tabix=tabix
    )
    
    # Detect LD source type for each study
    ld_source_types = _detect_ld_source_types(
        n_studies=n_studies,
        ld_paths_list=ld_paths_list,
        ld_maps_list=ld_maps_list,
        bfiles=bfiles,
        pfiles=pfiles,
        vcfs=vcfs,
        log=log,
        verbose=verbose
    )
    
    ############################################################################################
    # STEP 1: Filter and prepare sumstats
    ############################################################################################
    required_cols = ["SNPID", "CHR", "POS", "EA", "NEA"]
    sumstats = _filter_and_prepare_sumstats(
        sumstats=sumstats,
        required_cols=required_cols,
        log=log,
        verbose=verbose
    )
    
    ############################################################################################
    # STEP 2: Extract lead variants
    ############################################################################################
    sig_df = _extract_lead_variants(
        sumstats=sumstats,
        loci=loci,
        loci_chrpos=loci_chrpos,
        locus=locus,
        suffixes=suffixes,
        getlead_kwargs=getlead_kwargs,
        exclude_hla=exclude_hla,
        log=log,
        verbose=verbose
    )
    
    if len(sig_df) == 0:
        log.write(" -No lead variants found.", verbose=verbose)
        output_file_list = pd.DataFrame(columns=["SNPID","SNPID_LIST","LD_R_MATRIX","LOCUS_SUMSTATS"])
        return None, output_file_list, ""
    
    log.write(" -Processing {} lead variant(s)...".format(len(sig_df)), verbose=verbose)
    
    # Initialize output file list
    output_file_list = pd.DataFrame(columns=["SNPID","SNPID_LIST","LD_R_MATRIX","LOCUS_SUMSTATS"])
    plink_log = ""
    
    ############################################################################################
    # STEP 3: Process each lead variant
    ############################################################################################
    for locus_idx, (_, row) in enumerate(sig_df.iterrows()):
        log.write(" -Processing lead variant {} of {}: {} (chr{}:{})...".format(
            locus_idx + 1, len(sig_df), row["SNPID"], row["CHR"], row["POS"]), verbose=verbose)
        
        # STEP 3.1: Extract variants in locus
        locus_sumstats = _extract_variants_in_locus(sumstats, windowsizekb, locus=(row["CHR"], row["POS"]))
        locus_sumstats = locus_sumstats.dropna(subset=required_cols).copy()
        log.write(" -Variants in locus after NA filtering: {}".format(len(locus_sumstats)), verbose=verbose)
        
        if len(locus_sumstats) == 0:
            log.write(" -No variants in locus for {} (chr{}:{}). Skipping...".format(
                row["SNPID"], row["CHR"], row["POS"]), verbose=verbose)
            continue
        
        # STEP 3.2: Get LD maps for all studies
        all_ld_maps = _get_all_ld_maps(
            n_studies=n_studies,
            studies=studies,
            ld_source_types=ld_source_types,
            row=row,
            locus_sumstats=locus_sumstats,
            ld_maps_list=ld_maps_list,
            ld_map_dics_list=ld_map_dics_list,
            ld_map_kwargs=ld_map_kwargs,
            bfiles=bfiles,
            pfiles=pfiles,
            vcfs=vcfs,
            windowsizekb=windowsizekb,
            plink=plink,
            plink2=plink2,
            threads=threads,
            memory=memory,
            overwrite=overwrite,
            tabix_list=tabix_list,
            log=log,
            verbose=verbose,
            **kwargs
        )
        
        # STEP 3.3: Find common variants across all studies
        common_variants = _find_common_variants_across_studies(
            n_studies=n_studies,
            studies=studies,
            locus_sumstats=locus_sumstats,
            all_ld_maps=all_ld_maps,
            suffixes=suffixes,
            row=row,
            log=log,
            verbose=verbose
        )
        
        if len(common_variants) == 0:
            log.write(" -No common variants found across all studies for locus {} (chr{}:{}). Skipping to next locus...".format(
                row["SNPID"], row["CHR"], row["POS"]), verbose=verbose)
            continue
        
        log.write(" -Found {} common variants across all studies".format(len(common_variants)), verbose=verbose)
        
        # STEP 3.4: Process LD matrices for all studies
        all_ld_matrices, study_plink_log, locus_file_list = _process_ld_matrices_for_all_studies(
            n_studies=n_studies,
            studies=studies,
            ld_source_types=ld_source_types,
            common_variants=common_variants,
            row=row,
            ld_paths_list=ld_paths_list,
            ld_maps_list=ld_maps_list,
            ld_types_list=ld_types_list,
            ld_if_square=ld_if_square,
            ld_if_add_T=ld_if_add_T,
            bfiles=bfiles,
            pfiles=pfiles,
            vcfs=vcfs,
            windowsizekb=windowsizekb,
            out=out,
            group=group,
            plink=plink,
            plink2=plink2,
            threads=threads,
            mode=mode,
            memory=memory,
            overwrite=overwrite,
            extra_plink_option=extra_plink_option,
            tabix_list=tabix_list,
            suffixes=suffixes,
            log=log,
            verbose=verbose,
            **kwargs
        )
        
        plink_log += study_plink_log
        
        # Export common SNP list (same for all studies)
        matched_snp_list_path = "{}/{}_{}_{}.snplist.raw".format(out.rstrip("/"), group, row["SNPID"], windowsizekb)
        common_variants["SNPID"].to_csv(matched_snp_list_path, index=None, header=None)
        log.write(" -Exported common SNP list ({} variants) to: {}".format(len(common_variants), matched_snp_list_path), verbose=verbose)
        
        # Update SNPID_LIST in all rows for this locus to point to the common list
        locus_rows = locus_file_list["SNPID"] == row["SNPID"]
        if locus_rows.any():
            locus_file_list.loc[locus_rows, "SNPID_LIST"] = matched_snp_list_path
        
        # Add locus file list to main output
        output_file_list = pd.concat([output_file_list, locus_file_list], ignore_index=True)
    
    # After processing all loci, finalize output
    if len(output_file_list)>0:
        output_file_list["GROUP"] = group
        nloci = len(output_file_list) // n_studies  # Number of unique loci
        # Use the last processed locus SNPID for filename (or could use a generic name)
        if len(sig_df) > 0:
            last_snpid = sig_df.iloc[-1]["SNPID"]
            output_file_list_path =  "{}/{}_{}study_{}_{}kb.filelist".format(out.rstrip("/"), group, n_studies, last_snpid, windowsizekb)
        else:
            output_file_list_path =  "{}/{}_{}study_{}kb.filelist".format(out.rstrip("/"), group, n_studies, windowsizekb)
        output_file_list.to_csv(output_file_list_path,index=None,sep="\t")
        log.write(" -File list is saved to: {}".format(output_file_list_path),verbose=verbose)
        log.write(" -Finished LD matrix calculation for {} loci.".format(nloci),verbose=verbose)
    else:
        output_file_list_path=None
        log.write(" -No available lead variants.",verbose=verbose)
        log.write(" -Stopped LD matrix calculation.",verbose=verbose)

    return output_file_list_path, output_file_list,  plink_log


###################################################################################################################################################################
####################################################################################################
def _extract_variants_from_ld_matrix_m(
    merged_sumstats: pd.DataFrame, 
    r_matrix: np.ndarray, 
    out: str, 
    group: str, 
    row: pd.Series, 
    windowsizekb: int, 
    log: Log, 
    verbose: bool, 
    index: int
) -> str:
    """
    Extract sub-matrix from full LD matrix for matched variants.
    
    This function ensures the extracted LD matrix:
    1. Contains only variants in merged_sumstats (using _INDEX_BIM_{index+1})
    2. Has the same order as merged_sumstats
    3. Has allele flips corrected (negate rows/cols for flipped variants)
    
    The extracted matrix will have shape (len(merged_sumstats), len(merged_sumstats))
    and matches the variant order in merged_sumstats exactly.
    """
    # study suffixes starting from 1
    index_bim_header = "_INDEX_BIM_{}".format(index + 1) 
    flipped_header = "_FLIPPED_{}".format(index + 1) 
    
    # Get indices: these point to rows/columns in the full LD matrix
    # Each value in avaiable_index corresponds to a variant in merged_sumstats
    avaiable_index = merged_sumstats[index_bim_header].values 
    
    # Validate indices before extraction
    if len(avaiable_index) == 0:
        raise ValueError("No valid indices found for study {}. Cannot extract LD matrix.".format(index + 1))
    
    # Check if indices are within LD matrix bounds
    max_index = avaiable_index.max() if len(avaiable_index) > 0 else -1
    if max_index >= r_matrix.shape[0] or max_index < 0:
        raise ValueError("Invalid indices for study {}: max_index={}, matrix_shape={}".format(
            index + 1, max_index, r_matrix.shape))
    
    if r_matrix.size == 0:
        raise ValueError("LD matrix is empty for study {}. Cannot extract sub-matrix.".format(index + 1))
    
    # Get flip flags: True for variants where alleles are flipped
    flipped = merged_sumstats[flipped_header].values 
    
    # Extract sub-matrix using fancy indexing
    # np.ix_ creates index arrays for both dimensions, extracting the sub-matrix
    # The result has shape (len(merged_sumstats), len(merged_sumstats))
    # and the order matches merged_sumstats
    reduced_r_matrix = r_matrix[np.ix_(avaiable_index, avaiable_index)]
    
    # Correct for allele flips: if a variant's alleles are flipped, we negate
    # its row and column in the LD matrix. This is because LD is calculated
    # based on allele coding, and flipping alleles changes the sign of correlations.
    n_flipped = sum(flipped)
    if n_flipped > 0:
        log.write(" -Flipping LD matrix for {} variants...".format(n_flipped), verbose=verbose)
    reduced_r_matrix[flipped,:] = -1 * reduced_r_matrix[flipped,:]
    reduced_r_matrix[:,flipped] = -1 * reduced_r_matrix[:,flipped]

    output_prefix =  "{}/{}_{}_{}_{}".format(out.rstrip("/"),group,row["SNPID"],windowsizekb, index + 1)
    output_path = "{}.ld.gz".format(output_prefix)
    
    pd.DataFrame(reduced_r_matrix).to_csv(output_path,sep="\t",index=None,header=None)
    #reduced_r_matrix.to_csv("{}.ld.gz".format(output_prefix),se="\t")
    return output_path

