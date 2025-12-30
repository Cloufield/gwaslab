"""
Utilities for extracting and processing GWAS Catalog associations.

This module provides functions to retrieve associations from the GWAS Catalog API v2
and align them with summary statistics data.
"""

import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Tuple, Optional, Any, Union, List, Dict

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
from gwaslab.info.g_Log import Log


# ============================================================================
# Constants
# ============================================================================

# Mapping from API v2 columns to GWASLab standard columns with _GCV2 suffix
# This maps API response fields to GWASLab reserved headers
API_TO_GWASLAB_MAPPING = {
    # Basic identifiers
    'variant_rsID': 'rsID_GCV2',
    'association_id': 'associationId_GCV2',
    'accession_id': 'accession_id_GCV2',
    
    # Genomic coordinates (will be extracted from locations)
    # 'CHR_GCV2': extracted from locations
    # 'POS_GCV2': extracted from locations
    
    # Alleles
    'risk_allele_base': 'EA_GCV2',  # Effect allele from GWAS Catalog
    'snp_effect_allele': 'snp_effect_allele_GCV2',  # Keep original for reference
    
    # Statistics
    'beta': 'BETA_GCV2',  # Will be converted to numeric
    'p_value': 'P_GCV2',
    'risk_frequency': 'EAF_GCV2',  # Risk allele frequency
    'ci_lower': 'BETA_95L_GCV2',
    'ci_upper': 'BETA_95U_GCV2',
    'range': 'range_GCV2',
    
    # Traits and genes
    'efo_traits': 'efo_traits_GCV2',  # Will be processed
    'reported_trait': 'reported_trait_GCV2',
    'mapped_genes': 'GENENAME_GCV2',  # Will be joined
    'locations': 'locations_GCV2',  # Will be processed
    
    # Study information
    'pubmed_id': 'KNOWN_PUBMED_ID_GCV2',
    'first_author': 'KNOWN_AUTHOR_GCV2',
    'pvalue_description': 'pvalue_description_GCV2',
    'bg_efo_traits': 'bg_efo_traits_GCV2',
}

# Default columns to include (all available by default)
DEFAULT_GCV2_COLUMNS = [
    'rsID_GCV2', 'CHR_GCV2', 'POS_GCV2', 'EA_GCV2', 'BETA_GCV2', 'P_GCV2',
    'EAF_GCV2', 'BETA_95L_GCV2', 'BETA_95U_GCV2', 'associationId_GCV2',
    'efo_traits_GCV2', 'reported_trait_GCV2', 'GENENAME_GCV2',
    'KNOWN_PUBMED_ID_GCV2', 'KNOWN_AUTHOR_GCV2', 'accession_id_GCV2'
]

# Legacy column mappings (for backward compatibility)
COLUMN_MAPPING = {
    'association_id': 'associationId',
    'p_value': 'pvalue',
    'beta': 'betaNum',
    'or_per_copy_num': 'orPerCopyNum',
    'risk_frequency': 'riskFrequency',
    'standard_error': 'standardError',
    'functional_class': 'functionalClass',
    'accession_id': 'accession_id',
}

# Summary columns for final output (legacy)
SUMMARY_COLUMNS = [
    'GWASCATALOG_TRAIT', 'associationId', 'rsID', 'geneName',
    'RA', 'RAF', 'Beta', 'Unit', 'P-value', 'cohort',
    'initialSampleSize', 'publicationInfo.pubmedId',
    'functionalClass', 'gene.geneName'
]


# ============================================================================
# Helper Functions
# ============================================================================

def _extract_scalar_value(value: Any) -> Optional[Any]:
    """
    Extract scalar value from pandas Series, list, or array.
    
    Parameters
    ----------
    value : Any
        Value that might be a Series, list, array, or scalar
        
    Returns
    -------
    Optional[Any]
        Scalar value or None if extraction fails
    """
    if value is None:
        return None
    
    try:
        # Handle pandas Series
        if hasattr(value, 'iloc'):
            return value.iloc[0] if len(value) > 0 else None
        # Handle lists/tuples
        elif isinstance(value, (list, tuple)):
            return value[0] if len(value) > 0 else None
        # Handle numpy arrays
        elif hasattr(value, 'item'):
            return value.item()
        # Already a scalar
        else:
            return value
    except (TypeError, ValueError, AttributeError, IndexError):
        return None


def _is_valid_value(value: Any) -> bool:
    """
    Check if a value is valid (not None, not NaN, not empty string).
    
    Parameters
    ----------
    value : Any
        Value to check
        
    Returns
    -------
    bool
        True if value is valid, False otherwise
    """
    if value is None:
        return False
    
    # Check for pandas NA
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    
    # Check for string 'nan'
    if isinstance(value, str) and value.lower() == 'nan':
        return False
    
    # Check for empty string
    if isinstance(value, str) and not value.strip():
        return False
    
    return True


def _extract_risk_allele(snp_effect_allele: Any) -> Optional[str]:
    """
    Extract risk allele base from snp_effect_allele field.
    
    Parameters
    ----------
    snp_effect_allele : Any
        snp_effect_allele value (usually a list like ['rs123-A'])
        
    Returns
    -------
    Optional[str]
        Risk allele base (e.g., 'A') or None
    """
    if not isinstance(snp_effect_allele, list) or len(snp_effect_allele) == 0:
        return None
    
    risk_allele_str = snp_effect_allele[0]
    if not isinstance(risk_allele_str, str):
        return None
    
    # Split on '-' and get the last part (allele base)
    parts = risk_allele_str.split('-')
    return parts[-1] if len(parts) > 1 else None


def _extract_genes(mapped_genes: Any) -> Optional[str]:
    """
    Extract and join gene names from mapped_genes field.
    
    Parameters
    ----------
    mapped_genes : Any
        mapped_genes value (list of strings, some may be comma-separated)
        
    Returns
    -------
    Optional[str]
        Comma-separated gene names or None
    """
    if not isinstance(mapped_genes, list):
        return None
    
    all_genes = []
    for gene_str in mapped_genes:
        if isinstance(gene_str, str):
            genes = [g.strip() for g in gene_str.split(',') if g.strip()]
            all_genes.extend(genes)
    
    return ','.join(all_genes) if all_genes else None


def _extract_chr_pos_from_locations(locations: Any) -> Tuple[Optional[str], Optional[int]]:
    """
    Extract chromosome and position from locations field.
    
    Parameters
    ----------
    locations : Any
        Locations value (usually a list like ['12:111803962'])
        
    Returns
    -------
    Tuple[Optional[str], Optional[int]]
        (chromosome, position) or (None, None) if extraction fails
    """
    if not isinstance(locations, list) or len(locations) == 0:
        return None, None
    
    location_str = locations[0]
    if not isinstance(location_str, str):
        return None, None
    
    # Format: "12:111803962" or "chr12:111803962"
    parts = location_str.replace('chr', '').split(':')
    if len(parts) == 2:
        try:
            chr_val = parts[0]
            pos_val = int(parts[1])
            return chr_val, pos_val
        except (ValueError, IndexError):
            pass
    
    return None, None


def _extract_trait_names(efo_traits: Any) -> str:
    """
    Extract trait names from efo_traits field.
    
    Parameters
    ----------
    efo_traits : Any
        efo_traits value (list of dicts)
        
    Returns
    -------
    str
        Comma-separated trait names
    """
    if not isinstance(efo_traits, list):
        return ""
    
    trait_names = []
    for trait_info in efo_traits:
        if isinstance(trait_info, dict):
            trait_name = trait_info.get('efo_trait', trait_info.get('trait', ''))
            if trait_name:
                trait_names.append(trait_name)
    
    return ', '.join(trait_names) if trait_names else ""


def _transform_to_gcv2_format(df: pd.DataFrame, 
                             columns: Optional[list] = None) -> pd.DataFrame:
    """
    Transform API response to GWASLab format with _GCV2 suffix.
    
    Parameters
    ----------
    df : pd.DataFrame
        Association DataFrame from API
    columns : list, optional
        List of _GCV2 columns to include. If None, includes all available.
        
    Returns
    -------
    pd.DataFrame
        Transformed DataFrame with _GCV2 columns
    """
    if len(df) == 0:
        return df
    
    result_df = pd.DataFrame()
    
    # Map basic columns
    for api_col, gwaslab_col in API_TO_GWASLAB_MAPPING.items():
        if api_col in df.columns:
            result_df[gwaslab_col] = df[api_col]
    
    # Extract CHR and POS from locations
    if 'locations' in df.columns:
        chr_pos = df['locations'].apply(_extract_chr_pos_from_locations)
        result_df['CHR_GCV2'] = chr_pos.apply(lambda x: x[0] if x[0] else None)
        result_df['POS_GCV2'] = chr_pos.apply(lambda x: x[1] if x[1] else None)
    
    # Extract trait names from efo_traits
    if 'efo_traits' in df.columns:
        result_df['efo_traits_GCV2'] = df['efo_traits'].apply(_extract_trait_names)
    
    # Extract reported traits (join list)
    if 'reported_trait' in df.columns:
        result_df['reported_trait_GCV2'] = df['reported_trait'].apply(
            lambda x: ', '.join(x) if isinstance(x, list) else str(x) if x else ""
        )
    
    # Extract and join mapped genes
    if 'mapped_genes' in df.columns:
        result_df['GENENAME_GCV2'] = df['mapped_genes'].apply(
            lambda x: ', '.join(x) if isinstance(x, list) else str(x) if x else ""
        )
    
    # Convert BETA_GCV2 to numeric (extract from string if needed)
    if 'BETA_GCV2' in result_df.columns:
        result_df['BETA_GCV2'] = result_df['BETA_GCV2'].apply(_extract_numeric_from_string)
        result_df['BETA_GCV2'] = pd.to_numeric(result_df['BETA_GCV2'], errors='coerce')
    
    # Convert EAF_GCV2 to numeric (handle 'NR' and other non-numeric values)
    if 'EAF_GCV2' in result_df.columns:
        result_df['EAF_GCV2'] = pd.to_numeric(result_df['EAF_GCV2'], errors='coerce')
    
    # Select columns if specified
    if columns is not None:
        available_cols = [col for col in columns if col in result_df.columns]
        if available_cols:
            result_df = result_df[available_cols]
        else:
            # If none of the requested columns exist, return empty with those columns
            result_df = pd.DataFrame(columns=columns)
    
    return result_df


def _normalize_association_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize column names in association DataFrame.
    
    Parameters
    ----------
    df : pd.DataFrame
        Association DataFrame with API v2 column names
        
    Returns
    -------
    pd.DataFrame
        DataFrame with normalized column names
    """
    df = df.copy()
    
    # Apply column mapping
    for old_col, new_col in COLUMN_MAPPING.items():
        if old_col in df.columns:
            df = df.rename(columns={old_col: new_col})
    
    # Ensure associationId column exists
    if 'associationId' not in df.columns:
        if 'association_id' in df.columns:
            df['associationId'] = df['association_id']
        else:
            # Try to find any column with 'association' in the name
            assoc_cols = [col for col in df.columns if 'association' in col.lower()]
            if assoc_cols:
                df['associationId'] = df[assoc_cols[0]]
    
    # Ensure betaNum is numeric (handle cases where API returns strings)
    # Check if betaNum exists, if not, try 'beta' field
    beta_col = None
    if 'betaNum' in df.columns:
        beta_col = 'betaNum'
    elif 'beta' in df.columns:
        beta_col = 'beta'
    
    if beta_col:
        # Convert to numeric, handling string values
        df[beta_col] = pd.to_numeric(df[beta_col], errors='coerce')
        
        # If still has non-numeric values (object dtype), try extracting from strings
        if df[beta_col].dtype == 'object':
            df[beta_col] = df[beta_col].apply(_extract_numeric_from_string)
            df[beta_col] = pd.to_numeric(df[beta_col], errors='coerce')
        
        # Ensure we have betaNum column (rename if needed)
        if beta_col != 'betaNum':
            df['betaNum'] = df[beta_col]
    
    return df


# ============================================================================
# Main Functions
# ============================================================================

def _extract_associations(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame], 
    rsid: str = "rsID", 
                         log: Log = Log(), verbose: bool = True,
                         fetch_metadata: bool = False,
                         fetch_traits: Optional[bool] = None,
                         fetch_studies: Optional[bool] = None,
                         fetch_variants: Optional[bool] = None,
                         gcv2_columns: Optional[list] = None,
                         use_gcv2_format: bool = True) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Extract and process GWAS Catalog associations for variants in sumstats.
    
    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame with rsID column. Limited to 100 unique variants.
        If more than 100 unique variants are provided, only the first 100 will be processed.
    rsid : str, optional
        Name of the rsID column (default: "rsID")
    log : Log, optional
        Logging object
    verbose : bool, optional
        Whether to print log messages
    fetch_metadata : bool, optional
        If True, fetch additional metadata (traits, studies, variants).
        If False, only fetch associations (faster, fewer API calls).
        Default: True
    fetch_traits : bool, optional
        If True, fetch traits. If False, skip traits.
        If None, uses fetch_metadata value. Default: None
    fetch_studies : bool, optional
        If True, fetch studies. If False, skip studies.
        If None, uses fetch_metadata value. Default: None
    fetch_variants : bool, optional
        If True, fetch variants. If False, skip variants.
        If None, uses fetch_metadata value. Default: None
        
    Returns
    -------
    Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]
        (associations_full, associations_summary) or (None, None) if no associations
        
    Note
    ----
    The input is limited to 100 unique variants to prevent excessive API calls.
    If more variants are provided, only the first 100 will be processed and a warning will be issued.
    """
    # Set individual flags based on fetch_metadata if not explicitly set
    if fetch_traits is None:
        fetch_traits = fetch_metadata
    if fetch_studies is None:
        fetch_studies = fetch_metadata
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    if fetch_variants is None:
        fetch_variants = fetch_metadata
    
    assoc, traits, studies, variants = get_associations_from_gwascatalog(
        sumstats, rsid=rsid, log=log, verbose=verbose,
        fetch_traits=fetch_traits, fetch_studies=fetch_studies, fetch_variants=fetch_variants
    )
    
    if len(assoc) == 0:
        log.write("No associations!", verbose=verbose)
        return None, None
    
    # Transform to GCV2 format if requested
    if use_gcv2_format:
        log.write("Transforming associations to GWASLab GCV2 format...", verbose=verbose)
        
        # Get raw associations (before normalization) for proper GCV2 transformation
        # Fetch raw associations without normalization
        assoc_raw, _, _, _ = get_associations_from_gwascatalog(
            sumstats, rsid=rsid, log=log, verbose=verbose,
            fetch_traits=False, fetch_studies=False, fetch_variants=False,
            normalize_columns=False  # Skip normalization to get raw API format
        )
        
        if len(assoc_raw) > 0:
            assoc_gcv2 = _transform_to_gcv2_format(assoc_raw, columns=gcv2_columns)
            
            # Merge with sumstats on rsID
            if rsid in sumstats.columns and 'rsID_GCV2' in assoc_gcv2.columns:
                # Group by rsID to handle multiple associations per variant (take first)
                assoc_grouped = assoc_gcv2.groupby('rsID_GCV2').first().reset_index()
                
                # Merge GCV2 columns with sumstats
                merged = pd.merge(
                    sumstats,
                    assoc_grouped,
                    left_on=rsid,
                    right_on='rsID_GCV2',
                    how='left',
                    suffixes=('', '_dup')
                )
                # Remove duplicate columns if any
                merged = merged.loc[:, ~merged.columns.str.endswith('_dup')]
                log.write(f"Merged {len(merged)} variants with GWAS Catalog associations", verbose=verbose)
                return merged, merged  # Return same for both full and summary
            else:
                log.warning("Cannot merge: rsID column mismatch", verbose=verbose)
                return assoc_gcv2, assoc_gcv2
        else:
            log.write("No associations to transform", verbose=verbose)
            return None, None
    
    # Legacy format processing
    # Fix beta values (convert OR to beta if needed)
    assoc = _fix_beta(assoc)

    # Aggregate traits by associationId (if fetched)
    if len(traits) > 0:
        traits_agg = traits.groupby("associationId")[["trait", "shortForm"]].agg(
            lambda x: ",".join(x)
        ).reset_index()
        assoc = pd.merge(assoc, traits_agg, on="associationId", how="left")
    else:
        # Add empty columns if traits were not fetched
        if "trait" not in assoc.columns:
            assoc["trait"] = ""
        if "shortForm" not in assoc.columns:
            assoc["shortForm"] = ""
    
    # Merge with studies (if fetched)
    if len(studies) > 0:
        assoc = pd.merge(assoc, studies, on="associationId", how="left")
    
    # Merge with variants (if fetched)
    if len(variants) > 0:
        assoc = pd.merge(assoc, variants, on="associationId", how="left")
    
    # Rename columns for consistency
    # Only rename betaNum to Beta (not betaUnit, which contains strings)
    rename_dict = {}
    if "trait" in assoc.columns:
        rename_dict["trait"] = "GWASCATALOG_TRAIT"
    if "riskFrequency" in assoc.columns:
        rename_dict["riskFrequency"] = "RAF"
    # Only rename betaNum if it exists (don't use betaUnit which has strings)
    if "betaNum" in assoc.columns:
        rename_dict["betaNum"] = "Beta"
    elif "beta" in assoc.columns:
        rename_dict["beta"] = "Beta"
    if "pvalue" in assoc.columns:
        rename_dict["pvalue"] = "P-value"
    if "betaUnit" in assoc.columns:
        rename_dict["betaUnit"] = "Unit"
    
    if rename_dict:
        assoc = assoc.rename(columns=rename_dict)
    
    # Ensure Beta column exists and is numeric (handle cases where it might be a string with units)
    if "Beta" not in assoc.columns:
        # Try to find numeric beta column (not betaUnit which is a string)
        beta_cols = [col for col in assoc.columns 
                    if 'beta' in col.lower() and 'unit' not in col.lower() and 'direction' not in col.lower()]
        if beta_cols:
            assoc["Beta"] = assoc[beta_cols[0]]
        else:
            # Create empty Beta column if not found
            assoc["Beta"] = np.nan
    
    # Ensure Beta column is numeric (extract numbers from strings if needed)
    assoc = _ensure_numeric_beta(assoc, beta_col="Beta")
    
    # Create summary DataFrame
    available_cols = [col for col in SUMMARY_COLUMNS if col in assoc.columns]
    assoc_summary = assoc[available_cols].copy()
    
    # Align beta values with sumstats
    assoc_summary = _align_beta(sumstats, assoc_summary, log, verbose)
    assoc = _align_beta(sumstats, assoc, log, verbose)
    
    return assoc, assoc_summary


def get_associations_from_gwascatalog(sumstats: pd.DataFrame, rsid: str = "rsID",
                                     log: Log = Log(), verbose: bool = True,
                                     fetch_metadata: bool = True,
                                     fetch_traits: Optional[bool] = None,
                                     fetch_studies: Optional[bool] = None,
                                     fetch_variants: Optional[bool] = None,
                                     normalize_columns: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Get associations from GWAS Catalog API v2 for variants in sumstats.
    
    This function retrieves associations from the GWAS Catalog API v2 for each unique
    variant in the sumstats DataFrame. Optionally, it can also fetch additional metadata
    (traits, studies, variants) for each association.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame with rsID column. Limited to 100 unique variants.
        If more than 100 unique variants are provided, only the first 100 will be processed.
    rsid : str, optional
        Name of the rsID column (default: "rsID")
    log : Log, optional
        Logging object
    verbose : bool, optional
        Whether to print log messages
    fetch_metadata : bool, optional
        If True, fetch additional metadata (traits, studies, variants) for each association.
        If False, only fetch associations (faster, fewer API calls).
        Default: True
    fetch_traits : bool, optional
        If True, fetch traits. If False, skip traits.
        If None, uses fetch_metadata value. Default: None
    fetch_studies : bool, optional
        If True, fetch studies. If False, skip studies.
        If None, uses fetch_metadata value. Default: None
    fetch_variants : bool, optional
        If True, fetch variants. If False, skip variants.
        If None, uses fetch_metadata value. Default: None
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]
        (associations, traits, studies, variants) DataFrames.
        If metadata fetching is disabled, corresponding DataFrames will be empty.
        
    Note
    ----
    The input is limited to 100 unique variants to prevent excessive API calls.
    If more variants are provided, only the first 100 will be processed and a warning will be issued.
        
    Examples
    --------
    >>> # Fetch only associations (fastest)
    >>> assoc, traits, studies, variants = get_associations_from_gwascatalog(
    ...     sumstats, fetch_metadata=False
    ... )
    >>> 
    >>> # Fetch associations and traits only
    >>> assoc, traits, studies, variants = get_associations_from_gwascatalog(
    ...     sumstats, fetch_traits=True, fetch_studies=False, fetch_variants=False
    ... )
    """
    from gwaslab.extension.gwascatalog import GWASCatalogClient
    
    # Set individual flags based on fetch_metadata if not explicitly set
    if fetch_traits is None:
        fetch_traits = fetch_metadata
    if fetch_studies is None:
        fetch_studies = fetch_metadata
    if fetch_variants is None:
        fetch_variants = fetch_metadata
    
    client = GWASCatalogClient(verbose=verbose, log=log)
    
    # Get unique variants
    unique_sumstats = sumstats.dropna(subset=[rsid]).drop_duplicates(subset=[rsid])
    
    if len(unique_sumstats) == 0:
        log.write("No valid rsIDs found in sumstats", verbose=verbose)
        return (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    
    # Limit to 100 variants
    MAX_VARIANTS = 100
    if len(unique_sumstats) > MAX_VARIANTS:
        log.warning(f"Input contains {len(unique_sumstats)} unique variants. Limiting to first {MAX_VARIANTS} variants for API query.", verbose=verbose)
        unique_sumstats = unique_sumstats.head(MAX_VARIANTS)
    
    # Step 1: Get associations for each variant
    association = _fetch_associations(client, unique_sumstats, rsid, log, verbose, normalize_columns=normalize_columns)
    
    if len(association) == 0:
        log.write("No associations found", verbose=verbose)
        return (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    
    # Step 2: Get traits, studies, and variants for each unique association (if requested)
    if fetch_traits or fetch_studies or fetch_variants:
        traits, studies, variants = _fetch_metadata(
            client, association, rsid, log, verbose,
            fetch_traits=fetch_traits, fetch_studies=fetch_studies, fetch_variants=fetch_variants
        )
    else:
        log.write("Skipping metadata fetch (traits, studies, variants) as requested", verbose=verbose)
        traits = pd.DataFrame()
        studies = pd.DataFrame()
        variants = pd.DataFrame()
    
    return association, traits, studies, variants


def _fetch_associations(
    client: Any, 
    unique_sumstats: pd.DataFrame, 
    rsid: str,
    log: Log, 
    verbose: bool, 
    normalize_columns: bool = True
) -> pd.DataFrame:
    """
    Fetch associations from GWAS Catalog API for each variant.
    
    Parameters
    ----------
    client : GWASCatalogClient
        GWAS Catalog API client
    unique_sumstats : pd.DataFrame
        DataFrame with unique variants
    rsid : str
        Name of the rsID column
    log : Log
        Logging object
    verbose : bool
        Whether to print log messages
        
    Returns
    -------
    pd.DataFrame
        DataFrame with all associations
    """
    association = pd.DataFrame()
    empty_variants = []
    
    for index, row in unique_sumstats.iterrows():
        rs_id = row[rsid]
        log.write(f"Getting associations from GWAS Catalog for {rs_id}...", verbose=verbose)
        
        try:
            associations_df = client.get_associations(rs_id=rs_id, get_all=True)
            
            if not isinstance(associations_df, pd.DataFrame) or len(associations_df) == 0:
                empty_variants.append(rs_id)
                continue
            
            # Add rsID column
            associations_df[rsid] = rs_id
            
            # Extract risk allele
            if 'snp_effect_allele' in associations_df.columns:
                associations_df['RA'] = associations_df['snp_effect_allele'].apply(_extract_risk_allele)
            
            # Normalize column names (skip if normalize_columns=False for GCV2 format)
            if normalize_columns:
                associations_df = _normalize_association_columns(associations_df)
            
            # Extract gene names (only if normalizing)
            if normalize_columns and 'mapped_genes' in associations_df.columns:
                associations_df['geneName'] = associations_df['mapped_genes'].apply(_extract_genes)
            
            association = pd.concat([association, associations_df], ignore_index=True)
            log.write("", show_time=False, verbose=verbose)
            
        except Exception as e:
            log.warning(f"Error fetching associations for {rs_id}: {str(e)}", verbose=verbose)
            empty_variants.append(rs_id)
    
    if empty_variants:
        log.write(f"No associations found for {len(empty_variants)} variants: {empty_variants[:10]}...", verbose=verbose)
    
    log.write(f"Retrieved {len(association)} associations from GWAS Catalog", verbose=verbose)
    return association


def _fetch_metadata(
    client: Any, 
    association: pd.DataFrame, 
    rsid: str,
    log: Log, 
    verbose: bool,
    fetch_traits: bool = True,
    fetch_studies: bool = True,
    fetch_variants: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Fetch traits, studies, and variants metadata for each unique association.
    
    Parameters
    ----------
    client : GWASCatalogClient
        GWAS Catalog API client
    association : pd.DataFrame
        DataFrame with associations
    rsid : str
        Name of the rsID column
    log : Log
        Logging object
    verbose : bool
        Whether to print log messages
    fetch_traits : bool, optional
        If True, fetch traits. Default: True
    fetch_studies : bool, optional
        If True, fetch studies. Default: True
    fetch_variants : bool, optional
        If True, fetch variants. Default: True
        
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        (traits, studies, variants) DataFrames
    """
    traits = pd.DataFrame()
    studies = pd.DataFrame()
    variants = pd.DataFrame()

    if len(association) == 0:
        return traits, studies, variants
    
    unique_associations = association.drop_duplicates(subset=["associationId"])
    
    # Build list of what to fetch for logging
    fetch_list = []
    if fetch_traits:
        fetch_list.append("traits")
    if fetch_studies:
        fetch_list.append("studies")
    if fetch_variants:
        fetch_list.append("variants")
    
    if fetch_list:
        log.write(f'Fetching metadata: {", ".join(fetch_list)}', verbose=verbose)
    
    for index, row in unique_associations.iterrows():
        association_id = row["associationId"]
        
        if fetch_traits or fetch_studies or fetch_variants:
            log.write(f'Getting metadata for associationId: {association_id}...', verbose=verbose)
        
        # Extract traits (if requested)
        if fetch_traits:
            trait_rows = _extract_traits(row, association_id)
            if trait_rows:
                traits = pd.concat([traits, pd.DataFrame(trait_rows)], ignore_index=True)
        
        # Extract study information (if requested)
        if fetch_studies:
            study_row = _extract_study(client, row, association_id, log, verbose)
            if study_row:
                studies = pd.concat([studies, pd.DataFrame([study_row])], ignore_index=True)
        
        # Extract variant information (if requested)
        if fetch_variants:
            variant_row = _extract_variant(client, row, association_id, rsid, log, verbose)
            if variant_row:
                variants = pd.concat([variants, pd.DataFrame([variant_row])], ignore_index=True)
    
    return traits, studies, variants


def _extract_traits(row: pd.Series, association_id: str) -> list:
    """
    Extract trait information from association row.
    
    Parameters
    ----------
    row : pd.Series
        Association row
    association_id : str
        Association ID
        
    Returns
    -------
    list
        List of trait dictionaries
    """
    efo_traits = _extract_scalar_value(row.get('efo_traits'))
    
    if not isinstance(efo_traits, list):
        return []
    
    trait_rows = []
    for trait_info in efo_traits:
        if isinstance(trait_info, dict):
            trait_row = {
                'associationId': association_id,
                'trait': trait_info.get('trait', trait_info.get('efo_trait', '')),
                'shortForm': trait_info.get('shortForm', trait_info.get('short_form', '')),
            }
            trait_rows.append(trait_row)
    
    return trait_rows


def _extract_study(client, row: pd.Series, association_id: str,
                   log: Log, verbose: bool) -> Optional[dict]:
    """
    Extract study information for an association.
    
    Parameters
    ----------
    client : GWASCatalogClient
        GWAS Catalog API client
    row : pd.Series
        Association row
    association_id : str
        Association ID
    log : Log
        Logging object
    verbose : bool
        Whether to print log messages
        
    Returns
    -------
    Optional[dict]
        Study row dictionary or None
    """
    study_id = _extract_scalar_value(row.get('study_id') or row.get('accession_id'))
    
    if not _is_valid_value(study_id):
        return None
    
    try:
        study_data = client.get_studies(accession_id=str(study_id))
        
        if not isinstance(study_data, dict) or not study_data:
            return None
        
        # Handle cohort as list
        cohort = study_data.get('cohort', '')
        if isinstance(cohort, list):
            cohort = ', '.join(cohort) if cohort else ''
        
        return {
            'associationId': association_id,
            'accession_id': study_data.get('accession_id', study_id),
            'pubmed_id': study_data.get('pubmed_id', ''),
            'initialSampleSize': study_data.get('initial_sample_size', ''),
            'cohort': cohort,
        }
    except Exception as e:
        log.warning(f"Error fetching study for {study_id}: {str(e)}", verbose=verbose)
        return None


def _extract_variant(client, row: pd.Series, association_id: str, rsid: str,
                    log: Log, verbose: bool) -> Optional[dict]:
    """
    Extract variant information for an association.
    
    Parameters
    ----------
    client : GWASCatalogClient
        GWAS Catalog API client
    row : pd.Series
        Association row
    association_id : str
        Association ID
    rsid : str
        Name of the rsID column
    log : Log
        Logging object
    verbose : bool
        Whether to print log messages
        
    Returns
    -------
    Optional[dict]
        Variant row dictionary or None
    """
    rs_id_for_variant = _extract_scalar_value(row.get('variant_rsID') or row.get(rsid))
    
    if not _is_valid_value(rs_id_for_variant):
        return None
    
    try:
        variant_data = client.get_variants(rs_id=str(rs_id_for_variant))
        
        if not isinstance(variant_data, dict) or not variant_data:
            return None
        
        return {
            'associationId': association_id,
            'functionalClass': variant_data.get('functional_class', ''),
            'gene.geneName': row.get('geneName', ''),
        }
    except Exception as e:
        log.warning(f"Error fetching variant for {rs_id_for_variant}: {str(e)}", verbose=verbose)
        return None


# ============================================================================
# Beta Processing Functions
# ============================================================================

def _fix_beta(association: pd.DataFrame) -> pd.DataFrame:
    """
    Fix beta values by converting OR to beta or using range if beta is missing.
    
    Parameters
    ----------
    association : pd.DataFrame
        Association DataFrame
        
    Returns
    -------
    pd.DataFrame
        Association DataFrame with fixed beta values
    """
    association = association.copy()
    
    # Initialize missing columns
    for col in ["betaNum", "orPerCopyNum", "range"]:
        if col not in association.columns:
            association[col] = pd.NA
    
    # Convert OR to beta where beta is missing
    is_or_available = association["betaNum"].isna() & association["orPerCopyNum"].notna()
    association.loc[is_or_available, "betaNum"] = np.log(association.loc[is_or_available, "orPerCopyNum"])
    
    # Use range midpoint where both beta and OR are missing
    is_range_available = (
        association["betaNum"].isna() & 
        association["orPerCopyNum"].isna() & 
        association["range"].notna()
    )
    association.loc[is_range_available, "betaNum"] = association.loc[is_range_available, "range"].apply(_parse_range)
    
    return association


def _parse_range(x: str) -> float:
    """
    Parse range string and return midpoint of log-transformed values.
    
    Parameters
    ----------
    x : str
        Range string like "[0.5-2.0]"
        
    Returns
    -------
    float
        Midpoint of log-transformed range
    """
    try:
        range_list = x.strip("[|]").split("-")
        if len(range_list) != 2:
            return np.nan
        
        high = np.log(float(range_list[1]))
        low = np.log(float(range_list[0]))
        return (high + low) / 2
    except (ValueError, IndexError, AttributeError):
        return np.nan


def _extract_numeric_from_string(value: Any) -> Optional[float]:
    """
    Extract numeric value from string that may contain units.
    
    Parameters
    ----------
    value : Any
        Value that might be a string like "0.0441029 unit increase" or a number
        
    Returns
    -------
    Optional[float]
        Numeric value or None if extraction fails
    """
    if pd.isna(value) or value is None:
        return None
    
    # If already numeric, return as is
    if isinstance(value, (int, float, np.number)):
        return float(value)
    
    # If string, try to extract number
    if isinstance(value, str):
        # Remove common unit words and extract first number
        import re
        # Match pattern: optional sign, digits, optional decimal point, digits
        match = re.search(r'([+-]?\d+\.?\d*)', value)
        if match:
            try:
                return float(match.group(1))
            except (ValueError, AttributeError):
                pass
    
    return None


def _ensure_numeric_beta(df: pd.DataFrame, beta_col: str = "Beta") -> pd.DataFrame:
    """
    Ensure Beta column is numeric, extracting numbers from strings if needed.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with Beta column
    beta_col : str, optional
        Name of the beta column (default: "Beta")
        
    Returns
    -------
    pd.DataFrame
        DataFrame with numeric Beta column
    """
    if beta_col not in df.columns:
        return df
    
    df = df.copy()
    
    # Check if column contains strings
    if df[beta_col].dtype == 'object':
        # Try to convert to numeric first
        df[beta_col] = pd.to_numeric(df[beta_col], errors='coerce')
        
        # For remaining non-numeric values, try to extract numbers
        mask = df[beta_col].isna()
        if mask.any():
            df.loc[mask, beta_col] = df.loc[mask, beta_col].apply(_extract_numeric_from_string)
            # Convert to numeric again
            df[beta_col] = pd.to_numeric(df[beta_col], errors='coerce')
    
    return df


def _align_beta(sumstats: pd.DataFrame, assoc_df: pd.DataFrame, 
               log: Log, verbose: bool) -> pd.DataFrame:
    """
    Align beta values between GWAS Catalog and sumstats.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics DataFrame
    assoc_df : pd.DataFrame
        Association DataFrame from GWAS Catalog
    log : Log
        Logging object
    verbose : bool
        Whether to print log messages
        
    Returns
    -------
    pd.DataFrame
        Association DataFrame with aligned beta values
    """
    # Merge with sumstats
    log.write("Merging GWAS Catalog associations with Sumstats...", verbose=verbose)
    
    merge_cols = ["rsID", "CHR", "POS", "EA", "NEA", "BETA"]
    available_cols = [col for col in merge_cols if col in sumstats.columns]
    
    if "rsID" not in available_cols:
        log.warning("rsID column not found in sumstats, cannot merge", verbose=verbose)
        return assoc_df
    
    assoc_df = pd.merge(
        sumstats[available_cols], 
        assoc_df, 
        on="rsID", 
        how="inner"
    )
    
    # Ensure Beta column is numeric before calculations
    assoc_df = _ensure_numeric_beta(assoc_df, beta_col="Beta")
    
    # Initialize aligned beta columns (ensure they're float type)
    assoc_df["GC_Beta_GL_RA"] = np.nan
    assoc_df["GC_Beta_GL_EA"] = np.nan
    assoc_df["GC_Beta_GL_RA"] = assoc_df["GC_Beta_GL_RA"].astype(float)
    assoc_df["GC_Beta_GL_EA"] = assoc_df["GC_Beta_GL_EA"].astype(float)
    
    # Determine risk allele from sumstats (allele with positive beta)
    assoc_df["_GL_RA_BETA"] = assoc_df["BETA"].abs()
    assoc_df["_GL_RA"] = assoc_df["EA"]
    assoc_df["_GL_NRA"] = assoc_df["NEA"]
    
    # Flip if beta is negative
    is_flipped = assoc_df["BETA"] < 0
    assoc_df.loc[is_flipped, "_GL_RA"] = assoc_df.loc[is_flipped, "NEA"]
    assoc_df.loc[is_flipped, "_GL_NRA"] = assoc_df.loc[is_flipped, "EA"]
    
    # Align with risk allele (RA)
    is_matched = assoc_df["_GL_RA"] == assoc_df["RA"]
    log.write(f" -GL_RA matched RA for {sum(is_matched)} variants", verbose=verbose)
    # Ensure we're assigning numeric values
    beta_values = pd.to_numeric(assoc_df.loc[is_matched, "Beta"], errors='coerce')
    assoc_df.loc[is_matched, "GC_Beta_GL_RA"] = beta_values.values
    
    is_flipped_ra = assoc_df["_GL_NRA"] == assoc_df["RA"]
    log.write(f" -GL_NRA matched RA for {sum(is_flipped_ra)} variants. Aligning...", verbose=verbose)
    # Negate numeric values
    beta_values_flipped = pd.to_numeric(assoc_df.loc[is_flipped_ra, "Beta"], errors='coerce')
    assoc_df.loc[is_flipped_ra, "GC_Beta_GL_RA"] = -beta_values_flipped.values
    
    # Align with effect allele (EA)
    is_matched_ea = assoc_df["EA"] == assoc_df["RA"]
    log.write(f" -EA matched RA for {sum(is_matched_ea)} variants", verbose=verbose)
    beta_values_ea = pd.to_numeric(assoc_df.loc[is_matched_ea, "Beta"], errors='coerce')
    assoc_df.loc[is_matched_ea, "GC_Beta_GL_EA"] = beta_values_ea.values
    
    is_flipped_ea = assoc_df["NEA"] == assoc_df["RA"]
    log.write(f" -NEA matched RA for {sum(is_flipped_ea)} variants. Aligning...", verbose=verbose)
    beta_values_ea_flipped = pd.to_numeric(assoc_df.loc[is_flipped_ea, "Beta"], errors='coerce')
    assoc_df.loc[is_flipped_ea, "GC_Beta_GL_EA"] = -beta_values_ea_flipped.values
    
    # Report unmatched
    not_matched = (assoc_df["_GL_RA"] != assoc_df["RA"]) & (assoc_df["_GL_NRA"] != assoc_df["RA"])
    log.write(f" -NEA and EA not matching RA for {sum(not_matched)} variants", verbose=verbose)
    
    return assoc_df
