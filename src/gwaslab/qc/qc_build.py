from typing import TYPE_CHECKING, Optional, Union, List, Tuple
import re
import gc
import pandas as pd
import numpy as np
from itertools import repeat
from multiprocessing import  Pool
from functools import partial
from gwaslab.info.g_vchange_status import vchange_status
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_build import BUILD_MAPPINGS, BUILD_DISPLAY_NAMES

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def check_species_compatibility(species: Optional[str], required_species: Union[str, List[str]], operation_name: str, log: Optional[Log] = None, verbose: bool = True) -> bool:
    """
    Check if the current species is compatible with a species-specific operation.
    
    Parameters
    ----------
    species : str
        Current species name (case-insensitive)
    required_species : str or list
        Required species name(s) for the operation. Can be a single string or list of strings.
    operation_name : str
        Name of the operation (for error messages)
    log : Log, optional
        Logging object. If None, errors will be raised without logging.
    verbose : bool, default True
        Whether to print log messages
    
    Returns
    -------
    bool
        True if species is compatible, False otherwise
    
    Raises
    ------
    ValueError
        If species is not compatible with the required species
    """
    if species is None:
        species = "homo sapiens"
    
    species_lower = species.lower()
    
    # Normalize required_species to a list
    if isinstance(required_species, str):
        required_species_list = [required_species.lower()]
    else:
        required_species_list = [s.lower() for s in required_species]
    
    # Check if species matches any of the required species (direct match)
    if species_lower in required_species_list:
        return True
    
    # Also check for common aliases
    species_aliases = {
        "homo sapiens": ["human"],
        "human": ["homo sapiens"],
        "mus musculus": ["mouse"],
        "mouse": ["mus musculus"],
        "rattus norvegicus": ["rat"],
        "rat": ["rattus norvegicus"],
        "gallus gallus": ["chicken"],
        "chicken": ["gallus gallus"],
        "danio rerio": ["zebrafish"],
        "zebrafish": ["danio rerio"],
    }
    
    # Check if any alias matches
    compatible = False
    for req_species in required_species_list:
        # Check aliases of required species
        aliases = species_aliases.get(req_species, [])
        if species_lower in aliases:
            compatible = True
            break
        # Check reverse (if current species has aliases that match required)
        current_aliases = species_aliases.get(species_lower, [])
        if req_species in current_aliases:
            compatible = True
            break
    
    if not compatible:
        if isinstance(required_species, str):
            required_str = required_species
        else:
            required_str = " or ".join(required_species)
        
        error_msg = (
            f"Operation '{operation_name}' is only available for {required_str}. "
            f"Current species is '{species}'. Please use a {required_str} dataset or "
            f"use a species-appropriate alternative."
        )
        
        if log is not None:
            log.error(error_msg, verbose=verbose)
        raise ValueError(error_msg)
    
    return True

def _process_build(build: str,  
                   log: Log = Log(), 
                   verbose: bool = True,
                   species: str = "homo sapiens") -> str:
    """
    Process and normalize genome build identifier based on species.
    
    Parameters
    ----------
    build : str
        Genome build identifier (can be various formats)
    log : Log, optional
        Logging object
    verbose : bool, default True
        Whether to print log messages
    species : str, default "homo sapiens"
        Species name (case-insensitive)
    
    Returns
    -------
    str
        Normalized build identifier, or "99" if unknown
    """
    species_lower = species.lower()
    build_str = str(build).lower()
    
    # Get build mappings for this species
    if species_lower in BUILD_MAPPINGS:
        mappings = BUILD_MAPPINGS[species_lower]
        
        # Check each mapping group
        for aliases, normalized_build in mappings.items():
            if build_str in [str(a).lower() for a in aliases]:
                # Get species-specific display name
                display_name = BUILD_DISPLAY_NAMES.get((species_lower, normalized_build), normalized_build)
                log.write(f" -Genomic coordinates are based on {display_name}...", verbose=verbose)
                return normalized_build
        
        # If no match found, warn and return unknown
        log.warning(f"Version of genomic coordinates is unknown for {species} (build: {build})...", verbose=verbose)
        return "99"
    else:
        # Species not in mappings - return unknown
        log.warning(f"Species '{species}' not recognized for build processing. Using unknown build.", verbose=verbose)
        return "99"

def _set_build(
    sumstats_or_dataframe: Union['Sumstats', pd.DataFrame],
    build: str = "99",
    status: str = "STATUS",
    verbose: bool = True,
    log: Log = Log(),
    species: Optional[str] = None
) -> Tuple[pd.DataFrame, str]:
    """
    Set genome build in sumstats status column.
    
    Parameters
    ----------
    sumstats_or_dataframe : pd.DataFrame or Sumstats object
        Sumstats data or object
    build : str, default "99"
        Genome build identifier
    status : str, default "STATUS"
        Status column name
    verbose : bool, default True
        Whether to print log messages
    log : Log, optional
        Logging object
    species : str, optional
        Species name. If None, will try to extract from Sumstats object metadata.
    
    Returns
    -------
    tuple
        (sumstats DataFrame, processed_build)
    """
    import pandas as pd
    
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
        # If species not provided and we have a DataFrame, default to "homo sapiens"
        if species is None:
            species = "homo sapiens"
    else:
        sumstats = sumstats_or_dataframe.data
        # Try to extract species from Sumstats object metadata
        if species is None:
            try:
                species = sumstats_or_dataframe.meta.get("gwaslab", {}).get("species", "homo sapiens")
            except (AttributeError, KeyError):
                species = "homo sapiens"
    
    # Process build with species information
    processed_build = _process_build(build, log=log, verbose=verbose, species=species)
    
    # Update status column (all builds are now two-digit numeric strings)
    if len(processed_build) == 2 and processed_build != "99":
        # All builds are now two-digit numeric strings (e.g., "19", "10", "06", "09")
        if processed_build[0].isdigit() and processed_build[1].isdigit():
            sumstats[status] = vchange_status(sumstats[status], 1, "139", int(processed_build[0]) * 3)
            sumstats[status] = vchange_status(sumstats[status], 2, "89", int(processed_build[1]) * 3)
    
    return sumstats, processed_build

def _check_build(
    target_build: str,
    build: str = "99",
    status: str = "STATUS",
    verbose: bool = True,
    log: Log = Log(),
    species: Optional[str] = None
) -> bool:
    """
    Check if sumstats build matches target build.
    
    Parameters
    ----------
    target_build : str
        Target genome build identifier
    build : str, default "99"
        Current sumstats build identifier
    status : str, default "STATUS"
        Status column name (unused, kept for compatibility)
    verbose : bool, default True
        Whether to print log messages
    log : Log, optional
        Logging object
    species : str, optional
        Species name. If None, defaults to "homo sapiens".
    
    Returns
    -------
    bool
        True if builds match
    
    Raises
    ------
    ValueError
        If builds are unknown or don't match
    """
    if species is None:
        species = "homo sapiens"
    
    target_build = _process_build(target_build, log=log, verbose=verbose, species=species)
    processed_build = _process_build(build, log=log, verbose=verbose, species=species)
    
    if processed_build == "99":
        raise ValueError("Sumstats build is unknown. Please run infer_build() or set_build()")
    
    if target_build == "99": 
        raise ValueError("Target build is unknown.")
    
    if processed_build != target_build:
        raise ValueError(f"Please make sure sumstats build is {target_build} (current: {processed_build})")
    else:
        log.write(" -Sumstats build matches target build", verbose=verbose)
    
    return True

    
