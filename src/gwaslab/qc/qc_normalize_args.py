from typing import TYPE_CHECKING, Optional, Union, List, Tuple, Dict, Any
import pandas as pd
import re
from gwaslab.info.g_Log import Log
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.io.io_process_kwargs import normalize_series_inputs

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _parse_flanking(flanking_str: str) -> int:
    '''
    Parse flanking value, supporting both base pairs and kilobases.
    If 'kb' suffix is present (case-insensitive), converts to base pairs.
    
    Parameters:
    -----------
    flanking_str : str
        Flanking value, e.g., '500', '500kb', '500KB', or '500 kb'
    
    Returns:
    --------
    int
        Flanking size in base pairs
    '''
    flanking_str = flanking_str.strip()
    # Check for kb suffix (case-insensitive, with optional space)
    kb_match = re.match(r'^(.+?)\s*kb\s*$', flanking_str, re.IGNORECASE)
    if kb_match:
        # Has 'kb' suffix, convert to base pairs
        value = float(kb_match.group(1).strip())
        return int(value * 1000)
    else:
        # Assume base pairs
        return int(float(flanking_str))

def _normalize_region(
    region: Union[str, Tuple[Union[str, int], int, int], List[Union[str, int, float]]],
    chr_dict: Optional[Dict[str, int]] = None,
    sumstats: Optional[pd.DataFrame] = None,
    snpid: str = "SNPID",
    rsid: str = "rsID",
    chrom_col: str = "CHR",
    pos_col: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    log: Log = Log(),
    verbose: bool = True
) -> Optional[Tuple[int, int, int]]:
    '''
    Normalize a region input to a tuple (chr, start, end) with integer
    positions and standardized chromosome notation.
    Accepts tuples/lists or strings in multiple formats:
    - 'chr:start-end' (e.g., 'chr1:12345-67890')
    - 'chr:pos:flanking' (e.g., 'chr1:1500:500' or 'chr1:1500:500kb' -> chr1:1000-2000)
    - 'snpid:flanking' (e.g., 'rs123:500' or 'rs123:500kb' or '1:12345:A:T:500kb' -> requires sumstats parameter)
    '''
    if region is None:
        return None
    if chr_dict is None:
        chr_dict = get_chr_to_number()
    if isinstance(region, str):
        s = region.strip()
        parts = s.split(":")
        
        # Format: chr:start-end
        if ":" in s and "-" in s and len(parts) == 2:
            left, right = s.split(":", 1)
            start_str, end_str = right.split("-", 1)
            chrom = left.strip().upper().lstrip("CHR")
            start = int(float(start_str.strip()))
            end = int(float(end_str.strip()))
        
        # Format: chr:pos:flanking (3 parts, no "-", last part is numeric or has kb suffix)
        elif len(parts) == 3 and "-" not in s:
            try:
                chrom_str = parts[0].strip().upper().lstrip("CHR")
                pos_str = parts[1].strip()
                flanking_str = parts[2].strip()
                
                chrom = chrom_str
                center_pos = int(float(pos_str))
                flanking = _parse_flanking(flanking_str)
                start = center_pos - flanking
                end = center_pos + flanking
            except (ValueError, IndexError):
                # Not chr:pos:flanking format, might be snpid:flanking with colons in snpid
                raise ValueError("Region string must be in one of these formats: 'chr:start-end', 'chr:pos:flanking', or 'snpid:flanking'")
        
        # Format: snpid:flanking
        # Can be: "rs123:500" (2 parts) or "1:2:A:T:500" (5 parts with chr:pos:ea:nea:flanking)
        elif (len(parts) == 2 or len(parts) == 5) and "-" not in s:
            if sumstats is None:
                raise ValueError("Region format 'snpid:flanking' requires sumstats parameter to look up SNP position")
            
            # Handle 5-part format: chr:pos:ea:nea:flanking
            if len(parts) == 5:
                try:
                    flanking_str = parts[4].strip()
                    flanking = _parse_flanking(flanking_str)
                    # First 4 parts form the SNPID: chr:pos:ea:nea
                    snpid_value = ":".join(parts[0:4])
                except (ValueError, IndexError):
                    raise ValueError("Region string must be in one of these formats: 'chr:start-end', 'chr:pos:flanking', or 'snpid:flanking'")
            else:
                # 2-part format: snpid:flanking
                snpid_value = parts[0].strip()
                flanking_str = parts[1].strip()
                flanking = _parse_flanking(flanking_str)
            
            # Search for the SNP in sumstats
            found = None
            
            # Check if snpid_value is in chr:pos:ea:nea format (e.g., "1:12345:A:T")
            snpid_match = re.match(r'^(chr|Chr|CHR)?(\d+)[:_-](\d+)([:_-]([ATCG]+)[:_-]([ATCG]+))?$', snpid_value, flags=0)
            
            if snpid_match is not None:
                # Parse chr:pos:ea:nea format
                single_chrom = int(snpid_match.group(2))
                single_pos = int(snpid_match.group(3))
                
                if snpid_match.group(4) is not None:
                    # Has alleles: search by chr, pos, and alleles
                    single_ea = snpid_match.group(5)
                    single_nea = snpid_match.group(6)
                    
                    # Match alleles (can be swapped)
                    if ea in sumstats.columns and nea in sumstats.columns:
                        allele_match = ((sumstats[nea] == single_nea) & (sumstats[ea] == single_ea)) | \
                                      ((sumstats[nea] == single_ea) & (sumstats[ea] == single_nea))
                        mask = (sumstats[pos_col] == single_pos) & (sumstats[chrom_col] == single_chrom) & allele_match
                    else:
                        mask = (sumstats[pos_col] == single_pos) & (sumstats[chrom_col] == single_chrom)
                else:
                    # No alleles: search by chr and pos only
                    mask = (sumstats[pos_col] == single_pos) & (sumstats[chrom_col] == single_chrom)
                
                if mask.any():
                    found = sumstats.loc[mask].iloc[0]
            
            # If not found by coordinate format, try rsID
            if found is None and rsid in sumstats.columns:
                mask = sumstats[rsid] == snpid_value
                if mask.any():
                    found = sumstats.loc[mask].iloc[0]
            
            # If still not found, try SNPID column
            if found is None and snpid in sumstats.columns:
                mask = sumstats[snpid] == snpid_value
                if mask.any():
                    found = sumstats.loc[mask].iloc[0]
            
            if found is None:
                raise ValueError("SNP '{}' not found in sumstats".format(snpid_value))
            
            chrom = str(found[chrom_col]).strip().upper().lstrip("CHR")
            center_pos = int(float(found[pos_col]))
            start = center_pos - flanking
            end = center_pos + flanking
        
        else:
            raise ValueError("Region string must be in one of these formats: 'chr:start-end', 'chr:pos:flanking', or 'snpid:flanking'")
    else:
        if len(region) != 3:
            raise ValueError("Region must be a tuple/list of (chr, start, end)")
        chrom_raw, start_raw, end_raw = region
        chrom = str(chrom_raw).strip().upper().lstrip("CHR")
        start = int(float(start_raw))
        end = int(float(end_raw))
    if chrom in chr_dict:
        chrom = chr_dict[chrom]
    else:
        try:
            chrom = int(chrom)
        except Exception:
            raise ValueError("Chromosome '{}' is not recognized".format(chrom))
    if start > end:
        start, end = end, start
    log.write(" -Normalized region: (CHR={}, START={}, END={})".format(chrom, start, end), verbose=verbose)
    return (chrom, start, end)

@normalize_series_inputs(keys=["highlight","pinpoint","anno_set"])
def _normalize_group(
    items: Union[str, List[str], List[List[str]]],
    colors: Union[str, List[str]],
    item_name: str,
    log: Log,
    verbose: bool,
    is_chrpos_mode: bool = False
) -> Tuple[List[Union[str, List[str]]], Union[str, List[str]]]:
    """
    Helper to normalize and log highlight/pinpoint groups.
    """
    # 1. Single string -> list
    if isinstance(items, str):
        items = [items]

    # 2. Check if grouped (list of lists/tuples) or flat (list of strings)
    #    Note: If empty list is passed, access to items[0] would fail, 
    #    but the calling code usually checks len(items) > 0.
    is_grouped = pd.api.types.is_list_like(items[0])

    if is_grouped and (not is_chrpos_mode):
        # If grouped, ensure we have enough colors
        if not isinstance(colors, list):
             colors = [colors] * len(items)
        
        if len(items) != len(colors):
            log.warning(f"Number of {item_name} groups does not match number of provided colors.")
        
        # Log each group with its assigned color
        for i, item_set in enumerate(items):
            color_i = colors[i % len(colors)]
            log.write(f" -Set {i + 1} {item_name}s ({color_i}) : " + ",".join(item_set), verbose=verbose)
    else:
        # Flat list or chr:pos mode
        # If is_grouped is True here, it implies is_chrpos_mode is True (for highlight)
        content_str = items if is_grouped else ",".join(items)
        log.write(f" -{item_name.capitalize()}s ({colors}): {content_str}", verbose=verbose)
    
    return items, colors
