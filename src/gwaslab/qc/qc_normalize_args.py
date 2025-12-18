import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.io.io_process_kwargs import normalize_series_inputs

def _normalize_region(region, chr_dict=None, log=Log(), verbose=True):
    '''
    Normalize a region input to a tuple (chr, start, end) with integer
    positions and standardized chromosome notation.
    Accepts tuples/lists or strings like 'chr1:12345-67890'.
    '''
    if region is None:
        return None
    if chr_dict is None:
        chr_dict = get_chr_to_number()
    if isinstance(region, str):
        s = region.strip()
        if ":" in s and "-" in s:
            left, right = s.split(":", 1)
            start_str, end_str = right.split("-", 1)
            chrom = left.strip().upper().lstrip("CHR")
            start = int(float(start_str.strip()))
            end = int(float(end_str.strip()))
        else:
            raise ValueError("Region string must be in 'chr:start-end' format")
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
    items, 
    colors, 
    item_name, 
    log, 
    verbose, 
    is_chrpos_mode=False
):
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
