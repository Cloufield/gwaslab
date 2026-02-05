from typing import TYPE_CHECKING, Union, List, Dict, Set
import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
import re

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

@with_logging(
        start_to_msg="initialize gl.SumstatsSet",
        finished_msg="initializing gl.SumstatsSet"
)
def _extract_variant(variant_set: Union[List[str], List[List[Union[int, float]]], Set[str]], sumstats_dic: Dict[str, 'Sumstats'], log: Log = Log(), verbose: bool = True) -> pd.DataFrame:
    """
    Extract specified variants from multiple Sumstats objects and combine them into a single DataFrame.

    This function searches for variants across multiple GWAS summary statistics datasets
    and returns a combined DataFrame containing the matched variants with a "STUDY" column
    indicating the source dataset.

    Parameters
    ----------
    variant_set : Union[List[str], List[List[Union[int, float]]], Set[str]]
        A collection of variants to extract. Can be specified as:
        - List of SNPID strings (e.g., ["rs12345", "rs67890"])
        - List of [chromosome, position] pairs (e.g., [[1, 12345], [2, 67890]])
        - Set of SNPID strings
        - Variant strings in format "chr1:12345:A:T" (parsed for CHR:POS matching)
    sumstats_dic : Dict[str, Sumstats]
        A dictionary mapping study names (keys) to Sumstats objects (values).
    log : Log, optional
        A Log object for recording messages. Defaults to a new Log instance.
    verbose : bool, optional
        Whether to print log messages. Defaults to True.

    Returns
    -------
    pd.DataFrame
        A combined DataFrame containing the extracted variants from all studies.
        Includes columns:
        - "STUDY": The study name (key from sumstats_dic)
        - Standard columns if available: "SNPID", "EA", "NEA", "CHR", "POS",
          "BETA", "SE", "P", "MLOG10P", "EAF", "MAF", "STATUS"

    Notes
    -----
    Variant matching is performed using the following logic:
    - If variant is a list-like (e.g., [chrom, pos]): matches by CHR and POS columns
    - If variant is a string: matches by exact SNPID match
    - String variants matching the pattern "chr1:12345:A:T" (with optional "chr" prefix,
      using ":", "_", or "-" as separators) are also matched by CHR and POS

    Examples
    --------
    >>> # Extract variants by SNPID
    >>> result = _extract_variant(["rs12345", "rs67890"], {"study1": sumstats1, "study2": sumstats2})
    
    >>> # Extract variants by chromosome and position
    >>> result = _extract_variant([[1, 12345], [2, 67890]], {"study1": sumstats1})
    
    >>> # Extract variants using variant ID format
    >>> result = _extract_variant(["chr1:12345:A:T"], {"study1": sumstats1})
    """
    combined = pd.DataFrame()
    for key, sumstats_gls in sumstats_dic.items():
        log.write(" -{} : {}".format(key, sumstats_gls), verbose=verbose)

    for key, sumstats_gls in sumstats_dic.items():
        
        sumstats_single = sumstats_gls.data
        
        # create a boolean col with FALSE 
        is_extract = sumstats_single["SNPID"]!=sumstats_single["SNPID"]
        
        for variant in variant_set:
            
            if pd.api.types.is_list_like(variant):

                chrom=variant[0]
                pos=variant[1]

                is_extract = is_extract | ((sumstats_single["POS"] == pos ) &(sumstats_single["CHR"] == chrom))
            elif pd.api.types.is_string_dtype(type(variant)):
                
                is_extract = is_extract | (sumstats_single["SNPID"] == variant)

                a= re.search(r'^(chr|Chr|CHR)?(\d+)[:_-](\d+)[:_-][ATCG]+[:_-][ATCG]+$', variant, flags=0)
                if a is not None:
                    chrom=int(a[2])
                    pos=int(a[3])
                    is_extract = is_extract | ((sumstats_single["POS"] == pos ) &(sumstats_single["CHR"] == chrom))

        to_extract =  sumstats_single.loc[is_extract,:].copy()
        log.write(" -Extracted {} variants from {}".format(len(to_extract), key),verbose=verbose)
        to_extract["STUDY"] = key
        
        to_extract_cols=["STUDY"]

        default_cols=["SNPID","EA","NEA","CHR","POS","BETA","SE","P","MLOG10P","EAF","MAF","STATUS"]
        
        for i in default_cols:
            if i in sumstats_single.columns:
                to_extract_cols.append(i)

        combined = pd.concat([combined, to_extract[to_extract_cols]], ignore_index=True)
    return combined
