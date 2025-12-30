from shutil import which
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import scipy as sp
from typing import TYPE_CHECKING, Optional, List, Dict, Any, Union, Tuple
from allel import GenotypeArray
from allel import read_vcf
from allel import rogers_huff_r_between
import matplotlib as mpl
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from adjustText import adjust_text
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_recombination_rate
from gwaslab.util.util_in_filter_value import _get_flanking
from gwaslab.io.io_vcf import auto_check_vcf_chr_dict
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.qc.qc_normalize_args import _normalize_region
# unmatched SNP list 1 

# for each SNP in unmatched SNP list 1:
    # calculate LD for the variants
    
    # if the variant not in reference : return
    
    # if no variants reached threshold : return

    # if there are vailable possible proxies : 
        # from high to low: 
            #check if in outcome and exposure snp list
            #replace

def _extract_with_ld_proxy(
    snplist: Optional[List[str]] = None,
    common_sumstats: Optional[pd.DataFrame] = None,
    sumstats1: Optional[pd.DataFrame] = None,
    vcf_path: Optional[str] = None, 
    vcf_chr_dict: Optional[Dict[str, str]] = None,
    tabix: Optional[str] = None,
    log: Log = Log(), 
    verbose: bool = True, 
    windowsizekb: int = 100,
    ld_threshold: float = 0.8
) -> pd.DataFrame:
    """
    Find LD proxies for SNPs in the list using a VCF reference.

    Parameters
    ----------
    snplist : list or None
        List of SNP IDs to find proxies for.
    common_sumstats : pandas.DataFrame or None
        DataFrame with summary statistics for existing SNPs.
    sumstats1 : pandas.DataFrame or None
        DataFrame with summary statistics for SNPs to check.
    vcf_path : str or None
        Path to the VCF file for reference genotypes.
    vcf_chr_dict : dict or None
        Dictionary mapping chromosomes to VCF region strings.
    tabix : str or None
        Path to tabix executable.
    log : gwaslab.g_Log.Log
        Logging object.
    verbose : bool
        If True, write detailed logs.
    windowsizekb : int
        Size in kb for the flanking region around each SNP.
    ld_threshold : float
        Minimum R^2 value to consider a proxy valid.

    Returns
    -------
    pandas.DataFrame
        Extracted summary statistics including matched proxies.
    """
    
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)
    if tabix is None:
        tabix = which("tabix")
    vcf_chr_dict = auto_check_vcf_chr_dict(vcf_path=vcf_path, vcf_chr_dict=vcf_chr_dict, verbose=verbose, log=log)

    is_needed=[]
    no_need  =[]
    
    print(common_sumstats.head())
    for i in snplist:
        if i in common_sumstats["SNPID"].values:
            no_need.append(i)
        else:
            is_needed.append(i)

    extracted_sumstats  = common_sumstats.loc[common_sumstats["SNPID"].isin(no_need),:]
    
    ld_proxies  = pd.DataFrame()
    in_sumstats = sumstats1.loc[sumstats1["SNPID"].isin(is_needed),:]
    
    if len(in_sumstats)==0:
        log.write(" -No available variants for LD proxy checking...Skipping... ", verbose=verbose)
    else:
        log.write(" -{} available variants for LD proxy checking... ".format(len(in_sumstats)), verbose=verbose)

    for index,row in in_sumstats.iterrows():
        # determine SNP and select region
        snpid = row["SNPID"]
        chrom= int(row["CHR"])
        start= int(row["POS"]-windowsizekb*1000)
        end=   int(row["POS"]+windowsizekb*1000)

        region = (chrom, start, end)
        
        ###  #######################################################################################
        #is_flanking = common_sumstats["CHR"] == chrom & common_sumstats["CHR"]>start & common_sumstats["CHR"]<end
        #flanking_sumstats = common_sumstats.loc[is_flanking,:]
        flanking_sumstats = common_sumstats.query('CHR == @chrom and @start < POS < @end',engine='python').copy()
        
        log.write(" -Extract {} variants in flanking region of {} for checking: {}:{}-{}".format(len(flanking_sumstats), snpid, chrom, start, end), verbose=verbose)

        if len(flanking_sumstats)==0:
            log.write("  -No available variants in the region...Skipping!", verbose=verbose)
            continue
        
        _get_rsq_single(in_sumstats.loc[index,["POS","NEA_1","EA_1"]], 
                        row_pos=row["POS"], 
                        vcf_path=vcf_path, 
                        region=region,
                        log=log, 
                        verbose=verbose, 
                        vcf_chr_dict=vcf_chr_dict, 
                        tabix=tabix)


        flanking_sumstats = _get_rsq(row =in_sumstats.loc[index,["POS","NEA_1","EA_1"]],
                                     sumstats = flanking_sumstats, 
                                     row_pos=row["POS"], 
                                     vcf_path=vcf_path, 
                                     region=region,
                                     log=log, 
                                     verbose=verbose, 
                                     vcf_chr_dict=vcf_chr_dict, 
                                     tabix=tabix)
        if flanking_sumstats is None:
            log.write("  -{} is not found in the vcf...Skipping!".format(snpid))
            continue
        flanking_sumstats = flanking_sumstats.loc[flanking_sumstats["RSQ"]>ld_threshold,:]
        
        log.write("  -Variants in LD with {} (RSQ > {}): {}".format(snpid, ld_threshold,len(flanking_sumstats)), verbose=verbose)
        
        if len(flanking_sumstats)>0:
            flanking_sumstats["LD_REF_VARIANT"]= snpid
            for i,row_with_rsq in flanking_sumstats.iterrows():
                if row_with_rsq["SNPID"] in common_sumstats["SNPID"].values:
                    log.write("  -Proxy for {} is found: {} (LD RSQ= {})".format(snpid, row_with_rsq["SNPID"], row_with_rsq["RSQ"]))
                    row_with_rsq = pd.DataFrame(row_with_rsq)
                    ld_proxies = pd.concat([ld_proxies, row_with_rsq.T], ignore_index=True)
                    break
            
    
    extracted_sumstats = pd.concat([extracted_sumstats, ld_proxies],ignore_index=True)

    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return extracted_sumstats


def _extract_vcf_proxies_not_in_sumstats(
    ref_genotype: Optional[Dict[str, Any]], 
    vcf_proxies_df: pd.DataFrame, 
    chrom: int, 
    snpid: str, 
    ld_threshold: float, 
    common_sumstats: pd.DataFrame, 
    log: Log, 
    verbose: bool
) -> pd.DataFrame:
    """
    Extract proxy variants from VCF that are not in common_sumstats.
    
    Parameters
    ----------
    ref_genotype : dict
        VCF genotype data loaded from allel.read_vcf
    vcf_proxies_df : pandas.DataFrame
        DataFrame with SNPID and RSQ columns for VCF variants
    chrom : int
        Chromosome number
    snpid : str
        Reference SNP ID
    ld_threshold : float
        LD threshold
    common_sumstats : pandas.DataFrame
        Summary statistics
    log : Log
        Logging object
    verbose : bool
        Verbose flag
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with proxy variants from VCF not in sumstats
    """
    if ref_genotype is None or len(vcf_proxies_df) == 0:
        return pd.DataFrame()
    
    # Filter to variants not in common_sumstats and above threshold
    vcf_proxies_df = vcf_proxies_df.loc[
        (vcf_proxies_df["RSQ"] > ld_threshold) & 
        (~vcf_proxies_df["SNPID"].isin(common_sumstats["SNPID"].values)), :
    ]
    
    if len(vcf_proxies_df) == 0:
        return pd.DataFrame()
    
    # Create mapping from SNPID to index for efficient lookup
    vcf_id_to_idx = {ref_genotype["variants/ID"][i]: i 
                     for i in range(len(ref_genotype["variants/ID"]))}
    
    # Extract variant information from VCF
    vcf_variants = []
    for vcf_snpid in vcf_proxies_df["SNPID"].values:
        if vcf_snpid in vcf_id_to_idx:
            vcf_idx = vcf_id_to_idx[vcf_snpid]
            rsq_value = vcf_proxies_df.loc[vcf_proxies_df["SNPID"] == vcf_snpid, "RSQ"].values[0]
            
            # Extract ALT allele (first ALT if available, otherwise use REF)
            alt_alleles = ref_genotype["variants/ALT"][vcf_idx]
            ea = alt_alleles[0] if len(alt_alleles) > 0 else ref_genotype["variants/REF"][vcf_idx]
            
            vcf_variants.append({
                "SNPID": ref_genotype["variants/ID"][vcf_idx],
                "CHR": chrom,
                "POS": ref_genotype["variants/POS"][vcf_idx],
                "EA": ea,
                "NEA": ref_genotype["variants/REF"][vcf_idx],
                "RSQ": rsq_value,
                "LD_REF_VARIANT": snpid
            })
    
    if len(vcf_variants) > 0:
        result = pd.DataFrame(vcf_variants)
        log.write("  -Found {} proxy variants from VCF not in sumstats (RSQ > {})".format(
            len(result), ld_threshold), verbose=verbose)
        return result
    
    return pd.DataFrame()


@with_logging(
        start_to_msg="find LD proxies for variants",
        finished_msg="finding LD proxies for variants",
        start_cols=["SNPID","CHR","POS","EA","NEA"],
        start_function="_extract_ld_proxy"
)
def _extract_ld_proxy(
    snplist: Optional[List[str]] = None,
    common_sumstats: Optional[Union['Sumstats', pd.DataFrame]] = None,
    vcf_path: Optional[str] = None, 
    vcf_chr_dict: Optional[Dict[str, str]] = None,
    tabix: Optional[str] = None,
    log: Log = Log(), 
    verbose: bool = True, 
    windowsizekb: int = 500,
    ld_threshold: float = 0.8,
    include_all: bool = False
) -> pd.DataFrame:
    """
    Find LD proxies within the sumstats for SNPs in the list using a VCF reference.

    Parameters
    ----------
    snplist : list or None
        List of full SNPIDs to find proxies for.
    common_sumstats : Sumstats or pd.DataFrame or None
        Sumstats object or DataFrame with summary statistics for existing SNPs.
    vcf_path : str or None
        Path to the VCF file for reference genotypes.
    verbose : bool
        If True, write detailed logs.
    windowsizekb : int
        Size in kb for the flanking region around each SNP.
    ld_threshold : float
        Minimum R^2 value to consider a proxy valid.
    include_all : bool, default=False
        If True, include proxy variants from VCF that are not in common_sumstats.
        When False, only returns proxies that are already in common_sumstats.

    Returns
    -------
    pandas.DataFrame
        Extracted summary statistics including matched proxies sorted by LD strength.

    Less used parametrs
    -------
    vcf_chr_dict : dict or None
        Dictionary mapping chromosomes to VCF region strings.
    log : gwaslab.g_Log.Log
        Logging object.
    tabix : str or None
        Path to tabix executable.
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if common_sumstats is not None:
        if isinstance(common_sumstats, pd.DataFrame):
            common_sumstats = common_sumstats
        else:
            common_sumstats = common_sumstats.data
    
    # Initialize VCF setup
    log.write("Start to load reference genotype...", verbose=verbose)
    log.write(" -reference vcf path : "+ vcf_path, verbose=verbose)
    if tabix is None:
        tabix = which("tabix")
    vcf_chr_dict = auto_check_vcf_chr_dict(vcf_path=vcf_path, vcf_chr_dict=vcf_chr_dict, verbose=verbose, log=log)
    
    # Filter input SNPs
    in_sumstats = common_sumstats.loc[common_sumstats["SNPID"].isin(snplist), :].copy()
    
    if len(in_sumstats) == 0:
        log.write(" -No available variants for LD proxy checking...Skipping... ", verbose=verbose)
        return pd.DataFrame()
    
    log.write(" -{} available variants for LD proxy checking... ".format(len(in_sumstats)), verbose=verbose)
    
    # Pre-compute common_sumstats SNPID set for faster lookups
    common_snpid_set = set(common_sumstats["SNPID"].values)
    all_ld_proxies = []

    # Process each SNP
    for index, row in in_sumstats.iterrows():
        snpid = row["SNPID"]
        chrom = int(row["CHR"])
        pos = int(row["POS"])
        start = pos - windowsizekb * 1000
        end = pos + windowsizekb * 1000
        
        # Normalize region (ensures start < end, proper chromosome format)
        region = _normalize_region((chrom, start, end), log=log, verbose=verbose)
        chrom, start, end = region
        
        # Extract flanking variants from sumstats
        flanking_sumstats = common_sumstats.query(
            'CHR == @chrom and @start < POS < @end', engine='python'
        ).copy()
        
        log.write(" -Extract {} variants in flanking region of {} for checking: {}:{}-{}".format(
            len(flanking_sumstats), snpid, chrom, start, end), verbose=verbose)

        if len(flanking_sumstats) == 0:
            log.write("  -No available variants in the region...Skipping!", verbose=verbose)
            continue

        # Calculate LD for variants in sumstats
        flanking_sumstats = _get_rsq(
            row=row[["POS", "NEA", "EA"]],
            sumstats=flanking_sumstats, 
            row_pos=pos, 
            vcf_path=vcf_path, 
            region=region,
            log=log, 
            verbose=verbose, 
            vcf_chr_dict=vcf_chr_dict, 
            tabix=tabix
        )
        
        if flanking_sumstats is None:
            log.write("  -{} is not found in the vcf...Skipping!".format(snpid), verbose=verbose)
            continue
        
        # Filter by LD threshold
        flanking_sumstats = flanking_sumstats.loc[flanking_sumstats["RSQ"] > ld_threshold, :].copy()
        log.write("  -Variants in LD with {} (RSQ > {}): {}".format(
            snpid, ld_threshold, len(flanking_sumstats)), verbose=verbose)
        
        # Process sumstats proxies
        if len(flanking_sumstats) > 0:
            flanking_sumstats["LD_REF_VARIANT"] = snpid
            
            # Log top proxy if found in common_sumstats
            for _, proxy_row in flanking_sumstats.iterrows():
                if proxy_row["SNPID"] in common_snpid_set:
                    log.write("  -Top Proxy for {} is found: {} (LD RSQ= {})".format(
                        snpid, proxy_row["SNPID"], proxy_row["RSQ"]), verbose=verbose)
                    break
            
            all_ld_proxies.append(flanking_sumstats)
        
        # Get variants from VCF that are not in sumstats if option is enabled
        if include_all:
            try:
                vcf_proxies_df = _get_rsq_single(
                    row=row[["POS", "NEA", "EA"]],
                    row_pos=pos, 
                    vcf_path=vcf_path, 
                    region=region,
                    log=log, 
                    verbose=verbose, 
                    vcf_chr_dict=vcf_chr_dict, 
                    tabix=tabix
                )
                
                if vcf_proxies_df is not None and len(vcf_proxies_df) > 0:
                    # Load VCF data once to extract variant information
                    ref_genotype = read_vcf(
                        vcf_path,
                        region=vcf_chr_dict[chrom] + ":" + str(start) + "-" + str(end),
                        tabix=tabix
                    )
                    
                    vcf_proxies = _extract_vcf_proxies_not_in_sumstats(
                        ref_genotype, vcf_proxies_df, chrom, snpid, ld_threshold,
                        common_sumstats, log, verbose
                    )
                    
                    if len(vcf_proxies) > 0:
                        all_ld_proxies.append(vcf_proxies)
                        
            except Exception as e:
                log.write("  -Error calculating LD for VCF variants: {}".format(str(e)), verbose=verbose)
    
    # Combine all proxies
    if len(all_ld_proxies) == 0:
        log.write("No proxy variants were found at given LD rsq threshold.", verbose=verbose)
        return pd.DataFrame()
    
    ld_proxies = pd.concat(all_ld_proxies, ignore_index=True)
    
    # Drop REFINDEX column if it exists
    if "REFINDEX" in ld_proxies.columns:
        ld_proxies = ld_proxies.drop(columns=["REFINDEX"])
    
    log.write("Finished loading reference genotype successfully!", verbose=verbose)
    return ld_proxies.sort_values(by="RSQ", ascending=False)

        


def _get_rsq(
    row: pd.Series,
    sumstats: pd.DataFrame,
    row_pos: int,
    vcf_path: str, 
    region: Tuple[int, int, int],
    log: Log, 
    verbose: bool, 
    vcf_chr_dict: Dict[str, str],
    tabix: Optional[str]
) -> pd.DataFrame:
        #load genotype data of the targeted region
        ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
        
        if ref_genotype is None:
            log.warning("No data was retrieved. Skipping ...", verbose=verbose)
            ref_genotype=dict()
            ref_genotype["variants/POS"]=np.array([],dtype="int64")
            return None
        
        log.write("  -Retrieving index...", verbose=verbose)
        log.write("  -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
        #  match sumstats pos and ref pos: 
        # get ref index for its first appearance of sumstats pos
        #######################################################################################
        def match_varaint(x):
            # x: "POS,NEA,EA"
            if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
                if len(np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0])>1:  
                # multiple position matches
                    for j in np.where(ref_genotype["variants/POS"] == x.iloc[0])[0]:
                    # for each possible match, compare ref and alt
                        if x.iloc[1] == ref_genotype["variants/REF"][j]:
                            if x.iloc[2] in ref_genotype["variants/ALT"][j]:
                                return j
                        elif x.iloc[1] in ref_genotype["variants/ALT"][j]:
                            if x.iloc[2] == ref_genotype["variants/REF"][j]:
                                return j
                        else:
                            return None
                else: 
                    # single match
                    return np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0][0] 
            else:
                # no position match
                return None
        log.write("  -Matching variants using POS, NEA, EA ...", verbose=verbose)

        sumstats["REFINDEX"] = sumstats.loc[:,["POS","NEA","EA"]].apply(lambda x: match_varaint(x), axis=1)
        log.write("  -Matched variants in sumstats and vcf:{} ".format(sum(~sumstats["REFINDEX"].isna())))
        #############################################################################################
        lead_pos = row_pos

        # if lead pos is available: 
        if lead_pos in ref_genotype["variants/POS"]:
            
            # get ref index for lead snp
            lead_snp_ref_index = match_varaint(row)
            #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

            # non-na other snp index
            other_snps_ref_index = sumstats["REFINDEX"].dropna().astype("int").values
            # get genotype 
            lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
            other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
            
            log.write("  -Calculating Rsq...", verbose=verbose)
            
            if len(other_snp_genotype)>1:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
            else:
                valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)
            sumstats.loc[~sumstats["REFINDEX"].isna(),"RSQ"] = valid_r2
        else:
            log.write("  -Lead SNP not found in reference...", verbose=verbose)
            sumstats["RSQ"]=None
        
        sumstats["RSQ"] = sumstats["RSQ"].astype("float")
        return sumstats

def _check_if_in_sumstats2(row: pd.Series, sumstast: pd.DataFrame) -> None:
    pass

def _get_rsq_single(
    row: pd.Series,
    row_pos: int,
    vcf_path: str, 
    region: Tuple[int, int, int],
    log: Log, 
    verbose: bool, 
    vcf_chr_dict: Dict[str, str],
    tabix: Optional[str]
) -> Optional[pd.DataFrame]:
    #load genotype data of the targeted region
    ref_genotype = read_vcf(vcf_path,region=vcf_chr_dict[region[0]]+":"+str(region[1])+"-"+str(region[2]),tabix=tabix)
    
    if ref_genotype is None:
        log.warning("No data was retrieved. Skipping ...", verbose=verbose)
        ref_genotype=dict()
        ref_genotype["variants/POS"]=np.array([],dtype="int64")
        return None
    
    log.write("  -Retrieving index...", verbose=verbose)
    log.write("  -Ref variants in the region: {}".format(len(ref_genotype["variants/POS"])), verbose=verbose)
    #  match sumstats pos and ref pos: 
    # get ref index for its first appearance of sumstats pos
    #######################################################################################
    def match_varaint(x):
        # x: "POS,NEA,EA"
        if np.any(ref_genotype["variants/POS"] == x.iloc[0]):
            if len(np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0])>1:  
            # multiple position matches
                for j in np.where(ref_genotype["variants/POS"] == x.iloc[0])[0]:
                # for each possible match, compare ref and alt
                    if x.iloc[1] == ref_genotype["variants/REF"][j]:
                        if x.iloc[2] in ref_genotype["variants/ALT"][j]:
                            return j
                    elif x.iloc[1] in ref_genotype["variants/ALT"][j]:
                        if x.iloc[2] == ref_genotype["variants/REF"][j]:
                            return j
                    else:
                        return None
            else: 
                # single match
                return np.where(ref_genotype["variants/POS"] == x.iloc[0] )[0][0] 
        else:
            # no position match
            return None

    #############################################################################################
    lead_pos = row_pos

    # if lead pos is available: 
    if lead_pos in ref_genotype["variants/POS"]:
        
        # get ref index for lead snp
        lead_snp_ref_index = match_varaint(row)
        #lead_snp_ref_index = np.where(ref_genotype["variants/POS"] == lead_pos)[0][0]

        # non-na other snp index
        other_snps_ref_index = list(range(len(ref_genotype["calldata/GT"])))
        other_snps_ref_index.remove(lead_snp_ref_index)

        # get genotype 
        lead_snp_genotype = GenotypeArray([ref_genotype["calldata/GT"][lead_snp_ref_index]]).to_n_alt()
        other_snp_genotype = GenotypeArray(ref_genotype["calldata/GT"][other_snps_ref_index]).to_n_alt()
        
        log.write("  -Calculating Rsq...", verbose=verbose)
        
        if len(other_snp_genotype)>1:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype)[0],2)
        else:
            valid_r2= np.power(rogers_huff_r_between(lead_snp_genotype,other_snp_genotype),2)

        ld_proxy = pd.DataFrame( {"SNPID":ref_genotype["variants/ID"][other_snps_ref_index],"RSQ":valid_r2 })

    return  ld_proxy.sort_values(by="RSQ",ascending=False)
