from typing import TYPE_CHECKING, Optional, Tuple, Any
import subprocess
import numpy as np
import os
import pandas as pd
import tempfile
import shutil
from gwaslab.info.g_Log import Log
from gwaslab.io.io_plink import _process_plink_input_files
from gwaslab.extension import _checking_plink_version
from gwaslab.qc.qc_decorator import with_logging

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _match_sumstats_with_ref_bim(sumstats: pd.DataFrame, 
                                 ref_bim_all: pd.DataFrame, 
                                 has_allele_info: bool, 
                                 log: Log, 
                                 verbose: bool = True) -> int:
    """
    Match sumstats variants with reference BIM using CHR, POS, and optionally EA, NEA.
    
    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe with CHR, POS, and optionally EA, NEA columns
    ref_bim_all : pd.DataFrame
        Reference BIM dataframe with CHR_bim, POS_bim, SNPID, EA_bim, NEA_bim columns
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
    
    # Filter ref_bim to only include CHR/POS that exist in sumstats for efficiency
    # Convert to int for consistent comparison
    ref_bim_all["CHR_bim"] = ref_bim_all["CHR_bim"].astype(int)
    ref_bim_all["POS_bim"] = ref_bim_all["POS_bim"].astype(int)
    
    # Create set of (CHR, POS) tuples from sumstats for fast lookup
    sumstats_chr_pos = set(zip(sumstats["CHR"].astype(int), sumstats["POS"].astype(int)))
    
    # Create a boolean mask for filtering - more efficient than apply
    ref_bim_all["_chr_pos_tuple"] = list(zip(ref_bim_all["CHR_bim"], ref_bim_all["POS_bim"]))
    mask = ref_bim_all["_chr_pos_tuple"].isin(sumstats_chr_pos)
    ref_bim_filtered = ref_bim_all[mask].drop(columns=["_chr_pos_tuple"]).copy()
    
    log.write(" -Filtered reference BIM from {} to {} variants matching sumstats CHR/POS...".format(
        len(ref_bim_all), len(ref_bim_filtered)),verbose=verbose)
    
    if len(ref_bim_filtered) == 0:
        log.write(" -No matching CHR/POS found in reference BIM.",verbose=verbose)
        return 0
    
    # Create mapping dictionary for efficient lookup
    # Map (CHR, POS) -> list of variant dicts with SNPID, EA, NEA
    pos_to_variants = {}
    for _, row in ref_bim_filtered.iterrows():
        key = (int(row["CHR_bim"]), int(row["POS_bim"]))
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
            chr_pos = (int(sumstats.loc[idx, "CHR"]), int(sumstats.loc[idx, "POS"]))
            
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
            chr_pos = (int(sumstats.loc[idx, "CHR"]), int(sumstats.loc[idx, "POS"]))
            
            if chr_pos in pos_to_variants:
                # Use first variant at this position
                sumstats.loc[idx, "SNPID_bim"] = pos_to_variants[chr_pos][0]["SNPID"]
                n_matched += 1
        
        log.write(" -Matched {} variants using CHR, POS only...".format(n_matched),verbose=verbose)
    
    return n_matched

@with_logging(
        start_to_msg="perfrom clumping",
        finished_msg="clumping",
        start_cols=["SNPID","CHR","POS"],
        start_function=".clump()"
)
def _clump(gls: 'Sumstats', 
           vcf: Optional[str] = None, 
           scaled: bool = False, 
           out: Optional[str] = "clumping_plink2", 
           p: str = "P",
           mlog10p: str = "MLOG10P", 
           overwrite: bool = False, 
           study: Optional[str] = None, 
           bfile: Optional[str] = None, 
           pfile: Optional[str] = None,
           threads: int = 1, 
           memory: Optional[int] = None, 
           chrom: Optional[Any] = None, 
           clump_p1: float = 5e-8, 
           clump_p2: float = 5e-8, 
           clump_r2: float = 0.01, 
           clump_kb: int = 250,
           log: Log = Log(),
           verbose: bool = True,
           plink: str = "plink",
           plink2: str = "plink2") -> Tuple[pd.DataFrame, pd.DataFrame, str]:
    """
    Perform LD clumping of GWAS summary statistics using PLINK2.

    Parameters
    ----------
    vcf : str or None, optional
        Path or prefix to reference VCF or genotype data compatible with PLINK2.
        Used when deriving `--pfile` inputs.
    bfile : str or None, optional
        Prefix to PLINK binary files (`.bed/.bim/.fam`). May include "@"
        as a chromosome placeholder.
    pfile : str or None, optional
        Prefix to PLINK2 files (`.pgen/.pvar/.psam`). May include "@"
        as a chromosome placeholder.
    scaled : bool, optional
        If True, clump on `mlog10p` using PLINK2 `--clump-log10`. If False,
        clump on `p`.
    out : str or None, optional
        Output prefix. If None, uses "./{study}_clumpping".
    p : str, optional
        Column name of p-values in `gls.data`.
    mlog10p : str, optional
        Column name of -log10(p) in `gls.data`.
    overwrite : bool, optional
        Whether to overwrite any intermediate reference files produced while
        preparing inputs.
    study : str or None, optional
        Study name used when `out` is None.
    threads : int, optional
        Number of threads to pass to PLINK2 via `--threads`.
    memory : int or None, optional
        Memory limit (MB) for PLINK2 via `--memory`.
    chrom : any, optional
        Unused parameter kept for API compatibility.
    clump_p1 : float, optional
        Primary p-value threshold (`--clump-p1` or `--clump-log10-p1`).
    clump_p2 : float, optional
        Secondary p-value threshold (`--clump-p2` or `--clump-log10-p2`).
    clump_r2 : float, optional
        LD threshold (`--clump-r2`).
    clump_kb : int, optional
        Window size in kilobases (`--clump-kb`).
    log : gwaslab.g_Log.Log, optional
        Logger instance used for progress reporting.
    verbose : bool, optional
        Whether to emit verbose log messages.
    plink : str, optional
        Path to PLINK (v1). Not used directly in clumping.
    plink2 : str, optional
        Path to PLINK2 binary.

    Returns
    -------
    results_sumstats : pandas.DataFrame
        Subset of input summary statistics for clumped lead variants.
    results : pandas.DataFrame
        Concatenated PLINK2 `.clumps` output across processed chromosomes.
    plink_log : str
        Combined PLINK2 log output captured during execution.

    Workflow
    --------
    The clumping process follows these steps:
    
    1. **Filter significant variants**: Extract variants below the p-value threshold
       (clump_p1 or clump_p2) from the input sumstats.
    
    2. **Process reference files**: Convert VCF/BGEN to PLINK format (bfile/pfile) if needed,
       and load BIM/PVAR variant information for matching.
    
    3. **Match variants with reference**: Match sumstats variants with reference BIM using
       CHR, POS, and optionally EA/NEA to assign reference SNPIDs. This ensures PLINK
       uses consistent IDs that match the reference panel.
    
    4. **Create temporary input files**: For each chromosome, create a temporary SNPIDP file
       containing variant IDs and p-values in a temporary directory.
    
    5. **Run PLINK2 clumping**: Execute PLINK2 clumping for each chromosome separately,
       using the reference panel and temporary input files. PLINK2 identifies lead variants
       and their clumped variants based on LD (r²) within the specified window.
    
    6. **Process results**: Read and concatenate clumping results from all chromosomes,
       map BIM SNPIDs back to original sumstats SNPIDs, and filter sumstats to include
       only clumped lead variants.
    
    7. **Cleanup**: Delete temporary files and intermediate clumps output files after
       successful data reload.

    Notes
    -----
    - Writes temporary files in a temporary directory, which are automatically removed.
    - Produces per-chromosome output files "{out}.{chr}.clumps" which are deleted after
      successful reload if delete_files option is used.
    - Variant matching uses CHR, POS, EA, NEA to ensure ID consistency between sumstats
      and reference panel, preventing missing matches due to ID mismatches.

    Examples
    --------
    >>> results_sumstats, results, logstr = _clump(
    ...     bfile="ref/chr@",
    ...     clump_p1=5e-8,
    ...     clump_p2=1e-5,
    ...     clump_r2=0.1,
    ...     clump_kb=250,
    ...     threads=4
    ... )
    """
    ##start function with col checking##########################################################
    
    if out is None:
        out = f"./{study}_clumpping".lstrip('/')
    else:
        out = out.lstrip('/')
    sumstats_id = gls.id
    sumstats = gls.data
    gls.offload()

    ############################################################################################
    ## process reference
    log.write("Start to perform clumping...",verbose=verbose)
    log.write(" -Clumping parameters for PLINK2:",verbose=verbose)
    log.write("  -clump_p1 : {}...".format(clump_p1),verbose=verbose)
    log.write("  -clump_p2 : {}...".format(clump_p2),verbose=verbose)
    log.write("  -clump_kb : {}...".format(clump_kb),verbose=verbose)
    log.write("  -clump_r2 : {}...".format(clump_r2),verbose=verbose)
    if scaled == True:
        log.write(" -Clumping will be performed using {}".format(mlog10p),verbose=verbose)
        clump_log10_p1=-np.log10(clump_p1)
        clump_log10_p2=-np.log10(clump_p2)
        log.write("  -clump_log10_p1 : {}...".format(clump_log10_p1),verbose=verbose)
        log.write("  -clump_log10_p2 : {}...".format(clump_log10_p2),verbose=verbose)
        sumstats = sumstats.loc[sumstats[mlog10p]>min(clump_log10_p1,clump_log10_p2),:].copy()
    # extract lead variants
    else:
        log.write(" -Clumping will be performed using {}".format(p),verbose=verbose)
        sumstats = sumstats.loc[sumstats[p]<max(clump_p1,clump_p2),:].copy()

    if len(sumstats)==0:
        log.write(" -No significant variants after filtering.")
        gls.reload()
        return pd.DataFrame(), pd.DataFrame(), ""
    
    log.write(" -Significant variants on CHR: ",list(sumstats["CHR"].unique()),verbose=verbose)
    
    plink_log=""

    # process reference file
    bfile, plink_log, ref_bim,filetype = _process_plink_input_files(chrlist=sumstats["CHR"].unique(), 
                                                       bfile=bfile, 
                                                       pfile=pfile,
                                                       vcf=vcf, 
                                                       threads=threads,
                                                       plink_log=plink_log, 
                                                       log=log,
                                                       load_bim=True,
                                                       overwrite=overwrite)           
    
    # Concatenate all ref_bim dataframes into one for matching
    if len(ref_bim) > 0:
        ref_bim_all = pd.concat(ref_bim, ignore_index=True)
        # Convert CHR_bim to match sumstats CHR type (might be int or category)
        ref_bim_all["CHR_bim"] = ref_bim_all["CHR_bim"].astype(str).astype(int)
        log.write(" -Total variants in reference BIM: {}...".format(len(ref_bim_all)),verbose=verbose)
    else:
        ref_bim_all = pd.DataFrame()
        log.write(" -Warning: ref_bim is empty. Falling back to SNPID-only matching.",verbose=verbose)
    
    # Match sumstats with ref_bim using helper function
    # Check if EA and NEA columns exist in sumstats for allele-based matching
    has_allele_info = "EA" in sumstats.columns and "NEA" in sumstats.columns
    n_matched = _match_sumstats_with_ref_bim(sumstats, ref_bim_all, has_allele_info, log, verbose)
    
    # Create a temporary directory for all temp files
    temp_dir = tempfile.mkdtemp(prefix="gwaslab_clump_", suffix=f"_{sumstats_id}_")
    log.write(" -Created temporary directory: {}".format(temp_dir),verbose=verbose)
    
    ## process sumstats by CHR
    for i in sumstats["CHR"].unique():
        log.write(" -Processing sumstats for CHR {}...".format(i),verbose=verbose)
        
        if "@" in bfile:
            bfile_to_use = bfile.replace("@",str(i))
        else:
            bfile_to_use = bfile
        
        # checking # variants
        try:
            if filetype=="bfile":
                bim = pd.read_csv(bfile_to_use + ".bim",usecols=[1],header=None,sep="\s+")[1]
            else:
                bim = pd.read_csv(bfile_to_use + ".pvar",usecols=[2],header=None,comment="#",sep="\s+")[2]
            
            snplist = sumstats.loc[sumstats["CHR"]==i,"SNPID"]
            
            # Use SNPID_bim for matching instead of original SNPID
            is_on_both = sumstats.loc[sumstats["CHR"]==i, "SNPID_bim"].isin(bim)

            log.write(" -Variants in reference file: {}...".format(len(bim)),verbose=verbose)
            log.write(" -Variants in sumstats: {}...".format(len(snplist)),verbose=verbose)
            log.write(" -Variants available in both reference and sumstats: {}...".format(sum(is_on_both)),verbose=verbose)

            is_avaialable_variant = (sumstats["CHR"]==i) & (sumstats["SNPID_bim"].isin(bim))

            # Create temp file in temp directory using tempfile
            temp_file = os.path.join(temp_dir, "chr{}.SNPIDP".format(i))
            
            # Use SNPID_bim in the temp file for PLINK clump
            try:
                if scaled == True:
                    sumstats.loc[is_avaialable_variant,["SNPID_bim",mlog10p]].rename(columns={"SNPID_bim":"SNPID"}).to_csv(temp_file,index=False,sep="\t")
                else:
                    sumstats.loc[is_avaialable_variant,["SNPID_bim",p]].rename(columns={"SNPID_bim":"SNPID"}).to_csv(temp_file,index=False,sep="\t")
            except Exception as e:
                log.write(" -Error creating temp file for CHR {}: {}".format(i, str(e)),verbose=verbose)
        except Exception as e:
            log.write(" -Not available for: {}... Error: {}".format(i, str(e)),verbose=verbose)
        
    # create a empty dataframe for combining results from each CHR 
    results = pd.DataFrame()
    
    # Track clumps files for cleanup after successful reload
    clumps_files = []

    # clumping using plink
    try:
        for i in sumstats["CHR"].unique():
            chrom = i
            # temp file in temp directory
            clump = os.path.join(temp_dir, "chr{}.SNPIDP".format(chrom))
            # output prefix
            out_single_chr= out + ".{}".format(chrom)
            # Track clumps file
            clump_result_file = "{}.clumps".format(out_single_chr)
            clumps_files.append(clump_result_file)
            
            if "@" in bfile:
                bfile_to_use = bfile.replace("@",str(i))
            else:
                bfile_to_use = bfile

            # Check if temp file exists before proceeding
            if not os.path.exists(clump):
                log.write(" -Skipping clumping for CHR {}: temp file not found: {}".format(i, clump),verbose=verbose)
                continue
            
            log.write(" -Performing clumping for CHR {}...".format(i),verbose=verbose)
            log = _checking_plink_version(plink2=plink2, log=log)
            if memory is not None:
                memory_flag = "--memory {}".format(memory)
            else:
                memory_flag = ""
            
            if filetype=="bfile":
                file_flag = "--bfile {}".format(bfile_to_use) 
            else:
                file_flag = "--pfile {}".format(bfile_to_use) 
        
            if scaled == True:
                # clumping using LOG10P
                script = """
                {} \
                    {}\
                    --chr {} \
                    --clump {} \
                    --clump-log10 \
                    --clump-field {} \
                    --clump-snp-field SNPID \
                    --clump-log10-p1 {} \
                    --clump-log10-p2 {} \
                    --clump-r2 {} \
                    --clump-kb {} \
                    --threads {} {}\
                    --out {}
                """.format(plink2, file_flag, chrom, clump, mlog10p,clump_log10_p1, clump_log10_p2, clump_r2, clump_kb, threads, memory_flag, out_single_chr)    
            else:
                # clumping using P
                script = """
                {} \
                    {}\
                    --chr {} \
                    --clump {} \
                    --clump-field {} \
                    --clump-snp-field SNPID \
                    --clump-p1 {} \
                    --clump-p2 {} \
                    --clump-r2 {} \
                    --clump-kb {} \
                    --threads {} {}\
                    --out {}
                """.format(plink2,file_flag, chrom, clump, p, clump_p1, clump_p2, clump_r2, clump_kb, threads, memory_flag, out_single_chr)
            
            try:
                output = subprocess.check_output(script, stderr=subprocess.STDOUT, shell=True,text=True)
                log.write(" -Saved results for CHR {} to : {}".format(i,"{}.clumps".format(out_single_chr)),verbose=verbose)
                plink_log +=output + "\n"
            except subprocess.CalledProcessError as e:
                log.write(" -Error during clumping for CHR {}: {}".format(i, e.output),verbose=verbose)
                plink_log += e.output + "\n"
            
            # Try to read clumping results
            try:
                if os.path.exists(clump_result_file):
                    clumped = pd.read_csv(clump_result_file,sep="\s+")
                    results = pd.concat([results,clumped],ignore_index=True)
                else:
                    log.write(" -Clumping result file not found for CHR {}: {}".format(i, clump_result_file),verbose=verbose)
            except Exception as e:
                log.write(" -Failed to read clumping results for CHR {}: {}".format(i, str(e)),verbose=verbose)
            
    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
                log.write(" -Cleaned up temporary directory: {}".format(temp_dir),verbose=verbose)
            except Exception as e:
                log.write(" -Warning: Could not remove temporary directory {}: {}".format(temp_dir, str(e)),verbose=verbose)
    
    results = results.sort_values(by=["#CHROM","POS"]).rename(columns={"#CHROM":"CHR","ID":"SNPID"})
    
    # Map BIM SNPIDs back to original SNPIDs in results
    # Create a mapping from SNPID_bim to original SNPID
    if "SNPID_bim" in sumstats.columns:
        snpid_mapping = sumstats.set_index("SNPID_bim")["SNPID"].to_dict()
        # Update results SNPID column: if BIM SNPID exists in mapping, use original SNPID
        results["SNPID"] = results["SNPID"].map(snpid_mapping).fillna(results["SNPID"])
        log.write(" -Mapped BIM SNPIDs back to original SNPIDs in clumping results...",verbose=verbose)
    
    # Filter sumstats using the mapped SNPIDs
    results_sumstats = sumstats.loc[sumstats["SNPID"].isin(results["SNPID"]),:].copy()
    
    # Reload gls and clean up clumps files if successful
    # Pass clumps_files to reload() to delete them after successful reload
    gls.reload(delete_files=clumps_files)
    
    return results_sumstats, results, plink_log



       
