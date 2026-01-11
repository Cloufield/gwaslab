from typing import TYPE_CHECKING, Optional, List, Union
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.util.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.util.rwrapper.util_ex_r_runner import RScriptRunner
from gwaslab.util.general.util_ex_result_manager import ResultManager

if TYPE_CHECKING:
    from gwaslab.g_SumstatsMulti import SumstatsMulti

def _run_hyprcoloc(
    sumstats_multi: 'SumstatsMulti',
    r: str = "Rscript",
    study: str = "Group1",
    traits: Optional[List[str]] = None,
    types: Optional[List[str]] = None,
    loci: Optional[List[str]] = None,
    nstudy: int = 2,
    windowsizekb: int = 1000,
    build: str = "99",
    log: Log = Log(),
    verbose: bool = True,
    timeout: Optional[float] = 1800,  # 30 minutes default timeout
    stop_on_error: bool = False
) -> pd.DataFrame:
    """
    Run hyprcoloc analysis using the new R script execution framework.
    
    Args:
        sumstats_multi: SumstatsMulti object
        r: Path to Rscript executable
        study: Study name/identifier
        traits: List of trait names
        types: List of trait types (not currently used but kept for compatibility)
        loci: List of loci to analyze (if None, extracts significant loci)
        nstudy: Number of studies
        windowsizekb: Window size in kb for locus extraction
        build: Build version
        log: Log instance
        verbose: Whether to log verbosely
        timeout: Timeout for R script execution in seconds (default: 1800)
        stop_on_error: Whether to stop on first error (default: False)
    
    Returns:
        Combined DataFrame with hyprcoloc results
    """
    log.write(" Start to run hyprcoloc from command line:", verbose=verbose)
    
    # Check R version
    log = _checking_r_version(r, log)
    
    # Initialize R script runner and result manager
    runner = RScriptRunner(r=r, log=log, timeout=timeout, cleanup=True)
    result_manager = ResultManager(log=log)
    
    if traits is None:
        traits_to_form_string = ['"trait_{}"'.format(i+1) for i in range(nstudy)]
    else:
        traits_to_form_string = ['"{}"'.format(i) for i in traits]
    
    hyprcoloc_res_combined = pd.DataFrame()

    # Get loci to process
    if loci is None:
        log.write(" -Loci were not provided. All significant loci will be automatically extracted...", verbose=verbose)
        sig_df = _get_sig(sumstats_multi, variant_id="SNPID", chrom="CHR", pos="POS", p="P_MIN", build=build)
    else:
        sig_df = sumstats_multi.loc[sumstats_multi["SNPID"].isin(loci), :]

    # Process each locus
    for index, row in sig_df.iterrows():
        gc.collect()  # Garbage collection for memory management
        
        # Extract locus
        locus = row["SNPID"]
        log.write(" -Running hyprcoloc for locus : {}...".format(locus), verbose=verbose)

        # Prepare input files
        output_beta_cols = []
        output_se_cols = []

        for i in range(nstudy):
            output_beta_cols.append("BETA_{}".format(i+1))
            output_se_cols.append("SE_{}".format(i+1))

        matched_sumstats = _extract_variants_in_locus(
            sumstats_multi, 
            windowsizekb, 
            locus=(row["CHR"], row["POS"])
        )
        
        to_export = matched_sumstats[["SNPID"] + output_se_cols + output_beta_cols].dropna()
        
        if len(to_export) == 0:
            log.write(" -No shared variants in locus {}...skipping".format(locus), verbose=verbose)
            continue
        
        log.write(" -Number of shared variants in locus {} : {}...".format(locus, len(to_export)), verbose=verbose)
        
        # Prepare output file paths
        beta_file = "{}_{}_beta_cols.tsv.gz".format(study, locus)
        se_file = "{}_{}_se_cols.tsv.gz".format(study, locus)
        output_path = "{}_{}_{}studies.res".format(study, locus, nstudy)
        output_prefix = "{}_{}_{}studies".format(study, locus, nstudy)
        
        # Write input files
        to_export[["SNPID"] + output_beta_cols].to_csv(beta_file, index=None, sep="\t")
        to_export[["SNPID"] + output_se_cols].to_csv(se_file, index=None, sep="\t")
        
        # Generate R script
        rscript = '''
library(hyprcoloc)

betas<-read.csv("{beta_file}",row.names = 1,sep="\\t")
ses  <-read.csv("{se_file}",row.names = 1,sep="\\t")

betas <- as.matrix(betas)
ses <- as.matrix(ses)

traits <- c({traits_string})
snpid <- rownames(betas)

res <- hyprcoloc(betas, 
          ses, 
          trait.names = traits, 
          snp.id = snpid,
          snpscore=TRUE)

write.csv(res[[1]], "{output_path}",row.names = FALSE) 
        '''.format(
            beta_file=beta_file,
            se_file=se_file,
            output_path=output_path,
            traits_string=','.join(traits_to_form_string)
        )
        
        # Execute R script using new framework
        result = runner.execute(
            script_content=rscript,
            expected_outputs=[output_path],
            temp_prefix=f"hyprcoloc_{study}_{locus}",
            timeout=timeout,
            verbose=verbose,
            working_dir="."  # Use current directory for input/output files
        )
        
        # Trace the result using ResultManager (automatically tracks success/failure)
        record = result_manager.trace(
            result=result,
            identifier=locus,
            parameters={
                "study": study,
                "locus": locus,
                "nstudy": nstudy,
                "windowsizekb": windowsizekb
            },
            read_outputs=True,
            expected_files=[output_path]
        )
        
        # Process successful results
        if result.success:
            # Get the output from cached data (read by trace with read_outputs=True)
            if output_path in record.output_data and record.output_data[output_path] is not None:
                hyprcoloc_res = record.output_data[output_path]
                hyprcoloc_res["PREFIX"] = output_prefix
                hyprcoloc_res_combined = pd.concat(
                    [hyprcoloc_res_combined, hyprcoloc_res],
                    ignore_index=True
                )
                log.write(" -Successfully processed results for locus {}".format(locus), verbose=verbose)
                
                # Show preview of results
                result_manager.preview_dataframe(output_path, result=result, n_rows=3)
            else:
                log.warning("  -Could not read output file for locus {}".format(locus), verbose=verbose)
        else:
            # Error already logged by trace(), stop if requested
            if stop_on_error:
                log.write("  -Stopping execution due to error (stop_on_error=True)", verbose=verbose)
                break
        
        log.write(" -Finishing hyprcoloc for locus : {}...".format(locus), verbose=verbose)
    
    # Log summary statistics
    stats = result_manager.get_statistics()
    log.write(" -Execution summary: {} successful, {} failed out of {} total".format(
        stats["successful_executions"],
        stats["failed_executions"],
        stats["total_executions"]
    ), verbose=verbose)
    
    log.write("Finished clocalization using hyprcoloc.", verbose=verbose)
    return hyprcoloc_res_combined
