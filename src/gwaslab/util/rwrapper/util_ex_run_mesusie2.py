from typing import TYPE_CHECKING, Optional, Tuple, List, Dict
import os
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.util.rwrapper.util_ex_r_runner import RScriptRunner, RExecutionResult
from gwaslab.util.general.util_ex_result_manager import ResultManager
from gwaslab.util.general.util_path_manager import _path

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats


def _prepare_paths(
    filelist: pd.DataFrame,
    output_dir: Optional[str],
    log: Log,
    verbose: bool
) -> Dict[str, str]:
    """
    Prepare all file paths needed for MESuSiE execution.
    
    Args:
        filelist: DataFrame with study information
        output_dir: Output directory (if None, uses first sumstats directory)
        log: Log instance
        verbose: Verbose logging
    
    Returns:
        Dict with keys: output_prefix, pipcs_file, rds_file, diagnostic_files,
                       stacked_plot_file, working_dir, and basenames
    """
    # Determine output directory from first sumstats file
    first_sumstats = filelist.iloc[0]["LOCUS_SUMSTATS"]
    working_dir = output_dir if output_dir is not None else (
        os.path.dirname(first_sumstats) if os.path.dirname(first_sumstats) else "./"
    )
    
    # Get group and locus from first row
    group = filelist.iloc[0]["GROUP"]
    locus = filelist.iloc[0]["LOCUS"]
    
    # Generate output prefix using out parameter (since _path doesn't have 'group' parameter)
    output_prefix = _path(
        out=f"{group}_{locus}_mesusie",
        directory=working_dir,
        log=log,
        verbose=False
    )
    
    # Generate full paths for output files
    pipcs_file = _path(
        out=output_prefix,
        suffix="pipcs",
        log=log,
        verbose=False
    )
    rds_file = _path(
        out=output_prefix,
        suffix="rds",
        log=log,
        verbose=False
    )
    stacked_plot_file = _path(
        out=output_prefix,
        suffix="stacked_regions.png",
        log=log,
        verbose=False
    )
    cscs_index_file = _path(
        out=output_prefix,
        suffix="cscs_index",
        log=log,
        verbose=False
    )
    cspurity_file = _path(
        out=output_prefix,
        suffix="cspurity",
        log=log,
        verbose=False
    )
    cscs_category_file = _path(
        out=output_prefix,
        suffix="cscs_category",
        log=log,
        verbose=False
    )
    cscs_i_file = _path(
        out=output_prefix,
        suffix="cscs_i",
        log=log,
        verbose=False
    )
    cscs_snpid_file = _path(
        out=output_prefix,
        suffix="cscs_snpid",
        log=log,
        verbose=False
    )
    
    # Diagnostic files (one per study)
    diagnostic_files = []
    for index, row in filelist.iterrows():
        diagnostic_file = _path(
            out=f"{group}_{locus}_mesusie_diagnostic_{index}",
            suffix="png",
            directory=working_dir,
            log=log,
            verbose=False
        )
        diagnostic_files.append(diagnostic_file)
    
    # Get basenames for R script (files written relative to working_dir)
    pipcs_basename = os.path.basename(pipcs_file)
    rds_basename = os.path.basename(rds_file)
    stacked_plot_basename = os.path.basename(stacked_plot_file)
    cscs_index_basename = os.path.basename(cscs_index_file)
    cspurity_basename = os.path.basename(cspurity_file)
    cscs_category_basename = os.path.basename(cscs_category_file)
    cscs_i_basename = os.path.basename(cscs_i_file)
    cscs_snpid_basename = os.path.basename(cscs_snpid_file)
    diagnostic_basenames = [os.path.basename(f) for f in diagnostic_files]
    
    return {
        "output_prefix": output_prefix,
        "pipcs_file": pipcs_file,
        "rds_file": rds_file,
        "stacked_plot_file": stacked_plot_file,
        "cscs_index_file": cscs_index_file,
        "cspurity_file": cspurity_file,
        "cscs_category_file": cscs_category_file,
        "cscs_i_file": cscs_i_file,
        "cscs_snpid_file": cscs_snpid_file,
        "diagnostic_files": diagnostic_files,
        "pipcs_basename": pipcs_basename,
        "rds_basename": rds_basename,
        "stacked_plot_basename": stacked_plot_basename,
        "cscs_index_basename": cscs_index_basename,
        "cspurity_basename": cspurity_basename,
        "cscs_category_basename": cscs_category_basename,
        "cscs_i_basename": cscs_i_basename,
        "cscs_snpid_basename": cscs_snpid_basename,
        "diagnostic_basenames": diagnostic_basenames,
        "working_dir": working_dir,
        "group": group,
        "locus": locus
    }


def _build_r_script(
    filelist: pd.DataFrame,
    ns: List[int],
    fillldna: bool,
    paths: Dict[str, str],
    L: int = 10,
    mesusie_kwargs: str = ""
) -> Tuple[str, str]:
    """
    Build R script content for MESuSiE execution.
    
    Args:
        filelist: DataFrame with study information
        ns: List of sample sizes for each study
        fillldna: Fill NA values in LD matrix with 0
        paths: Path dictionary from _prepare_paths
        L: Number of effects (default: 10)
        mesusie_kwargs: Additional kwargs for meSuSie_core
    
    Returns:
        Tuple of (rscript_content, mesusie_call_string)
    """
    group = paths["group"]
    locus = paths["locus"]
    pipcs_basename = paths["pipcs_basename"]
    rds_basename = paths["rds_basename"]
    stacked_plot_basename = paths["stacked_plot_basename"]
    cscs_index_basename = paths["cscs_index_basename"]
    cspurity_basename = paths["cspurity_basename"]
    cscs_category_basename = paths["cscs_category_basename"]
    cscs_i_basename = paths["cscs_i_basename"]
    cscs_snpid_basename = paths["cscs_snpid_basename"]
    diagnostic_basenames = paths["diagnostic_basenames"]
    
    # Build MESuSiE call string for logging
    mesusie_call = "meSuSie_core(ld_list, summ_stat_list, L = {}{})".format(
        L,
        mesusie_kwargs
    )
    
    # Initialize R script
    r_script_init = '''
library(MESuSiE)
ld_list <- list()
summ_stat_list <- list()
    '''
    
    r_scripts_for_loading = [r_script_init]
    
    # Build loading section for each study
    study0 = None
    for index, row in filelist.iterrows():
        if index == 0:
            study0 = row["STUDY"]
        
        study = row["STUDY"]
        ld_r_matrix = row["LD_R_MATRIX"]
        sumstats = row["LOCUS_SUMSTATS"]
        n = ns[index] if index < len(ns) else ns[0]  # Fallback to first n if index out of range
        diagnostic_basename = diagnostic_basenames[index] if index < len(diagnostic_basenames) else f"diagnostic_{group}_{locus}_{index}.png"
        
        rscript = '''
sum{index} <- read.csv("{sumstats}", sep="\\t")
sum{index}$Z <- sum{index}$Beta/sum{index}$Se
sum{index}$N <- {n}
ld{index} <- read.csv("{ld_r_matrix}", sep="\\t", header=FALSE)
{fill_na}
names(ld{index}) <- sum{index}$SNP
ld_list${study} <- as.matrix(ld{index})
summ_stat_list${study} <- sum{index}

png(filename="{diagnostic_basename}")
diagnostic <- kriging_rss(summ_stat_list${study}$Z, ld_list${study})
diagnostic$plot
dev.off()
        '''.format(
            index=index,
            study=study,
            n=n,
            sumstats=sumstats,
            ld_r_matrix=ld_r_matrix,
            fill_na="ld{index}[is.na(ld{index})] <- 0".format(index=index) if fillldna else "",
            diagnostic_basename=diagnostic_basename
        )
        r_scripts_for_loading.append(rscript)
    
    rscript_loading = "".join(r_scripts_for_loading)
    
    # Build computing section
    rscript_computing = '''
MESuSiE_res <- {mesusie_call}
    '''.format(mesusie_call=mesusie_call)
    
    # Build output section
    rscript_output = '''
saveRDS(MESuSiE_res, file = "{rds_basename}")
pips <- cbind(summ_stat_list${study0}$SNP, summ_stat_list${study0}$CHR, summ_stat_list${study0}$POS, MESuSiE_res$pip_config)
colnames(pips)[1] <- "SNPID"
colnames(pips)[2] <- "CHR"
colnames(pips)[3] <- "POS"
pips <- data.frame(pips)
pips[c("CREDIBLE_SET_INDEX")] <- 0 
pips[c("CS_CATEGORY")] <- NA
for (i in 1:length(MESuSiE_res$cs$cs)) {{
    pips[MESuSiE_res$cs$cs[[i]], c("CREDIBLE_SET_INDEX")] <- i
    pips[MESuSiE_res$cs$cs[[i]], c("CS_CATEGORY")] <- MESuSiE_res$cs$cs_category[[i]]
}}
write.csv(pips, "{pipcs_basename}", row.names = FALSE)

write.csv(MESuSiE_res$cs$cs_index, "{cscs_index_basename}", row.names = FALSE)
write.csv(MESuSiE_res$cs$purity, "{cspurity_basename}", row.names = FALSE)
write.csv(MESuSiE_res$cs$cs_category, "{cscs_category_basename}", row.names = FALSE)

for (p in MESuSiE_res$cs$cs) {{
  write(p, "{cscs_i_basename}", append=TRUE, sep="\\t", ncolumns=10000000)
  write(summ_stat_list${study0}$SNP[p], "{cscs_snpid_basename}", append=TRUE, sep="\\t", ncolumns=10000000)
}}
    '''.format(
        rds_basename=rds_basename,
        study0=study0,
        pipcs_basename=pipcs_basename,
        cscs_index_basename=cscs_index_basename,
        cspurity_basename=cspurity_basename,
        cscs_category_basename=cscs_category_basename,
        cscs_i_basename=cscs_i_basename,
        cscs_snpid_basename=cscs_snpid_basename
    )
    
    # Build plotting section
    rscript_plotting = '''
png(filename="{stacked_plot_basename}")
MESuSiE_Plot(MESuSiE_res, ld_list, summ_stat_list)
dev.off()
    '''.format(stacked_plot_basename=stacked_plot_basename)
    
    rscript = rscript_loading + rscript_computing + rscript_output + rscript_plotting
    
    return rscript, mesusie_call


def _process_mesusie_result(
    result: RExecutionResult,
    result_manager: ResultManager,
    filelist: pd.DataFrame,
    paths: Dict[str, str],
    delete: bool,
    show_diagnostic: bool,
    log: Log,
    verbose: bool
) -> Optional[pd.DataFrame]:
    """
    Process successful MESuSiE execution result.
    
    Args:
        result: Execution result from R script
        result_manager: ResultManager instance
        filelist: DataFrame with study information
        paths: Path dictionary from _prepare_paths
        delete: Whether to delete output files after reading
        show_diagnostic: Whether to display diagnostic images
        log: Log instance
        verbose: Verbose logging
    
    Returns:
        DataFrame with results, or None if processing failed
    """
    if not result.success:
        return None
    
    pipcs_file = paths["pipcs_file"]
    
    # Read and process output file
    if pipcs_file not in result.output_files:
        log.warning("  -Expected output file not found: {}".format(pipcs_file), verbose=verbose)
        return None
    
    pipcs_path = result.output_files[pipcs_file]
    pip_cs = pd.read_csv(pipcs_path)
    
    # Add locus information from first row
    pip_cs["LOCUS"] = filelist.iloc[0]["LOCUS"]
    pip_cs["GROUP"] = filelist.iloc[0]["GROUP"]
    
    # Show preview of results
    result_manager.preview_dataframe(pipcs_file, result=result, n_rows=3)
    
    # Handle file deletion
    if delete:
        os.remove(pipcs_file)
        log.write("  -Removed output file: {}".format(pipcs_file), verbose=verbose)
    else:
        log.write("  -MESuSiE result summary to: {}".format(pipcs_file), verbose=verbose)
    
    # Display diagnostic images if requested
    if show_diagnostic:
        for i, diagnostic_file in enumerate(paths["diagnostic_files"]):
            if diagnostic_file in result.output_files:
                result_manager.show_image(
                    diagnostic_file,
                    result=result,
                    title=f"MESuSiE Diagnostic Plot - {paths['locus']} (Study {i+1})"
                )
        
        # Also show stacked plot if available
        if paths["stacked_plot_file"] in result.output_files:
            result_manager.show_image(
                paths["stacked_plot_file"],
                result=result,
                title=f"MESuSiE Stacked Regions - {paths['locus']}"
            )
    
    return pip_cs


def _create_execution_log(
    result_manager: ResultManager,
    result: RExecutionResult,
    rscript: str,
    filelist: pd.DataFrame,
    paths: Dict[str, str],
    working_dir: str,
    r_path: str,
    log: Log,
    verbose: bool
) -> None:
    """
    Create log file for R script execution.
    
    Args:
        result_manager: ResultManager instance
        result: Execution result
        rscript: R script content
        filelist: DataFrame with study information
        paths: Path dictionary from _prepare_paths
        working_dir: Working directory
        r_path: Path to Rscript executable
        log: Log instance
        verbose: Verbose logging
    """
    r_log_file = _path(
        out=f"{paths['group']}_{paths['locus']}_mesusie",
        result_type="log",
        directory=working_dir,
        tmp=True,
        suffix="log",
        log=log,
        verbose=verbose
    )
    
    try:
        metadata = {
            "Group": paths["group"],
            "Locus": paths["locus"],
            "N_studies": len(filelist),
            "Studies": ", ".join(filelist["STUDY"].unique().tolist()),
            "Output prefix": paths["output_prefix"]
        }
        result_manager.create_r_log(
            result=result,
            script_content=rscript,
            log_file_path=r_log_file,
            working_dir=working_dir,
            metadata=metadata,
            r_path=r_path,
            package_name="MESuSiE",
            verbose=verbose
        )
    except Exception as e:
        log.warning("  -Could not save R log to file: {}".format(str(e)), verbose=verbose)


@with_logging(
    start_to_msg="run cross-ancestry finemapping using MESuSiE from command line",
    finished_msg="running cross-ancestry finemapping using MESuSiE from command line",
    start_cols=["SNPID", "CHR", "POS"],
    start_function=".run_mesusie()",
    check_tools=["r", "r:MESuSiE"]
)
def _run_mesusie(
    filepath: Optional[str],
    r: str = "Rscript",
    ns: Optional[List[int]] = None,
    fillldna: bool = True,
    delete: bool = False,
    ncols: Optional[List[int]] = None,
    log: Log = Log(),
    verbose: bool = True,
    timeout: Optional[float] = None,
    show_diagnostic: bool = True,
    L: int = 10,
    mesusie_kwargs: str = "",
    out: Optional[str] = None,
    reference_data: Optional[pd.DataFrame] = None,
    group_name: Optional[str] = None,
    study_name: Optional[str] = None
) -> pd.DataFrame:
    """
    Run cross-ancestry finemapping using MESuSiE from command line using the RScriptRunner framework.
    
    Args:
        filepath: Path to filelist containing study information
        r: Path to Rscript executable (default: "Rscript")
        ns: List of sample sizes for each study (default: None, uses ncols if provided)
        fillldna: Fill NA values in LD matrix with 0 (default: True)
        delete: Delete output files after reading (default: False)
        ncols: List of sample sizes (alternative to ns, kept for compatibility)
        log: Log instance (default: new Log)
        verbose: Verbose logging (default: True)
        timeout: Timeout for R script execution in seconds (default: None)
        show_diagnostic: Display diagnostic images using matplotlib if found (default: True)
        L: Number of effects (default: 10)
        mesusie_kwargs: Additional kwargs for meSuSie_core (default: "")
        out: Output directory (default: None, uses sumstats directory)
        reference_data: Optional DataFrame with SNPID, CHR, POS for merging (default: None)
        group_name: Optional group name to add to results (default: None)
        study_name: Optional study name to add to results (default: None)
    
    Returns:
        DataFrame with finemapping results (processed with CHR/POS merged and GROUP/STUDY set if provided)
    """
    # Log citation
    log.write(" -Citation: Gao, B., & Zhou, X. (2024). MESuSiE enables scalable and powerful multi-ancestry fine-mapping of causal variants in genome-wide association studies. Nature genetics, 56(1), 170-179.", verbose=verbose)
    
    # Early return if no filepath
    if filepath is None:
        log.write(" -File path is None.", verbose=verbose)
        log.write("Finished finemapping using MESuSiE.", verbose=verbose)
        return pd.DataFrame()
    
    filelist = pd.read_csv(filepath, sep="\t")
    
    if len(filelist) == 0:
        log.write(" -Filelist is empty", verbose=verbose)
        return pd.DataFrame()
    
    # Get unique loci
    unique_loci = filelist["LOCUS"].unique() if "LOCUS" in filelist.columns else filelist["SNPID"].unique()
    log.write(" -Found {} unique locus/loci".format(len(unique_loci)), verbose=verbose)
    
    # Get number of studies
    if "STUDY" in filelist.columns:
        nstudy = len(filelist["STUDY"].unique())
    elif "SUBSTUDY" in filelist.columns:
        nstudy = filelist["SUBSTUDY"].max()
    else:
        raise ValueError("Cannot determine number of studies from filelist")
    
    log.write(" -Number of studies: {}".format(nstudy), verbose=verbose)
    
    # Initialize runners and managers (reused for all loci)
    runner = RScriptRunner(
        r=r,
        log=log,
        timeout=timeout,
        temp_dir=out,
        cleanup=True
    )
    result_manager = ResultManager(log=log)
    
    all_results = []
    
    # Process each locus
    for locus_idx, locus in enumerate(unique_loci):
        log.write(" -Processing locus {} ({}/{})...".format(locus, locus_idx + 1, len(unique_loci)), 
                 verbose=verbose)
        
        # Get rows for this locus
        if "LOCUS" in filelist.columns:
            locus_rows = filelist[filelist["LOCUS"] == locus]
        else:
            locus_rows = filelist[filelist["SNPID"] == locus]
        
        if len(locus_rows) != nstudy:
            log.warning("  -Locus {} has {} rows, expected {}. Skipping...".format(
                locus, len(locus_rows), nstudy
            ), verbose=verbose)
            continue
        
        # Sort by SUBSTUDY or STUDY to ensure consistent order
        if "SUBSTUDY" in locus_rows.columns:
            locus_rows = locus_rows.sort_values("SUBSTUDY")
        elif "STUDY" in locus_rows.columns:
            locus_rows = locus_rows.sort_values("STUDY")
        
        # Reset index for consistent indexing
        locus_rows = locus_rows.reset_index(drop=True)
        
        # Determine sample sizes for this locus
        locus_ns = ns
        if locus_ns is None:
            if ncols is not None:
                locus_ns = ncols
            else:
                log.warning("  -No sample sizes provided (ns or ncols), using default values", verbose=verbose)
                locus_ns = [10000] * len(locus_rows)  # Default fallback
        
        # Ensure ns matches number of studies
        if len(locus_ns) != len(locus_rows):
            if len(locus_ns) == 1:
                locus_ns = locus_ns * len(locus_rows)
            else:
                log.warning("  -Sample sizes length ({}) doesn't match number of studies ({}). Using first {} values...".format(
                    len(locus_ns), len(locus_rows), len(locus_rows)
                ), verbose=verbose)
                locus_ns = locus_ns[:len(locus_rows)] if len(locus_ns) >= len(locus_rows) else locus_ns + [locus_ns[0]] * (len(locus_rows) - len(locus_ns))
        
        log.write("  -Ns: {}".format(", ".join(map(str, locus_ns))), verbose=verbose)
        
        # Prepare paths for this locus
        paths = _prepare_paths(locus_rows, out, log, verbose)
        
        log.write("  -Running for group: {} - locus: {}".format(paths["group"], paths["locus"]), verbose=verbose)
        log.write("  -Number of studies: {}".format(len(locus_rows)), verbose=verbose)
        for index, row in locus_rows.iterrows():
            log.write("  -Study {}: {} - {}".format(index, row["STUDY"], row["SNPID"]), verbose=verbose)
            log.write("    -Locus sumstats: {}".format(row["LOCUS_SUMSTATS"]), verbose=verbose)
            log.write("    -LD r matrix: {}".format(row["LD_R_MATRIX"]), verbose=verbose)
        
        log.write("  -Output prefix: {}".format(paths["output_prefix"]), verbose=verbose)
        
        # Build R script for this locus
        rscript, mesusie_call = _build_r_script(
            filelist=locus_rows,
            ns=locus_ns,
            fillldna=fillldna,
            paths=paths,
            L=L,
            mesusie_kwargs=mesusie_kwargs
        )
        
        log.write("  -MESuSiE script: {}".format(mesusie_call), verbose=verbose)
        
        # Collect all expected output files
        expected_outputs = [
            paths["pipcs_file"],
            paths["rds_file"],
            paths["stacked_plot_file"],
            paths["cscs_index_file"],
            paths["cspurity_file"],
            paths["cscs_category_file"],
            paths["cscs_i_file"],
            paths["cscs_snpid_file"]
        ] + paths["diagnostic_files"]
        
        # Execute R script
        log.write("  -Running MESuSiE from command line...", verbose=verbose)
        result = runner.execute(
            script_content=rscript,
            expected_outputs=expected_outputs,
            temp_prefix=f"mesusie_{paths['group']}_{paths['locus']}",
            temp_suffix=".R",
            timeout=timeout,
            verbose=False,
            working_dir=paths["working_dir"]
        )
        
        # Trace result
        result_manager.trace(
            result=result,
            identifier=f"{paths['group']}_{paths['locus']}",
            parameters={
                "group": paths["group"],
                "locus": paths["locus"],
                "n_studies": len(locus_rows),
                "output_prefix": paths["output_prefix"]
            }
        )
        
        # Create execution log
        _create_execution_log(
            result_manager=result_manager,
            result=result,
            rscript=rscript,
            filelist=locus_rows,
            paths=paths,
            working_dir=paths["working_dir"],
            r_path=r,
            log=log,
            verbose=verbose
        )
        
        # Process results
        pip_cs = _process_mesusie_result(
            result=result,
            result_manager=result_manager,
            filelist=locus_rows,
            paths=paths,
            delete=delete,
            show_diagnostic=show_diagnostic,
            log=log,
            verbose=verbose
        )
        
        if pip_cs is None:
            log.warning("  -Failed to process MESuSiE results for locus {}".format(locus), verbose=verbose)
            continue
        
        all_results.append(pip_cs)
        log.write("  -Completed locus {}: {} variants".format(locus, len(pip_cs)), verbose=verbose)
    
    # Combine all results
    if len(all_results) == 0:
        log.write(" -No results generated", verbose=verbose)
        return pd.DataFrame()
    
    combined_results = pd.concat(all_results, ignore_index=True)
    log.write(" -Combined results from {} loci: {} total variants".format(
        len(all_results), len(combined_results)
    ), verbose=verbose)
    
    # Post-process results: merge CHR/POS, set GROUP/STUDY, rename columns
    if not combined_results.empty:
        # Merge CHR and POS from reference data if not present
        if reference_data is not None:
            if "CHR" not in combined_results.columns or "POS" not in combined_results.columns:
                log.write(" -Merging CHR and POS from reference data...", verbose=verbose)
                combined_results = pd.merge(
                    combined_results,
                    reference_data[["SNPID", "CHR", "POS"]],
                    on="SNPID",
                    how="left"
                )
        
        # Set GROUP and STUDY columns if provided (overrides filelist values if provided)
        if group_name is not None:
            combined_results["GROUP"] = group_name
        if study_name is not None:
            combined_results["STUDY"] = study_name
        
        # Rename columns if needed (for compatibility)
        rename_map = {
            "cs": "CREDIBLE_SET_INDEX",
            "variable_prob": "PIP",
            "variable": "N_SNP"
        }
        for old_col, new_col in rename_map.items():
            if old_col in combined_results.columns and new_col not in combined_results.columns:
                combined_results = combined_results.rename(columns={old_col: new_col})
    
    return combined_results
