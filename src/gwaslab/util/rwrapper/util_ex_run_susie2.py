from typing import TYPE_CHECKING, Optional, Union, Dict, Tuple
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_susie_version
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.util.rwrapper.util_ex_r_runner import RScriptRunner, RExecutionResult
from gwaslab.util.general.util_ex_result_manager import ResultManager
from gwaslab.util.general.util_path_manager import _path

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats


def _prepare_paths(
    row: pd.Series,
    sumstats: str,
    study: str,
    output_dir: Optional[str],
    log: Log,
    verbose: bool
) -> Dict[str, str]:
    """
    Prepare all file paths needed for SuSieR execution.
    
    Args:
        row: Row from filelist DataFrame
        sumstats: Path to locus sumstats file
        study: Study name
        output_dir: Output directory (if None, uses sumstats directory)
        log: Log instance
        verbose: Verbose logging
    
    Returns:
        Dict with keys: output_prefix, pipcs_file, diagnostic_file, 
                       pipcs_basename, diagnostic_basename, working_dir
    """
    # Determine output directory
    working_dir = output_dir if output_dir is not None else (
        os.path.dirname(sumstats) if os.path.dirname(sumstats) else "./"
    )
    
    # Generate output prefix (with or without base_name from sumstats)
    base_name = os.path.basename(sumstats).replace(".sumstats.gz", "").replace(".sumstats", "")
    if base_name:
        # base_name already contains study and locus information, so use it directly
        # and only add the analysis suffix
        output_prefix = _path(
            out=f"{base_name}_susie",
            directory=working_dir,
            log=log,
            verbose=False
        )
    else:
        output_prefix = _path(
            study=study,
            snpid=row["SNPID"],
            analysis="susie",
            directory=working_dir,
            log=log,
            verbose=False
        )
    
    # Generate full paths for output files
    pipcs_file_full = _path(
        out=output_prefix,
        suffix="pipcs",
        log=log,
        verbose=False
    )
    diagnostic_file_full = _path(
        out=output_prefix,
        suffix="png",
        log=log,
        verbose=False
    )
    
    # Get basenames for R script (files written relative to working_dir)
    pipcs_basename = os.path.basename(pipcs_file_full)
    diagnostic_basename = os.path.basename(diagnostic_file_full)
    
    return {
        "output_prefix": output_prefix,
        "pipcs_file": pipcs_file_full,
        "diagnostic_file": diagnostic_file_full,
        "pipcs_basename": pipcs_basename,
        "diagnostic_basename": diagnostic_basename,
        "working_dir": working_dir
    }


def _build_r_script(
    sumstats: str,
    ld_r_matrix: str,
    row: pd.Series,
    mode: str,
    n: Optional[Union[int, str]],
    max_iter: int,
    min_abs_corr: float,
    refine: str,
    L: int,
    susie_kwargs: str,
    fillldna: bool,
    pipcs_basename: str,
    diagnostic_basename: str
) -> Tuple[str, str]:
    """
    Build R script content for SuSieR execution.
    
    Args:
        sumstats: Path to sumstats file
        ld_r_matrix: Path to LD matrix file
        row: Row from filelist DataFrame
        mode: Mode for susie_rss ("z" or "bs")
        n: Sample size (or None to use mean from sumstats)
        max_iter: Maximum iterations
        min_abs_corr: Minimum absolute correlation
        refine: Refine parameter
        L: Number of effects
        susie_kwargs: Additional kwargs for susie_rss
        fillldna: Fill NA values in LD matrix with 0
        pipcs_basename: Basename for PIPCS output file
        diagnostic_basename: Basename for diagnostic image
    
    Returns:
        Tuple of (rscript_content, susie_rss_call_string)
    """
    # Determine input format based on mode
    if mode == "z":
        input_params = "z= sumstats$Z,"
        kriging_input = "sumstats$Z"
    else:
        input_params = "bhat = sumstats$BETA,shat = sumstats$SE"
        kriging_input = "sumstats$BETA/sumstats$SE"
    
    # Build susie_rss call string for logging
    susie_rss_call = "susie_rss({}, n = {}, R = R, max_iter = {}, min_abs_corr={}, refine = {}, L = {}{})".format(
        input_params,
        n if n is not None else "n",
        max_iter,
        min_abs_corr,
        refine,
        L,
        susie_kwargs
    )
    
    # Build full R script
    rscript = '''
library(susieR)

sumstats <- read.csv("{}", sep="\t")

R <- as.matrix(read.csv("{}", sep="\t", header=FALSE))
{}

n <- floor(mean(sumstats$N))

fitted_rss1 <- susie_rss({}, n = {}, R = R, max_iter = {}, min_abs_corr={}, refine = {}, L = {}{})

susie_fitted_summary <- summary(fitted_rss1)

output <- susie_fitted_summary$vars
output$SNPID <- sumstats$SNPID[susie_fitted_summary$vars$variable]
output$LOCUS <- "{}"
output$STUDY <- "{}"

write.csv(output, "{}", row.names = FALSE)

png(filename="{}")
diagnostic <- kriging_rss({}, R, n=n)
diagnostic$plot
dev.off()
    '''.format(
        sumstats,
        ld_r_matrix,
        "R[is.na(R)] <- 0" if fillldna else "",
        input_params,
        n if n is not None else "n",
        max_iter,
        min_abs_corr,
        refine,
        L,
        susie_kwargs,
        row["SNPID"],
        row["STUDY"],
        pipcs_basename,
        diagnostic_basename,
        kriging_input
    )
    
    return rscript, susie_rss_call


def _process_susie_result(
    result: RExecutionResult,
    result_manager: ResultManager,
    row: pd.Series,
    paths: Dict[str, str],
    delete: bool,
    show_diagnostic: bool,
    log: Log,
    verbose: bool
) -> Optional[pd.DataFrame]:
    """
    Process successful SuSieR execution result.
    
    Args:
        result: Execution result from R script
        result_manager: ResultManager instance
        row: Row from filelist DataFrame
        paths: Path dictionary from _prepare_paths
        delete: Whether to delete output files after reading
        show_diagnostic: Whether to display diagnostic image
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
    pip_cs["LOCUS"] = row["SNPID"]
    pip_cs["STUDY"] = row["STUDY"]
    
    # Show preview of results
    result_manager.preview_dataframe(pipcs_file, result=result, n_rows=3)
    
    # Handle file deletion
    if delete:
        os.remove(pipcs_file)
        log.write("  -Removed output file: {}".format(pipcs_file))
    else:
        log.write("  -SuSieR result summary to: {}".format(pipcs_file))
    
    # Display diagnostic image if requested
    if show_diagnostic:
        result_manager.show_image(
            paths["diagnostic_file"],
            result=result,
            title=f"SuSieR Diagnostic Plot - {row['SNPID']} ({row['STUDY']})"
        )
    
    return pip_cs


def _create_execution_log(
    result_manager: ResultManager,
    result: RExecutionResult,
    rscript: str,
    row: pd.Series,
    study: str,
    output_prefix: str,
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
        row: Row from filelist DataFrame
        study: Study name
        output_prefix: Output prefix
        working_dir: Working directory
        r_path: Path to Rscript executable
        log: Log instance
        verbose: Verbose logging
    """
    r_log_file = _path(
        study=study,
        snpid=row['SNPID'],
        analysis="susie",
        result_type="log",
        directory=working_dir,
        tmp=True,
        suffix="log",
        log=log,
        verbose=verbose
    )
    
    try:
        metadata = {
            "Study": study,
            "SNPID": row['SNPID'],
            "Output prefix": output_prefix
        }
        result_manager.create_r_log(
            result=result,
            script_content=rscript,
            log_file_path=r_log_file,
            working_dir=working_dir,
            metadata=metadata,
            r_path=r_path,
            package_name="susieR",
            verbose=verbose
        )
    except Exception as e:
        log.warning("  -Could not save R log to file: {}".format(str(e)), verbose=verbose)


@with_logging(
    start_to_msg="run finemapping using SuSieR from command line",
    finished_msg="running finemapping using SuSieR from command line",
    start_cols=["SNPID", "CHR", "POS"],
    start_function=".run_susie_rss()"
)
def _run_susie_rss(
    gls: 'Sumstats',
    filepath: Optional[str],
    r: str = "Rscript",
    mode: str = "bs",
    out: Optional[str] = None,
    max_iter: int = 100,
    min_abs_corr: float = 0.5,
    refine: str = "FALSE",
    L: int = 10,
    fillldna: bool = True,
    n: Optional[Union[int, str]] = None,
    delete: bool = False,
    susie_kwargs: str = "",
    log: Log = Log(),
    verbose: bool = True,
    timeout: Optional[float] = None,
    show_diagnostic: bool = True
) -> pd.DataFrame:
    """
    Run finemapping using SuSieR from command line using the RScriptRunner framework.
    
    Args:
        gls: Sumstats object
        filepath: Path to filelist containing study information
        r: Path to Rscript executable (default: "Rscript")
        mode: Mode for susie_rss ("z" or "bs", default: "bs")
        out: Output directory (default: None, uses sumstats directory)
        max_iter: Maximum iterations for susie_rss (default: 100)
        min_abs_corr: Minimum absolute correlation (default: 0.5)
        refine: Refine parameter (default: "FALSE")
        L: Number of effects (default: 10)
        fillldna: Fill NA values in LD matrix with 0 (default: True)
        n: Sample size (default: None, uses mean from sumstats)
        delete: Delete output files after reading (default: False)
        susie_kwargs: Additional kwargs for susie_rss (default: "")
        log: Log instance (default: new Log)
        verbose: Verbose logging (default: True)
        timeout: Timeout for R script execution in seconds (default: None)
        show_diagnostic: Display diagnostic image using matplotlib if found (default: True)
    
    Returns:
        DataFrame with finemapping results
    """
    # Early return if no filepath
    if filepath is None:
        log.write(" -File path is None.")
        log.write("Finished finemapping using SuSieR.")
        return pd.DataFrame()
    
    gls.offload()
    
    filelist = pd.read_csv(filepath, sep="\t")
    locus_pip_cs = pd.DataFrame()
    
    # Check R and SuSieR versions
    log = _checking_r_version(r, log)
    log = _check_susie_version(r, log)
    
    # Initialize runners and managers
    runner = RScriptRunner(
        r=r,
        log=log,
        timeout=timeout,
        temp_dir=out,
        cleanup=True
    )
    result_manager = ResultManager(log=log)
    
    # Process each locus
    for index, row in filelist.iterrows():
        gc.collect()
        
        study = row["STUDY"]
        ld_r_matrix = row["LD_R_MATRIX"]
        sumstats = row["LOCUS_SUMSTATS"]
        
        # Prepare paths
        paths = _prepare_paths(row, sumstats, study, out, log, verbose)
        
        log.write(" -Running for: {} - {}".format(row["SNPID"], row["STUDY"]))
        log.write("  -Locus sumstats:{}".format(sumstats))
        log.write("  -LD r matrix:{}".format(ld_r_matrix))
        log.write("  -output_prefix:{}".format(paths["output_prefix"]))
        
        # Build R script
        rscript, susie_rss_call = _build_r_script(
            sumstats=sumstats,
            ld_r_matrix=ld_r_matrix,
            row=row,
            mode=mode,
            n=n,
            max_iter=max_iter,
            min_abs_corr=min_abs_corr,
            refine=refine,
            L=L,
            susie_kwargs=susie_kwargs,
            fillldna=fillldna,
            pipcs_basename=paths["pipcs_basename"],
            diagnostic_basename=paths["diagnostic_basename"]
        )
        
        log.write("  -SuSieR script: {}".format(susie_rss_call))
        
        # Execute R script
        log.write("  -Running SuSieR from command line...")
        result = runner.execute(
            script_content=rscript,
            expected_outputs=[paths["pipcs_file"], paths["diagnostic_file"]],
            temp_prefix=f"susie_{study}_{row['SNPID']}",
            temp_suffix=".R",
            timeout=timeout,
            verbose=False,
            working_dir=paths["working_dir"]
        )
        
        # Trace result
        result_manager.trace(
            result=result,
            identifier=row['SNPID'],
            parameters={
                "study": study,
                "locus": row["SNPID"],
                "output_prefix": paths["output_prefix"]
            }
        )
        
        # Create execution log
        _create_execution_log(
            result_manager=result_manager,
            result=result,
            rscript=rscript,
            row=row,
            study=study,
            output_prefix=paths["output_prefix"],
            working_dir=paths["working_dir"],
            r_path=r,
            log=log,
            verbose=verbose
        )
        
        # Process results
        pip_cs = _process_susie_result(
            result=result,
            result_manager=result_manager,
            row=row,
            paths=paths,
            delete=delete,
            show_diagnostic=show_diagnostic,
            log=log,
            verbose=verbose
        )
        
        if pip_cs is not None:
            locus_pip_cs = pd.concat([locus_pip_cs, pip_cs], ignore_index=True)
    
    # Finalize results
    gls.reload()
    
    if not locus_pip_cs.empty:
        locus_pip_cs = locus_pip_cs.rename(columns={
            "variable": "N_SNP",
            "variable_prob": "PIP",
            "cs": "CREDIBLE_SET_INDEX"
        })
        locus_pip_cs = pd.merge(
            locus_pip_cs,
            gls.data[["SNPID", "CHR", "POS"]],
            on="SNPID",
            how="left"
        )
    
    return locus_pip_cs


def _get_cs_lead(pipcs: pd.DataFrame) -> pd.DataFrame:
    """
    Extract lead variants from credible sets.
    
    Args:
        pipcs: DataFrame with PIPCS results
    
    Returns:
        DataFrame with lead variants (highest PIP per credible set)
    """
    leads = pipcs.loc[pipcs["CREDIBLE_SET_INDEX"] > 0, :]
    leads = leads.sort_values(by="PIP", ascending=False).drop_duplicates(
        subset=["STUDY", "LOCUS", "CREDIBLE_SET_INDEX"]
    )
    return leads
