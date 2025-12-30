"""
Decorator for downstream analysis functions with logging, validation, and tool checking.

This module provides a decorator that adds standardized logging, argument validation,
column checking, and tool availability/version checking to downstream analysis functions.
"""

from typing import TYPE_CHECKING, Optional, List, Tuple, Dict, Any, Callable, Union
import gc
import os
import pandas as pd
import sys
import shutil
import subprocess
import numpy as np
from functools import wraps
import inspect

from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import _get_version
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build
from gwaslab.extension import _checking_r_version
from gwaslab.extension import _check_tool_availability

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

# Constants
TOOL_CHECK_TIMEOUT = 5
SUPPORTED_SPECIAL_TOOLS = ["python", "r", "plink", "plink2", "tabix", "bcftools", "scdrs"]


###############################################################################################################
# Main Decorator
###############################################################################################################

def with_analysis_logging(start_to_msg: str, 
                          finished_msg: str,
                          start_function: Optional[str] = None,
                          start_cols: Optional[List[Union[str, Tuple[str, ...], List[str]]]] = None,
                          must_kwargs: Optional[List[str]] = None,
                          show_shape: bool = True,
                          check_tools: Optional[List[str]] = None
                          ) -> Callable:
    """
    Decorator to add standardized logging, argument checks, column validation,
    and tool availability checks to downstream analysis functions.

    Parameters
    ----------
    start_to_msg : str
        Message describing the operation start.
    finished_msg : str
        Message describing the operation completion.
    start_function : str or None, default None
        Function label used in logs when reporting missing columns or args.
    start_cols : list or None, default None
        Required columns to check in the input DataFrame prior to execution.
    must_kwargs : list or None, default None
        Argument names that must be provided (non-None) for the wrapped func.
    show_shape : bool, default True
        Log DataFrame shape before and after the function call when available.
    check_tools : list or None, default None
        List of tools to check for availability and version. Supported: "python", "r", "plink", "plink2",
        "tabix", "bcftools", "scdrs", or any other command-line tool name. 
        - For "r": checks the "r" argument from function call, or searches for "R"/"Rscript" in PATH
        - For R packages: use format "r:PackageName" (e.g., "r:TwoSampleMR", "r:susieR") to check if 
          the R package is installed and get its version
        - For "python": checks the "python" argument, or uses sys.executable, or searches for "python3"/"python" in PATH
        - For "plink"/"plink2": checks the "plink"/"plink2" argument, or searches PATH
        - For "tabix", "bcftools", "scdrs": checks the corresponding argument, or searches PATH
        - Other tools: searches PATH and checks version with --version flag
        All tools are checked for both availability and version information.

    Behavior
    --------
    - Logs references (threads, VCF/FASTA/TSV) and start/end messages.
    - Validates required args and columns; skips execution if columns missing.
    - Checks DataFrame shapes before and after execution.
    - Checks tool availability (Python, R, and other tools) before execution.
    - Preserves original function metadata via `functools.wraps`.
    """
    # Normalize parameters
    start_cols = start_cols or []
    must_kwargs = must_kwargs or []
    check_tools = check_tools or []
    start_function = start_function or "this step"

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            # Setup: bind arguments and extract common parameters
            bound_kwargs = _setup_function_call(func, args, kwargs)
            log = bound_kwargs.arguments.get('log', Log())
            verbose = bound_kwargs.arguments.get('verbose', True)
            
            # Extract sumstats from various possible sources
            sumstats = _extract_sumstats(bound_kwargs)
            
            # Log references and process build
            _log_references(bound_kwargs, log, verbose)
            _process_build_argument(bound_kwargs, log, verbose)
            
            # Log operation start
            log.log_operation_start(start_to_msg, version=_get_version(), verbose=verbose)
            
            # Validate tools, arguments, and columns
            if check_tools:
                check_tools_availability(check_tools, bound_kwargs, log, verbose)
            
            _validate_required_kwargs(must_kwargs, bound_kwargs, log, verbose, start_function)
            
            # Find Sumstats object for shape tracking
            sumstats_obj = _find_sumstats_object(args, bound_kwargs, log)
            if sumstats is None and sumstats_obj is not None:
                sumstats = _get_sumstats_dataframe(sumstats_obj)
            
            # Execute function with validation
            result = _execute_with_validation(
                func, args, kwargs, sumstats, sumstats_obj, 
                start_cols, start_function, show_shape, log, verbose
            )
            
            # Log operation finish
            log.log_operation_finish(finished_msg, verbose=verbose)
            return result
        
        return wrapper
    return decorator


###############################################################################################################
# Setup and Extraction Functions
###############################################################################################################

def _setup_function_call(func: Callable, args: Tuple[Any, ...], kwargs: Dict[str, Any]) -> inspect.BoundArguments:
    """Bind function arguments and apply defaults."""
    sig = inspect.signature(func)
    bound_kwargs = sig.bind(*args, **kwargs)
    bound_kwargs.apply_defaults()
    return bound_kwargs


def _extract_sumstats(bound_kwargs: inspect.BoundArguments) -> Optional[pd.DataFrame]:
    """Extract sumstats DataFrame from bound arguments."""
    insumstats = bound_kwargs.arguments.get('insumstats', None)
    if insumstats is not None:
        return insumstats
    return bound_kwargs.arguments.get('sumstats', None)


def _log_references(bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> None:
    """Log thread and reference file information."""
    threads = bound_kwargs.arguments.get('threads', None)
    if threads is not None:
        log.log_threads(threads, verbose=verbose)
    
    ref_vcf = bound_kwargs.arguments.get('ref_vcf', None)
    if ref_vcf is not None:
        log.log_reference_path("VCF", ref_vcf, verbose=verbose)
    
    ref_fasta = bound_kwargs.arguments.get('ref_fasta', None)
    if ref_fasta is not None:
        log.log_reference_path("FASTA", ref_fasta, verbose=verbose)
    
    ref_tsv = bound_kwargs.arguments.get('ref_tsv', None)
    if ref_tsv is not None:
        log.log_reference_path("TSV", ref_tsv, verbose=verbose)


def _process_build_argument(bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> None:
    """Process and normalize build argument if present."""
    if "build" in bound_kwargs.arguments:
        bound_kwargs.arguments["build"] = _process_build(
            bound_kwargs.arguments["build"], 
            log=log, 
            verbose=verbose
        )
        # Update args and kwargs after build processing
        bound_kwargs.args = bound_kwargs.args
        bound_kwargs.kwargs = bound_kwargs.kwargs


def _find_sumstats_object(args: Tuple[Any, ...], bound_kwargs: inspect.BoundArguments, log: Log) -> Optional['Sumstats']:
    """Find Sumstats object instance from various sources."""
    sumstats_obj = None
    try:
        from gwaslab.g_Sumstats import Sumstats
        
        # Check if first argument is a Sumstats instance
        if args and isinstance(args[0], Sumstats):
            return args[0]
        
        # Check if log has a reference to Sumstats object
        sumstats_obj = getattr(log, '_sumstats_obj', None)
        if sumstats_obj is not None:
            return sumstats_obj
        
        # Check if any kwarg is a Sumstats instance
        for value in bound_kwargs.arguments.values():
            if isinstance(value, Sumstats):
                return value
    except ImportError:
        # If import fails, try log reference as fallback
        sumstats_obj = getattr(log, '_sumstats_obj', None)
    
    return sumstats_obj


def _get_sumstats_dataframe(sumstats_obj: 'Sumstats') -> Optional[pd.DataFrame]:
    """Extract DataFrame from Sumstats object."""
    try:
        return sumstats_obj.data
    except AttributeError:
        return None


###############################################################################################################
# Validation Functions
###############################################################################################################

def _validate_required_kwargs(must_kwargs: List[str], bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool, start_function: str) -> None:
    """Validate that all required keyword arguments are provided and non-empty."""
    is_valid = True
    for key in must_kwargs:
        value = bound_kwargs.arguments.get(key, None)
        is_valid = is_valid and check_arg(log, verbose, key, value, start_function)
    
    if not is_valid:
        raise ValueError(f"{must_kwargs} must be provided.")


def _execute_with_validation(func: Callable, args: Tuple[Any, ...], kwargs: Dict[str, Any], 
                             sumstats: Optional[pd.DataFrame], sumstats_obj: Optional['Sumstats'], 
                             start_cols: List[Union[str, Tuple[str, ...], List[str]]], 
                             start_function: str, show_shape: bool, log: Log, verbose: bool) -> Any:
    """Execute function with column and shape validation."""
    initial_shape = None
    
    if sumstats is not None:
        # Check shape before execution
        if show_shape:
            initial_shape = (len(sumstats), len(sumstats.columns))
            check_dataframe_shape(
                sumstats=sumstats, 
                log=log, 
                verbose=verbose,
                sumstats_obj=sumstats_obj
            )
        
        # Check required columns
        is_enough_col = check_col(
            sumstats.columns, 
            verbose=verbose, 
            log=log, 
            cols=start_cols, 
            function=start_function
        )
        
        if not is_enough_col:
            return sumstats  # Return early if columns are missing
    
    # Execute the function
    result = func(*args, **kwargs)
    
    # Check shape after execution if it changed
    if sumstats is not None and show_shape:
        final_shape = (len(sumstats), len(sumstats.columns))
        if initial_shape != final_shape:
            check_dataframe_shape(
                sumstats=sumstats, 
                log=log, 
                verbose=verbose, 
                sumstats_obj=sumstats_obj
            )
    
    return result


###############################################################################################################
# Column and Argument Checking Functions
###############################################################################################################

def check_col(df_col_names: Any, verbose: bool = True, log: Log = Log(), 
             cols: Optional[List[Union[str, Tuple[str, ...], List[str]]]] = None, 
             function: Optional[str] = None) -> bool:
    """
    Verify presence of required columns prior to executing a processing step.

    Parameters
    ----------
    df_col_names : Iterable[str]
        Column names of the input DataFrame.
    verbose : bool, default True
        Whether to print log messages to stdout.
    log : gwaslab.g_Log.Log
        Logger used for messages and warnings.
    cols : list
        Required columns (strings) or tuples/lists to check as pairs/groups.
    function : str or None
        Function label used for contextual logging.

    Returns
    -------
    bool
        True when all required columns are present, False otherwise.
    """
    if not cols:
        return True
    
    not_in_df = []
    for col_spec in cols:
        if isinstance(col_spec, str):
            # Single column check
            if col_spec not in df_col_names:
                not_in_df.append(col_spec)
        else:
            # Paired/group check - all columns in group must be present
            for col in col_spec:
                if col not in df_col_names:
                    not_in_df.append(col)

    if not_in_df:
        function_label = f" for {function} " if function else " "
        log.warning(
            f"Necessary columns{function_label}were not detected: {','.join(not_in_df)}"
        )
        return False
    
    return True


def check_arg(log: Log, verbose: bool, key: str, value: Any, start_function: str) -> bool:
    """Check if a required argument is valid (non-None and non-empty)."""
    if value is None:
        log.warning(
            f"{start_function} requires non-None argument: {key}", 
            verbose=verbose
        )
        return False
    
    if isinstance(value, (str, bytes)) and str(value).strip() == "":
        log.warning(
            f"{start_function} requires non-empty argument: {key}", 
            verbose=verbose
        )
        return False
    
    return True


###############################################################################################################
# Tool Checking Functions
###############################################################################################################

def check_tools_availability(check_tools: List[str], bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Check availability and version of required tools for downstream analysis.
    
    Parameters
    ----------
    check_tools : list
        List of tool names to check (e.g., ["python", "r", "plink", "tabix", "scdrs", "r:TwoSampleMR"]).
        Tool names can match parameter names (e.g., "r", "python", "plink", "plink2", "tabix", "scdrs").
        For R packages, use format "r:PackageName" (e.g., "r:TwoSampleMR", "r:susieR").
    bound_kwargs : BoundArguments
        Bound arguments from function signature inspection.
    log : Log
        Logger instance.
    verbose : bool
        Whether to print log messages.
    
    Returns
    -------
    dict
        Dictionary with tool names as keys and dicts with "available", "path", and "version" as values.
    """
    results = {}
    generic_tools = []
    
    for tool in check_tools:
        # Check if this is an R package specification (format: "r:PackageName")
        if tool.lower().startswith("r:") and ":" in tool:
            package_name = tool.split(":", 1)[1]
            result = _check_r_package(package_name, bound_kwargs, log, verbose)
            results[tool] = result
        else:
            tool_lower = tool.lower()
            
            if tool_lower in SUPPORTED_SPECIAL_TOOLS:
                # Handle special tools with custom logic
                result = _check_special_tool(tool_lower, bound_kwargs, log, verbose)
                results[tool] = result
            else:
                # Generic tools - check in batch later
                generic_tools.append(tool)
    
    # Check generic tools using existing function
    if generic_tools:
        generic_results = _check_generic_tools(generic_tools, log, verbose)
        results.update(generic_results)
    
    return results


def _check_special_tool(tool_name: str, bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Optional[str]]:
    """Check a special tool (python, r, plink, etc.) with custom logic."""
    tool_checkers = {
        "python": _check_python_tool,
        "r": _check_r_tool,
        "plink": lambda kwargs, log, verbose: _check_plink_tool("plink", kwargs, log, verbose),
        "plink2": lambda kwargs, log, verbose: _check_plink_tool("plink2", kwargs, log, verbose),
        "tabix": lambda kwargs, log, verbose: _check_generic_path_tool("tabix", kwargs, log, verbose),
        "bcftools": lambda kwargs, log, verbose: _check_generic_path_tool("bcftools", kwargs, log, verbose),
        "scdrs": lambda kwargs, log, verbose: _check_generic_path_tool("scdrs", kwargs, log, verbose),
    }
    
    checker = tool_checkers.get(tool_name)
    if checker:
        return checker(bound_kwargs, log, verbose)
    
    # Fallback to generic check
    return _check_generic_path_tool(tool_name, bound_kwargs, log, verbose)


def _check_python_tool(bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Optional[str]]:
    """Check Python availability and version."""
    result = {"available": False, "path": None, "version": None}
    
    # Try to get Python path from arguments, sys.executable, or PATH
    python_path = bound_kwargs.arguments.get("python", None)
    if python_path is None:
        python_path = sys.executable
    if python_path is None:
        python_path = shutil.which("python3") or shutil.which("python")
    
    if not python_path:
        log.warning(" -Python not found in PATH and 'python' argument not provided", verbose=verbose)
        return result
    
    result["path"] = python_path
    result["available"] = True
    
    # Get version
    version = _get_tool_version(python_path, log, verbose, "Python")
    result["version"] = version
    
    return result


def _check_r_tool(bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Optional[str]]:
    """Check R availability and version."""
    result = {"available": False, "path": None, "version": None}
    
    # Try to get R path from arguments or PATH
    r_path = bound_kwargs.arguments.get("r", None)
    if r_path is None:
        r_path = shutil.which("R") or shutil.which("Rscript")
    
    if not r_path:
        log.warning(" -R not found in PATH and 'r' argument not provided", verbose=verbose)
        return result
    
    result["path"] = r_path
    result["available"] = True
    
    # Use existing R version checking function
    try:
        _checking_r_version(r_path, log=log, verbose=verbose)
        # Also extract version for result dict
        version = _get_tool_version(r_path, log, verbose, "R", silent=True)
        result["version"] = version
    except Exception as e:
        log.warning(f" -R found at {r_path} but version check failed: {e}", verbose=verbose)
    
    return result


def _check_r_package(package_name: str, bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Optional[str]]:
    """
    Check if an R package is installed and get its version.
    
    Parameters
    ----------
    package_name : str
        Name of the R package to check (e.g., "TwoSampleMR", "susieR").
    bound_kwargs : BoundArguments
        Bound arguments from function signature inspection.
    log : Log
        Logger instance.
    verbose : bool
        Whether to print log messages.
    
    Returns
    -------
    dict
        Dictionary with "available", "path" (R path), and "version" keys.
    """
    result = {"available": False, "path": None, "version": None}
    
    # First, get R path (required for checking packages)
    r_path = bound_kwargs.arguments.get("r", None)
    if r_path is None:
        r_path = shutil.which("R") or shutil.which("Rscript")
    
    if not r_path:
        log.warning(
            f" -R package '{package_name}' cannot be checked: R not found in PATH and 'r' argument not provided",
            verbose=verbose
        )
        return result
    
    result["path"] = r_path
    
    # Create temporary R script to check package version
    temp_r = None
    try:
        # Create a unique temporary file
        temp_r = f"_gwaslab_r_pkg_check_{package_name}_{np.random.randint(1, 99999999)}.R"
        
        # Write R script to check package version
        rscript = f'''if (!requireNamespace("{package_name}", quietly = TRUE)) {{
    cat("PACKAGE_NOT_INSTALLED")
    quit(status=1)
}} else {{
    cat(as.character(packageVersion("{package_name}")))
}}'''
        
        with open(temp_r, "w") as f:
            f.write(rscript)
        
        # Run R script
        try:
            output = subprocess.check_output(
                f"{r_path} {temp_r}",
                stderr=subprocess.STDOUT,
                shell=True,
                text=True,
                timeout=TOOL_CHECK_TIMEOUT
            )
            version = output.strip()
            
            if version and version != "PACKAGE_NOT_INSTALLED":
                result["available"] = True
                result["version"] = version
                log.write(f" -R package '{package_name}' version: {version}", verbose=verbose)
            else:
                log.warning(f" -R package '{package_name}' is not installed", verbose=verbose)
        except subprocess.CalledProcessError:
            log.warning(f" -R package '{package_name}' is not installed", verbose=verbose)
        except subprocess.TimeoutExpired:
            log.warning(f" -R package '{package_name}' version check timed out", verbose=verbose)
        except Exception as e:
            log.warning(f" -R package '{package_name}' check failed: {e}", verbose=verbose)
    
    finally:
        # Clean up temporary file
        if temp_r and os.path.exists(temp_r):
            try:
                os.remove(temp_r)
            except Exception:
                pass  # Ignore cleanup errors
    
    return result


def _check_plink_tool(tool_name: str, bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Optional[str]]:
    """Check PLINK or PLINK2 availability and version."""
    result = {"available": False, "path": None, "version": None}
    
    # Get PLINK path from arguments or PATH
    plink_path = bound_kwargs.arguments.get(tool_name, None)
    if plink_path is None:
        plink_path = shutil.which(tool_name)
    
    if not plink_path:
        log.warning(
            f" -{tool_name.upper()} not found in PATH and '{tool_name}' argument not provided", 
            verbose=verbose
        )
        return result
    
    result["path"] = plink_path
    result["available"] = True
    
    # Use existing PLINK version checking function
    try:
        if tool_name == "plink":
            _checking_plink_version(plink=plink_path, log=log, verbose=verbose)
        else:
            _checking_plink_version(plink2=plink_path, log=log, verbose=verbose)
        # Also extract version for result dict
        version = _get_tool_version(plink_path, log, verbose, tool_name.upper(), silent=True)
        result["version"] = version
    except Exception as e:
        log.warning(
            f" -{tool_name.upper()} found at {plink_path} but version check failed: {e}", 
            verbose=verbose
        )
    
    return result


def _check_generic_path_tool(tool_name: str, bound_kwargs: inspect.BoundArguments, log: Log, verbose: bool) -> Dict[str, Optional[str]]:
    """Check a generic tool that may have a parameter or fall back to PATH."""
    result = {"available": False, "path": None, "version": None}
    
    # Try to get tool path from arguments or PATH
    tool_path = bound_kwargs.arguments.get(tool_name, None)
    if tool_path is None:
        tool_path = shutil.which(tool_name)
    
    if not tool_path:
        log.warning(
            f" -{tool_name} not found in PATH and '{tool_name}' argument not provided", 
            verbose=verbose
        )
        return result
    
    result["path"] = tool_path
    result["available"] = True
    
    # Get version
    version = _get_tool_version(tool_path, log, verbose, tool_name)
    result["version"] = version
    
    return result


def _check_generic_tools(tool_names: List[str], log: Log, verbose: bool) -> Dict[str, Dict[str, Optional[str]]]:
    """Check generic tools using the existing batch checking function."""
    try:
        return _check_tool_availability(tools=tuple(tool_names), log=log, verbose=verbose)
    except Exception as e:
        log.warning(f" -Tool availability check failed: {e}", verbose=verbose)
        # Return unavailable results for all tools
        return {
            tool: {"available": False, "path": None, "version": None}
            for tool in tool_names
        }


def _get_tool_version(tool_path: str, log: Log, verbose: bool, tool_display_name: str, silent: bool = False) -> Optional[str]:
    """
    Get version information for a tool by running --version command.
    
    Parameters
    ----------
    tool_path : str
        Path to the tool executable.
    log : Log
        Logger instance.
    verbose : bool
        Whether to print log messages.
    tool_display_name : str
        Display name for the tool in log messages.
    silent : bool, default False
        If True, don't log messages (useful when another function already logged).
    
    Returns
    -------
    str or None
        Version string or None if version check failed.
    """
    try:
        output = subprocess.check_output(
            f"{tool_path} --version",
            stderr=subprocess.STDOUT,
            shell=True,
            text=True,
            timeout=TOOL_CHECK_TIMEOUT
        )
        version_line = output.strip().splitlines()[0] if output else ""
        
        if not silent:
            log.write(f" -{tool_display_name} version: {version_line}", verbose=verbose)
        
        return version_line
    except subprocess.TimeoutExpired:
        if not silent:
            log.warning(
                f" -{tool_display_name} version check timed out", 
                verbose=verbose
            )
        return None
    except Exception as e:
        if not silent:
            log.warning(
                f" -{tool_display_name} found at {tool_path} but version check failed: {e}", 
                verbose=verbose
            )
        return None
