import gc
import pandas as pd
from functools import partial
from functools import wraps

from gwaslab.g_vchange_status import vchange_status
from gwaslab.g_vchange_status import status_match
from gwaslab.g_vchange_status import change_status
from gwaslab.g_Log import Log
from gwaslab.g_version import _get_version

from gwaslab.qc.qc_check_datatype import check_datatype_for_cols
from gwaslab.qc.qc_check_datatype import check_dataframe_shape
from gwaslab.qc.qc_build import _process_build

def with_logging(start_to_msg, 
                 finished_msg,
                 start_function=None,
                 start_cols = None,
                 must_kwargs=None,
                 show_shape=True,
                 check_dtype=False,
                 fix=False
                 ):
    """
    Decorator to add standardized logging, argument checks, and optional dtype
    verification to QC/processing functions.

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
    check_dtype : bool, default False
        If True, verify dtypes for `start_cols` using `check_datatype_for_cols`.
    fix : bool, default True
        If True and `check_dtype` is enabled, attempt dtype fixes where possible.

    Behavior
    --------
    - Logs references (threads, VCF/FASTA/TSV) and start/end messages.
    - Validates required args and columns; skips execution if columns missing.
    - Optionally checks DataFrame dtypes and applies fixes.
    - Preserves original function metadata via `functools.wraps`.
    """
    if fix ==True:
        check_dtype=True
    if start_cols is None:
        start_cols = []
    if must_kwargs is None:
        must_kwargs = []
    if start_function is None:
        start_function = "this step"

    def decorator(func):
        @wraps(func)  # This preserves the original function's metadata including __doc__
        def wrapper(*args, **kwargs):
            import inspect
            sig = inspect.signature(func)
            bound_kwargs = sig.bind(*args, **kwargs)
            bound_kwargs.apply_defaults()
            log = bound_kwargs.arguments.get('log', Log())
            verbose = bound_kwargs.arguments.get('verbose', True)
            
            #############################################################################################

            insumstats = bound_kwargs.arguments.get('insumstats', None)
            if insumstats is None:
                sumstats = bound_kwargs.arguments.get('sumstats', None)
            else:
                sumstats = insumstats
            #############################################################################################
            n_cores = bound_kwargs.arguments.get('n_cores', None)
            ref_vcf = bound_kwargs.arguments.get('ref_vcf', None)
            ref_fasta = bound_kwargs.arguments.get('ref_fasta', None)
            ref_tsv = bound_kwargs.arguments.get('ref_tsv',None)
            if n_cores is not None:
                log.write(" -Number of threads/cores to use: {}".format(n_cores), verbose=verbose)
            if ref_vcf is not None:
                log.write(" -Reference VCF: {}".format(ref_vcf), verbose=verbose)
            if ref_fasta is not None:
                log.write(" -Reference FASTA: {}".format(ref_fasta), verbose=verbose)
            if ref_tsv is not None:
                log.write(" -Reference TSV: {}".format(ref_tsv), verbose=verbose)

            if "build" in bound_kwargs.arguments:
                bound_kwargs.arguments["build"] = _process_build(bound_kwargs.arguments["build"], log=log, verbose=verbose)
                args = bound_kwargs.args
                kwargs = bound_kwargs.kwargs

            # Log start message
            log.write(f"Start to {start_to_msg} ...({_get_version()})", verbose=verbose)
            
            # must_arg can not be None
            #############################################################################################
            is_kwargs_valid = True
            for key in must_kwargs:
                value = bound_kwargs.arguments.get(key,None)
                is_kwargs_valid = is_kwargs_valid & check_arg(log, verbose, key, value, start_function)
            if is_kwargs_valid==False:
                raise ValueError("{} must be provided.".format(must_kwargs))

            #############################################################################################
            if sumstats is not None:
                # check sumstats shape, columns
                initial_shape = None
                if show_shape:
                    initial_shape = (len(sumstats), len(sumstats.columns))
                    check_dataframe_shape(sumstats=sumstats, 
                                        log=log, 
                                        verbose=verbose)  
                
                is_enough_col = check_col(sumstats.columns, 
                                        verbose=verbose, 
                                        log=log, 
                                        cols=start_cols, 
                                        function=start_function)
                if check_dtype:
                    check_datatype_for_cols(sumstats=sumstats, cols=start_cols, fix=fix)

                if is_enough_col==True:
                    # Execute the original function
                    result = func(*args, **kwargs)
                else:
                    # not enough cols, return sumstats
                    return sumstats
            else:
                # Execute the original function
                result = func(*args, **kwargs)

            if sumstats is not None and show_shape:
                final_shape = (len(sumstats), len(sumstats.columns))
                if initial_shape != final_shape:
                    check_dataframe_shape(sumstats=sumstats, log=log, verbose=verbose)
            # Log finish message
            log.write(f"Finished {finished_msg}.", verbose=verbose)
            return result
        return wrapper
    return decorator


###############################################################################################################

def check_col(df_col_names, verbose=True, log=Log(), cols=None, function=None):
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
    not_in_df=[]
    for i in cols:
        if type(i) is str:
            # single check
            if i in df_col_names:
                continue
            else:
                not_in_df.append(i)
        else:
            # paried check
            count=0
            for j in i:
                if j not in df_col_names:
                    not_in_df.append(j)
                    count+=1

    if len(not_in_df)>0:
        if function is None:
            to_show_title=" "
        else:
            to_show_title = " for {} ".format(function)
        log.warning("Necessary columns{}were not detected:{}".format(to_show_title, ",".join(not_in_df)))
        return False
    
    return True

def check_arg(log, verbose, key, value, start_function):
    if value is None:
        log.warning("{} requires non-None argument: {}".format(start_function, key), verbose=verbose)
        return False
    if isinstance(value, (str, bytes)) and str(value).strip() == "":
        log.warning("{} requires non-empty argument: {}".format(start_function, key), verbose=verbose)
        return False
    return True
