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

def with_logging(start_to_msg, 
                 finished_msg,
                 start_function=None,
                 start_cols = None,
                 must_args=None,
                 show_shape=True,
                 check_dtype=False
                 ):
    if start_cols is None:
        start_cols = []
    if must_args is None:
        must_args = []
    if start_function is None:
        start_function = "this step"

    def decorator(func):
        @wraps(func)  # This preserves the original function's metadata including __doc__
        def wrapper(*args, **kwargs):
            import inspect
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            log = bound_args.arguments.get('log', Log())
            verbose = bound_args.arguments.get('verbose', True)
            
            #############################################################################################

            insumstats = bound_args.arguments.get('insumstats', None)
            if insumstats is None:
                sumstats = bound_args.arguments.get('sumstats', None)
            else:
                sumstats = insumstats
            #############################################################################################
            n_cores = bound_args.arguments.get('n_cores', None)
            ref_vcf = bound_args.arguments.get('ref_vcf', None)
            ref_fasta = bound_args.arguments.get('ref_fasta', None)
            ref_tsv = bound_args.arguments.get('ref_tsv',None)
            if n_cores is not None:
                log.write(" -Number of threads/cores to use: {}".format(n_cores))
            if ref_vcf is not None:
                log.write(" -Reference VCF: {}".format(ref_vcf))
            if ref_fasta is not None:
                log.write(" -Reference FASTA: {}".format(ref_fasta))
            if ref_tsv is not None:
                log.write(" -Reference TSV: {}".format(ref_tsv))

            # Log start message
            log.write(f"Start to {start_to_msg} ...({_get_version()})", verbose=verbose)
            
            # must_arg can not be None
            #############################################################################################
            is_args_valid = True
            for key in must_args:
                value = bound_args.arguments.get(key,None)
                is_args_valid = is_args_valid & check_arg(log, verbose, key, value, start_function)
            if is_args_valid==False:
                raise ValueError("{} must be provided.".format(must_args))

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
                    check_datatype_for_cols(sumstats=sumstats, cols=start_cols)

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
def start_to(sumstats, 
             log, 
             verbose, 
             start_line,
             end_line,
             start_cols,
             start_function,
             ref_vcf=None,
             ref_fasta=None,
             n_cores=None,
             ref_tsv=None,
             **kwargs
             ):
    
    log.write("Start to {}...{}".format(start_line,_get_version()), verbose=verbose)
    
    if sumstats is not None:
        check_dataframe_shape(sumstats=sumstats, 
                            log=log, 
                            verbose=verbose)  
        is_enough_col = check_col(sumstats.columns, 
                                verbose=verbose, 
                                log=log, 
                                cols=start_cols, 
                                function=start_function)
        if is_enough_col==True:
            if n_cores is not None:
                log.write(" -Number of threads/cores to use: {}".format(n_cores))
            if ref_vcf is not None:
                log.write(" -Reference VCF: {}".format(ref_vcf))
            if ref_fasta is not None:
                log.write(" -Reference FASTA: {}".format(ref_fasta))
            if ref_tsv is not None:
                log.write(" -Reference TSV: {}".format(ref_tsv))
            
            is_args_valid = True
            for key, value in kwargs.items():
                is_args_valid = is_args_valid & check_arg(log, verbose, key, value, start_function)
            is_enough_col = is_args_valid & is_enough_col

        if  is_enough_col == False:
            skipped(log, verbose, end_line)
        return is_enough_col
    else:
        return None

def finished(log, verbose, end_line):
    log.write("Finished {}.".format(end_line), verbose=verbose)
    gc.collect()

def skipped(log, verbose, end_line):
    log.write("Skipped {}.".format(end_line), verbose=verbose)
    gc.collect()

def check_arg(log, verbose, key, value, function):
    if value is None:
        log.warning("Necessary argument {} for {} is not provided!".format(key, function))
        return False
    return True

def check_col(df_col_names, verbose=True, log=Log(), cols=None, function=None):
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
        skipped(log, verbose, end_line=function)
        return False
    
    return True
