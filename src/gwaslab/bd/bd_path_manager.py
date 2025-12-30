from typing import Optional, Any
import os
from gwaslab.info.g_Log import Log


def _path(*args: Any,
                 out: Optional[str] = None,
                 directory: Optional[str] = None,
                 tmp: bool = False,
                 prefix: Optional[str] = None,
                 study: Optional[str] = None,
                 nstudy: Optional[str] = None,
                 trait: Optional[str] = None,
                 exposure: Optional[str] = None,
                 outcome: Optional[str] = None,
                 chrom: Optional[str] = None,
                 rsid: Optional[str] = None,
                 snpid: Optional[str] = None,
                 locus: Optional[str] = None,
                 loci: Optional[str] = None, 
                 analysis: Optional[str] = None,
                 mode: Optional[str] = None,
                 method: Optional[str] = None,
                 ancestry: Optional[str] = None,
                 population: Optional[str] = None,
                 sample_size: Optional[str] = None,
                 genotyping: Optional[str] = None, 
                 pid: Optional[str] = None,
                 build: Optional[str] = None,
                 suffix: Optional[str] = None,
                 result_type: Optional[str] = None,
                 subdirectory: Optional[str] = None,
                 log: Log = Log(),
                 verbose: bool = True
                 ) -> str:
    """
    Create a file path for gwaslab-generated files based on various components.
    Supports both general file paths and downstream analysis result files.

    Parameters:
    out (str): Full path to output file or directory. If provided, it takes precedence over other path components.
    directory (str): Directory path where the file should be saved. Overrides any directory in 'out'.
    tmp (bool): If True, adds '_gwaslab' prefix and uses temporary directory logic.
    prefix (str): Prefix to be added at the beginning of the filename.
    study (str): Study identifier.
    nstudy (str): Numbered study identifier.
    trait (str): Trait name.
    exposure (str): Exposure variable name.
    outcome (str): Outcome variable name.
    chrom (str): Chromosome number.
    rsid (str): SNP rsID.
    snpid (str): SNP identifier.
    locus (str): Locus name.
    loci (str): Multiple loci names.
    analysis (str): Analysis type (e.g., 'mtag', 'clumping', 'prs').
    mode (str): Mode of operation.
    method (str): Method used.
    ancestry (str): Ancestry information.
    population (str): Population identifier.
    sample_size (str): Sample size.
    genotyping (str): Genotyping platform.
    pid (str): Process ID.
    build (str): Genome build version.
    suffix (str): File extension to be appended (e.g., 'tsv', 'csv', 'png').
    result_type (str): Type of result file (e.g., 'summary', 'plot', 'table', 'log', 'report', 'intermediate').
                      Useful for organizing downstream analysis results.
    subdirectory (str): Subdirectory name to create nested directory structures for results.
                       If provided, creates a subdirectory within the main directory.
    log (Log): Logging object for progress messages.
    verbose (bool): If True, logs progress messages.

    Returns:
    str: Constructed file path based on provided components, with appropriate directory and suffix.
    """
    out_basename = ""

    if out is not None:
        if os.path.isdir(out):
            directory = out
            log.write( "Directory detected: {}".format(directory), verbose=verbose)
        else:
            directory = os.path.dirname(out)
            log.write( "Directory detected: {}".format(directory), verbose=verbose)
            out_basename = os.path.basename(out)
            log.write( "Basename detected: {}".format(out_basename), verbose=verbose)

    if out_basename == "":
        # create default path
        ###############################################################################################################################################
        path_list = []

        # ordered components excluding directory, suffix, result_type, subdirectory
        path_order = [
            "prefix", "study", "nstudy", "trait", "exposure", "outcome", 
            "chrom", "rsid", "snpid", "locus", "loci", "analysis", "mode", "method", 
            "ancestry", "population","sample_size", "genotyping", "pid", "build", "result_type"
        ]
        
        # create default path
        ###############################################################################################################################################
        if tmp == True:
            pid = id(path_list)
        
        ###############################################################################################################################################
        all_kwargs = locals()
        for key,value in all_kwargs.items():
            if key in path_order:
                if value is not None:
                    if value !=False:
                        path_list.append(value)

        # merge path components
        path_list = list(map(str, path_list))  + list(map(str, args))
        path_list = list(map(lambda x: x.replace("_","-") , path_list))

        if tmp == True:
            path_list.insert(0, "_gwaslab")

        log.write( "Path component detected: {}".format(path_list), verbose=verbose)

        path = "_".join(path_list)
    else:
        # use user-provided path
        path = out_basename
    ###############################################################################################################################################

    # add directory and subdirectory for result files
    if directory is not None:
        # Create subdirectory structure for downstream analysis results if specified
        if subdirectory is not None:
            directory = os.path.join(directory, subdirectory)
            # Create subdirectory if it doesn't exist (for result file organization)
            if not os.path.exists(directory):
                try:
                    os.makedirs(directory, exist_ok=True)
                    log.write("Created subdirectory for results: {}".format(directory), verbose=verbose)
                except OSError:
                    log.write("Warning: Could not create subdirectory: {}".format(directory), verbose=verbose)
        path = os.path.join(directory, path)
    else:
        # Handle subdirectory even when main directory is not specified
        if subdirectory is not None:
            subdir_path = os.path.join("./", subdirectory)
            if not os.path.exists(subdir_path):
                try:
                    os.makedirs(subdir_path, exist_ok=True)
                    log.write("Created subdirectory for results: {}".format(subdir_path), verbose=verbose)
                except OSError:
                    log.write("Warning: Could not create subdirectory: {}".format(subdir_path), verbose=verbose)
            path = os.path.join(subdir_path, path)
        else:
            path = os.path.join("./", path)
    ###############################################################################################################################################

    # add file extension
    if suffix is not None:
        path = ".".join([path, suffix])
    ###############################################################################################################################################

    log.write( "Creating path: {}".format(path), verbose=verbose)

    return path

def _process_out(out: str) -> tuple:
    out_dirname = os.path.dirname(out)
    out_basename = os.path.basename(out)
    return out_dirname, out_basename
