import os
from gwaslab.g_Log import Log


def _path( *args,
                 out = None,
                 directory = None,
                 tmp = False,
                 prefix=None,
                 study = None,
                 nstudy=None,
                 trait = None,
                 exposure = None,
                 outcome = None,
                 chrom = None,
                 rsid = None,
                 snpid = None,
                 locus = None,
                 loci = None, 
                 analysis = None,
                 mode = None,
                 method = None,
                 ancestry = None,
                 population = None,
                 sample_size = None,
                 genotyping = None, 
                 pid = None,
                 build = None,
                 suffix = None,
                 log = Log(),
                 verbose=True
                 ):
    """
    Create a file path for gwaslab-generated files based on various components.

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
    analysis (str): Analysis type.
    mode (str): Mode of operation.
    method (str): Method used.
    ancestry (str): Ancestry information.
    population (str): Population identifier.
    sample_size (str): Sample size.
    genotyping (str): Genotyping platform.
    pid (str): Process ID.
    build (str): Genome build version.
    suffix (str): File extension to be appended.
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

        # ordered components excluding directory, suffix
        path_order = [
            "prefix", "study", "nstudy", "trait", "exposure", "outcome", 
            "chrom", "rsid", "snpid", "locus", "loci", "analysis", "mode_methods", 
            "ancestry", "population","sample_size", "genotyping", "pid", "build"
        ]
        
        # create default path
        ###############################################################################################################################################
        if tmp == True:
            pid = id(path_list)
        
        ###############################################################################################################################################
        all_args = locals()
        for key,value in all_args.items():
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

    # add directory
    if directory is not None:
        path = os.path.join(directory, path)
    else:
        path = os.path.join("./", path)
    ###############################################################################################################################################

    # add file extension
    if suffix is not None:
        path = ".".join([path, suffix])
    ###############################################################################################################################################

    log.write( "Creating path: {}".format(path), verbose=verbose)

    return path

def _process_out(out):
    out_dirname = os.path.dirname(out)
    out_basename = os.path.base_name(out)
