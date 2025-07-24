import os
from gwaslab.g_Log import Log


def _path( *args,
                 directory = None,
                 tmp = False,
                 prefix=None,
                 study = None,
                 chrom = None,
                 rsid = None,
                 snpid = None,
                 locus = None,
                 loci = None, 
                 analysis = None,
                 mode = None,
                 method = None,
                 pid = None,
                 suffix = None,
                 build = None,
                 log = Log(),
                 verbose=True
                 ):
    
    path_list = []

    # get file name
    if tmp == True:
        path_list.append("_gwaslab")
        pid = id(path_list)
        
    if prefix is not None:
        path_list.append(prefix)
    
    if study is not None:
        path_list.append(study)
    
    if chrom is not None:
        path_list.append(chrom)
    
    if snpid is not None:
        path_list.append(snpid)

    if rsid is not None:
        path_list.append(rsid)

    if locus is not None:
        path_list.append(locus)
    
    if loci is not None:
        path_list.append(loci)
    
    if analysis is not None:
        path_list.append(analysis)
    
    if mode is not None:
        path_list.append(mode)
    
    if method is not None:
        path_list.append(method)
    
    if build is not None:
        path_list.append(build)

    if pid is not None:
        path_list.append(pid)
    
    # merge path components
    path_list = list(map(str, path_list))  + list(map(str, args))
    path = "_".join(path_list)
    
    # add directory
    if directory is not None:
        path = os.path.join(directory, path)
    
    # add file extension
    if suffix is not None:
        path = ".".join([path, suffix])
    
    log.write( " -Creating path: {}".format(path), verbose=verbose)

    return path