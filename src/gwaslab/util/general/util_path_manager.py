from typing import Optional, Any
import os
import re
import time
import uuid
from gwaslab.info.g_Log import Log


def _sanitize_path_component(component: str) -> str:
    """
    Sanitize a path component by removing or replacing invalid filesystem characters.
    
    Parameters:
    component (str): Path component to sanitize
    
    Returns:
    str: Sanitized path component safe for filesystem use
    """
    if not isinstance(component, str):
        component = str(component)
    
    # Replace invalid filesystem characters (Windows and Unix)
    # Invalid chars: < > : " | ? * \ / and control characters
    invalid_chars = r'[<>:"|?*\\/\x00-\x1f]'
    component = re.sub(invalid_chars, '_', component)
    
    # Replace multiple consecutive dashes/underscores with single underscore
    component = re.sub(r'[-_]+', '_', component)
    
    # Remove leading/trailing dots, dashes, and spaces (problematic on some filesystems)
    component = component.strip('. -_')
    
    # Replace spaces with underscores
    component = component.replace(' ', '_')
    
    # Limit component length to prevent issues
    max_component_length = 255
    if len(component) > max_component_length:
        component = component[:max_component_length]
    
    # Ensure component is not empty
    if not component:
        component = "unnamed"
    
    return component


def _generate_unique_id(use_uuid: bool = True) -> str:
    """
    Generate a unique identifier for file paths.
    
    Parameters:
    use_uuid (bool): If True, use UUID. If False, use timestamp-based ID.
    
    Returns:
    str: Unique identifier string
    """
    if use_uuid:
        # Use short UUID (first 8 characters) for brevity
        return uuid.uuid4().hex[:8]
    else:
        # Use timestamp with process ID for uniqueness
        timestamp = int(time.time() * 1000000) % 1000000000  # microseconds, last 9 digits
        pid = os.getpid() % 10000  # last 4 digits of PID
        return f"{timestamp:09d}_{pid:04d}"


def _handle_path_collision(path: str, max_attempts: int = 1000) -> str:
    """
    Handle path collisions by appending a unique suffix if the path already exists.
    
    Parameters:
    path (str): Original file path
    max_attempts (int): Maximum number of attempts to find a unique path
    
    Returns:
    str: Unique file path that doesn't exist
    """
    if not os.path.exists(path):
        return path
    
    # Split path into directory, basename, and extension
    directory = os.path.dirname(path)
    basename = os.path.basename(path)
    
    # Try to split extension
    if '.' in basename:
        name_parts = basename.rsplit('.', 1)
        base_name = name_parts[0]
        extension = '.' + name_parts[1]
    else:
        base_name = basename
        extension = ''
    
    # Try appending unique IDs until we find a non-existent path
    for attempt in range(1, max_attempts + 1):
        unique_id = _generate_unique_id(use_uuid=True)
        new_basename = f"{base_name}_{unique_id}{extension}"
        new_path = os.path.join(directory, new_basename)
        
        if not os.path.exists(new_path):
            return new_path
    
    # Fallback: use timestamp if all attempts fail
    timestamp = int(time.time() * 1000) % 1000000000
    new_basename = f"{base_name}_{timestamp}{extension}"
    return os.path.join(directory, new_basename)


def _validate_path_length(path: str, max_length: int = 4096) -> str:
    """
    Validate and truncate path if it exceeds filesystem limits.
    
    Parameters:
    path (str): File path to validate
    max_length (int): Maximum path length (default 4096 for most filesystems)
    
    Returns:
    str: Validated path (truncated if necessary)
    """
    if len(path) <= max_length:
        return path
    
    # Split path components
    directory = os.path.dirname(path)
    basename = os.path.basename(path)
    
    # Calculate available length for basename
    # Reserve space for directory separator and safety margin
    available_length = max_length - len(directory) - 10  # 10 char safety margin
    
    if available_length < 10:  # Minimum basename length
        # If directory is too long, truncate it
        directory = directory[:max_length - 200]  # Reserve 200 for basename
        available_length = max_length - len(directory) - 10
    
    # Truncate basename if needed
    if len(basename) > available_length:
        if '.' in basename:
            name_parts = basename.rsplit('.', 1)
            base_name = name_parts[0]
            extension = '.' + name_parts[1]
            # Reserve space for extension
            base_name = base_name[:available_length - len(extension) - 1]
            basename = base_name + extension
        else:
            basename = basename[:available_length]
    
    return os.path.join(directory, basename)


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
    tmp (bool): If True, adds '_gwaslab' prefix, uses temporary directory logic, and adds unique ID.
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
        # Generate unique ID for temporary files if needed and pid not provided
        if tmp == True and pid is None:
            pid = _generate_unique_id(use_uuid=True)
            log.write("Generated unique ID for temporary file: {}".format(pid), verbose=verbose)
        
        ###############################################################################################################################################
        all_kwargs = locals()
        for key,value in all_kwargs.items():
            if key in path_order:
                if value is not None:
                    if value !=False:
                        path_list.append(value)

        # merge path components
        path_list = list(map(str, path_list))  + list(map(str, args))
        
        # Sanitize all path components to ensure filesystem compatibility
        path_list = [_sanitize_path_component(component) for component in path_list]

        if tmp == True:
            path_list.insert(0, "_gwaslab")

        log.write( "Path component detected: {}".format(path_list), verbose=verbose)

        path = "_".join(path_list)
    else:
        # use user-provided path, but sanitize it for safety
        # Preserve extension if present
        if '.' in out_basename:
            name_parts = out_basename.rsplit('.', 1)
            base_name = _sanitize_path_component(name_parts[0])
            extension = '.' + _sanitize_path_component(name_parts[1])
            path = base_name + extension
        else:
            path = _sanitize_path_component(out_basename)
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
        # Sanitize suffix (remove any invalid characters)
        suffix = _sanitize_path_component(suffix)
        path = ".".join([path, suffix])
    ###############################################################################################################################################

    # Validate path length
    path = _validate_path_length(path)
    
    # Handle path collisions (only for files, not directories)
    # Only check if the path has a suffix (likely a file) and tmp is False
    # (tmp files should always be unique due to unique ID)
    if suffix is not None and not tmp:
        original_path = path
        path = _handle_path_collision(path)
        if path != original_path:
            log.write("Path collision detected, using unique path: {}".format(path), verbose=verbose)

    log.write( "Creating path: {}".format(path), verbose=verbose)

    return path

def _process_out(out: str) -> tuple:
    out_dirname = os.path.dirname(out)
    out_basename = os.path.basename(out)
    return out_dirname, out_basename
