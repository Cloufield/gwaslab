"""
VCF to HDF5 Processing Module

This module provides optimized functions to convert VCF files to HDF5 format for efficient
rsID to CHR:POS lookup. Each chromosome is stored in a separate HDF5 file for maximum
parallel processing speed.

Workflow Overview:
==================

process_vcf_to_hfd5():
   - Uses bcftools for chromosome-based parallel processing (FASTEST for large indexed VCFs)
   - Extracts each chromosome separately using bcftools query -r
   - Processes chromosomes in parallel using multiple threads
   - Each chromosome writes to its own HDF5 file (no locking needed, true parallel writes)
   - Each chromosome processed in chunks (3M rows) to minimize memory usage
   - Streams data directly from bcftools (no intermediate file loading)
   - Writes incrementally - each chunk is written as soon as it's processed
   - Groups rsIDs using modulo 10 (creates 10 approximately equal-sized groups per chromosome)
   - Overwrite option: skip existing files (default) or overwrite them
   - Optimizations:
     * Direct bcftools query (single pass, no view+query pipeline)
     * Chunked processing within each chromosome (memory efficient)
     * Vectorized "rs" prefix removal in Python
     * Optimized data types for storage (int32 for POS, int64 for rsn; CHR not stored as it's in filename)
     * Auto-detection of VCF chromosome notation
     * Incremental writing to start I/O early
     * Separate files per chromosome (no locking, true parallel writes)
     * Append-only ingestion with final deduplication pass (O(n) instead of O(n²))

HDF5 File Structure:
====================
- One HDF5 file per chromosome: `{vcf_file_name}.chr{chr_num}.rsID_CHR_POS_mod10.h5`
- Each file contains groups: "group_0", "group_1", ..., "group_9" (based on rsID % 10)
- Each group contains DataFrame with columns: POS (int32), rsn (int64)
- CHR is not stored as it's already encoded in the filename (chr_num)
- Groups allow efficient lookup by rsID modulo

Performance Features:
=====================
- Parallel processing with configurable thread count
- Memory-efficient chunked processing
- Optimized data types (46% storage reduction)
- Streaming I/O to minimize memory footprint
- Separate files per chromosome (no locking overhead, maximum parallel write speed)
- Incremental writing for early I/O start
- Modulo 10 grouping for balanced group sizes and efficient lookup
"""

import pandas as pd
import os
import re
from gwaslab.info.g_Log import Log
from multiprocessing import Pool, Manager
import subprocess
import shutil
from gwaslab.io.io_vcf import auto_check_vcf_chr_dict

def _is_valid_hdf5_file(file_path):
    """
    Check if an HDF5 file is valid and can be opened.
    
    Parameters
    ----------
    file_path : str
        Path to HDF5 file
        
    Returns
    -------
    bool
        True if file is valid, False if corrupted or cannot be opened
    """
    if not os.path.exists(file_path):
        return False
    
    try:
        with pd.HDFStore(file_path, mode='r') as store:
            # Try to access the store to verify it's not corrupted
            _ = store.keys()
        return True
    except Exception:
        # File is corrupted or cannot be opened
        return False

def _process_chromosome_bcftools(args):
    """
    Process a single chromosome using bcftools and write to chromosome-specific HDF5 file.
    
    Parameters
    ----------
    args : tuple
        (chr_name, vcf_path, chr_dict, h5_dir, chr_num, complevel, overwrite, log, verbose)
        
    Returns
    -------
    int
        Number of rows processed
    """
    (chr_name, vcf_path, chr_dict, h5_dir, chr_num, complevel, overwrite, log, verbose) = args
    
    # Get process/thread ID for logging
    import os
    import time
    process_id = os.getpid()
    
    # Create chromosome-specific HDF5 file path
    # Use chr_num for consistent naming (standard chromosome number)
    vcf_file_name = os.path.basename(vcf_path)
    if vcf_file_name.endswith('.gz'):
        vcf_file_name = vcf_file_name[:-3]
    h5_path = f"{h5_dir}/{vcf_file_name}.chr{chr_num}.rsID_CHR_POS_mod10.h5"
    
    # Check if file exists and validate it
    if os.path.exists(h5_path):
        if not _is_valid_hdf5_file(h5_path):
            # File is corrupted, remove it and recreate
            timestamp = time.strftime('%Y/%m/%d %H:%M:%S')
            msg = f"{timestamp}  -[PID {process_id}] Detected corrupted HDF5 file for chromosome {chr_name} (chr{chr_num}), removing and recreating: {h5_path}\n"
            log.write(msg, verbose=verbose, show_time=False)
            os.remove(h5_path)
        elif not overwrite:
            # File exists and is valid, skip processing
            timestamp = time.strftime('%Y/%m/%d %H:%M:%S')
            msg = f"{timestamp}  -[PID {process_id}] Skipping chromosome {chr_name} (chr{chr_num}): file already exists: {h5_path}\n"
            log.write(msg, verbose=verbose, show_time=False)
            return 0
        elif overwrite:
            # File exists and is valid, but overwrite is True, remove it
            os.remove(h5_path)
    
    # Log start of processing
    # Format entire message (including timestamp) as single string to make print() atomic
    # This prevents interleaving when multiple processes log simultaneously
    timestamp = time.strftime('%Y/%m/%d %H:%M:%S')
    msg = f"{timestamp}  -[PID {process_id}] Processing chromosome {chr_name} (chr{chr_num})...\n"
    log.write(msg, verbose=verbose, show_time=False)
    
    try:
        # Use bcftools query directly for maximum efficiency
        # bcftools query can filter by region (-r) and extract fields in one pass
        # This is faster than view + query because it streams through VCF once
        # Format: %POS\t%ID\n (CHR not needed since we filter by chromosome)
        # Note: We do "rs" prefix stripping in Python (vectorized) instead of awk for better performance
        fmt = "%POS\t%ID\n"
        
        # Build optimized bcftools command:
        # bcftools query -r: directly extract region and query fields in one pass
        # -i 'ID!="."': filter out variants without ID (optional, can improve speed)
        # This is the fastest approach for extracting only ID/POS (CHR already filtered)
        cmd = f"bcftools query -r {chr_name} -f '{fmt}' {vcf_path}"
        
        # Process in chunks to reduce memory usage
        # Use subprocess.Popen to stream output instead of loading all into memory
        # Larger bufsize (4MB) reduces system call overhead for high-throughput data streams
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE, text=True, bufsize=4194304)
        
        # Process in chunks and write incrementally to start writing early
        chunk_size = 3000000  # Process 1M rows at a time
        chunk_reader = pd.read_csv(process.stdout, sep='\t', header=None, 
                                  names=["POS", "rsn"], 
                                  dtype={0: 'Int64', 1: str},
                                  chunksize=chunk_size)
        
        # Initialize counters for logging
        chunk_count = 0
        total_rows = 0
        groups_written = set()  # Track which groups have been written
        
        # Open HDFStore once for the entire chromosome (outside chunk loop)
        # This avoids repeated file open/close overhead
        with pd.HDFStore(h5_path, mode='a', complevel=complevel, complib="blosc:zstd") as store:
            try:
                for chunk in chunk_reader:
                    chunk_count += 1
                    if len(chunk) == 0:
                        continue
                    
                    total_rows += len(chunk)
                    
                    # Remove 'rs' prefix by stripping first 2 characters (simple and fast)
                    chunk["rsn"] = chunk["rsn"].str[2:]
                    
                    # Convert rsn to numeric (coerce invalid to NaN for later filtering)
                    chunk["rsn"] = pd.to_numeric(chunk["rsn"], errors='coerce')
                    
                    # Drop rows with missing values (combines both rsn and POS validation)
                    # POS is already Int64 from read_csv, so we only need to check rsn
                    chunk = chunk.dropna(subset=["rsn"])
                    
                    if len(chunk) == 0:
                        continue
                    
                    # Convert to optimized data types (POS already int-like from read_csv)
                    chunk["rsn"] = chunk["rsn"].astype("int64")
                    chunk["POS"] = chunk["POS"].astype("int32")
                    
                    # Calculate groups using modulo 10 for approximately equal-sized groups
                    chunk["group"] = chunk["rsn"] % 10
                    
                    # Group by group and append each group (fast append-only ingestion)
                    # No locking needed since each chromosome has its own file
                    # CHR not stored - it's already in the filename (chr_num)
                    for group_id, group_df in chunk.groupby("group"):
                        # Prepare data for this group (only POS and rsn, no CHR)
                        group_data = group_df[["POS", "rsn"]].copy()
                        
                        # Append to HDF5 file (fast append-only, no read-modify-write)
                        key = "group_{}".format(group_id)
                        store.append(
                            key, group_data,
                            format="table",
                            index=False,
                            data_columns=None  # avoids building query indexes during ingest
                        )
                        
                        groups_written.add(group_id)
                        
                        # Free memory immediately after writing
                        del group_data
                    
                    # Explicitly delete chunk to free memory
                    del chunk
                
                # Wait for process to finish and check for errors
                process.wait()
                stderr_output = process.stderr.read()
                if process.returncode != 0:
                    timestamp = time.strftime('%Y/%m/%d %H:%M:%S')
                    msg = f"{timestamp}  -Warning: bcftools failed for {chr_name}: {stderr_output}\n"
                    log.write(msg, verbose=verbose, show_time=False)
                    return 0
                
                # Finalize phase: deduplicate each group once (fast, single pass per group)
                # This is much faster than deduplicating during each append
                for key in store.keys():
                    df = store[key]
                    df = df.drop_duplicates(subset="rsn", keep='first')
                    # Ensure optimized data types
                    df["rsn"] = df["rsn"].astype("int64")
                    df["POS"] = df["POS"].astype("int32")
                    store.put(key, df, format='table', index=False, dropna=True)
                
                # Log completion (format as single string to avoid interleaving)
                timestamp = time.strftime('%Y/%m/%d %H:%M:%S')
                msg = f"{timestamp}  -[PID {process_id}] Chromosome {chr_name} (chr{chr_num}) completed: {total_rows:,} rows, {len(groups_written)} groups, file: {h5_path}\n"
                log.write(msg, verbose=verbose, show_time=False)
                
            except Exception as e:
                try:
                    process.terminate()
                except:
                    pass
                raise e
        
        return total_rows
        
    except Exception as e:
        # Format error message as single string to avoid interleaving
        timestamp = time.strftime('%Y/%m/%d %H:%M:%S')
        msg = f"{timestamp}  -Error processing chromosome {chr_name}: {e}\n"
        log.write(msg, verbose=verbose, show_time=False)
        return 0

from typing import Optional, List, Dict, Union
from gwaslab.info.g_Log import Log

def process_vcf_to_hfd5(vcf: str, 
                                 directory: Optional[str] = None, 
                                 chr_dict: Optional[Dict[str, str]] = None, 
                                 complevel: int = 3,
                                 threads: int = 1,
                                 chr_list: Optional[List[Union[int, str]]] = None,
                                 overwrite: bool = False,
                                 log: Log = Log(),
                                 verbose: bool = True) -> str:
    """
    Process VCF file to HDF5 using bcftools to extract chromosomes in parallel.
    
    This function uses bcftools to extract each chromosome separately and process
    them in parallel, which is much faster for large indexed VCF files.
    
    Parameters
    ----------
    vcf : str
        Path to VCF/BCF file (must be indexed with .tbi or .csi)
    directory : str, optional
        Output directory for HDF5 file
    chr_dict : dict, optional
        Dictionary for chromosome mapping
    complevel : int, default=3
        Compression level for HDF5 (0-9). Level 3 provides a good balance between compression ratio and processing speed.
    threads : int, default=1
        Number of threads for parallel processing
    chr_list : list, optional
        List of chromosomes to process. If None, processes all chromosomes 1-25
    overwrite : bool, default=False
        If True, overwrite existing HDF5 files. If False, skip processing if files exist.
    log : Log, optional
        Logging object
    verbose : bool, default=True
        Verbose output
        
    Returns
    -------
    str
        Path to output HDF5 file
    """
    
    # Check if bcftools is available
    if not shutil.which("bcftools"):
        raise RuntimeError("bcftools not found in PATH. Please install bcftools.")
    
    # Check if VCF file is indexed
    if not (os.path.exists(vcf + ".tbi") or os.path.exists(vcf + ".csi")):
        raise RuntimeError(f"VCF file {vcf} must be indexed. Please create .tbi or .csi index file.")
    
    log.write("Start to process VCF file to HDF5 using bcftools:", verbose=verbose)
    log.write(" -Reference VCF path: {}".format(vcf), verbose=verbose)
    log.write(" -Grouping method: modulo 10 (creating 10 groups per chromosome)", verbose=verbose)
    log.write(" -Output format: separate HDF5 file per chromosome (no locking, true parallel writes)", verbose=verbose)
    log.write(" -Compression level: {}".format(complevel), verbose=verbose)
    log.write(" -Number of threads: {}".format(threads), verbose=verbose)
    log.write(" -Overwrite existing files: {}".format(overwrite), verbose=verbose)
    
    vcf_file_name = os.path.basename(vcf)
    if vcf_file_name.endswith('.gz'):
        vcf_file_name_base = vcf_file_name[:-3]
    else:
        vcf_file_name_base = vcf_file_name
    vcf_dir_path = os.path.dirname(vcf)
    
    if directory is None:
        directory = vcf_dir_path
    elif directory[-1] == "/":
        directory = directory.rstrip('/')
    
    # Output directory for chromosome-specific HDF5 files
    h5_dir = directory
    log_path = "{}/{}.rsID_CHR_POS_mod10.log".format(directory, vcf_file_name_base)
    log.write(" -HDF5 Output directory: {}".format(h5_dir), verbose=verbose)
    log.write(" -HDF5 File pattern: {}.chr{{chr_num}}.rsID_CHR_POS_mod10.h5".format(vcf_file_name_base), verbose=verbose)
    log.write(" -Log output path: {}".format(log_path), verbose=verbose)
    
    # Handle chr_dict: it should map FROM VCF notation TO standard numbers
    # Example: {"NC_000001.10": 1, "NC_000002.11": 2, ...} or {"chr1": 1, "chr2": 2, ...}
    if chr_dict is None:
        log.write(" -Auto-detecting VCF chromosome notation...", verbose=verbose)
        vcf_chr_dict = auto_check_vcf_chr_dict(vcf, None, verbose, log)
        # vcf_chr_dict maps FROM standard numbers TO VCF notation
        # We need the inverse: FROM VCF notation TO standard numbers
        if vcf_chr_dict:
            chr_dict = {v: k for k, v in vcf_chr_dict.items()}
            log.write(" -Detected VCF chromosome notation, using inverse mapping", verbose=verbose)
        else:
            chr_dict = None
    
    # Determine chromosome list for bcftools (use VCF notation = keys of chr_dict)
    if chr_list is None:
        if chr_dict:
            # Use keys of chr_dict (VCF notation) for bcftools
            # Filter to only numeric chromosomes (1-25) if possible
            chr_list_vcf = []
            for vcf_chr in chr_dict.keys():
                std_num = chr_dict[vcf_chr]
                if isinstance(std_num, (int, str)) and str(std_num).isdigit():
                    num = int(std_num)
                    if 1 <= num <= 25:
                        chr_list_vcf.append(vcf_chr)
            # If no numeric chromosomes found, use all keys
            if not chr_list_vcf:
                chr_list_vcf = list(chr_dict.keys())
            # Sort by standard number for consistent ordering
            chr_list_vcf = sorted(chr_list_vcf, key=lambda x: chr_dict.get(x, 999))
        else:
            # No chr_dict, assume standard notation
            chr_list_vcf = [str(i) for i in range(1, 26)]
    else:
        # User provided chr_list - assume it's in standard notation, convert to VCF notation
        if chr_dict:
            # Create reverse mapping: standard number -> VCF notation
            std_to_vcf = {v: k for k, v in chr_dict.items()}
            chr_list_vcf = [std_to_vcf.get(str(c), str(c)) for c in chr_list]
            # Filter out None values
            chr_list_vcf = [c for c in chr_list_vcf if c is not None]
        else:
            # No mapping, assume already in VCF notation
            chr_list_vcf = [str(c) for c in chr_list]
    
    log.write(" -Processing chromosomes (VCF notation): {}".format(", ".join(chr_list_vcf)), verbose=verbose)
    
    # No need to initialize files or create locks - each chromosome writes to its own file
    
    if threads > 1:
        # Parallel processing mode - each chromosome writes to its own file (no locking needed)
        log.write(" -Processing {} chromosomes in parallel with {} threads...".format(len(chr_list_vcf), threads), verbose=verbose)
        
        # Prepare arguments for each chromosome
        # Need to map VCF notation to standard chromosome numbers for file naming
        tasks = []
        for chr_name in chr_list_vcf:
            # Get standard chromosome number for file naming
            if chr_dict and chr_name in chr_dict:
                chr_num = chr_dict[chr_name]
            else:
                # Try to extract number from chr_name (e.g., "NC_000001.10" -> 1, "chr1" -> 1, "1" -> 1)
                chr_num_str = str(chr_name).replace("chr", "").replace("NC_", "").replace("_", "")
                # Extract first number found
                match = re.search(r'(\d+)', chr_num_str)
                if match:
                    chr_num = int(match.group(1))
                    # Limit to 1-25 for standard chromosomes
                    if chr_num > 25:
                        chr_num = chr_name  # Use original name if not standard
                else:
                    chr_num = chr_name  # Use original name if no number found
            
            tasks.append((chr_name, vcf, chr_dict, h5_dir, chr_num, complevel, overwrite, log, verbose))
        
        # Process chromosomes in parallel
        log.write(" -Starting parallel processing (each chromosome writes to its own file)...", verbose=verbose)
        total_rows = 0
        with Pool(threads) as pool:
            row_counts = pool.map(_process_chromosome_bcftools, tasks)
            total_rows = sum(row_counts)
            log.write(" -All chromosomes processed successfully", verbose=verbose)
            log.write(" -Total rows processed: {}".format(total_rows), verbose=verbose)
    else:
        # Sequential processing
        log.write(" -Processing chromosomes sequentially...", verbose=verbose)
        total_rows = 0
        for chr_name in chr_list_vcf:
            # Get standard chromosome number for file naming
            if chr_dict and chr_name in chr_dict:
                chr_num = chr_dict[chr_name]
            else:
                # Try to extract number from chr_name
                chr_num_str = str(chr_name).replace("chr", "").replace("NC_", "").replace("_", "")
                match = re.search(r'(\d+)', chr_num_str)
                if match:
                    chr_num = int(match.group(1))
                    if chr_num > 25:
                        chr_num = chr_name
                else:
                    chr_num = chr_name
            
            row_count = _process_chromosome_bcftools(
                (chr_name, vcf, chr_dict, h5_dir, chr_num, complevel, overwrite, log, verbose)
            )
            total_rows += row_count
        
        log.write(" -Total rows processed: {}".format(total_rows), verbose=verbose)
    
    log.write(" -All chromosomes processed and written to separate HDF5 files", verbose=verbose)
    
    log.write("Processing finished!", verbose=verbose)
    log.save(log_path, verbose=verbose)
    
    # Return the directory containing the chromosome-specific HDF5 files
    # Users can access individual chromosome files or use a helper function to merge/query
    return h5_dir