"""
GWASLab Standard Format (GSF) - Efficient storage format for GWAS sumstats.
GSF uses Parquet format internally for optimal compression and performance.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import pyarrow.compute as pc
import json
import os
from pathlib import Path
from typing import Optional, Dict, List, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_reserved_headers import _reserved_headers_dict, _get_headers
from gwaslab.info.g_meta import _init_meta


def _get_optimal_column_order(sumstats: pd.DataFrame) -> List[str]:
    """Get optimal column order: info -> stats -> others."""
    info_cols = _get_headers(mode="info")
    stats_cols = _get_headers(mode="stats")
    others_cols = [h for h in _get_headers() if h not in info_cols and h not in stats_cols]
    
    existing_cols = set(sumstats.columns)
    ordered_cols = []
    
    for col_list in [info_cols, stats_cols, others_cols]:
        for col in col_list:
            if col in existing_cols:
                ordered_cols.append(col)
    
    remaining_cols = [col for col in sumstats.columns if col not in ordered_cols]
    ordered_cols.extend(remaining_cols)
    
    return ordered_cols


def _optimize_dtypes_for_parquet(sumstats: pd.DataFrame) -> pd.DataFrame:
    """Optimize data types for GSF (Parquet-based): category -> string, preserve int/float types."""
    sumstats = sumstats.copy()
    
    for col in sumstats.columns:
        if col in _reserved_headers_dict:
            target_dtype = _reserved_headers_dict[col].get("Datatype", "")
            
            if sumstats[col].dtype.name == "category":
                sumstats[col] = sumstats[col].astype("string")
            elif target_dtype == "string" and sumstats[col].dtype.name not in ["string", "object"]:
                sumstats[col] = sumstats[col].astype("string")
            elif target_dtype == "Int64" and sumstats[col].dtype.name not in ["Int64", "int64"]:
                if sumstats[col].dtype.name.startswith("int"):
                    sumstats[col] = sumstats[col].astype("Int64")
            elif target_dtype == "float64" and sumstats[col].dtype.name not in ["float64", "Float64"]:
                if sumstats[col].dtype.name.startswith("float"):
                    sumstats[col] = sumstats[col].astype("float64")
            elif target_dtype == "float32" and sumstats[col].dtype.name not in ["float32", "Float32"]:
                if sumstats[col].dtype.name.startswith("float"):
                    sumstats[col] = sumstats[col].astype("float32")
    
    return sumstats


def _prepare_metadata(meta: Optional[Dict] = None) -> Dict[str, str]:
    """Prepare metadata for GSF file metadata."""
    parquet_meta = {
        "gwaslab_format": "gsf",
        "gwaslab_version": "1.0",
    }
    
    if meta is not None:
        if "gwaslab" in meta:
            gwaslab_meta = meta["gwaslab"]
            parquet_meta["gwaslab_study_name"] = str(gwaslab_meta.get("study_name", "Unknown"))
            parquet_meta["gwaslab_genome_build"] = str(gwaslab_meta.get("genome_build", "Unknown"))
            parquet_meta["gwaslab_species"] = str(gwaslab_meta.get("species", "Unknown"))
        
        try:
            parquet_meta["gwaslab_full_meta"] = json.dumps(meta)
        except (TypeError, ValueError):
            parquet_meta["gwaslab_full_meta"] = json.dumps({"error": "metadata_not_serializable"})
    
    return parquet_meta


def write_gsf(
    sumstats_or_dataframe,
    path: Union[str, Path],
    meta: Optional[Dict] = None,
    partition_cols: Optional[List[str]] = None,
    compression: str = "zstd",
    log: Optional[Log] = None,
    verbose: bool = True
) -> str:
    """
    Save sumstats to GSF (GWASLab Standard Format) file.
    
    Parameters
    ----------
    sumstats_or_dataframe : Sumstats or pd.DataFrame
        Sumstats object or DataFrame to process.
    path : str or Path
        Output path (.gsf file or directory for partitioned)
    meta : dict, optional
        GWASLab metadata dictionary
    partition_cols : list of str, optional
        Columns to partition by (e.g., ["CHR"])
    compression : str, default "zstd"
        Compression codec: "snappy", "gzip", "brotli", "zstd", "lz4"
    log : Log, optional
        Log object
    verbose : bool, default True
        Print progress messages
        
    Examples
    --------
    >>> write_gsf(sumstats, "sumstats.gsf")
    >>> write_gsf(sumstats, "sumstats_partitioned", partition_cols=["CHR"])
        
    Returns
    -------
    str
        Path to written file
    """
    import pandas as pd
    # Handle both DataFrame and Sumstats object
    if isinstance(sumstats_or_dataframe, pd.DataFrame):
        sumstats = sumstats_or_dataframe
    else:
        sumstats = sumstats_or_dataframe.data
    
    if log is None:
        log = Log()
    
    if verbose:
        log.write("Writing GSF format (Parquet-based)...", verbose=verbose)
        log.write(f"  -Input variants: {len(sumstats):,}", verbose=verbose)
        log.write(f"  -Input columns: {len(sumstats.columns)}", verbose=verbose)
    
    # Check for empty DataFrame
    if sumstats.empty:
        log.write("Warning: DataFrame is empty", verbose=verbose)
    
    if verbose:
        log.write("  -Optimizing data types for Parquet...", verbose=verbose)
    sumstats = _optimize_dtypes_for_parquet(sumstats)
    ordered_cols = _get_optimal_column_order(sumstats)
    sumstats = sumstats[ordered_cols]
    
    if verbose:
        log.write(f"  -Column order optimized: {len(ordered_cols)} columns", verbose=verbose)
        log.write(f"  -Compression: {compression}", verbose=verbose)
    
    parquet_meta = _prepare_metadata(meta)
    if verbose and meta:
        if "gwaslab" in meta:
            gwaslab_meta = meta.get("gwaslab", {})
            if "genome_build" in gwaslab_meta:
                log.write(f"  -Genome build: {gwaslab_meta['genome_build']}", verbose=verbose)
            if "study_name" in gwaslab_meta:
                log.write(f"  -Study: {gwaslab_meta['study_name']}", verbose=verbose)
        log.write("  -Metadata will be embedded in Parquet file", verbose=verbose)
    
    table = pa.Table.from_pandas(sumstats, preserve_index=False)
    
    if parquet_meta:
        existing_meta = table.schema.metadata or {}
        combined_meta = {**existing_meta}
        for key, value in parquet_meta.items():
            combined_meta[key.encode()] = str(value).encode()
        table = table.replace_schema_metadata(combined_meta)
    
    path = Path(path)
    if partition_cols:
        if verbose:
            log.write(f"  -Partitioning by: {', '.join(partition_cols)}", verbose=verbose)
        if path.suffix == ".gsf":
            path = path.parent / path.stem
        path.mkdir(parents=True, exist_ok=True)
        output_path = str(path)
        if verbose:
            log.write(f"  -Writing partitioned Parquet dataset to: {output_path}", verbose=verbose)
        pq.write_to_dataset(
            table,
            root_path=output_path,
            compression=compression,
            write_statistics=True,
            partition_cols=partition_cols
        )
    else:
        if path.suffix != ".gsf":
            path = path.with_suffix(".gsf")
        # Ensure parent directory exists (skip if parent is current dir or root)
        parent = path.parent
        if parent != Path('.') and str(parent) != '/':
            try:
                parent.mkdir(parents=True, exist_ok=True)
            except (OSError, PermissionError):
                # Parent might already exist or be root, which is fine
                pass
        output_path = str(path)
        if verbose:
            log.write(f"  -Writing Parquet file to: {output_path}", verbose=verbose)
        pq.write_table(
            table,
            output_path,
            compression=compression,
            write_statistics=True
        )
    
    if verbose:
        # Calculate file size (handle both files and directories)
        if os.path.isdir(output_path):
            # For partitioned datasets, calculate total size
            parquet_files = []
            for dirpath, dirnames, filenames in os.walk(output_path):
                for filename in filenames:
                    if filename.endswith('.parquet'):
                        parquet_files.append(os.path.join(dirpath, filename))
            total_size = sum(os.path.getsize(f) for f in parquet_files) / (1024 * 1024)
            log.write(f"  -Partitioned dataset written: {len(parquet_files)} Parquet file(s)", verbose=verbose)
            log.write(f"Successfully wrote GSF partitioned dataset ({total_size:.2f} MB total)", verbose=verbose)
        else:
            file_size = os.path.getsize(output_path) / (1024 * 1024)
            log.write(f"Successfully wrote GSF file ({file_size:.2f} MB)", verbose=verbose)
    
    return output_path


def _parse_filter_string(filter_str: str, schema: pa.Schema = None):
    """
    Parse filter string to PyArrow Expression.
    
    Supports operators: ==, !=, <, <=, >, >=, in
    Supports logical operators: & (AND), | (OR)
    
    Examples:
        "P < 5e-8"
        "CHR == 1"
        "P < 5e-8 & CHR == 1"
        "CHR in [1, 2, 3]"
    """
    import re
    
    # Split by & or | while handling brackets in 'in' operator
    parts = []
    current = ""
    bracket_depth = 0
    
    i = 0
    while i < len(filter_str):
        char = filter_str[i]
        
        if char == '[':
            bracket_depth += 1
            current += char
        elif char == ']':
            bracket_depth -= 1
            current += char
        elif bracket_depth == 0:
            if char == '&' and (i == 0 or filter_str[i-1] != '&'):
                if current.strip():
                    parts.append(('&', current.strip()))
                current = ""
            elif char == '|' and (i == 0 or filter_str[i-1] != '|'):
                if current.strip():
                    parts.append(('|', current.strip()))
                current = ""
            else:
                current += char
        else:
            current += char
        
        i += 1
    
    if current.strip():
        parts.append(('&', current.strip()))
    
    
    # Track if we need type casting (string fields compared to numeric values)
    needs_type_cast = False
    cast_info = {}  # {col: target_type}
    
    # Parse each condition
    expressions = []
    for op, condition in parts:
        condition = condition.strip()
        if not condition:
            continue
        
        # Try 'in' operator first (more specific pattern)
        in_match = re.match(r'(\w+)\s+in\s+\[(.*?)\]', condition)
        if in_match:
            col = in_match.group(1)
            value_str = in_match.group(2)
            # Parse list values
            values = [v.strip().strip('"\'') for v in value_str.split(',') if v.strip()]
            # Try to convert to numbers
            try:
                values = [float(v) if '.' in v or 'e' in v.lower() else int(v) for v in values]
            except ValueError:
                pass  # Keep as strings
            
            field = ds.field(col)
            # Check if we need type casting
            if schema and col in schema.names:
                field_type = schema.field(col).type
                if all(isinstance(v, (int, float)) for v in values):
                    if pa.types.is_string(field_type) or pa.types.is_large_string(field_type):
                        needs_type_cast = True
                        if all(isinstance(v, int) or (isinstance(v, float) and v.is_integer()) for v in values):
                            cast_info[col] = pa.int64()
                        else:
                            cast_info[col] = pa.float64()
                        # For now, use string comparison (will be fixed after casting)
                        expr = field.isin([str(v) for v in values])
                    else:
                        expr = field.isin(values)
                else:
                    expr = field.isin(values)
            else:
                expr = field.isin(values)
            
            expressions.append((op, expr))
            continue
        
        # Try comparison operators
        compare_match = re.match(r'(\w+)\s*(==|!=|<=|>=|<|>)\s*(.+)', condition)
        if compare_match:
            col = compare_match.group(1)
            op_str = compare_match.group(2)
            value_str = compare_match.group(3).strip().strip('"\'')
            
            # Try to convert value to number
            is_numeric = False
            try:
                if 'e' in value_str.lower() or 'E' in value_str:
                    value = float(value_str)
                    is_numeric = True
                elif '.' in value_str:
                    value = float(value_str)
                    is_numeric = True
                else:
                    value = int(value_str)
                    is_numeric = True
            except ValueError:
                value = value_str  # Keep as string
            
            field = ds.field(col)
            
            # Check if we need type casting
            if schema and col in schema.names:
                field_type = schema.field(col).type
                if is_numeric and (pa.types.is_string(field_type) or pa.types.is_large_string(field_type)):
                    needs_type_cast = True
                    if isinstance(value, int):
                        cast_info[col] = pa.int64()
                    elif isinstance(value, float):
                        cast_info[col] = pa.float64()
                    # For now, use string comparison (will be fixed after casting)
                    if op_str == "==":
                        expr = field == str(value)
                    elif op_str == "!=":
                        expr = field != str(value)
                    else:
                        # For <, <=, >, >=, we need numeric comparison, so mark for casting
                        expr = None  # Will be handled after casting
                else:
                    # No casting needed, use normal comparison
                    if op_str == "==":
                        expr = field == value
                    elif op_str == "!=":
                        expr = field != value
                    elif op_str == "<":
                        expr = field < value
                    elif op_str == "<=":
                        expr = field <= value
                    elif op_str == ">":
                        expr = field > value
                    elif op_str == ">=":
                        expr = field >= value
                    else:
                        raise ValueError(f"Unsupported operator: {op_str}")
            else:
                # No schema info, assume types match
                if op_str == "==":
                    expr = field == value
                elif op_str == "!=":
                    expr = field != value
                elif op_str == "<":
                    expr = field < value
                elif op_str == "<=":
                    expr = field <= value
                elif op_str == ">":
                    expr = field > value
                elif op_str == ">=":
                    expr = field >= value
                else:
                    raise ValueError(f"Unsupported operator: {op_str}")
            
            if expr is not None:
                expressions.append((op, expr))
            else:
                # Store condition info for later (after casting)
                expressions.append((op, (col, op_str, value)))
            continue
        
        # If no pattern matched
        raise ValueError(f"Could not parse filter condition: {condition}")
    
    # If type casting is needed, return None for filter_expr and the cast info
    if needs_type_cast:
        return None, True, cast_info, filter_str
    
    # Combine expressions with logical operators
    if not expressions:
        raise ValueError("No valid filter conditions found")
    
    result = expressions[0][1]
    for i in range(1, len(expressions)):
        op, expr = expressions[i]
        if op == '&':
            result = result & expr
        elif op == '|':
            result = result | expr
    
    return result, False, {}, None


def load_gsf(
    path: Union[str, Path],
    columns: Optional[List[str]] = None,
    filters: Optional[str] = None,
    verbose: bool = True
) -> 'Sumstats':
    """
    Load GWAS sumstats from GSF (GWASLab Standard Format) file.
    
    Parameters
    ----------
    path : str or Path
        Path to .gsf file or directory (for partitioned)
    columns : list of str, optional
        Columns to read (None = all columns)
    filters : str, optional
        Filter string for predicate pushdown. Supports:
        - Operators: ==, !=, <, <=, >, >=, in
        - Logical operators: & (AND), | (OR)
        Examples:
          "P < 5e-8"
          "CHR == 1"
          "P < 5e-8 & CHR == 1"
          "CHR in [1, 2, 3]"
          "P < 5e-8 | P > 0.99"
    verbose : bool, default True
        Print progress messages
        
    Returns
    -------
    Sumstats
        GWASLab Sumstats object
        
    Examples
    --------
    >>> mysumstats = gl.load_gsf("sumstats.gsf")
    >>> mysumstats = gl.load_gsf("sumstats.gsf", columns=["CHR", "POS", "BETA", "P"])
    >>> mysumstats = gl.load_gsf("sumstats.gsf", filters="CHR == 1")
    >>> mysumstats = gl.load_gsf("sumstats.gsf", filters="P < 5e-8 & CHR == 1")
    >>> mysumstats = gl.load_gsf("sumstats.gsf", filters="CHR in [1, 2, 3]")
    """
    log = Log()
    path = Path(path)
    
    # Check if path exists
    if not path.exists():
        raise FileNotFoundError(f"GSF file or directory not found: {path}")
    
    if verbose:
        if path.is_dir():
            log.write(f"Loading GSF partitioned dataset: {path}", verbose=verbose)
        else:
            log.write(f"Loading GSF file (Parquet-based): {path}", verbose=verbose)
            file_size = os.path.getsize(path) / (1024 * 1024)
            log.write(f"  -File size: {file_size:.2f} MB", verbose=verbose)
    
    try:
        # Use PyArrow dataset API for better filtering and performance
        # This supports both files and directories, and provides efficient predicate pushdown
        dataset = ds.dataset(str(path), format="parquet")
        
        # Get schema for type-aware filter parsing
        schema = dataset.schema
        
        if verbose:
            log.write(f"  -Parquet schema: {len(schema)} columns", verbose=verbose)
            if path.is_dir():
                # Count parquet files in directory
                parquet_files = list(path.rglob("*.parquet"))
                if parquet_files:
                    log.write(f"  -Partitioned dataset: {len(parquet_files)} Parquet file(s)", verbose=verbose)
        
        # Parse filter string to PyArrow Expression
        filter_expr = None
        needs_type_cast = False
        cast_info = {}
        filter_str = None
        if filters is not None:
            if verbose:
                log.write(f"  -Filter expression: {filters}", verbose=verbose)
            filter_expr, needs_type_cast, cast_info, filter_str = _parse_filter_string(filters, schema=schema)
            
            if needs_type_cast and verbose:
                log.write(f"  -Type casting required for: {', '.join(cast_info.keys())}", verbose=verbose)
                for col, target_type in cast_info.items():
                    log.write(f"    * {col}: {schema.field(col).type} -> {target_type}", verbose=verbose)
        
        if columns is not None:
            if verbose:
                log.write(f"  -Column selection: {len(columns)} column(s) specified", verbose=verbose)
        
        # If type casting is needed (string fields compared to numeric values),
        # we need to read first, cast, then filter
        if needs_type_cast:
            if verbose:
                log.write("  -Reading data with type casting (predicate pushdown limited)...", verbose=verbose)
            # Read all data first (or with column selection)
            table = dataset.to_table(columns=columns)
            # Cast string columns to numeric types
            for col, target_type in cast_info.items():
                if col in table.column_names:
                    table = table.set_column(
                        table.schema.get_field_index(col),
                        col,
                        pc.cast(table[col], target_type)
                    )
            # Convert to pandas, apply filter, convert back
            df = table.to_pandas()
            # Apply filter using pandas query (now types match)
            try:
                df = df.query(filter_str, engine='python')
                table = pa.Table.from_pandas(df)
                if verbose:
                    log.write(f"  -Filter applied after type casting", verbose=verbose)
            except Exception as e:
                if verbose:
                    log.write(f"Warning: Could not apply filter after type casting: {e}. Loading all data.", verbose=verbose)
        else:
            if verbose:
                if filter_expr is not None:
                    log.write("  -Using predicate pushdown for filtering...", verbose=verbose)
                elif columns is not None:
                    log.write("  -Using column pruning...", verbose=verbose)
            # Read with filters and column selection (predicate pushdown)
            table = dataset.to_table(columns=columns, filter=filter_expr)
    except Exception as e:
        raise IOError(f"Failed to read GSF file: {e}") from e
    
    df = table.to_pandas()
    
    if verbose:
        log.write(f"  -Loaded {len(df):,} variants", verbose=verbose)
        log.write(f"  -Loaded {len(df.columns)} columns", verbose=verbose)
    
    # Extract metadata
    meta = None
    schema_meta = table.schema.metadata
    if schema_meta:
        for key, value in schema_meta.items():
            key_str = key.decode() if isinstance(key, bytes) else key
            if key_str == "gwaslab_full_meta":
                value_str = value.decode() if isinstance(value, bytes) else value
                try:
                    meta = json.loads(value_str)
                    if verbose:
                        log.write("  -Metadata found in GSF file", verbose=verbose)
                        if "gwaslab" in meta:
                            gwaslab_meta = meta.get("gwaslab", {})
                            if "genome_build" in gwaslab_meta:
                                log.write(f"    * Genome build: {gwaslab_meta['genome_build']}", verbose=verbose)
                            if "study_name" in gwaslab_meta:
                                log.write(f"    * Study: {gwaslab_meta['study_name']}", verbose=verbose)
                            if "gwaslab_version" in gwaslab_meta:
                                log.write(f"    * GWASLab version: {gwaslab_meta['gwaslab_version']}", verbose=verbose)
                except (json.JSONDecodeError, TypeError):
                    if verbose:
                        log.write("  -Warning: Could not parse metadata from GSF file", verbose=verbose)
                    pass
                break
    else:
        if verbose:
            log.write("  -No metadata found in GSF file", verbose=verbose)
    
    if verbose:
        log.write(f"Successfully loaded GSF file: {len(df):,} variants, {len(df.columns)} columns", verbose=verbose)
    
    # Reconstruct Sumstats object exactly as it would be normally
    # Bypass _preformat since data is already in GWASLab format
    # Import here to avoid circular import
    from gwaslab.g_Sumstats import Sumstats
    sumstats_obj = object.__new__(Sumstats)
    
    # Initialize log first (required for other initializations)
    sumstats_obj.log = Log()
    sumstats_obj.log._sumstats_obj = sumstats_obj
    
    # Set data directly (already in correct format)
    sumstats_obj.data = df
    sumstats_obj._last_shape = None
    
    # Set metadata (use saved metadata or initialize new)
    if meta is not None:
        sumstats_obj.meta = meta
    else:
        sumstats_obj.meta = _init_meta()
    
    # Initialize other required attributes
    from gwaslab.bd.bd_path_manager import _path
    from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config
    from gwaslab.info.g_version import _show_version
    
    # Set build from metadata if available
    if "gwaslab" in sumstats_obj.meta and "genome_build" in sumstats_obj.meta["gwaslab"]:
        sumstats_obj._build = sumstats_obj.meta["gwaslab"]["genome_build"]
    else:
        sumstats_obj._build = "99"
    
    # Initialize visualization parameters
    sumstats_obj.viz_params = VizParamsManager()
    load_viz_config(sumstats_obj.viz_params)
    
    # Set object ID and temp path
    sumstats_obj.id = id(sumstats_obj)
    sumstats_obj.tmp_path = _path(pid=sumstats_obj.id, log=sumstats_obj.log, verbose=False)
    
    # Initialize downstream analysis result manager
    from gwaslab.downstream.ds_result_manager import DownstreamResultManager
    sumstats_obj.downstream = DownstreamResultManager()
    
    # Show version info if verbose
    if verbose:
        _show_version(sumstats_obj.log, verbose=verbose)
        sumstats_obj.log.write("Loaded GSF file (bypassed _preformat - data already in GWASLab format)", verbose=verbose)
    
    return sumstats_obj
