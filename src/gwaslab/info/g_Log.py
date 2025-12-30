import time
from typing import TYPE_CHECKING, Optional, Tuple, Any

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

class Log():
    def __init__(self) -> None:
        self.log_text=str(time.strftime('%Y/%m/%d %H:%M:%S'))+ " " + "Sumstats Object created."+ "\n"
    
    def write(self, *message: Any, end: str = "\n", show_time: bool = True, verbose: bool = True) -> None:
        if show_time is True:
            if verbose: print(str(time.strftime('%Y/%m/%d %H:%M:%S')),*message,end=end)
            self.log_text = self.log_text + str(time.strftime('%Y/%m/%d %H:%M:%S')) + " " + " ".join(map(str,message)) + end
        else:
            if verbose: print(*message,end=end)
            self.log_text = self.log_text + " ".join(map(str,message)) + end
    
    def warning(self, *message: Any, end: str = "\n", show_time: bool = True, verbose: bool = True) -> None:
        self.write("#WARNING! {}".format(" ".join(map(str,message))), 
                   end=end, 
                   show_time=show_time,
                   verbose=verbose)

    def show(self) -> None:
        print(self.log_text)

    def save(self, path: str, verbose: bool = True) -> None:
        with open(path,"w") as f:
            if verbose: print(str(time.strftime('%Y/%m/%d %H:%M:%S')) + " " + " -Save log file to : ", path)
            f.write(self.log_text)

    def log(self, *message: Any, end: str = "\n", show_time: bool = True, verbose: bool = True) -> None:
        if show_time is True:
            if verbose: print(str(time.strftime('%Y/%m/%d %H:%M:%S')),*message,end=end)
            self.log_text = self.log_text + str(time.strftime('%Y/%m/%d %H:%M:%S')) + " " + " ".join(map(str,message)) + end
        else:
            if verbose: print(*message,end=end)
            self.log_text = self.log_text + " ".join(map(str,message)) + end

    def get_log_for_last_operation(self) -> str:
        last_log = []
        rows = self.log_text.strip().split("\n")

        # Iterate backwards to find the last operation
        for line in reversed(rows):
            last_log.append(line)
            if "Start to" in line:
                break

        return "".join(reversed(last_log))
    
    def combine(self, log: 'Log', pre: bool = True) -> None:
        if pre ==True:
            self.log_text = "{}\n{}".format(log.log_text, self.log_text)
        else:
            self.log_text = "{}\n{}".format(self.log_text, log.log_text)
    
    # ============================================================================
    # Standardized logging methods for common operations
    # ============================================================================
    
    def _get_indent(self, indent: int = 0, indent_size: int = 2) -> str:
        """
        Generate indentation string.
        
        Parameters
        ----------
        indent : int, default 0
            Indentation level (0 = no indent, 1 = 2 spaces, 2 = 4 spaces, etc.)
        indent_size : int, default 2
            Number of spaces per indent level
        
        Returns
        -------
        str
            Indentation string
        """
        return " " * (indent * indent_size)
    
    def log_variants_filtered(self, count: int, reason: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log filtered variants count with optional reason.
        
        Parameters
        ----------
        count : int
            Number of variants filtered out
        reason : str, optional
            Reason for filtering (e.g., "with P > 0.05", "based on threshold")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if reason:
            self.write(f"{indent_str} -Filtered out variants {reason}: {count}", verbose=verbose)
        else:
            self.write(f"{indent_str} -Filtered out variants: {count}", verbose=verbose)
    
    def log_variants_removed(self, count: int, reason: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log removed variants count with optional reason.
        
        Parameters
        ----------
        count : int
            Number of variants removed
        reason : str, optional
            Reason for removal (e.g., "with nan in CHR column", "duplicate variants")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if reason:
            self.write(f"{indent_str} -Removed variants {reason}: {count}", verbose=verbose)
        else:
            self.write(f"{indent_str} -Removed variants: {count}", verbose=verbose)
    
    def log_variants_kept(self, count: int, reason: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log kept variants count with optional reason.
        
        Parameters
        ----------
        count : int
            Number of variants kept
        reason : str, optional
            Reason/condition for keeping (e.g., "with P < 0.05", "in specified regions")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if reason:
            self.write(f"{indent_str} -Keeping variants {reason}: {count}", verbose=verbose)
        else:
            self.write(f"{indent_str} -Keeping variants: {count}", verbose=verbose)
    
    def log_column_added(self, column_name: str, verbose: bool = True, indent: int = 0) -> None:
        """
        Log addition of a column.
        
        Parameters
        ----------
        column_name : str
            Name of the added column
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str} -Added column: {column_name}", verbose=verbose)
    
    def log_column_dropped(self, column_name: str, reason: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log removal of a column.
        
        Parameters
        ----------
        column_name : str
            Name of the dropped column
        reason : str, optional
            Reason for dropping (e.g., "all values are invalid", "duplicate")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if reason:
            self.write(f"{indent_str} -Dropped column {column_name} ({reason})", verbose=verbose)
        else:
            self.write(f"{indent_str} -Dropped column: {column_name}", verbose=verbose)
    
    def log_shape_change(self, before_shape: Tuple[int, int], after_shape: Tuple[int, int], verbose: bool = True, indent: int = 0) -> None:
        """
        Log DataFrame shape change.
        
        Parameters
        ----------
        before_shape : tuple
            Shape before operation (rows, cols)
        after_shape : tuple
            Shape after operation (rows, cols)
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        rows_before, cols_before = before_shape
        rows_after, cols_after = after_shape
        
        row_change = rows_after - rows_before
        col_change = cols_after - cols_before
        
        if row_change != 0 or col_change != 0:
            changes = []
            if row_change != 0:
                changes.append(f"{row_change:+d} rows")
            if col_change != 0:
                changes.append(f"{col_change:+d} columns")
            change_str = ", ".join(changes)
            self.write(f"{indent_str} -Shape changed: {rows_before} x {cols_before} -> {rows_after} x {cols_after} ({change_str})", verbose=verbose)
    
    def log_dataframe_shape(self, sumstats: Any, verbose: bool = True, indent: int = 0, sumstats_obj: Optional['Sumstats'] = None) -> None:
        """
        Log current DataFrame shape and memory usage.
        Only logs if shape or memory has changed from the last check.
        
        Parameters
        ----------
        sumstats : pandas.DataFrame
            DataFrame to log shape for
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        sumstats_obj : Sumstats, optional
            Sumstats object instance. If provided, checks _last_shape and _last_memory
            attributes to skip logging if both shape and memory are unchanged.
            If None, will try to get from self._sumstats_obj.
        """
        indent_str = self._get_indent(indent)
        try:
            import pandas as pd
            current_shape = (len(sumstats), len(sumstats.columns))
            memory_in_mb = sumstats.memory_usage().sum() / 1024 / 1024
            
            # If sumstats_obj not provided, try to get it from log object
            if sumstats_obj is None:
                sumstats_obj = getattr(self, '_sumstats_obj', None)
            
            # Check if shape or memory changed (if tracking is enabled)
            if sumstats_obj is not None and hasattr(sumstats_obj, '_last_shape'):
                # Check both shape and memory
                shape_unchanged = sumstats_obj._last_shape == current_shape
                
                # If shape unchanged, check memory (if tracking is enabled)
                if shape_unchanged:
                    if hasattr(sumstats_obj, '_last_memory'):
                        # Compare memory with small tolerance (0.01 MB) to account for floating point precision
                        memory_unchanged = abs(sumstats_obj._last_memory - memory_in_mb) < 0.01
                        if memory_unchanged:
                            # Shape and memory unchanged, skip logging
                            return
                    else:
                        # Shape unchanged and memory tracking not enabled, skip logging
                        return
                
                # Update tracked shape and memory
                sumstats_obj._last_shape = current_shape
                if hasattr(sumstats_obj, '_last_memory'):
                    sumstats_obj._last_memory = memory_in_mb
            
            self.write(f"{indent_str} -Current Dataframe shape : {len(sumstats)} x {len(sumstats.columns)} ; Memory usage: {memory_in_mb:.2f} MB", verbose=verbose)
        except Exception:
            self.warning("Error: cannot get Dataframe shape...", verbose=verbose)
    
    def log_operation_start(self, operation_name: str, version: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log start of an operation with optional version.
        
        Parameters
        ----------
        operation_name : str
            Description of the operation (e.g., "filter variants by condition")
        version : str, optional
            Version string to append
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if version:
            self.write(f"{indent_str}Start to {operation_name} ...({version})", verbose=verbose)
        else:
            self.write(f"{indent_str}Start to {operation_name} ...", verbose=verbose)
    
    def log_operation_finish(self, operation_name: str, verbose: bool = True, indent: int = 0) -> None:
        """
        Log completion of an operation.
        
        Parameters
        ----------
        operation_name : str
            Description of the operation (e.g., "filtering variants")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str}Finished {operation_name}.", verbose=verbose)
    
    def log_filtering_condition(self, column: str, operator: str, threshold: Any, count: int, action: str = "Removing", verbose: bool = True, indent: int = 0) -> None:
        """
        Log filtering condition with count of affected variants.
        
        Parameters
        ----------
        column : str
            Column name being filtered
        operator : str
            Comparison operator (">", "<", "=", ">=", "<=")
        threshold : str or numeric
            Threshold value
        count : int
            Number of variants affected
        action : str, default "Removing"
            Action being taken ("Removing", "Keeping", etc.)
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str} -{action} variants with {column} {operator} {threshold}: {count}", verbose=verbose)
    
    def log_operation(self, message: str, prefix: str = " -", verbose: bool = True, indent: int = 0) -> None:
        """
        Log a general operation message with optional prefix.
        
        Parameters
        ----------
        message : str
            Operation message to log
        prefix : str, default " -"
            Prefix to add before the message
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str}{prefix}{message}", verbose=verbose)
    
    def log_status_change(self, digit: int, before: Any, after: Any, count: Optional[int] = None, reason: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log STATUS column changes (digit position updates).
        
        Parameters
        ----------
        digit : int
            Digit position being changed (1-indexed, 1=leftmost, 7=rightmost)
        before : str or list
            Original digit value(s) (e.g., "9" or ["4", "5"])
        after : str or list
            New digit value(s) (e.g., "1" or ["1", "2"])
        count : int, optional
            Number of variants affected
        reason : str, optional
            Reason for the status change (e.g., "genome build update")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if isinstance(before, list):
            before_str = "/".join(str(b) for b in before)
        else:
            before_str = str(before)
        
        if isinstance(after, list):
            after_str = "/".join(str(a) for a in after)
        else:
            after_str = str(after)
        
        msg_parts = [f"{indent_str} -Updated STATUS digit {digit}: {before_str} -> {after_str}"]
        
        if count is not None:
            msg_parts.append(f"({count} variants)")
        
        if reason:
            msg_parts.append(f"- {reason}")
        
        self.write(" ".join(msg_parts), verbose=verbose)
    
    def log_datatype_change(self, column: str, from_dtype: str, to_dtype: str, success: bool = True, status: Optional[str] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log datatype conversion for a column.
        
        This function handles both conversion attempts (before conversion) and 
        conversion results (after conversion, with success/failure status).
        
        Parameters
        ----------
        column : str
            Column name being converted
        from_dtype : str
            Original datatype
        to_dtype : str
            Target datatype
        success : bool, default True
            Whether the conversion was successful (used when status is None)
        status : str, optional
            Status of conversion: "attempt" (before conversion), "success" (after successful conversion),
            or "failed" (after failed conversion). If None, uses success parameter to determine status.
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        
        # Determine status
        if status is None:
            status = "success" if success else "failed"
        
        # Log based on status
        if status == "attempt":
            self.write(f"{indent_str} -Trying to convert datatype for {column}: {from_dtype} -> {to_dtype}...", verbose=verbose)
        elif status == "success":
            self.write(f"{indent_str} -Converted datatype for {column}: {from_dtype} -> {to_dtype}", verbose=verbose)
        elif status == "failed":
            self.write(f"{indent_str} -Failed to convert datatype for {column}: {from_dtype} -> {to_dtype}", verbose=verbose)
        else:
            # Fallback for unknown status
            self.write(f"{indent_str} -Datatype conversion for {column}: {from_dtype} -> {to_dtype} ({status})", verbose=verbose)
    
    def log_formula(self, target_column: str, formula: str, source_columns: Optional[list] = None, verbose: bool = True, indent: int = 0) -> None:
        """
        Log formula/calculation used to fill or compute a column.
        
        Parameters
        ----------
        target_column : str
            Column being filled/computed
        formula : str
            Formula description (e.g., "from P", "from Z", "from BETA/SE")
        source_columns : list, optional
            Source columns used in the calculation
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        if source_columns:
            sources = ", ".join(source_columns)
            self.write(f"{indent_str}    Filling {target_column} {formula} (using {sources})...", verbose=verbose)
        else:
            self.write(f"{indent_str}    Filling {target_column} {formula}...", verbose=verbose)
    
    def log_reference_path(self, ref_type: str, path: str, verbose: bool = True, indent: int = 0) -> None:
        """
        Log reference file path (VCF, FASTA, TSV, etc.).
        
        Parameters
        ----------
        ref_type : str
            Type of reference ("VCF", "FASTA", "TSV", "BED", etc.)
        path : str
            Path to the reference file
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str} -Reference {ref_type}: {path}", verbose=verbose)
    
    def log_threads(self, threads: int, verbose: bool = True, indent: int = 0) -> None:
        """
        Log number of threads/cores being used.
        
        Parameters
        ----------
        threads : int
            Number of threads/cores
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str} -Number of threads/cores to use: {threads}", verbose=verbose)
    
    def log_variants_with_condition(self, condition: str, count: int, verbose: bool = True, indent: int = 0) -> None:
        """
        Log count of variants matching a specific condition.
        
        Parameters
        ----------
        condition : str
            Condition description (e.g., "standardized chromosome notation", "fixable chromosome notations")
        count : int
            Number of variants matching the condition
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str} -Variants with {condition}: {count}", verbose=verbose)
    
    def log_variants_count(self, count: int, description: str = "variants", verbose: bool = True, indent: int = 0) -> None:
        """
        Log a count of variants with optional description.
        
        Parameters
        ----------
        count : int
            Number of variants
        description : str, default "variants"
            Description of what is being counted (e.g., "variants could be fixed", "variants to check")
        verbose : bool, default True
            Whether to print to stdout
        indent : int, default 0
            Indentation level for nested operations (0 = no indent, 1 = 2 spaces, etc.)
        """
        indent_str = self._get_indent(indent)
        self.write(f"{indent_str} -Number of {description}: {count}", verbose=verbose)