"""
Simplified Result Management Framework

Provides result tracking, preview, and log file management for script executions.
"""

import os
import pandas as pd
from typing import Optional, Dict, List, Any, Union
from dataclasses import dataclass, field
from datetime import datetime

from gwaslab.util.rwrapper.util_ex_r_runner import RExecutionResult
from gwaslab.util.pwrapper.util_ex_python_runner import PythonExecutionResult
from gwaslab.util.cwrapper.util_ex_command_runner import CommandExecutionResult
from gwaslab.info.g_Log import Log
from gwaslab.util.general.util_wrapper_log import create_r_log, create_python_log, create_command_log
from gwaslab.util.general.util_path_manager import _path

try:
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Union type for all execution results
ExecutionResult = Union[RExecutionResult, PythonExecutionResult, CommandExecutionResult]


@dataclass
class ExecutionRecord:
    """Record of a single execution with its result and metadata."""
    identifier: Optional[str]
    result: ExecutionResult
    timestamp: str
    parameters: Dict[str, Any] = field(default_factory=dict)
    output_data: Dict[str, Any] = field(default_factory=dict)  # Cached output files


class ResultManager:
    """
    Simplified result manager for tracking script executions, showing previews, and managing logs.
    
    Key Features:
    - Trace execution results and track success/failure
    - Show previews of DataFrames and images
    - Create log files for R, Python, and command executions
    - Manage file workflow
    
    Typical Usage:
    -------------
    >>> manager = ResultManager(log=my_log)
    >>> result = runner.execute(...)
    >>> manager.trace(result, identifier="locus_1", parameters={...})
    >>> manager.create_r_log(result, script_content, log_file_path, ...)
    >>> manager.preview_dataframe("output.csv", n_rows=5)
    >>> manager.show_image("diagnostic.png", title="Diagnostic Plot")
    """
    
    def __init__(self, log: Optional[Log] = None):
        """
        Initialize ResultManager.
        
        Args:
            log: Log instance for logging (default: None, creates new Log)
        """
        self.log = log if log is not None else Log()
        self.records: List[ExecutionRecord] = []
    
    def trace(
        self,
        result: ExecutionResult,
        identifier: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
        read_outputs: bool = False,
        expected_files: Optional[List[str]] = None
    ) -> ExecutionRecord:
        """
        Trace an execution result (record it and optionally read outputs).
        
        Args:
            result: ExecutionResult from script execution
            identifier: Optional identifier (e.g., locus ID, study name)
            parameters: Optional parameters used for this execution
            read_outputs: Whether to read output files immediately (default: False)
            expected_files: List of files to read if read_outputs=True
        
        Returns:
            ExecutionRecord containing the result and metadata
        """
        record = ExecutionRecord(
            identifier=identifier,
            result=result,
            timestamp=datetime.now().isoformat(),
            parameters=parameters or {}
        )
        
        # Read outputs if requested
        if read_outputs and expected_files:
            record.output_data = self._read_outputs(result, expected_files)
        
        self.records.append(record)
        
        # Log success/failure
        if result.success:
            self.log.write(
                f"  -Execution successful: {identifier or 'unnamed'}",
                verbose=True
            )
        else:
            self.log.warning(
                f"  -Execution failed: {identifier or 'unnamed'}",
                verbose=True
            )
            if result.errors:
                for error in result.errors:
                    self.log.warning(f"    Error: {error}", verbose=True)
        
        return record
    
    def _read_outputs(
        self,
        result: ExecutionResult,
        expected_files: List[str]
    ) -> Dict[str, Any]:
        """Read output files from a result."""
        outputs = {}
        
        for file_name in expected_files:
            if file_name in result.output_files:
                file_path = result.output_files[file_name]
                try:
                    # Try to read as CSV/TSV
                    if file_path.endswith('.csv') or file_path.endswith('.res'):
                        outputs[file_name] = pd.read_csv(file_path)
                    elif file_path.endswith('.tsv') or file_path.endswith('.tsv.gz'):
                        outputs[file_name] = pd.read_csv(file_path, sep='\t')
                    else:
                        # Read as text
                        with open(file_path, 'r') as f:
                            outputs[file_name] = f.read()
                except Exception as e:
                    outputs[file_name] = None
                    self.log.warning(f"Error reading output file {file_name}: {str(e)}")
            else:
                outputs[file_name] = None
        
        return outputs
    
    def preview_dataframe(
        self,
        file_name: str,
        result: Optional[ExecutionResult] = None,
        n_rows: int = 5,
        show_info: bool = True
    ) -> Optional[pd.DataFrame]:
        """
        Show a preview of a DataFrame from an output file.
        
        Args:
            file_name: Name of the output file (as in result.output_files)
            result: ExecutionResult to read from (if None, uses last successful result)
            n_rows: Number of rows to show in preview (default: 5)
            show_info: Whether to show DataFrame info (default: True)
        
        Returns:
            DataFrame if found, None otherwise
        """
        # Find result to use
        if result is None:
            # Use last successful result
            for record in reversed(self.records):
                if record.result.success and file_name in record.result.output_files:
                    result = record.result
                    break
        
        if result is None or file_name not in result.output_files:
            self.log.warning(f"  -Output file not found: {file_name}", verbose=True)
            return None
        
        # Use the path from output_files directly (already resolved and validated)
        file_path = result.output_files[file_name]
        
        # Verify file exists (should always be true since it's in output_files)
        if not os.path.exists(file_path):
            self.log.warning(f"  -Output file path does not exist: {file_path}", verbose=True)
            return None
        
        try:
            # Read the file - try CSV first, then TSV
            # .pipcs files are CSV files written by R
            if (file_path.endswith('.csv') or file_path.endswith('.res') or 
                file_path.endswith('.pipcs')):
                df = pd.read_csv(file_path)
            elif file_path.endswith('.tsv') or file_path.endswith('.tsv.gz'):
                df = pd.read_csv(file_path, sep='\t')
            else:
                # Try reading as CSV as fallback (many R outputs are CSV even without .csv extension)
                try:
                    df = pd.read_csv(file_path)
                except Exception:
                    self.log.warning(f"  -File {file_name} is not a CSV/TSV file", verbose=True)
                    return None
            
            # Show preview
            self.log.write(f"  -Preview of {file_name}:", verbose=True)
            self.log.write(f"    Shape: {df.shape[0]} rows × {df.shape[1]} columns", verbose=True)
            
            if show_info:
                self.log.write(f"    Columns: {', '.join(df.columns.tolist()[:10])}", verbose=True)
                if len(df.columns) > 10:
                    self.log.write(f"    ... and {len(df.columns) - 10} more columns", verbose=True)
            
            # Show head
            preview_df = df.head(n_rows)
            self.log.write(f"\n    First {n_rows} rows:", verbose=True)
            self.log.write(f"\n{preview_df.to_string()}\n", verbose=True)
            
            return df
            
        except Exception as e:
            self.log.warning(f"  -Error reading {file_name}: {str(e)}", verbose=True)
            return None
    
    def show_image(
        self,
        file_name: str,
        result: Optional[ExecutionResult] = None,
        title: Optional[str] = None,
        figsize: tuple = (10, 6)
    ) -> bool:
        """
        Display an image file using matplotlib.
        
        Args:
            file_name: Name of the image file (as in result.output_files)
            result: ExecutionResult to read from (if None, uses last successful result)
            title: Optional title for the plot
            figsize: Figure size tuple (default: (10, 6))
        
        Returns:
            True if image was displayed, False otherwise
        """
        if not MATPLOTLIB_AVAILABLE:
            self.log.warning("  -matplotlib not available, cannot display image", verbose=True)
            return False
        
        # Find result to use
        if result is None:
            # Use last successful result
            for record in reversed(self.records):
                if record.result.success and file_name in record.result.output_files:
                    result = record.result
                    break
        
        if result is None or file_name not in result.output_files:
            self.log.warning(f"  -Image file not found: {file_name}", verbose=True)
            return False
        
        # Use the path from output_files directly (already resolved and validated by RScriptRunner)
        # The path is already absolute and verified to exist (as shown in log file)
        file_path = result.output_files[file_name]
        
        # Simple existence check (should always pass since RScriptRunner validates it)
        if not os.path.exists(file_path):
            self.log.warning(f"  -Image file does not exist: {file_path}", verbose=True)
            return False
        
        try:
            img = mpimg.imread(file_path)
            plt.figure(figsize=figsize)
            plt.imshow(img)
            plt.axis('off')
            if title:
                plt.title(title)
            plt.tight_layout()
            plt.show()
            self.log.write(f"  -Displayed image: {file_name}", verbose=True)
            return True
        except Exception as e:
            self.log.warning(f"  -Could not display image: {str(e)}", verbose=True)
            return False
    
    def create_r_log(
        self,
        result: ExecutionResult,
        script_content: str,
        log_file_path: str,
        working_dir: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        r_path: str = "Rscript",
        package_name: Optional[str] = None,
        verbose: bool = True
    ) -> str:
        """
        Create a log file for R script execution.
        
        Args:
            result: ExecutionResult from R script execution
            script_content: Content of the R script
            log_file_path: Path to save the log file
            working_dir: Working directory (default: None)
            metadata: Additional metadata (e.g., {"Study": "X", "SNPID": "Y"})
            r_path: Path to Rscript executable (default: "Rscript")
            package_name: R package name for version info (e.g., "susieR")
            verbose: Whether to log verbosely (default: True)
        
        Returns:
            Path to the created log file
        """
        return create_r_log(
            result=result,
            script_content=script_content,
            log_file_path=log_file_path,
            working_dir=working_dir,
            metadata=metadata,
            r_path=r_path,
            package_name=package_name,
            log=self.log,
            verbose=verbose
        )
    
    def create_python_log(
        self,
        result: ExecutionResult,
        script_content: str,
        log_file_path: str,
        working_dir: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        python_path: str = "python",
        package_name: Optional[str] = None,
        verbose: bool = True
    ) -> str:
        """
        Create a log file for Python script execution.
        
        Args:
            result: ExecutionResult from Python script execution
            script_content: Content of the Python script
            log_file_path: Path to save the log file
            working_dir: Working directory (default: None)
            metadata: Additional metadata
            python_path: Path to Python executable (default: "python")
            package_name: Python package name for version info
            verbose: Whether to log verbosely (default: True)
        
        Returns:
            Path to the created log file
        """
        return create_python_log(
            result=result,
            script_content=script_content,
            log_file_path=log_file_path,
            working_dir=working_dir,
            metadata=metadata,
            python_path=python_path,
            package_name=package_name,
            log=self.log,
            verbose=verbose
        )
    
    def create_command_log(
        self,
        result: ExecutionResult,
        command: str,
        log_file_path: str,
        working_dir: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        tool_path: Optional[str] = None,
        verbose: bool = True
    ) -> str:
        """
        Create a log file for command-line execution.
        
        Args:
            result: ExecutionResult from command execution
            command: Command that was executed
            log_file_path: Path to save the log file
            working_dir: Working directory (default: None)
            metadata: Additional metadata
            tool_path: Path to the tool executable
            verbose: Whether to log verbosely (default: True)
        
        Returns:
            Path to the created log file
        """
        return create_command_log(
            result=result,
            command=command,
            log_file_path=log_file_path,
            working_dir=working_dir,
            metadata=metadata,
            tool_path=tool_path,
            log=self.log,
            verbose=verbose
        )
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get execution statistics.
        
        Returns:
            Dict with execution statistics
        """
        if not self.records:
            return {
                "total_executions": 0,
                "successful_executions": 0,
                "failed_executions": 0,
                "success_rate": 0.0,
                "total_execution_time": 0.0,
                "average_execution_time": 0.0
            }
        
        total = len(self.records)
        successful = sum(1 for r in self.records if r.result.success)
        failed = total - successful
        total_time = sum(r.result.execution_time for r in self.records)
        avg_time = total_time / total if total > 0 else 0.0
        
        return {
            "total_executions": total,
            "successful_executions": successful,
            "failed_executions": failed,
            "success_rate": successful / total if total > 0 else 0.0,
            "total_execution_time": total_time,
            "average_execution_time": avg_time
        }
    
    def get_failed_executions(self) -> List[ExecutionRecord]:
        """Get all failed execution records."""
        return [r for r in self.records if not r.result.success]
    
    def get_successful_executions(self) -> List[ExecutionRecord]:
        """Get all successful execution records."""
        return [r for r in self.records if r.result.success]
    
    def clear(self) -> None:
        """Clear all stored records."""
        self.records.clear()