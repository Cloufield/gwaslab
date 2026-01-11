"""
Python Script Execution Framework

Provides a centralized, robust framework for executing Python scripts with:
- Timeout support
- Proper temporary file management
- Structured error handling
- Safe command execution
- Result validation
- Conda environment support
"""

import subprocess
import tempfile
import os
import time
import shutil
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any, Tuple

from gwaslab.info.g_Log import Log


@dataclass
class PythonExecutionResult:
    """Standardized result structure for Python script execution."""
    success: bool
    output: str  # stdout/stderr combined
    exit_code: int
    output_files: Dict[str, str] = field(default_factory=dict)  # expected_file -> actual_path
    errors: List[str] = field(default_factory=list)
    execution_time: float = 0.0
    script_path: Optional[str] = None
    temp_dir: Optional[str] = None


class PythonScriptRunner:
    """
    Centralized Python script execution with timeout, proper temp file management,
    and structured error handling.
    """
    
    def __init__(
        self,
        python: str = "python",
        log: Optional[Log] = None,
        timeout: Optional[float] = None,
        temp_dir: Optional[str] = None,
        cleanup: bool = True,
        conda_env: Optional[str] = None
    ):
        """
        Initialize Python script runner.
        
        Args:
            python: Path to Python executable (default: "python")
            log: Log instance for logging (default: None, creates new Log)
            timeout: Default timeout in seconds (default: None, no timeout)
            temp_dir: Directory for temporary files (default: None, uses system temp)
            cleanup: Whether to clean up temp files after execution (default: True)
            conda_env: Optional conda environment name to activate before execution
        """
        self.python = python
        self.log = log if log is not None else Log()
        self.timeout = timeout
        self.temp_dir = temp_dir
        self.cleanup = cleanup
        self.conda_env = conda_env
    
    def execute(
        self,
        script_content: str,
        expected_outputs: Optional[List[str]] = None,
        temp_prefix: str = "gwaslab_python_script",
        temp_suffix: str = ".py",
        timeout: Optional[float] = None,
        verbose: bool = True,
        working_dir: Optional[str] = None,
        python_args: Optional[List[str]] = None
    ) -> PythonExecutionResult:
        """
        Execute a Python script with proper error handling and validation.
        
        Args:
            script_content: The Python script content as a string
            expected_outputs: List of expected output file names (relative paths)
            temp_prefix: Prefix for temporary script file
            temp_suffix: Suffix for temporary script file (default: ".py")
            timeout: Override default timeout (in seconds)
            verbose: Whether to log verbosely
            working_dir: Working directory for script execution (default: temp_dir or current dir)
            python_args: Additional arguments to pass to Python (e.g., ["-u"] for unbuffered)
        
        Returns:
            PythonExecutionResult with execution details
        """
        if expected_outputs is None:
            expected_outputs = []
        
        if python_args is None:
            python_args = []
        
        timeout = timeout if timeout is not None else self.timeout
        temp_dir = self.temp_dir
        temp_script_path = None
        temp_dir_created = False
        
        result = PythonExecutionResult(
            success=False,
            output="",
            exit_code=-1
        )
        
        start_time = time.time()
        
        try:
            # Create temporary directory if needed
            if temp_dir is None:
                temp_dir = tempfile.mkdtemp(prefix="gwaslab_python_")
                temp_dir_created = True
                result.temp_dir = temp_dir
            else:
                os.makedirs(temp_dir, exist_ok=True)
                result.temp_dir = temp_dir
            
            # Create temporary Python script file
            temp_script_path = self._create_temp_python_script(
                script_content,
                temp_dir,
                temp_prefix,
                temp_suffix
            )
            result.script_path = temp_script_path
            
            # Validate Python script if needed
            validation_errors = self._validate_python_script(script_content, expected_outputs)
            if validation_errors:
                result.errors.extend(validation_errors)
                self.log.warning(f"Python script validation warnings: {validation_errors}", verbose=verbose)
            
            # Log script content (for debugging)
            self.log.write(f" -Executing Python script: {temp_script_path}", verbose=verbose)
            if verbose:
                self.log.write(f" -Script content preview (first 500 chars):\n{script_content[:500]}...", verbose=verbose)
            
            # Execute Python script
            working_dir = working_dir if working_dir is not None else temp_dir
            output, exit_code = self._run_python_script(
                temp_script_path,
                timeout=timeout,
                working_dir=working_dir,
                python_args=python_args,
                verbose=verbose
            )
            
            result.output = output
            result.exit_code = exit_code
            result.execution_time = time.time() - start_time
            
            # Check if execution was successful
            if exit_code == 0:
                # Validate and collect output files
                result.output_files = self._collect_output_files(
                    expected_outputs,
                    working_dir,
                    verbose=verbose
                )
                
                # Check if all expected outputs exist
                missing_outputs = [f for f in expected_outputs if f not in result.output_files]
                if missing_outputs:
                    result.errors.append(f"Missing expected output files: {missing_outputs}")
                    result.success = False
                else:
                    result.success = True
                    self.log.write(f" -Python script executed successfully in {result.execution_time:.2f}s", verbose=verbose)
            else:
                result.success = False
                result.errors.append(f"Python script exited with code {exit_code}")
                self.log.warning(f"Python script execution failed with exit code {exit_code}", verbose=verbose)
                if output:
                    self.log.write(f" -Error output:\n{output}", verbose=verbose)
        
        except subprocess.TimeoutExpired as e:
            result.success = False
            result.execution_time = time.time() - start_time
            result.errors.append(f"Python script execution timed out after {timeout}s")
            result.output = str(e)
            self.log.warning(f"Python script execution timed out after {timeout}s", verbose=verbose)
        
        except FileNotFoundError:
            result.success = False
            result.errors.append(f"Python executable not found: {self.python}")
            self.log.warning(f"Python executable not found: {self.python}", verbose=verbose)
        
        except Exception as e:
            result.success = False
            result.execution_time = time.time() - start_time
            result.errors.append(f"Unexpected error during Python script execution: {str(e)}")
            result.output = str(e)
            self.log.warning(f"Unexpected error during Python script execution: {str(e)}", verbose=verbose)
        
        finally:
            # Cleanup temporary files
            if self.cleanup:
                self._cleanup_temp_files(temp_script_path, temp_dir, temp_dir_created, verbose=verbose)
            else:
                if temp_script_path:
                    result.script_path = temp_script_path
                if temp_dir_created:
                    result.temp_dir = temp_dir
        
        return result
    
    def _create_temp_python_script(
        self,
        script_content: str,
        temp_dir: str,
        prefix: str,
        suffix: str
    ) -> str:
        """Create a temporary Python script file."""
        # Create unique filename with process ID and timestamp
        pid = os.getpid()
        timestamp = int(time.time() * 1000) % 1000000
        unique_name = f"{prefix}_{pid}_{timestamp}{suffix}"
        temp_script_path = os.path.join(temp_dir, unique_name)
        
        with open(temp_script_path, "w") as f:
            f.write(script_content)
        
        # Make script executable on Unix systems
        try:
            os.chmod(temp_script_path, 0o755)
        except Exception:
            pass  # Ignore on Windows
        
        return temp_script_path
    
    def _validate_python_script(
        self,
        script_content: str,
        required_files: List[str]
    ) -> List[str]:
        """
        Validate that Python script references expected files.
        Returns list of validation error messages (empty if valid).
        """
        errors = []
        
        # Basic validation: check if script content is not empty
        if not script_content or not script_content.strip():
            errors.append("Python script content is empty")
        
        # Check for common issues in file paths
        for file_path in required_files:
            if not file_path or not isinstance(file_path, str):
                errors.append(f"Invalid file path: {file_path}")
        
        return errors
    
    def _run_python_script(
        self,
        script_path: str,
        timeout: Optional[float] = None,
        working_dir: Optional[str] = None,
        python_args: Optional[List[str]] = None,
        verbose: bool = True
    ) -> Tuple[str, int]:
        """
        Execute Python script using subprocess.
        Returns (output, exit_code).
        """
        if python_args is None:
            python_args = []
        
        # Build command
        if self.conda_env:
            # Use conda run if conda environment is specified
            cmd = ["conda", "run", "-n", self.conda_env, self.python] + python_args + [script_path]
        else:
            # Use list form to avoid shell injection
            cmd = [self.python] + python_args + [script_path]
        
        try:
            process = subprocess.run(
                cmd,
                cwd=working_dir,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False  # Don't raise on non-zero exit
            )
            
            # Combine stdout and stderr
            output = ""
            if process.stdout:
                output += process.stdout
            if process.stderr:
                output += process.stderr
            
            return output, process.returncode
        
        except subprocess.TimeoutExpired:
            raise
        except FileNotFoundError:
            raise
    
    def _collect_output_files(
        self,
        expected_outputs: List[str],
        working_dir: str,
        verbose: bool = True
    ) -> Dict[str, str]:
        """
        Collect and validate output files.
        Returns dict mapping expected_file -> actual_path.
        """
        output_files = {}
        
        for expected_file in expected_outputs:
            # Try to find the file in working directory
            actual_path = os.path.join(working_dir, expected_file)
            
            # Also try just the expected_file name (in case it's an absolute path)
            if not os.path.exists(actual_path):
                if os.path.isabs(expected_file) and os.path.exists(expected_file):
                    actual_path = expected_file
                else:
                    # Try current directory
                    if os.path.exists(expected_file):
                        actual_path = expected_file
                    else:
                        continue  # File not found
            
            if os.path.exists(actual_path) and os.path.isfile(actual_path):
                output_files[expected_file] = actual_path
                self.log.write(f" -Found output file: {actual_path}", verbose=verbose)
            else:
                self.log.warning(f" -Expected output file not found: {expected_file}", verbose=verbose)
        
        return output_files
    
    def _cleanup_temp_files(
        self,
        temp_script_path: Optional[str],
        temp_dir: Optional[str],
        temp_dir_created: bool,
        verbose: bool = True
    ) -> None:
        """Clean up temporary files and directories."""
        try:
            if temp_script_path and os.path.exists(temp_script_path):
                os.remove(temp_script_path)
                self.log.write(f" -Removed temp Python script: {temp_script_path}", verbose=verbose)
            
            if temp_dir_created and temp_dir and os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir)
                    self.log.write(f" -Removed temp directory: {temp_dir}", verbose=verbose)
                except Exception as e:
                    self.log.warning(f" -Could not remove temp directory {temp_dir}: {e}", verbose=verbose)
        except Exception as e:
            self.log.warning(f" -Error during cleanup: {e}", verbose=verbose)


# Helper functions

def create_temp_python_script(
    script_content: str,
    prefix: str = "gwaslab_python_script",
    suffix: str = ".py",
    temp_dir: Optional[str] = None
) -> str:
    """
    Create a temporary Python script file.
    
    Args:
        script_content: The Python script content
        prefix: Prefix for the temporary file
        suffix: Suffix for the temporary file (default: ".py")
        temp_dir: Directory for temp file (default: system temp)
    
    Returns:
        Path to the created temporary script file
    """
    if temp_dir is None:
        temp_dir = tempfile.gettempdir()
    
    pid = os.getpid()
    timestamp = int(time.time() * 1000) % 1000000
    unique_name = f"{prefix}_{pid}_{timestamp}{suffix}"
    temp_script_path = os.path.join(temp_dir, unique_name)
    
    with open(temp_script_path, "w") as f:
        f.write(script_content)
    
    # Make script executable on Unix systems
    try:
        os.chmod(temp_script_path, 0o755)
    except Exception:
        pass  # Ignore on Windows
    
    return temp_script_path


def validate_python_script(
    script_content: str,
    required_files: List[str]
) -> List[str]:
    """
    Validate that Python script references expected files.
    
    Args:
        script_content: The Python script content
        required_files: List of required file paths
    
    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []
    
    if not script_content or not script_content.strip():
        errors.append("Python script content is empty")
    
    for file_path in required_files:
        if not file_path or not isinstance(file_path, str):
            errors.append(f"Invalid file path: {file_path}")
    
    return errors


def read_python_output_files(
    result: PythonExecutionResult,
    expected_files: List[str]
) -> Dict[str, Any]:
    """
    Read and return output files from Python execution result.
    
    Args:
        result: PythonExecutionResult from script execution
        expected_files: List of expected file names
    
    Returns:
        Dict mapping file names to their contents (as DataFrames for CSV, strings for text)
    """
    import pandas as pd
    
    outputs = {}
    
    for file_name in expected_files:
        if file_name in result.output_files:
            file_path = result.output_files[file_name]
            try:
                # Try to read as CSV first
                if file_path.endswith('.csv'):
                    outputs[file_name] = pd.read_csv(file_path)
                elif file_path.endswith('.tsv') or file_path.endswith('.tsv.gz'):
                    outputs[file_name] = pd.read_csv(file_path, sep='\t')
                else:
                    # Read as text
                    with open(file_path, 'r') as f:
                        outputs[file_name] = f.read()
            except Exception as e:
                outputs[file_name] = None
                result.errors.append(f"Error reading output file {file_name}: {str(e)}")
    
    return outputs
