"""
Wrapper Logging Framework

Provides a centralized logging framework for R, Python, and command-line script execution.
Creates comprehensive log files with version information, execution details, and outputs.
"""

import os
import time
import subprocess
from datetime import datetime
from typing import Optional, Dict, List, Any, Union
from dataclasses import dataclass

from gwaslab.info.g_Log import Log
from gwaslab.info.g_version import gwaslab_info


@dataclass
class ExecutionResult:
    """Generic execution result structure compatible with R, Python, and Command runners."""
    success: bool
    output: str
    exit_code: int
    output_files: Dict[str, str] = None
    errors: List[str] = None
    execution_time: float = 0.0
    script_path: Optional[str] = None
    temp_dir: Optional[str] = None
    command: Optional[str] = None  # For command-line execution


class WrapperLogger:
    """
    Centralized logging for script execution (R, Python, command-line).
    Creates comprehensive log files with version information, execution details, and outputs.
    """
    
    def __init__(
        self,
        log: Optional[Log] = None,
        verbose: bool = True
    ):
        """
        Initialize wrapper logger.
        
        Args:
            log: Log instance for logging (default: None, creates new Log)
            verbose: Whether to log verbosely (default: True)
        """
        self.log = log if log is not None else Log()
        self.verbose = verbose
    
    def create_log_file(
        self,
        result: ExecutionResult,
        script_type: str,  # "R", "Python", or "Command"
        script_content: Optional[str] = None,
        command: Optional[str] = None,
        log_file_path: str = None,
        working_dir: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        version_info: Optional[Dict[str, str]] = None,
        tool_path: Optional[str] = None,
        package_name: Optional[str] = None
    ) -> str:
        """
        Create a comprehensive log file for script execution.
        
        Args:
            result: Execution result (from RScriptRunner, PythonScriptRunner, or CommandRunner)
            script_type: Type of script ("R", "Python", or "Command")
            script_content: Content of the script (for R/Python scripts)
            command: Command that was executed (for command-line)
            log_file_path: Path to the log file (if None, will be generated)
            working_dir: Working directory where log should be saved
            metadata: Additional metadata to include (e.g., {"Study": "X", "SNPID": "Y"})
            version_info: Custom version information dict (if None, will be auto-detected)
            tool_path: Path to the tool executable (e.g., "Rscript", "python", "plink")
            package_name: Package name to check version (e.g., "susieR", "pandas")
        
        Returns:
            Path to the created log file
        """
        if log_file_path is None:
            # Generate unique log file name
            timestamp = int(time.time() * 1000) % 1000000000
            pid = os.getpid()
            log_file_path = os.path.join(
                working_dir or os.getcwd(),
                f"gwaslab_{script_type.lower()}_{pid}_{timestamp}.log"
            )
        
        # Ensure working directory exists
        if working_dir:
            os.makedirs(working_dir, exist_ok=True)
        
        # Get version information
        if version_info is None:
            version_info = self._get_version_info(
                script_type=script_type,
                tool_path=tool_path,
                package_name=package_name
            )
        
        # Convert result to ExecutionResult if needed
        exec_result = self._normalize_result(result)
        
        try:
            with open(log_file_path, "w") as f:
                # Write version information
                self._write_version_section(f, version_info)
                
                # Write execution time
                self._write_execution_time_section(f, exec_result)
                
                # Write metadata if provided
                if metadata:
                    self._write_metadata_section(f, metadata)
                
                # Write output files
                if exec_result.output_files:
                    self._write_output_files_section(f, exec_result.output_files, working_dir)
                
                # Write script/command content
                if script_type.upper() == "COMMAND":
                    self._write_command_section(f, command or exec_result.command or "N/A")
                else:
                    script_label = f"{script_type.upper()} SCRIPT"
                    self._write_script_section(f, script_content or "N/A", script_label)
                
                # Write execution output
                self._write_output_section(f, exec_result.output)
            
            self.log.write(f"  -Log file saved to: {log_file_path}", verbose=self.verbose)
            return log_file_path
        
        except Exception as e:
            self.log.warning(f"  -Could not save log file: {str(e)}", verbose=self.verbose)
            return log_file_path
    
    def _get_version_info(
        self,
        script_type: str,
        tool_path: Optional[str] = None,
        package_name: Optional[str] = None
    ) -> Dict[str, str]:
        """Get version information for the execution environment."""
        version_info = {}
        
        # Get GWASLab version (always)
        try:
            version_info["GWASLab"] = gwaslab_info()["version"]
        except Exception:
            version_info["GWASLab"] = "Unable to determine"
        
        # Get tool version based on script type
        if script_type.upper() == "R":
            if tool_path:
                try:
                    r_version_output = subprocess.check_output(
                        f"{tool_path} --version",
                        stderr=subprocess.STDOUT,
                        shell=True,
                        text=True,
                        timeout=5
                    )
                    version_info["R"] = r_version_output.strip()
                except Exception as e:
                    version_info["R"] = f"Unable to determine ({str(e)})"
            
            if package_name:
                try:
                    rscript_version = f'cat(as.character(packageVersion("{package_name}")))'
                    import tempfile
                    import random
                    temp_r_version = tempfile.NamedTemporaryFile(
                        mode='w',
                        suffix='.R',
                        delete=False,
                        prefix=f"_gwaslab_version_check_{random.randint(1, 99999999)}_"
                    )
                    temp_r_version.write(rscript_version)
                    temp_r_version.close()
                    
                    susie_output = subprocess.check_output(
                        f"{tool_path or 'Rscript'} {temp_r_version.name}",
                        stderr=subprocess.STDOUT,
                        shell=True,
                        text=True,
                        timeout=5
                    )
                    version_info[package_name] = susie_output.strip()
                    os.remove(temp_r_version.name)
                except Exception as e:
                    version_info[package_name] = f"Unable to determine ({str(e)})"
                    if os.path.exists(temp_r_version.name):
                        try:
                            os.remove(temp_r_version.name)
                        except:
                            pass
        
        elif script_type.upper() == "PYTHON":
            if tool_path:
                try:
                    python_version_output = subprocess.check_output(
                        [tool_path, "--version"],
                        stderr=subprocess.STDOUT,
                        text=True,
                        timeout=5
                    )
                    version_info["Python"] = python_version_output.strip()
                except Exception:
                    try:
                        python_version_output = subprocess.check_output(
                            f"{tool_path} --version",
                            stderr=subprocess.STDOUT,
                            shell=True,
                            text=True,
                            timeout=5
                        )
                        version_info["Python"] = python_version_output.strip()
                    except Exception as e:
                        version_info["Python"] = f"Unable to determine ({str(e)})"
            
            if package_name:
                try:
                    import importlib.metadata
                    version_info[package_name] = importlib.metadata.version(package_name)
                except Exception:
                    try:
                        import pkg_resources
                        version_info[package_name] = pkg_resources.get_distribution(package_name).version
                    except Exception as e:
                        version_info[package_name] = f"Unable to determine ({str(e)})"
        
        elif script_type.upper() == "COMMAND":
            if tool_path:
                try:
                    # Try common version flags
                    for flag in ["--version", "-v", "-V", "version"]:
                        try:
                            cmd_output = subprocess.check_output(
                                [tool_path, flag],
                                stderr=subprocess.STDOUT,
                                text=True,
                                timeout=5
                            )
                            version_info[tool_path] = cmd_output.strip()[:200]  # Limit length
                            break
                        except:
                            continue
                    else:
                        version_info[tool_path] = "Unable to determine version"
                except Exception as e:
                    version_info[tool_path] = f"Unable to determine ({str(e)})"
        
        return version_info
    
    def _normalize_result(self, result: Any) -> ExecutionResult:
        """Convert any execution result type to ExecutionResult."""
        if isinstance(result, ExecutionResult):
            return result
        
        # Handle dataclass results from runners
        return ExecutionResult(
            success=getattr(result, 'success', False),
            output=getattr(result, 'output', ''),
            exit_code=getattr(result, 'exit_code', -1),
            output_files=getattr(result, 'output_files', {}),
            errors=getattr(result, 'errors', []),
            execution_time=getattr(result, 'execution_time', 0.0),
            script_path=getattr(result, 'script_path', None),
            temp_dir=getattr(result, 'temp_dir', None),
            command=getattr(result, 'command', None)
        )
    
    def _write_version_section(self, f, version_info: Dict[str, str]):
        """Write version information section."""
        f.write("=" * 80 + "\n")
        f.write("VERSION INFORMATION\n")
        f.write("=" * 80 + "\n")
        for tool, version in version_info.items():
            f.write(f"{tool} version: {version}\n")
        f.write("=" * 80 + "\n")
        f.write("\n")
    
    def _write_execution_time_section(self, f, result: ExecutionResult):
        """Write execution time section."""
        execution_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write("=" * 80 + "\n")
        f.write("EXECUTION TIME\n")
        f.write("=" * 80 + "\n")
        f.write(f"Execution started: {execution_time}\n")
        if result.execution_time > 0:
            f.write(f"Execution duration: {result.execution_time:.2f} seconds\n")
        f.write(f"Exit code: {result.exit_code}\n")
        f.write(f"Success: {result.success}\n")
        if result.errors:
            f.write(f"Errors: {len(result.errors)}\n")
        f.write("=" * 80 + "\n")
        f.write("\n")
    
    def _write_metadata_section(self, f, metadata: Dict[str, Any]):
        """Write metadata section."""
        f.write("=" * 80 + "\n")
        f.write("METADATA\n")
        f.write("=" * 80 + "\n")
        for key, value in metadata.items():
            f.write(f"{key}: {value}\n")
        f.write("=" * 80 + "\n")
        f.write("\n")
    
    def _write_output_files_section(
        self,
        f,
        output_files: Dict[str, str],
        working_dir: Optional[str] = None
    ):
        """Write output files section."""
        f.write("=" * 80 + "\n")
        f.write("OUTPUT FILES\n")
        f.write("=" * 80 + "\n")
        if output_files:
            f.write("Generated files:\n")
            for expected_file, actual_path in output_files.items():
                try:
                    file_size = os.path.getsize(actual_path) if os.path.exists(actual_path) else 0
                    file_size_mb = file_size / (1024 * 1024) if file_size > 0 else 0
                    f.write(f"  - {expected_file}\n")
                    f.write(f"    Path: {os.path.abspath(actual_path)}\n")
                    f.write(f"    Size: {file_size_mb:.2f} MB ({file_size} bytes)\n")
                    f.write(f"    Exists: {os.path.exists(actual_path)}\n")
                except Exception as e:
                    f.write(f"  - {expected_file}\n")
                    f.write(f"    Path: {actual_path}\n")
                    f.write(f"    Error getting file info: {str(e)}\n")
        else:
            f.write("No output files specified.\n")
        f.write("=" * 80 + "\n")
        f.write("\n")
    
    def _write_script_section(self, f, script_content: str, script_label: str = "SCRIPT"):
        """Write script content section."""
        f.write("=" * 80 + "\n")
        f.write(f"{script_label}\n")
        f.write("=" * 80 + "\n")
        f.write(script_content)
        f.write("\n\n")
    
    def _write_command_section(self, f, command: str):
        """Write command section."""
        f.write("=" * 80 + "\n")
        f.write("COMMAND\n")
        f.write("=" * 80 + "\n")
        f.write(f"{command}\n")
        f.write("\n")
    
    def _write_output_section(self, f, output: str):
        """Write execution output section."""
        f.write("=" * 80 + "\n")
        f.write("EXECUTION OUTPUT\n")
        f.write("=" * 80 + "\n")
        if output:
            f.write(output)
        else:
            f.write("No output from execution.\n")
        f.write("\n")


# Convenience functions
def create_r_log(
    result: Any,
    script_content: str,
    log_file_path: str,
    working_dir: Optional[str] = None,
    metadata: Optional[Dict[str, Any]] = None,
    r_path: str = "Rscript",
    package_name: Optional[str] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> str:
    """Create a log file for R script execution."""
    logger = WrapperLogger(log=log, verbose=verbose)
    return logger.create_log_file(
        result=result,
        script_type="R",
        script_content=script_content,
        log_file_path=log_file_path,
        working_dir=working_dir,
        metadata=metadata,
        tool_path=r_path,
        package_name=package_name
    )


def create_python_log(
    result: Any,
    script_content: str,
    log_file_path: str,
    working_dir: Optional[str] = None,
    metadata: Optional[Dict[str, Any]] = None,
    python_path: str = "python",
    package_name: Optional[str] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> str:
    """Create a log file for Python script execution."""
    logger = WrapperLogger(log=log, verbose=verbose)
    return logger.create_log_file(
        result=result,
        script_type="Python",
        script_content=script_content,
        log_file_path=log_file_path,
        working_dir=working_dir,
        metadata=metadata,
        tool_path=python_path,
        package_name=package_name
    )


def create_command_log(
    result: Any,
    command: str,
    log_file_path: str,
    working_dir: Optional[str] = None,
    metadata: Optional[Dict[str, Any]] = None,
    tool_path: Optional[str] = None,
    log: Optional[Log] = None,
    verbose: bool = True
) -> str:
    """Create a log file for command-line execution."""
    logger = WrapperLogger(log=log, verbose=verbose)
    return logger.create_log_file(
        result=result,
        script_type="Command",
        command=command,
        log_file_path=log_file_path,
        working_dir=working_dir,
        metadata=metadata,
        tool_path=tool_path
    )
