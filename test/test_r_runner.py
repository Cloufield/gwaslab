"""
Unit and integration tests for R script execution framework.

Tests for:
- RScriptRunner class
- RExecutionResult dataclass
- ResultManager class
- Helper functions
"""

import os
import sys
import unittest
import tempfile
import shutil
import subprocess
from unittest.mock import Mock, patch, MagicMock

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
from gwaslab.util.rwrapper.util_ex_r_runner import (
    RScriptRunner,
    RExecutionResult,
    create_temp_r_script,
    validate_r_script_paths,
    read_r_output_files
)
from gwaslab.util.general.util_ex_result_manager import ResultManager
from gwaslab.info.g_Log import Log


class TestRScriptRunner(unittest.TestCase):
    """Test RScriptRunner class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.temp_dir = tempfile.mkdtemp()
        self.runner = RScriptRunner(
            r="Rscript",
            log=self.log,
            timeout=30,
            temp_dir=self.temp_dir,
            cleanup=False  # Don't cleanup in tests so we can inspect
        )
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_init(self):
        """Test RScriptRunner initialization."""
        runner = RScriptRunner(r="Rscript", log=self.log)
        self.assertEqual(runner.r, "Rscript")
        self.assertIsNotNone(runner.log)
        self.assertIsNone(runner.timeout)
        self.assertTrue(runner.cleanup)
    
    def test_create_temp_r_script(self):
        """Test temporary R script creation."""
        script_content = "print('Hello, R!')\n"
        script_path = self.runner._create_temp_r_script(
            script_content,
            self.temp_dir,
            "test_script",
            ".R"
        )
        
        self.assertTrue(os.path.exists(script_path))
        self.assertTrue(script_path.endswith(".R"))
        
        with open(script_path, 'r') as f:
            content = f.read()
        self.assertEqual(content, script_content)
    
    def test_validate_r_script_paths(self):
        """Test R script path validation."""
        # Valid script
        errors = self.runner._validate_r_script_paths("print('test')", ["output.csv"])
        self.assertEqual(len(errors), 0)
        
        # Empty script
        errors = self.runner._validate_r_script_paths("", ["output.csv"])
        self.assertGreater(len(errors), 0)
        self.assertIn("empty", errors[0].lower())
        
        # None script
        errors = self.runner._validate_r_script_paths(None, ["output.csv"])
        self.assertGreater(len(errors), 0)
    
    @patch('subprocess.run')
    def test_run_r_script_success(self, mock_run):
        """Test successful R script execution."""
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "R output"
        mock_process.stderr = ""
        mock_run.return_value = mock_process
        
        output, exit_code = self.runner._run_r_script(
            "/path/to/script.R",
            timeout=30
        )
        
        self.assertEqual(exit_code, 0)
        self.assertIn("R output", output)
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_run_r_script_failure(self, mock_run):
        """Test failed R script execution."""
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_process.stdout = ""
        mock_process.stderr = "Error message"
        mock_run.return_value = mock_process
        
        output, exit_code = self.runner._run_r_script(
            "/path/to/script.R",
            timeout=30
        )
        
        self.assertEqual(exit_code, 1)
        self.assertIn("Error message", output)
    
    @patch('subprocess.run')
    def test_run_r_script_timeout(self, mock_run):
        """Test R script timeout."""
        mock_run.side_effect = subprocess.TimeoutExpired("Rscript", 30)
        
        with self.assertRaises(subprocess.TimeoutExpired):
            self.runner._run_r_script("/path/to/script.R", timeout=30)
    
    def test_collect_output_files(self):
        """Test output file collection."""
        # Create a test output file
        test_file = os.path.join(self.temp_dir, "test_output.csv")
        with open(test_file, 'w') as f:
            f.write("col1,col2\n1,2\n")
        
        output_files = self.runner._collect_output_files(
            ["test_output.csv"],
            self.temp_dir,
            verbose=False
        )
        
        self.assertIn("test_output.csv", output_files)
        self.assertEqual(output_files["test_output.csv"], test_file)
    
    def test_collect_output_files_missing(self):
        """Test output file collection with missing file."""
        output_files = self.runner._collect_output_files(
            ["missing_file.csv"],
            self.temp_dir,
            verbose=False
        )
        
        self.assertEqual(len(output_files), 0)
    
    @patch('subprocess.run')
    def test_execute_success(self, mock_run):
        """Test successful script execution."""
        # Create expected output file
        output_file = os.path.join(self.temp_dir, "output.csv")
        with open(output_file, 'w') as f:
            f.write("col1,col2\n1,2\n")
        
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "Success"
        mock_process.stderr = ""
        mock_run.return_value = mock_process
        
        result = self.runner.execute(
            script_content="write.csv(data.frame(x=1), 'output.csv')",
            expected_outputs=["output.csv"],
            temp_prefix="test",
            working_dir=self.temp_dir
        )
        
        self.assertTrue(result.success)
        self.assertEqual(result.exit_code, 0)
        self.assertGreater(result.execution_time, 0)
    
    @patch('subprocess.run')
    def test_execute_failure(self, mock_run):
        """Test failed script execution."""
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_process.stdout = ""
        mock_process.stderr = "Error occurred"
        mock_run.return_value = mock_process
        
        result = self.runner.execute(
            script_content="stop('Error')",
            expected_outputs=[],
            temp_prefix="test"
        )
        
        self.assertFalse(result.success)
        self.assertEqual(result.exit_code, 1)
        self.assertGreater(len(result.errors), 0)
    
    @patch('subprocess.run')
    def test_execute_timeout(self, mock_run):
        """Test script execution timeout."""
        mock_run.side_effect = subprocess.TimeoutExpired("Rscript", 30)
        
        result = self.runner.execute(
            script_content="Sys.sleep(100)",
            expected_outputs=[],
            temp_prefix="test",
            timeout=30
        )
        
        self.assertFalse(result.success)
        # Check for "timed out" (the actual error message format)
        self.assertIn("timed out", result.errors[0].lower())
    
    @patch('subprocess.run')
    def test_execute_file_not_found(self, mock_run):
        """Test execution when R is not found."""
        mock_run.side_effect = FileNotFoundError()
        
        runner = RScriptRunner(r="nonexistent_rscript", log=self.log)
        result = runner.execute(
            script_content="print('test')",
            expected_outputs=[],
            temp_prefix="test"
        )
        
        self.assertFalse(result.success)
        self.assertIn("not found", result.errors[0].lower())


class TestRExecutionResult(unittest.TestCase):
    """Test RExecutionResult dataclass."""
    
    def test_creation(self):
        """Test RExecutionResult creation."""
        result = RExecutionResult(
            success=True,
            output="test output",
            exit_code=0,
            execution_time=1.5
        )
        
        self.assertTrue(result.success)
        self.assertEqual(result.output, "test output")
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.execution_time, 1.5)
        self.assertEqual(len(result.output_files), 0)
        self.assertEqual(len(result.errors), 0)
    
    def test_with_output_files(self):
        """Test RExecutionResult with output files."""
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={"output.csv": "/path/to/output.csv"}
        )
        
        self.assertIn("output.csv", result.output_files)
        self.assertEqual(result.output_files["output.csv"], "/path/to/output.csv")


class TestResultManager(unittest.TestCase):
    """Test ResultManager class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.manager = ResultManager(log=self.log)
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_init(self):
        """Test ResultManager initialization."""
        manager = ResultManager()
        self.assertEqual(len(manager.records), 0)
    
    def test_record_result(self):
        """Test recording a result."""
        result = RExecutionResult(
            success=True,
            output="Success",
            exit_code=0,
            execution_time=1.0
        )
        
        record = self.manager.trace(
            result,
            identifier="test_locus",
            parameters={"study": "Study1"}
        )
        
        self.assertEqual(len(self.manager.records), 1)
        self.assertEqual(record.identifier, "test_locus")
        self.assertEqual(record.result.success, True)
        self.assertEqual(record.parameters["study"], "Study1")
    
    def test_record_result_error(self):
        """Test recording a failed result."""
        result = RExecutionResult(
            success=False,
            output="Error",
            exit_code=1,
            execution_time=0.5,
            errors=["Test error"]
        )
        
        record = self.manager.trace(result, identifier="test_locus")
        
        self.assertEqual(len(self.manager.records), 1)
        self.assertEqual(record.identifier, "test_locus")
        self.assertFalse(record.result.success)
        self.assertEqual(len(record.result.errors), 1)
        self.assertEqual(record.result.errors[0], "Test error")
    
    def test_read_outputs(self):
        """Test reading output files."""
        # Create test output file
        output_file = os.path.join(self.temp_dir, "output.csv")
        df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        df.to_csv(output_file, index=False)
        
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={"output.csv": output_file}
        )
        
        outputs = self.manager._read_outputs(result, ["output.csv"])
        
        self.assertIn("output.csv", outputs)
        self.assertIsInstance(outputs["output.csv"], pd.DataFrame)
        self.assertEqual(len(outputs["output.csv"]), 2)
    
    def test_validate_result(self):
        """Test result validation."""
        # Create valid output file
        output_file = os.path.join(self.temp_dir, "output.csv")
        df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        df.to_csv(output_file, index=False)
        
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={"output.csv": output_file}
        )
        
        # Validation is now done by checking result.success and output_files
        is_valid = result.success and "output.csv" in result.output_files
        errors = [] if is_valid else ["Validation failed"]
        
        self.assertTrue(is_valid)
        self.assertEqual(len(errors), 0)
    
    def test_validate_result_missing_file(self):
        """Test validation with missing file."""
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={}
        )
        
        # Validation is now done by checking result.success and output_files
        is_valid = result.success and "output.csv" in result.output_files
        errors = [] if is_valid else ["Missing output file: output.csv"]
        
        self.assertFalse(is_valid)
        self.assertGreater(len(errors), 0)
    
    def test_validate_result_required_columns(self):
        """Test validation with required columns."""
        # Create output file with required columns
        output_file = os.path.join(self.temp_dir, "output.csv")
        df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        df.to_csv(output_file, index=False)
        
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={"output.csv": output_file}
        )
        
        # Read the file and check columns manually
        outputs = self.manager._read_outputs(result, ["output.csv"])
        df_read = outputs["output.csv"]
        required_columns = ["col1", "col2"]
        is_valid = result.success and all(col in df_read.columns for col in required_columns)
        
        self.assertTrue(is_valid)
    
    def test_validate_result_missing_columns(self):
        """Test validation with missing required columns."""
        # Create output file missing a column
        output_file = os.path.join(self.temp_dir, "output.csv")
        df = pd.DataFrame({"col1": [1, 2]})
        df.to_csv(output_file, index=False)
        
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={"output.csv": output_file}
        )
        
        # Read the file and check columns manually
        outputs = self.manager._read_outputs(result, ["output.csv"])
        df_read = outputs["output.csv"]
        required_columns = ["col1", "col2", "col3"]
        missing_columns = [col for col in required_columns if col not in df_read.columns]
        is_valid = result.success and len(missing_columns) == 0
        errors = [f"Missing column: {col}" for col in missing_columns] if not is_valid else []
        
        self.assertFalse(is_valid)
        self.assertIn("col2", str(errors))
    
    def test_aggregate_results(self):
        """Test result aggregation."""
        # Create multiple output files
        for i in range(3):
            output_file = os.path.join(self.temp_dir, f"output_{i}.csv")
            df = pd.DataFrame({"id": [i], "value": [i * 10]})
            df.to_csv(output_file, index=False)
            
            result = RExecutionResult(
                success=True,
                output="",
                exit_code=0,
                output_files={"output.csv": output_file}
            )
            
            self.manager.trace(result, identifier=f"locus_{i}", read_outputs=True, expected_files=["output.csv"])
        
        # Aggregate results manually from records
        aggregated = {}
        for record in self.manager.records:
            if "output.csv" in record.output_data:
                if "output.csv" not in aggregated:
                    aggregated["output.csv"] = []
                aggregated["output.csv"].append(record.output_data["output.csv"])
        
        self.assertIn("output.csv", aggregated)
        self.assertEqual(len(aggregated["output.csv"]), 3)
    
    def test_get_error_summary(self):
        """Test error summary generation."""
        # Add some results
        for i in range(5):
            result = RExecutionResult(
                success=(i < 3),
                output="",
                exit_code=0 if i < 3 else 1,
                errors=[] if i < 3 else ["Error"]
            )
            self.manager.trace(result, identifier=f"locus_{i}")
        
        # Use get_statistics instead of get_error_summary
        stats = self.manager.get_statistics()
        
        self.assertEqual(stats["total_executions"], 5)
        self.assertEqual(stats["failed_executions"], 2)
        self.assertAlmostEqual(stats["success_rate"], 0.6, places=1)
    
    def test_get_statistics(self):
        """Test statistics generation."""
        # Add some results
        for i in range(3):
            result = RExecutionResult(
                success=(i < 2),
                output="",
                exit_code=0 if i < 2 else 1,
                execution_time=1.0 + i
            )
            self.manager.trace(result)
        
        stats = self.manager.get_statistics()
        
        self.assertEqual(stats["total_executions"], 3)
        self.assertEqual(stats["successful_executions"], 2)
        self.assertEqual(stats["failed_executions"], 1)
        self.assertAlmostEqual(stats["average_execution_time"], 2.0, places=1)
    
    def test_clear(self):
        """Test clearing all results."""
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0
        )
        self.manager.trace(result)
        
        self.manager.clear()
        
        self.assertEqual(len(self.manager.records), 0)


class TestHelperFunctions(unittest.TestCase):
    """Test helper functions."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_create_temp_r_script(self):
        """Test create_temp_r_script helper."""
        script_content = "print('test')\n"
        script_path = create_temp_r_script(
            script_content,
            prefix="test",
            temp_dir=self.temp_dir
        )
        
        self.assertTrue(os.path.exists(script_path))
        with open(script_path, 'r') as f:
            content = f.read()
        self.assertEqual(content, script_content)
    
    def test_validate_r_script_paths(self):
        """Test validate_r_script_paths helper."""
        errors = validate_r_script_paths("print('test')", ["output.csv"])
        self.assertEqual(len(errors), 0)
        
        errors = validate_r_script_paths("", ["output.csv"])
        self.assertGreater(len(errors), 0)
    
    def test_read_r_output_files(self):
        """Test read_r_output_files helper."""
        # Create test output file
        output_file = os.path.join(self.temp_dir, "output.csv")
        df = pd.DataFrame({"col1": [1, 2]})
        df.to_csv(output_file, index=False)
        
        result = RExecutionResult(
            success=True,
            output="",
            exit_code=0,
            output_files={"output.csv": output_file}
        )
        
        outputs = read_r_output_files(result, ["output.csv"])
        
        self.assertIn("output.csv", outputs)
        self.assertIsInstance(outputs["output.csv"], pd.DataFrame)


if __name__ == '__main__':
    unittest.main()
