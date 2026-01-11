"""
Comprehensive test suite for _path function and related utilities in util_path_manager.

Tests cover:
- Basic path generation with various components
- Unique ID generation for temporary files
- Path sanitization (invalid characters, special cases)
- Collision detection and handling
- Path length validation
- Edge cases (empty paths, None values, user-provided paths)
- Directory and subdirectory handling
- File extension handling
"""

import os
import sys
import unittest
import tempfile
import shutil
from pathlib import Path

# Add parent directory to path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.general.util_path_manager import (
    _path,
    _sanitize_path_component,
    _generate_unique_id,
    _handle_path_collision,
    _validate_path_length
)
from gwaslab.info.g_Log import Log


class TestPathSanitization(unittest.TestCase):
    """Test path component sanitization."""
    
    def test_sanitize_invalid_characters(self):
        """Test that invalid filesystem characters are replaced."""
        # Test various invalid characters
        test_cases = [
            ("file<name>", "file_name_"),
            ('file"name"', "file_name_"),
            ("file|name", "file_name"),
            ("file?name", "file_name"),
            ("file*name", "file_name"),
            ("file\\name", "file_name"),
            ("file/name", "file_name"),
            ("file:name", "file_name"),
        ]
        
        for input_str, expected_pattern in test_cases:
            result = _sanitize_path_component(input_str)
            # Should not contain invalid characters
            self.assertNotIn("<", result)
            self.assertNotIn(">", result)
            self.assertNotIn('"', result)
            self.assertNotIn("|", result)
            self.assertNotIn("?", result)
            self.assertNotIn("*", result)
            self.assertNotIn("\\", result)
            self.assertNotIn("/", result)
            self.assertNotIn(":", result)
    
    def test_sanitize_consecutive_dashes(self):
        """Test that consecutive dashes/underscores are normalized to underscores."""
        result = _sanitize_path_component("file---name___test")
        # Should not have more than one consecutive underscore
        self.assertNotIn("---", result)
        self.assertNotIn("___", result)
        # Should use underscores
        self.assertIn("_", result)
    
    def test_sanitize_leading_trailing_chars(self):
        """Test that leading/trailing problematic characters are removed."""
        test_cases = [
            (".filename", "filename"),
            ("filename.", "filename"),
            ("-filename", "filename"),
            ("filename-", "filename"),
            (" filename ", "filename"),
        ]
        
        for input_str, expected in test_cases:
            result = _sanitize_path_component(input_str)
            self.assertFalse(result.startswith((".", "-", " ")))
            self.assertFalse(result.endswith((".", "-", " ")))
    
    def test_sanitize_spaces(self):
        """Test that spaces are replaced with underscores."""
        result = _sanitize_path_component("file name with spaces")
        self.assertNotIn(" ", result)
        self.assertIn("_", result)
    
    def test_sanitize_length_limit(self):
        """Test that long components are truncated."""
        long_string = "a" * 300
        result = _sanitize_path_component(long_string)
        self.assertLessEqual(len(result), 255)
    
    def test_sanitize_empty_string(self):
        """Test that empty strings become 'unnamed'."""
        result = _sanitize_path_component("")
        self.assertEqual(result, "unnamed")
    
    def test_sanitize_non_string(self):
        """Test that non-string types are converted."""
        result = _sanitize_path_component(12345)
        self.assertIsInstance(result, str)


class TestUniqueIdGeneration(unittest.TestCase):
    """Test unique ID generation."""
    
    def test_generate_uuid_id(self):
        """Test UUID-based ID generation."""
        id1 = _generate_unique_id(use_uuid=True)
        id2 = _generate_unique_id(use_uuid=True)
        
        # Should be 8 characters (hex)
        self.assertEqual(len(id1), 8)
        self.assertEqual(len(id2), 8)
        
        # Should be different (very high probability)
        self.assertNotEqual(id1, id2)
        
        # Should be hexadecimal
        self.assertTrue(all(c in '0123456789abcdef' for c in id1))
    
    def test_generate_timestamp_id(self):
        """Test timestamp-based ID generation."""
        id1 = _generate_unique_id(use_uuid=False)
        id2 = _generate_unique_id(use_uuid=False)
        
        # Should have format: 9 digits _ 4 digits
        parts = id1.split("_")
        self.assertEqual(len(parts), 2)
        self.assertEqual(len(parts[0]), 9)
        self.assertEqual(len(parts[1]), 4)
        
        # Should be numeric
        self.assertTrue(parts[0].isdigit())
        self.assertTrue(parts[1].isdigit())


class TestPathCollisionHandling(unittest.TestCase):
    """Test path collision detection and handling."""
    
    def setUp(self):
        """Set up temporary directory for tests."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_no_collision(self):
        """Test that non-existent paths are returned unchanged."""
        test_path = os.path.join(self.temp_dir, "nonexistent_file.txt")
        result = _handle_path_collision(test_path)
        self.assertEqual(result, test_path)
    
    def test_collision_handling(self):
        """Test that collisions are handled by appending unique ID."""
        # Create an existing file
        existing_path = os.path.join(self.temp_dir, "existing_file.txt")
        with open(existing_path, 'w') as f:
            f.write("test")
        
        # Should return a different path
        result = _handle_path_collision(existing_path)
        self.assertNotEqual(result, existing_path)
        self.assertFalse(os.path.exists(result))
        
        # Should preserve directory and extension
        self.assertEqual(os.path.dirname(result), os.path.dirname(existing_path))
        self.assertTrue(result.endswith(".txt"))
    
    def test_collision_with_extension(self):
        """Test collision handling preserves file extension."""
        existing_path = os.path.join(self.temp_dir, "test.csv")
        with open(existing_path, 'w') as f:
            f.write("test")
        
        result = _handle_path_collision(existing_path)
        self.assertTrue(result.endswith(".csv"))
        self.assertIn("test_", result)
    
    def test_collision_without_extension(self):
        """Test collision handling for files without extension."""
        existing_path = os.path.join(self.temp_dir, "testfile")
        with open(existing_path, 'w') as f:
            f.write("test")
        
        result = _handle_path_collision(existing_path)
        self.assertNotEqual(result, existing_path)
        self.assertIn("testfile_", result)


class TestPathLengthValidation(unittest.TestCase):
    """Test path length validation."""
    
    def test_normal_length_path(self):
        """Test that normal length paths are unchanged."""
        normal_path = "/path/to/normal/file.txt"
        result = _validate_path_length(normal_path)
        self.assertEqual(result, normal_path)
    
    def test_long_path_truncation(self):
        """Test that very long paths are truncated."""
        # Create a path longer than 4096 characters
        long_basename = "a" * 5000
        long_path = os.path.join("/tmp", long_basename + ".txt")
        
        result = _validate_path_length(long_path)
        self.assertLessEqual(len(result), 4096)
        # Should preserve extension
        self.assertTrue(result.endswith(".txt"))
    
    def test_long_path_preserves_extension(self):
        """Test that extension is preserved when truncating."""
        long_basename = "a" * 5000
        long_path = os.path.join("/tmp", long_basename + ".csv")
        
        result = _validate_path_length(long_path)
        self.assertTrue(result.endswith(".csv"))


class TestPathFunctionBasic(unittest.TestCase):
    """Test basic _path function functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_basic_path_generation(self):
        """Test basic path generation with study and analysis."""
        path = _path(
            study="test_study",
            analysis="clumping",
            directory=self.temp_dir,
            suffix="tsv",
            log=self.log,
            verbose=False
        )
        
        self.assertTrue(path.startswith(self.temp_dir))
        self.assertTrue(path.endswith(".tsv"))
        self.assertIn("test_study", path)
        self.assertIn("clumping", path)
    
    def test_path_with_multiple_components(self):
        """Test path generation with multiple components."""
        path = _path(
            study="study1",
            trait="trait1",
            chrom="1",
            analysis="mtag",
            directory=self.temp_dir,
            suffix="csv",
            log=self.log,
            verbose=False
        )
        
        self.assertIn("study1", path)
        self.assertIn("trait1", path)
        self.assertIn("1", path)
        self.assertIn("mtag", path)
        self.assertTrue(path.endswith(".csv"))
    
    def test_path_with_directory(self):
        """Test path generation with directory."""
        path = _path(
            study="test",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        self.assertTrue(path.startswith(self.temp_dir))
        self.assertTrue(os.path.dirname(path) == self.temp_dir)
    
    def test_path_with_subdirectory(self):
        """Test path generation with subdirectory."""
        path = _path(
            study="test",
            directory=self.temp_dir,
            subdirectory="results",
            log=self.log,
            verbose=False
        )
        
        self.assertIn("results", path)
        self.assertTrue(os.path.exists(os.path.dirname(path)))
    
    def test_path_with_suffix(self):
        """Test path generation with file extension."""
        path = _path(
            study="test",
            suffix="png",
            log=self.log,
            verbose=False
        )
        
        self.assertTrue(path.endswith(".png"))
    
    def test_path_with_user_provided_out(self):
        """Test path generation with user-provided output path."""
        user_path = os.path.join(self.temp_dir, "user_file.txt")
        path = _path(
            out=user_path,
            log=self.log,
            verbose=False
        )
        
        # Path should be sanitized (uses underscores)
        # But directory and extension should be preserved
        self.assertEqual(os.path.dirname(path), os.path.dirname(user_path))
        self.assertTrue(path.endswith(".txt"))
        self.assertIn("user", path.lower())
    
    def test_path_with_user_provided_directory(self):
        """Test path generation when out is a directory."""
        path = _path(
            out=self.temp_dir,
            study="test",
            suffix="tsv",
            log=self.log,
            verbose=False
        )
        
        self.assertTrue(path.startswith(self.temp_dir))
        self.assertTrue(path.endswith(".tsv"))


class TestPathFunctionTemporary(unittest.TestCase):
    """Test _path function with temporary file handling."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_tmp_flag_adds_prefix(self):
        """Test that tmp=True adds _gwaslab prefix."""
        path = _path(
            study="test",
            tmp=True,
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        basename = os.path.basename(path)
        self.assertTrue(basename.startswith("_gwaslab"))
    
    def test_tmp_flag_generates_unique_id(self):
        """Test that tmp=True generates unique ID when pid not provided."""
        path1 = _path(
            study="test",
            tmp=True,
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        path2 = _path(
            study="test",
            tmp=True,
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Paths should be different due to unique IDs
        self.assertNotEqual(path1, path2)
    
    def test_tmp_flag_with_provided_pid(self):
        """Test that provided pid is used instead of generating one."""
        path = _path(
            study="test",
            tmp=True,
            pid="custom123",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        self.assertIn("custom123", path)


class TestPathFunctionSanitization(unittest.TestCase):
    """Test _path function sanitization features."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_path_sanitizes_invalid_characters(self):
        """Test that invalid characters in components are sanitized."""
        path = _path(
            study="test<study>",
            analysis="analysis|test",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Should not contain invalid characters
        self.assertNotIn("<", path)
        self.assertNotIn(">", path)
        self.assertNotIn("|", path)
    
    def test_path_sanitizes_user_provided_path(self):
        """Test that user-provided paths are sanitized."""
        user_path = os.path.join(self.temp_dir, "file<name>.txt")
        path = _path(
            out=user_path,
            log=self.log,
            verbose=False
        )
        
        # Should be sanitized
        self.assertNotIn("<", path)
        self.assertNotIn(">", path)
    
    def test_path_sanitizes_suffix(self):
        """Test that file suffix is sanitized."""
        path = _path(
            study="test",
            suffix="png<invalid>",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Extension should be sanitized
        self.assertNotIn("<", path)
        self.assertNotIn(">", path)


class TestPathFunctionCollisionHandling(unittest.TestCase):
    """Test _path function collision handling."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_path_handles_collision(self):
        """Test that path collisions are handled automatically."""
        # Create first file
        path1 = _path(
            study="test",
            analysis="analysis",
            directory=self.temp_dir,
            suffix="tsv",
            log=self.log,
            verbose=False
        )
        
        # Create the file
        with open(path1, 'w') as f:
            f.write("test")
        
        # Generate same path again - should get different path
        path2 = _path(
            study="test",
            analysis="analysis",
            directory=self.temp_dir,
            suffix="tsv",
            log=self.log,
            verbose=False
        )
        
        # Should be different due to collision handling
        self.assertNotEqual(path1, path2)
        self.assertFalse(os.path.exists(path2))  # New path shouldn't exist yet
    
    def test_tmp_paths_dont_collide(self):
        """Test that tmp paths don't collide due to unique IDs."""
        path1 = _path(
            study="test",
            tmp=True,
            directory=self.temp_dir,
            suffix="tsv",
            log=self.log,
            verbose=False
        )
        
        path2 = _path(
            study="test",
            tmp=True,
            directory=self.temp_dir,
            suffix="tsv",
            log=self.log,
            verbose=False
        )
        
        # Should be different due to unique IDs
        self.assertNotEqual(path1, path2)


class TestPathFunctionEdgeCases(unittest.TestCase):
    """Test edge cases for _path function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_path_with_none_values(self):
        """Test that None values are handled correctly."""
        path = _path(
            study="test",
            trait=None,
            exposure=None,
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Should not contain None or error
        self.assertIsInstance(path, str)
        self.assertNotIn("None", path)
    
    def test_path_with_empty_strings(self):
        """Test that empty strings are handled."""
        path = _path(
            study="",
            analysis="test",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Should still generate valid path
        self.assertIsInstance(path, str)
        self.assertIn("test", path)
    
    def test_path_with_args(self):
        """Test path generation with *args."""
        path = _path(
            "arg1", "arg2",
            study="test",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        self.assertIn("arg1", path)
        self.assertIn("arg2", path)
    
    def test_path_with_long_components(self):
        """Test path with very long component names."""
        long_study = "a" * 1000
        path = _path(
            study=long_study,
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Should be truncated but still valid
        self.assertIsInstance(path, str)
        self.assertLess(len(path), 5000)  # Should be reasonable length
    
    def test_path_with_special_characters(self):
        """Test path with various special characters."""
        path = _path(
            study="test study with spaces",
            analysis="analysis-with-dashes",
            trait="trait_with_underscores",
            directory=self.temp_dir,
            log=self.log,
            verbose=False
        )
        
        # Should be sanitized but still readable
        self.assertIsInstance(path, str)
        # Spaces should be converted
        self.assertNotIn(" ", path)


class TestPathFunctionIntegration(unittest.TestCase):
    """Integration tests for _path function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def test_complete_workflow(self):
        """Test complete workflow with all features."""
        # Generate path with all components
        path = _path(
            study="my_study",
            trait="trait_name",
            chrom="1",
            analysis="clumping",
            method="plink",
            directory=self.temp_dir,
            subdirectory="results",
            suffix="tsv",
            result_type="summary",
            log=self.log,
            verbose=False
        )
        
        # Verify all components are present
        self.assertTrue(path.startswith(self.temp_dir))
        self.assertIn("results", path)
        self.assertTrue(path.endswith(".tsv"))
        self.assertIn("my_study", path.lower())
        self.assertIn("trait_name", path.lower())
        self.assertIn("clumping", path.lower())
        
        # Verify directory was created
        self.assertTrue(os.path.exists(os.path.dirname(path)))
    
    def test_tmp_file_workflow(self):
        """Test temporary file generation workflow."""
        path = _path(
            study="temp_study",
            analysis="test",
            tmp=True,
            directory=self.temp_dir,
            suffix="tmp",
            log=self.log,
            verbose=False
        )
        
        # Should have _gwaslab prefix
        basename = os.path.basename(path)
        self.assertTrue(basename.startswith("_gwaslab"))
        
        # Should have unique ID
        self.assertNotEqual(len(basename), len("_gwaslab_temp_study_test.tmp"))
        
        # Should be valid path
        self.assertTrue(path.startswith(self.temp_dir))
        self.assertTrue(path.endswith(".tmp"))


if __name__ == "__main__":
    unittest.main()
