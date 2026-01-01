"""
Comprehensive test suite for _filter_bed function.

Tests cover:
- Unit tests: Filtering in/out, different chromosome formats, edge cases
- Correctness: Compare results with expected outcomes
- Performance: Speed benchmarks with various data sizes
- Memory: Memory usage measurements
"""

import os
import sys
import unittest
import tempfile
import gzip
import time
import gc
import tracemalloc
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_filter_value import _filter_bed, _filter_region
from gwaslab.info.g_Log import Log


# ============================================================================
# Helper functions for generating test data
# ============================================================================

def create_test_sumstats():
    """Create test sumstats DataFrame with various chromosomes and positions."""
    rows = [
        # Chromosome 1 variants
        {"CHR": 1, "POS": 1000, "P": 5e-8, "EA": "A", "NEA": "G", "SNPID": "1:1000_G_A"},
        {"CHR": 1, "POS": 5000, "P": 1e-6, "EA": "C", "NEA": "T", "SNPID": "1:5000_T_C"},
        {"CHR": 1, "POS": 10000, "P": 0.05, "EA": "G", "NEA": "A", "SNPID": "1:10000_A_G"},
        {"CHR": 1, "POS": 20000, "P": 0.1, "EA": "T", "NEA": "C", "SNPID": "1:20000_C_T"},
        
        # Chromosome 2 variants
        {"CHR": 2, "POS": 1500, "P": 1e-3, "EA": "A", "NEA": "G", "SNPID": "2:1500_G_A"},
        {"CHR": 2, "POS": 8000, "P": 0.01, "EA": "C", "NEA": "T", "SNPID": "2:8000_T_C"},
        {"CHR": 2, "POS": 15000, "P": 0.2, "EA": "G", "NEA": "A", "SNPID": "2:15000_A_G"},
        
        # Chromosome X variants
        {"CHR": 23, "POS": 5000, "P": 1e-5, "EA": "A", "NEA": "G", "SNPID": "X:5000_G_A"},
        {"CHR": 23, "POS": 12000, "P": 0.03, "EA": "C", "NEA": "T", "SNPID": "X:12000_T_C"},
    ]
    return pd.DataFrame(rows)


def create_test_sumstats_string_chr():
    """Create test sumstats with string chromosome format."""
    df = create_test_sumstats()
    df["CHR"] = df["CHR"].replace({23: "X"})
    return df


def create_test_sumstats_chr_prefix():
    """Create test sumstats with chr-prefixed chromosome format."""
    df = create_test_sumstats()
    df["CHR"] = "chr" + df["CHR"].astype(str).replace({"23": "X"})
    return df


def generate_sumstats(n_variants: int, n_chromosomes: int = 5, seed: int = 42) -> pd.DataFrame:
    """Generate simulated sumstats DataFrame for performance testing."""
    np.random.seed(seed)
    
    # For very large datasets, use more efficient generation
    if n_variants > 100000:
        # Generate all data at once using numpy
        chr_nums = np.zeros(n_variants, dtype=int)
        positions_list = np.zeros(n_variants, dtype=int)
        
        total_generated = 0
        for chr_num in range(1, min(n_chromosomes + 1, 23)):
            if total_generated >= n_variants:
                break
            remaining = n_variants - total_generated
            n_chr_variants = min(remaining, max(1, n_variants // min(n_chromosomes, 22)))
            if n_chr_variants == 0:
                continue
            
            end_idx = min(total_generated + n_chr_variants, n_variants)
            actual_count = end_idx - total_generated
            
            positions = np.sort(np.random.randint(1, 50_000_000, size=actual_count))
            chr_nums[total_generated:end_idx] = chr_num
            positions_list[total_generated:end_idx] = positions
            total_generated = end_idx
        
        # Create DataFrame directly
        df = pd.DataFrame({
            "CHR": chr_nums,
            "POS": positions_list,
            "P": np.random.uniform(1e-8, 0.05, size=n_variants),
            "EA": np.random.choice(["A", "C", "G", "T"], size=n_variants),
            "NEA": np.random.choice(["A", "C", "G", "T"], size=n_variants),
            "BETA": np.random.normal(0, 0.1, size=n_variants),
            "SE": np.random.uniform(0.01, 0.1, size=n_variants)
        })
        # Generate SNPID more efficiently
        df["SNPID"] = df["CHR"].astype(str) + ":" + df["POS"].astype(str) + "_A_G"
        return df
    else:
        # Original method for smaller datasets
        variants = []
        for chr_num in range(1, min(n_chromosomes + 1, 6)):
            n_chr_variants = max(1, n_variants // min(n_chromosomes, 5))
            if n_chr_variants == 0:
                continue
            positions = np.sort(np.random.randint(1, 10_000_000, size=n_chr_variants))
            
            for pos in positions:
                variants.append({
                    "CHR": chr_num,
                    "POS": int(pos),
                    "P": np.random.uniform(1e-8, 0.05),
                    "EA": np.random.choice(["A", "C", "G", "T"]),
                    "NEA": np.random.choice(["A", "C", "G", "T"]),
                    "SNPID": f"{chr_num}:{pos}_A_G",
                    "BETA": np.random.normal(0, 0.1),
                    "SE": np.random.uniform(0.01, 0.1)
                })
        
        return pd.DataFrame(variants)


def generate_bed_file(n_regions: int, n_chromosomes: int = 5, 
                      region_size: int = 10000, seed: int = 42) -> str:
    """Generate a BED file with specified number of regions."""
    np.random.seed(seed)
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
    bed_path = temp_file.name
    
    regions_written = 0
    max_chr = min(n_chromosomes + 1, 23)
    max_pos = 50_000_000 if n_regions > 1000 else 5_000_000
    
    for chr_num in range(1, max_chr):
        if regions_written >= n_regions:
            break
        n_chr_regions = max(1, n_regions // min(n_chromosomes, 22))
        for _ in range(n_chr_regions):
            if regions_written >= n_regions:
                break
            start = np.random.randint(0, max_pos)
            end = start + region_size
            temp_file.write(f"{chr_num}\t{start}\t{end}\n")
            regions_written += 1
    
    temp_file.close()
    return bed_path


# ============================================================================
# Unit Tests
# ============================================================================

class TestFilterBedUnit(unittest.TestCase):
    """Unit test suite for _filter_bed function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.log = Log()
        self.sumstats_numeric = create_test_sumstats()
        self.sumstats_string = create_test_sumstats_string_chr()
        self.sumstats_chr = create_test_sumstats_chr_prefix()
    
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def _create_bed_file(self, content, filename="test.bed", compressed=False):
        """Create a BED file with given content."""
        bed_path = os.path.join(self.temp_dir, filename)
        if compressed:
            with gzip.open(bed_path, 'wt') as f:
                f.write(content)
        else:
            with open(bed_path, 'w') as f:
                f.write(content)
        return bed_path
    
    def test_filter_bed_keep_numeric_chr(self):
        """Test filtering in with numeric chromosomes."""
        # BED: chr1: 0-6000 (0-based) = positions 1-6000 (1-based)
        bed_content = "1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep variants at positions 1000 and 5000 on chr1
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["CHR"] == 1))
        self.assertTrue(all(result["POS"].isin([1000, 5000])))
    
    def test_filter_bed_remove_numeric_chr(self):
        """Test filtering out with numeric chromosomes."""
        bed_content = "1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=False,
            verbose=False,
            log=self.log
        )
        
        # Should remove variants at positions 1000 and 5000 on chr1
        self.assertEqual(len(result), 7)  # 9 total - 2 removed
        self.assertFalse(any((result["CHR"] == 1) & (result["POS"].isin([1000, 5000]))))
    
    def test_filter_bed_string_chr(self):
        """Test filtering with string chromosome format in sumstats."""
        # BED with numeric chr, sumstats with string chr
        bed_content = "1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_string,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should work correctly with format conversion
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["POS"].isin([1000, 5000])))
    
    def test_filter_bed_chr_prefix_bed(self):
        """Test filtering with chr-prefixed BED file."""
        bed_content = "chr1\t0\t6000\nchr2\t1000\t9000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep variants in both regions
        self.assertEqual(len(result), 4)  # 2 from chr1, 2 from chr2
        self.assertTrue(any(result["CHR"] == 1))
        self.assertTrue(any(result["CHR"] == 2))
    
    def test_filter_bed_chr_prefix_sumstats(self):
        """Test filtering with chr-prefixed sumstats."""
        bed_content = "1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_chr,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should work correctly with format conversion
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["CHR"] == "chr1"))
    
    def test_filter_bed_multiple_regions(self):
        """Test filtering with multiple BED regions."""
        bed_content = "1\t0\t6000\n2\t7000\t10000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep variants from both regions
        self.assertEqual(len(result), 3)  # 2 from chr1, 1 from chr2
        chr1_variants = result[result["CHR"] == 1]
        chr2_variants = result[result["CHR"] == 2]
        self.assertEqual(len(chr1_variants), 2)
        self.assertEqual(len(chr2_variants), 1)
        self.assertEqual(chr2_variants.iloc[0]["POS"], 8000)
    
    def test_filter_bed_chr_x(self):
        """Test filtering with X chromosome."""
        bed_content = "X\t0\t10000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep X chromosome variant at position 5000
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["CHR"], 23)
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_bed_boundary_conditions(self):
        """Test filtering with boundary positions."""
        # BED: chr1: 999-5001 (0-based) = positions 1000-5001 (1-based)
        bed_content = "1\t999\t5001\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should include positions 1000 and 5000 (boundaries)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["POS"].isin([1000, 5000])))
    
    def test_filter_bed_exact_boundary(self):
        """Test filtering with exact boundary match."""
        # BED: chr1: 999-1001 (0-based) = positions 1000-1001 (1-based)
        bed_content = "1\t999\t1001\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should include position 1000 (exact start boundary)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 1000)
    
    def test_filter_bed_compressed(self):
        """Test filtering with compressed BED file."""
        bed_content = "1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content, filename="test.bed.gz", compressed=True)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should work the same as uncompressed
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["CHR"] == 1))
    
    def test_filter_bed_empty_bed(self):
        """Test filtering with BED file that has no matching regions."""
        bed_content = "1\t100000\t200000\n"  # Region far from any variants
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should return empty DataFrame when no matches
        self.assertEqual(len(result), 0)
        self.assertTrue(isinstance(result, pd.DataFrame))
    
    def test_filter_bed_no_matches(self):
        """Test filtering when no variants match BED regions."""
        bed_content = "1\t50000\t60000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should return empty DataFrame
        self.assertEqual(len(result), 0)
    
    def test_filter_bed_all_variants_match(self):
        """Test filtering when all variants match."""
        bed_content = "1\t0\t30000\n2\t0\t20000\n23\t0\t20000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep all variants
        self.assertEqual(len(result), len(self.sumstats_numeric))
    
    def test_filter_bed_with_header(self):
        """Test filtering with BED file containing header lines."""
        bed_content = "# browser position chr1:0-10000\n# track name=test\n1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should skip header lines and process correctly
        self.assertEqual(len(result), 2)
    
    def test_filter_bed_multiple_chromosomes(self):
        """Test filtering across multiple chromosomes."""
        bed_content = "1\t0\t6000\n2\t7000\t10000\n23\t0\t10000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep variants from all three chromosomes
        self.assertEqual(len(result), 4)  # 2 from chr1, 1 from chr2, 1 from chr23
        self.assertTrue(any(result["CHR"] == 1))
        self.assertTrue(any(result["CHR"] == 2))
        self.assertTrue(any(result["CHR"] == 23))
    
    def test_filter_bed_overlapping_regions(self):
        """Test filtering with overlapping BED regions."""
        # Two overlapping regions on chr1
        bed_content = "1\t0\t6000\n1\t4000\t8000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should keep variants in either region (union)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["POS"].isin([1000, 5000])))
    
    def test_filter_bed_custom_chrom_pos_columns(self):
        """Test filtering with custom chromosome and position column names."""
        bed_content = "1\t0\t6000\n"
        bed_path = self._create_bed_file(bed_content)
        
        result = _filter_bed(
            self.sumstats_numeric,
            path=bed_path,
            chrom="CHR",
            pos="POS",
            keep=True,
            verbose=False,
            log=self.log
        )
        
        # Should work correctly
        self.assertEqual(len(result), 2)
        self.assertTrue(all(result["CHR"] == 1))
    
    def test_filter_bed_error_no_path(self):
        """Test that error is raised when path is None."""
        with self.assertRaises(ValueError):
            _filter_bed(
                self.sumstats_numeric,
                path=None,
                verbose=False,
                log=self.log
            )
    
    def test_filter_bed_error_file_not_found(self):
        """Test that error is raised when BED file doesn't exist."""
        with self.assertRaises(FileNotFoundError):
            _filter_bed(
                self.sumstats_numeric,
                path="/nonexistent/path.bed",
                verbose=False,
                log=self.log
            )


# ============================================================================
# Correctness Tests
# ============================================================================

class TestFilterBedCorrectness(unittest.TestCase):
    """Correctness tests comparing results with expected outcomes."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
    
    def test_correctness_simple_case(self):
        """Test simple case: single region."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        bed_content = "1\t0\t6000\n"
        
        bed_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        bed_path.write(bed_content)
        bed_path.close()
        bed_path_name = bed_path.name
        
        try:
            result_keep = _filter_bed(sumstats.copy(), path=bed_path_name, keep=True, verbose=False, log=self.log)
            result_remove = _filter_bed(sumstats.copy(), path=bed_path_name, keep=False, verbose=False, log=self.log)
            
            # Verify keep results
            kept_positions = set(result_keep["POS"].tolist())
            self.assertEqual(kept_positions, {1000, 5000})
            
            # Verify remove results
            removed_positions = set(result_remove["POS"].tolist())
            self.assertEqual(removed_positions, {10000, 2000})
        finally:
            os.unlink(bed_path_name)
    
    def test_correctness_multiple_regions(self):
        """Test multiple regions, multiple chromosomes."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 2, 2, 3],
            "POS": [1000, 5000, 1500, 8000, 3000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8, 1e-9]
        })
        bed_content = "1\t0\t6000\n2\t1000\t9000\n"
        
        bed_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        bed_path.write(bed_content)
        bed_path.close()
        bed_path_name = bed_path.name
        
        try:
            result_keep = _filter_bed(sumstats.copy(), path=bed_path_name, keep=True, verbose=False, log=self.log)
            result_remove = _filter_bed(sumstats.copy(), path=bed_path_name, keep=False, verbose=False, log=self.log)
            
            kept_positions = set(result_keep["POS"].tolist())
            self.assertEqual(kept_positions, {1000, 1500, 5000, 8000})
            
            removed_positions = set(result_remove["POS"].tolist())
            self.assertEqual(removed_positions, {3000})
        finally:
            os.unlink(bed_path_name)
    
    def test_correctness_boundary_conditions(self):
        """Test boundary conditions."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1],
            "POS": [1000, 2000, 3000],
            "P": [1e-6, 1e-7, 1e-5]
        })
        bed_content = "1\t999\t3001\n"  # BED: 999-3001 (0-based) = 1000-3001 (1-based)
        
        bed_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        bed_path.write(bed_content)
        bed_path.close()
        bed_path_name = bed_path.name
        
        try:
            result_keep = _filter_bed(sumstats.copy(), path=bed_path_name, keep=True, verbose=False, log=self.log)
            kept_positions = set(result_keep["POS"].tolist())
            self.assertEqual(kept_positions, {1000, 2000, 3000})
        finally:
            os.unlink(bed_path_name)
    
    def test_correctness_no_matches(self):
        """Test no matches."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1],
            "POS": [1000, 5000],
            "P": [1e-6, 1e-7]
        })
        bed_content = "1\t100000\t200000\n"  # Far from variants
        
        bed_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        bed_path.write(bed_content)
        bed_path.close()
        bed_path_name = bed_path.name
        
        try:
            result_keep = _filter_bed(sumstats.copy(), path=bed_path_name, keep=True, verbose=False, log=self.log)
            result_remove = _filter_bed(sumstats.copy(), path=bed_path_name, keep=False, verbose=False, log=self.log)
            
            self.assertEqual(len(result_keep), 0)
            self.assertEqual(len(result_remove), 2)
        finally:
            os.unlink(bed_path_name)
    
    def test_correctness_all_matches(self):
        """Test all matches."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1],
            "POS": [1000, 5000],
            "P": [1e-6, 1e-7]
        })
        bed_content = "1\t0\t100000\n"  # Covers all
        
        bed_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        bed_path.write(bed_content)
        bed_path.close()
        bed_path_name = bed_path.name
        
        try:
            result_keep = _filter_bed(sumstats.copy(), path=bed_path_name, keep=True, verbose=False, log=self.log)
            result_remove = _filter_bed(sumstats.copy(), path=bed_path_name, keep=False, verbose=False, log=self.log)
            
            self.assertEqual(len(result_keep), 2)
            self.assertEqual(len(result_remove), 0)
        finally:
            os.unlink(bed_path_name)


# ============================================================================
# Performance Tests
# ============================================================================

class TestFilterBedPerformance(unittest.TestCase):
    """Performance tests for speed and memory usage."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
    
    def test_performance_small_dataset(self):
        """Test performance with small dataset."""
        n_variants = 1000
        n_regions = 50
        
        sumstats = generate_sumstats(n_variants)
        bed_path = generate_bed_file(n_regions)
        
        try:
            # Warm-up
            _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            gc.collect()
            
            # Benchmark
            start = time.time()
            result = _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            elapsed = time.time() - start
            
            # Verify reasonable performance
            self.assertLess(elapsed, 5.0)  # Should complete in < 5 seconds
            self.assertGreater(len(result), 0)  # Should match some variants
        finally:
            os.unlink(bed_path)
    
    def test_performance_medium_dataset(self):
        """Test performance with medium dataset."""
        n_variants = 100000
        n_regions = 1000
        
        sumstats = generate_sumstats(n_variants)
        bed_path = generate_bed_file(n_regions)
        
        try:
            # Warm-up
            _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            gc.collect()
            
            # Benchmark
            start = time.time()
            result = _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            elapsed = time.time() - start
            
            # Verify reasonable performance
            self.assertLess(elapsed, 10.0)  # Should complete in < 10 seconds
            self.assertGreater(len(result), 0)
        finally:
            os.unlink(bed_path)
    
    def test_performance_large_dataset(self):
        """Test performance with large dataset (1M variants)."""
        n_variants = 1000000
        n_regions = 5000
        
        sumstats = generate_sumstats(n_variants)
        bed_path = generate_bed_file(n_regions)
        
        try:
            # Warm-up
            _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            gc.collect()
            
            # Benchmark
            start = time.time()
            result = _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            elapsed = time.time() - start
            
            # Verify reasonable performance (should be < 5 seconds for 1M variants)
            self.assertLess(elapsed, 5.0)
            self.assertGreater(len(result), 0)
            
            # Calculate throughput
            throughput = n_variants / elapsed
            self.assertGreater(throughput, 100000)  # Should process > 100K variants/sec
        finally:
            os.unlink(bed_path)
    
    def test_memory_usage_small(self):
        """Test memory usage with small dataset."""
        n_variants = 10000
        n_regions = 100
        
        sumstats = generate_sumstats(n_variants)
        bed_path = generate_bed_file(n_regions)
        
        try:
            tracemalloc.start()
            
            result = _filter_bed(sumstats.copy(), path=bed_path, keep=True, verbose=False, log=self.log)
            
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            
            # Memory should be reasonable (less than 100MB for 10K variants)
            peak_mb = peak / 1024 / 1024
            self.assertLess(peak_mb, 100)
            
            self.assertGreater(len(result), 0)
        finally:
            os.unlink(bed_path)


# ============================================================================
# Tests for _filter_region with special region options
# ============================================================================

class TestFilterRegionSpecial(unittest.TestCase):
    """Test suite for _filter_region with special region options (HLA, HIGH_LD, PAR)."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
    
    def test_filter_region_hla(self):
        """Test filtering with region='HLA' option."""
        # Create sumstats with variants in and outside HLA region (chr6: 29,500,000-33,500,000 for hg19)
        sumstats = pd.DataFrame({
            "CHR": [1, 6, 6, 6, 2],
            "POS": [1000, 30000000, 31000000, 32000000, 2000],  # 3 variants in HLA region
            "P": [1e-6, 1e-7, 1e-5, 1e-8, 1e-9]
        })
        
        result = _filter_region(sumstats.copy(), region="HLA", build="19", verbose=False, log=self.log)
        
        # Should exclude variants in HLA region (chr6: 29,500,000-33,500,000)
        # Variants at 30M, 31M, 32M should be removed
        self.assertEqual(len(result), 2)  # Only chr1 and chr2 variants remain
        self.assertFalse(any(result["CHR"] == 6))  # No chr6 variants should remain
    
    def test_filter_region_hla_hg38(self):
        """Test filtering with region='HLA' option for hg38."""
        # HLA region for hg38 is similar: chr6: 29,602,238-33,409,896
        sumstats = pd.DataFrame({
            "CHR": [1, 6, 6, 2],
            "POS": [1000, 30000000, 32000000, 2000],  # 2 variants in HLA region
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region="HLA", build="38", verbose=False, log=self.log)
        
        # Should exclude variants in HLA region
        self.assertEqual(len(result), 2)  # Only chr1 and chr2 variants remain
        self.assertFalse(any(result["CHR"] == 6))
    
    def test_filter_region_high_ld(self):
        """Test filtering with region='HIGH_LD' option."""
        # Create sumstats with variants in and outside high LD regions
        # High LD regions typically include HLA and other regions
        sumstats = pd.DataFrame({
            "CHR": [1, 6, 6, 23, 2],
            "POS": [1000, 30000000, 31000000, 5000, 2000],  # 2 variants in high LD (HLA)
            "P": [1e-6, 1e-7, 1e-5, 1e-8, 1e-9]
        })
        
        result = _filter_region(sumstats.copy(), region="HIGH_LD", build="19", verbose=False, log=self.log)
        
        # Should exclude variants in high LD regions
        # The exact number depends on what's in the high_ld_hla_hg19.bed.gz file
        # But at minimum, HLA variants should be removed
        self.assertLessEqual(len(result), len(sumstats))
        # Chr6 variants in HLA should be removed
        chr6_variants = result[result["CHR"] == 6]
        self.assertEqual(len(chr6_variants), 0)
    
    def test_filter_region_high_ld_hg38(self):
        """Test filtering with region='HIGH_LD' option for hg38."""
        sumstats = pd.DataFrame({
            "CHR": [1, 6, 6, 2],
            "POS": [1000, 30000000, 32000000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region="HIGH_LD", build="38", verbose=False, log=self.log)
        
        # Should exclude variants in high LD regions
        self.assertLessEqual(len(result), len(sumstats))
        chr6_variants = result[result["CHR"] == 6]
        self.assertEqual(len(chr6_variants), 0)
    
    def test_filter_region_par(self):
        """Test filtering with region='PAR' option."""
        # PAR regions for hg19: chr23: 60,001-2,699,520 (PAR1) and 154,931,044-155,260,560 (PAR2)
        sumstats = pd.DataFrame({
            "CHR": [1, 23, 23, 23, 2],
            "POS": [1000, 100000, 1000000, 155000000, 2000],  # 3 variants in PAR regions
            "P": [1e-6, 1e-7, 1e-5, 1e-8, 1e-9]
        })
        
        result = _filter_region(sumstats.copy(), region="PAR", build="19", verbose=False, log=self.log)
        
        # Should exclude variants in PAR regions
        # Variants at 100000 (in PAR1) and 155000000 (in PAR2) should be removed
        self.assertLessEqual(len(result), len(sumstats))
        chr23_variants = result[result["CHR"] == 23]
        # The variant at 1000000 might be in PAR1, so we check that at least some were removed
        self.assertLess(len(chr23_variants), 3)
    
    def test_filter_region_par_hg38(self):
        """Test filtering with region='PAR' option for hg38."""
        # PAR regions for hg38 are similar but coordinates differ slightly
        sumstats = pd.DataFrame({
            "CHR": [1, 23, 23, 2],
            "POS": [1000, 100000, 155000000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region="PAR", build="38", verbose=False, log=self.log)
        
        # Should exclude variants in PAR regions
        self.assertLessEqual(len(result), len(sumstats))
        chr23_variants = result[result["CHR"] == 23]
        # At least some PAR variants should be removed
        self.assertLess(len(chr23_variants), 2)
    
    def test_filter_region_case_insensitive(self):
        """Test that region string options are case-insensitive."""
        sumstats = pd.DataFrame({
            "CHR": [1, 6, 6, 2],
            "POS": [1000, 30000000, 31000000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        # Test lowercase
        result_lower = _filter_region(sumstats.copy(), region="hla", build="19", verbose=False, log=self.log)
        # Test uppercase
        result_upper = _filter_region(sumstats.copy(), region="HLA", build="19", verbose=False, log=self.log)
        # Test mixed case
        result_mixed = _filter_region(sumstats.copy(), region="HlA", build="19", verbose=False, log=self.log)
        
        # All should produce the same result
        self.assertEqual(len(result_lower), len(result_upper))
        self.assertEqual(len(result_upper), len(result_mixed))
        self.assertTrue(result_lower.equals(result_upper))
    
    def test_filter_region_invalid_string(self):
        """Test that invalid region string raises ValueError."""
        sumstats = pd.DataFrame({
            "CHR": [1, 2],
            "POS": [1000, 2000],
            "P": [1e-6, 1e-7]
        })
        
        with self.assertRaises(ValueError):
            _filter_region(sumstats.copy(), region="INVALID", build="19", verbose=False, log=self.log)
    
    def test_filter_region_regular_list(self):
        """Test that regular list format [chr, start, end] still works."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region=[1, 2000, 8000], verbose=False, log=self.log)
        
        # Should keep variants in region [2000, 8000] on chr1
        self.assertEqual(len(result), 1)  # Only position 5000 should be kept
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_region_string_format_chr_start_end(self):
        """Test region string format 'chr:start-end'."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region="chr1:2000-8000", verbose=False, log=self.log)
        
        # Should keep variants in region [2000, 8000] on chr1
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_region_string_format_chr_pos_flanking(self):
        """Test region string format 'chr:pos:flanking'."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region="1:5000:3000", verbose=False, log=self.log)
        
        # Should keep variants in region [2000, 8000] on chr1 (5000 ± 3000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_region_string_format_flanking_kb(self):
        """Test region string format 'chr:pos:flanking' with kb suffix."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        result = _filter_region(sumstats.copy(), region="1:5000:3kb", verbose=False, log=self.log)
        
        # Should keep variants in region [2000, 8000] on chr1 (5000 ± 3000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_region_string_format_snpid_flanking(self):
        """Test region string format 'chr:pos:ea:nea:flanking'."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8],
            "EA": ["A", "C", "G", "T"],
            "NEA": ["G", "T", "A", "C"]
        })
        
        result = _filter_region(sumstats.copy(), region="1:5000:C:T:3000", verbose=False, log=self.log)
        
        # Should keep variants in region [2000, 8000] on chr1 (5000 ± 3000)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_region_normalized_swapped_start_end(self):
        """Test that swapped start/end in list format is normalized."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        # Start > end should be normalized
        result = _filter_region(sumstats.copy(), region=[1, 8000, 2000], verbose=False, log=self.log)
        
        # Should keep variants in region [2000, 8000] on chr1 (normalized)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 5000)
    
    def test_filter_region_normalized_chromosome_format(self):
        """Test that different chromosome formats are normalized."""
        sumstats = pd.DataFrame({
            "CHR": [1, 1, 1, 2],
            "POS": [1000, 5000, 10000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        # Test with "chr1" prefix
        result1 = _filter_region(sumstats.copy(), region="chr1:2000-8000", verbose=False, log=self.log)
        # Test with "1" (no prefix)
        result2 = _filter_region(sumstats.copy(), region="1:2000-8000", verbose=False, log=self.log)
        
        # Both should produce the same result
        self.assertEqual(len(result1), len(result2))
        self.assertTrue(result1.equals(result2))
    
    def test_filter_region_none(self):
        """Test that region=None returns original sumstats."""
        sumstats = pd.DataFrame({
            "CHR": [1, 2],
            "POS": [1000, 2000],
            "P": [1e-6, 1e-7]
        })
        
        result = _filter_region(sumstats.copy(), region=None, verbose=False, log=self.log)
        
        # Should return original sumstats unchanged
        self.assertEqual(len(result), len(sumstats))
        self.assertTrue(result.equals(sumstats))
    
    def test_filter_region_build_parameter(self):
        """Test that build parameter is properly passed to filtering functions."""
        sumstats = pd.DataFrame({
            "CHR": [1, 6, 6, 2],
            "POS": [1000, 30000000, 31000000, 2000],
            "P": [1e-6, 1e-7, 1e-5, 1e-8]
        })
        
        # Test with build="19"
        result_19 = _filter_region(sumstats.copy(), region="HLA", build="19", verbose=False, log=self.log)
        # Test with build="38"
        result_38 = _filter_region(sumstats.copy(), region="HLA", build="38", verbose=False, log=self.log)
        
        # Both should work without errors
        self.assertIsInstance(result_19, pd.DataFrame)
        self.assertIsInstance(result_38, pd.DataFrame)
        # Both should exclude chr6 variants
        self.assertEqual(len(result_19[result_19["CHR"] == 6]), 0)
        self.assertEqual(len(result_38[result_38["CHR"] == 6]), 0)


if __name__ == "__main__":
    # Run all tests
    unittest.main(verbosity=2)
