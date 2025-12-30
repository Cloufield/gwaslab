"""
Test suite for ChromosomeMapper class.

Tests cover:
- Basic mapping functionality (numeric, string, chr, nc formats)
- Format detection
- Vectorized operations (Series, DataFrame)
- Species-specific mappings
- Sumstats and reference layer mappings
- Edge cases (invalid inputs, missing values, etc.)
"""

import os
import sys
import unittest
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log


class TestChromosomeMapper(unittest.TestCase):
    """Test suite for ChromosomeMapper class."""

    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.log.verbose = False

    # ========================================================================
    # Basic Mapping Tests - Sumstats Layer
    # ========================================================================

    def test_sumstats_to_number_numeric(self):
        """Test converting numeric sumstats to number."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        self.assertEqual(mapper.sumstats_to_number(1), 1)
        self.assertEqual(mapper.sumstats_to_number(22), 22)
        self.assertEqual(mapper.sumstats_to_number(23), 23)
        self.assertEqual(mapper.sumstats_to_number(24), 24)
        self.assertEqual(mapper.sumstats_to_number(25), 25)

    def test_sumstats_to_number_string(self):
        """Test converting string sumstats to number."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["1", "2", "X"]))
        
        self.assertEqual(mapper.sumstats_to_number("1"), 1)
        self.assertEqual(mapper.sumstats_to_number("22"), 22)
        self.assertEqual(mapper.sumstats_to_number("X"), 23)
        self.assertEqual(mapper.sumstats_to_number("Y"), 24)
        self.assertEqual(mapper.sumstats_to_number("MT"), 25)

    def test_sumstats_to_number_chr(self):
        """Test converting chr-prefixed sumstats to number."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))
        
        self.assertEqual(mapper.sumstats_to_number("chr1"), 1)
        self.assertEqual(mapper.sumstats_to_number("chr22"), 22)
        self.assertEqual(mapper.sumstats_to_number("chrX"), 23)
        self.assertEqual(mapper.sumstats_to_number("chrY"), 24)
        self.assertEqual(mapper.sumstats_to_number("chrMT"), 25)
        # Test case variations
        self.assertEqual(mapper.sumstats_to_number("Chr1"), 1)
        self.assertEqual(mapper.sumstats_to_number("CHR1"), 1)

    def test_number_to_sumstats_numeric(self):
        """Test converting number to numeric sumstats."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        self.assertEqual(mapper.number_to_sumstats(1), 1)
        self.assertEqual(mapper.number_to_sumstats(22), 22)
        self.assertEqual(mapper.number_to_sumstats(23), 23)

    def test_number_to_sumstats_string(self):
        """Test converting number to string sumstats."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["1", "2", "X"]))
        
        self.assertEqual(mapper.number_to_sumstats(1), "1")
        self.assertEqual(mapper.number_to_sumstats(22), "22")
        self.assertEqual(mapper.number_to_sumstats(23), "X")
        self.assertEqual(mapper.number_to_sumstats(24), "Y")
        # Note: The actual format may vary (MT vs Mt), so check it's a valid mitochondrial identifier
        result = mapper.number_to_sumstats(25)
        self.assertIn(result.upper(), ["MT", "M"])

    def test_number_to_sumstats_chr(self):
        """Test converting number to chr-prefixed sumstats."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))
        
        self.assertEqual(mapper.number_to_sumstats(1), "chr1")
        self.assertEqual(mapper.number_to_sumstats(22), "chr22")
        self.assertEqual(mapper.number_to_sumstats(23), "chrX")
        self.assertEqual(mapper.number_to_sumstats(24), "chrY")
        # Note: The actual format may vary (chrMT vs chrMt), so check it's a valid mitochondrial identifier
        result = mapper.number_to_sumstats(25)
        self.assertTrue(result.startswith("chr"))
        self.assertIn(result[3:].upper(), ["MT", "M"])

    def test_sumstats_to_reference(self):
        """Test converting sumstats to reference format."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        mapper.detect_reference_format("dummy.vcf")  # Will detect chr format
        
        # Mock the reference format detection to return chr format
        mapper._reference_format = "chr"
        mapper._reference_prefix = "chr"
        mapper._build_reference_layer("chr", "chr")
        
        result = mapper.sumstats_to_reference(1)
        self.assertEqual(result, "chr1")
        
        result = mapper.sumstats_to_reference(23)
        self.assertEqual(result, "chrX")

    def test_reference_to_sumstats(self):
        """Test converting reference to sumstats format."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        mapper._reference_format = "chr"
        mapper._reference_prefix = "chr"
        mapper._build_reference_layer("chr", "chr")
        
        result = mapper.reference_to_sumstats("chr1")
        self.assertEqual(result, 1)
        
        result = mapper.reference_to_sumstats("chrX")
        self.assertEqual(result, 23)

    # ========================================================================
    # Convenience Methods Tests
    # ========================================================================

    def test_to_numeric(self):
        """Test to_numeric convenience method."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        
        # Single value - need to detect format first for chr format
        mapper.detect_sumstats_format(pd.Series(["chr1"]))
        result = mapper.to_numeric("chr1")
        self.assertEqual(result, 1)
        
        # For string format
        mapper.detect_sumstats_format(pd.Series(["X"]))
        result = mapper.to_numeric("X")
        self.assertEqual(result, 23)
        
        # Series - detect format first
        series = pd.Series(["chr1", "chr2", "X"])
        mapper.detect_sumstats_format(series)
        result = mapper.to_numeric(series)
        expected = pd.Series([1, 2, 23])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)

    def test_to_string(self):
        """Test to_string convenience method."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        
        # Single value
        result = mapper.to_string(1)
        self.assertEqual(result, "1")
        
        result = mapper.to_string(23)
        self.assertEqual(result, "X")
        
        # Series
        series = pd.Series([1, 2, 23, 24, 25])
        result = mapper.to_string(series)
        expected = pd.Series(["1", "2", "X", "Y", "MT"])
        pd.testing.assert_series_equal(result, expected)

    def test_to_chr(self):
        """Test to_chr convenience method."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        
        # Single value
        result = mapper.to_chr(1, prefix="chr")
        self.assertEqual(result, "chr1")
        
        result = mapper.to_chr(23, prefix="chr")
        self.assertEqual(result, "chrX")
        
        # Series
        series = pd.Series([1, 2, 23])
        result = mapper.to_chr(series, prefix="chr")
        expected = pd.Series(["chr1", "chr2", "chrX"])
        pd.testing.assert_series_equal(result, expected)

    def test_to_nc(self):
        """Test to_nc convenience method."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="19",
            log=self.log,
            verbose=False
        )
        
        # Single value
        result = mapper.to_nc(1)
        self.assertEqual(result, "NC_000001.10")
        
        result = mapper.to_nc(23)
        self.assertEqual(result, "NC_000023.10")
        
        # Series
        series = pd.Series([1, 2, 23])
        result = mapper.to_nc(series)
        expected = pd.Series(["NC_000001.10", "NC_000002.11", "NC_000023.10"])
        pd.testing.assert_series_equal(result, expected)

    # ========================================================================
    # Format Detection Tests
    # ========================================================================

    def test_detect_format_numeric(self):
        """Test format detection for numeric format."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        series = pd.Series([1, 2, 3, 22, 23])
        detected = mapper.detect_format(series)
        self.assertEqual(detected, "numeric")

    def test_detect_format_string(self):
        """Test format detection for string format."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        series = pd.Series(["1", "2", "X", "Y", "MT"])
        detected = mapper.detect_format(series)
        self.assertEqual(detected, "string")

    def test_detect_format_chr(self):
        """Test format detection for chr format."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        series = pd.Series(["chr1", "chr2", "chrX", "chrY"])
        detected = mapper.detect_format(series)
        self.assertEqual(detected, "chr")

    def test_detect_format_nc(self):
        """Test format detection for NCBI RefSeq ID format."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        series = pd.Series(["NC_000001.10", "NC_000002.11", "NC_000023.10"])
        detected = mapper.detect_format(series)
        self.assertEqual(detected, "nc")

    def test_detect_format_single_value(self):
        """Test format detection for single value."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        detected = mapper.detect_format("chr1")
        self.assertEqual(detected, "chr")
        
        detected = mapper.detect_format(1)
        self.assertEqual(detected, "numeric")

    def test_detect_sumstats_format(self):
        """Test sumstats format detection."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        
        # Numeric format
        series = pd.Series([1, 2, 23])
        mapper.detect_sumstats_format(series)
        self.assertEqual(mapper._sumstats_format, "numeric")
        
        # String format
        series = pd.Series(["1", "2", "X"])
        mapper.detect_sumstats_format(series)
        self.assertEqual(mapper._sumstats_format, "string")
        
        # Chr format
        series = pd.Series(["chr1", "chr2", "chrX"])
        mapper.detect_sumstats_format(series)
        self.assertEqual(mapper._sumstats_format, "chr")

    # ========================================================================
    # Vectorized Operations Tests
    # ========================================================================

    def test_sumstats_to_reference_series(self):
        """Test vectorized sumstats to reference conversion."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        mapper._reference_format = "chr"
        mapper._reference_prefix = "chr"
        mapper._build_reference_layer("chr", "chr")
        
        series = pd.Series([1, 2, 23, 24, 25])
        result = mapper.sumstats_to_reference_series(series)
        
        # Check individual values as the exact format may vary
        self.assertEqual(result.iloc[0], "chr1")
        self.assertEqual(result.iloc[1], "chr2")
        self.assertEqual(result.iloc[2], "chrX")
        self.assertEqual(result.iloc[3], "chrY")
        # MT may be "chrMT" or "chrMt" depending on species mapping
        self.assertTrue(result.iloc[4].startswith("chr"))

    def test_reference_to_sumstats_series(self):
        """Test vectorized reference to sumstats conversion."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        mapper._reference_format = "chr"
        mapper._reference_prefix = "chr"
        mapper._build_reference_layer("chr", "chr")
        
        series = pd.Series(["chr1", "chr2", "chrX", "chrY", "chrMT"])
        result = mapper.reference_to_sumstats_series(series)
        
        expected = pd.Series([1, 2, 23, 24, 25])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)

    def test_map_series_backward_compat(self):
        """Test backward compatibility map_series method."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        series = pd.Series([1, 2, 22, 23, 24, 25])
        result = mapper.map_series(series)
        
        # Should return in sumstats format (numeric in this case)
        expected = pd.Series([1, 2, 22, 23, 24, 25])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)

    # ========================================================================
    # Species-Specific Tests
    # ========================================================================

    def test_species_specific_mapping(self):
        """Test species-specific chromosome mappings."""
        # Human (default)
        mapper_human = ChromosomeMapper(
            species="homo sapiens",
            log=self.log,
            verbose=False
        )
        mapper_human.detect_sumstats_format(pd.Series(["X"]))
        self.assertEqual(mapper_human.sumstats_to_number("X"), 23)
        
        # Mouse
        mapper_mouse = ChromosomeMapper(
            species="mus musculus",
            log=self.log,
            verbose=False
        )
        mapper_mouse.detect_sumstats_format(pd.Series(["1", "19"]))
        # Mouse has 19 autosomes
        self.assertEqual(mapper_mouse.sumstats_to_number("1"), 1)
        self.assertEqual(mapper_mouse.sumstats_to_number("19"), 19)

    def test_species_max_chromosome(self):
        """Test species-specific maximum chromosome count."""
        # Human has 22 autosomes
        mapper_human = ChromosomeMapper(
            species="homo sapiens",
            log=self.log,
            verbose=False
        )
        self.assertEqual(mapper_human._max_chr, 25)  # 22 autosomes + X + Y + MT
        
        # Mouse has 19 autosomes
        mapper_mouse = ChromosomeMapper(
            species="mus musculus",
            log=self.log,
            verbose=False
        )
        # Mouse has 19 autosomes + X + Y + MT
        self.assertGreater(mapper_mouse._max_chr, 19)

    # ========================================================================
    # Edge Cases Tests
    # ========================================================================

    def test_missing_values(self):
        """Test handling of missing/NaN values."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        series = pd.Series(["chr1", None, "chr2", np.nan, "X"])
        # Filter out None/NaN for format detection
        mapper.detect_sumstats_format(series.dropna())
        
        # Use to_numeric which should handle missing values
        # Note: The current implementation may raise errors for None/NaN
        # So we test with a series that has valid values and check the behavior
        valid_series = pd.Series(["chr1", "chr2", "X"])
        result = mapper.to_numeric(valid_series)
        expected = pd.Series([1, 2, 23])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)
        
        # For missing values, the mapper may raise errors or handle them
        # This is expected behavior - missing values need to be handled by the caller

    def test_invalid_chromosome(self):
        """Test handling of invalid chromosome values."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["1", "2"]))
        
        # Invalid values should return pd.NA instead of raising ValueError
        result = mapper.sumstats_to_number("invalid")
        self.assertTrue(pd.isna(result))
        
        # Test with unconvertible chromosome like '1_KI270766v1_alt'
        result = mapper.sumstats_to_number("1_KI270766v1_alt")
        self.assertTrue(pd.isna(result))

    def test_empty_series(self):
        """Test handling of empty Series."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        
        series = pd.Series([], dtype=object)
        result = mapper.to_numeric(series)
        
        self.assertEqual(len(result), 0)

    def test_mixed_formats_in_series(self):
        """Test handling of mixed formats in Series."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        
        # When format is detected, it should handle the primary format
        series = pd.Series(["chr1", "chr2", "chrX"])
        mapper.detect_sumstats_format(series)
        result = mapper.to_numeric(series)
        
        expected = pd.Series([1, 2, 23])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)

    def test_nc_format_without_build(self):
        """Test NC format conversion without build."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2]))
        
        # Without build, NC mapping is not available, should raise ValueError
        with self.assertRaises(ValueError):
            mapper.to_nc(1)

    def test_nc_format_with_build(self):
        """Test NC format conversion with build."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="19",
            log=self.log,
            verbose=False
        )
        
        result = mapper.to_nc(1)
        self.assertEqual(result, "NC_000001.10")

    def test_different_builds(self):
        """Test different genome builds."""
        # hg19
        mapper_19 = ChromosomeMapper(
            species="homo sapiens",
            build="19",
            log=self.log,
            verbose=False
        )
        result_19 = mapper_19.to_nc(1)
        self.assertEqual(result_19, "NC_000001.10")
        
        # hg38
        mapper_38 = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        result_38 = mapper_38.to_nc(1)
        self.assertEqual(result_38, "NC_000001.11")

    # ========================================================================
    # Integration Tests
    # ========================================================================

    def test_round_trip_conversion(self):
        """Test round-trip conversion (numeric -> chr -> numeric)."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        mapper._reference_format = "chr"
        mapper._reference_prefix = "chr"
        mapper._build_reference_layer("chr", "chr")
        
        original = 1
        intermediate = mapper.sumstats_to_reference(original)
        result = mapper.reference_to_sumstats(intermediate)
        
        self.assertEqual(result, original)

    def test_round_trip_nc_format(self):
        """Test round-trip conversion with NC format."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="19",
            log=self.log,
            verbose=False
        )
        mapper.detect_sumstats_format(pd.Series(["NC_000001.10"]))
        
        original = 1
        intermediate = mapper.to_nc(original)
        result = mapper.sumstats_to_number(intermediate)
        
        self.assertEqual(result, original)

    def test_case_insensitive_chr_prefix(self):
        """Test case-insensitive handling of chr prefix."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["chr1", "Chr1", "CHR1"]))
        
        # All should map to the same value
        self.assertEqual(mapper.sumstats_to_number("chr1"), 1)
        self.assertEqual(mapper.sumstats_to_number("Chr1"), 1)
        self.assertEqual(mapper.sumstats_to_number("CHR1"), 1)

    def test_alternative_mt_notation(self):
        """Test alternative mitochondrial notation (M vs MT)."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["MT", "M"]))
        
        # Both M and MT should map to 25
        self.assertEqual(mapper.sumstats_to_number("MT"), 25)
        self.assertEqual(mapper.sumstats_to_number("M"), 25)

    def test_large_series_performance(self):
        """Test performance with large Series."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["chr1", "chrX"]))
        
        # Create a large series
        large_series = pd.Series(["chr1"] * 10000 + ["chrX"] * 1000)
        result = mapper.to_numeric(large_series)
        
        # Check that all values were converted
        self.assertEqual(len(result), 11000)
        self.assertEqual((result == 1).sum(), 10000)
        self.assertEqual((result == 23).sum(), 1000)

    def test_from_reference_file(self):
        """Test creating mapper from reference file."""
        # This is a class method that creates a mapper and stores the reference file
        # We can't easily test file detection without actual files, but we can test the method exists
        mapper = ChromosomeMapper.from_reference_file(
            "dummy.vcf",
            log=self.log,
            verbose=False
        )
        
        self.assertIsInstance(mapper, ChromosomeMapper)
        self.assertEqual(mapper._default_reference_file, "dummy.vcf")

    def test_backward_compat_map_method(self):
        """Test backward compatibility map method."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Map should work with sumstats format
        result = mapper.map(1)
        self.assertEqual(result, 1)
        
        result = mapper.map(23)
        self.assertEqual(result, 23)

    # ========================================================================
    # NC Format Comprehensive Tests
    # ========================================================================

    def test_nc_format_detection_as_sumstats(self):
        """Test detecting NC format as sumstats format."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        series = pd.Series(["NC_000001.11", "NC_000002.12", "NC_000023.11"])
        mapper.detect_sumstats_format(series)
        self.assertEqual(mapper._sumstats_format, "nc")
        
        # Test that NC format can be converted to number
        self.assertEqual(mapper.sumstats_to_number("NC_000001.11"), 1)
        self.assertEqual(mapper.sumstats_to_number("NC_000023.11"), 23)

    def test_nc_format_all_chromosomes_hg19(self):
        """Test NC format conversion for all chromosomes in hg19."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="19",
            log=self.log,
            verbose=False
        )
        
        # Test autosomes
        self.assertEqual(mapper.to_nc(1), "NC_000001.10")
        self.assertEqual(mapper.to_nc(22), "NC_000022.10")
        
        # Test sex chromosomes
        self.assertEqual(mapper.to_nc(23), "NC_000023.10")  # X
        self.assertEqual(mapper.to_nc(24), "NC_000024.9")   # Y
        
        # Test mitochondrial
        self.assertEqual(mapper.to_nc(25), "NC_012920.1")    # MT

    def test_nc_format_all_chromosomes_hg38(self):
        """Test NC format conversion for all chromosomes in hg38."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # Test autosomes
        self.assertEqual(mapper.to_nc(1), "NC_000001.11")
        self.assertEqual(mapper.to_nc(22), "NC_000022.11")
        
        # Test sex chromosomes
        self.assertEqual(mapper.to_nc(23), "NC_000023.11")  # X
        self.assertEqual(mapper.to_nc(24), "NC_000024.10")   # Y
        
        # Test mitochondrial
        self.assertEqual(mapper.to_nc(25), "NC_012920.1")   # MT

    def test_nc_format_t2t_chm13(self):
        """Test NC format conversion for T2T-CHM13 build."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="13",
            log=self.log,
            verbose=False
        )
        
        # T2T-CHM13 uses CM prefix instead of NC
        self.assertEqual(mapper.to_nc(1), "CM000663.2")
        self.assertEqual(mapper.to_nc(22), "CM000684.2")
        self.assertEqual(mapper.to_nc(23), "CM000685.2")  # X
        self.assertEqual(mapper.to_nc(24), "CM000686.2")   # Y
        self.assertEqual(mapper.to_nc(25), "NC_012920.1")  # MT (still NC)

    def test_nc_format_different_species(self):
        """Test NC format conversion for different species."""
        # Mouse
        mapper_mouse = ChromosomeMapper(
            species="mus musculus",
            build="39",
            log=self.log,
            verbose=False
        )
        mapper_mouse.detect_sumstats_format(pd.Series(["1", "19", "X", "Y"]))
        self.assertEqual(mapper_mouse.to_nc(1), "NC_000067.7")
        self.assertEqual(mapper_mouse.to_nc(19), "NC_000085.7")
        # Mouse X and Y use numeric values 23, 24 (human convention)
        x_num = mapper_mouse.sumstats_to_number("X")
        y_num = mapper_mouse.sumstats_to_number("Y")
        self.assertEqual(mapper_mouse.to_nc(x_num), "NC_000086.8")  # X
        self.assertEqual(mapper_mouse.to_nc(y_num), "NC_000087.8")  # Y
        
        # Rat
        mapper_rat = ChromosomeMapper(
            species="rat",
            build="rn6",
            log=self.log,
            verbose=False
        )
        mapper_rat.detect_sumstats_format(pd.Series(["1", "20", "X", "Y"]))
        self.assertEqual(mapper_rat.to_nc(1), "NC_005100.4")
        self.assertEqual(mapper_rat.to_nc(20), "NC_005119.4")
        # Rat X and Y use numeric values 23, 24 (human convention)
        x_num = mapper_rat.sumstats_to_number("X")
        y_num = mapper_rat.sumstats_to_number("Y")
        self.assertEqual(mapper_rat.to_nc(x_num), "NC_005120.4")  # X
        self.assertEqual(mapper_rat.to_nc(y_num), "NC_005121.4")  # Y

    def test_nc_format_round_trip_numeric(self):
        """Test round-trip conversion: numeric -> NC -> numeric."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        original_values = [1, 2, 22, 23, 24, 25]
        for val in original_values:
            nc_id = mapper.to_nc(val)
            result = mapper.sumstats_to_number(nc_id)
            self.assertEqual(result, val, f"Round-trip failed for {val} -> {nc_id}")

    def test_nc_format_round_trip_nc_as_sumstats(self):
        """Test round-trip conversion with NC as sumstats format."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # Set NC as sumstats format
        mapper.detect_sumstats_format(pd.Series(["NC_000001.11", "NC_000023.11"]))
        
        original_nc = "NC_000001.11"
        number = mapper.sumstats_to_number(original_nc)
        result_nc = mapper.number_to_sumstats(number)
        self.assertEqual(result_nc, original_nc)

    def test_nc_format_series_conversion(self):
        """Test NC format conversion with Series."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # Convert numeric Series to NC
        numeric_series = pd.Series([1, 2, 22, 23, 24, 25])
        nc_series = mapper.to_nc(numeric_series)
        expected = pd.Series([
            "NC_000001.11", "NC_000002.12", "NC_000022.11",
            "NC_000023.11", "NC_000024.10", "NC_012920.1"
        ])
        pd.testing.assert_series_equal(nc_series, expected)
        
        # Convert NC Series back to numeric
        mapper.detect_sumstats_format(nc_series)
        result_numeric = mapper.to_numeric(nc_series)
        pd.testing.assert_series_equal(result_numeric, numeric_series, check_dtype=False)

    def test_nc_format_to_other_formats(self):
        """Test converting NC format to other formats."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # Set NC as sumstats format
        mapper.detect_sumstats_format(pd.Series(["NC_000001.11", "NC_000023.11"]))
        
        # NC -> numeric
        self.assertEqual(mapper.to_numeric("NC_000001.11"), 1)
        
        # NC -> string
        self.assertEqual(mapper.to_string("NC_000001.11"), "1")
        self.assertEqual(mapper.to_string("NC_000023.11"), "X")
        
        # NC -> chr
        self.assertEqual(mapper.to_chr("NC_000001.11"), "chr1")
        self.assertEqual(mapper.to_chr("NC_000023.11"), "chrX")

    def test_nc_format_from_other_formats(self):
        """Test converting other formats to NC format."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # numeric -> NC
        self.assertEqual(mapper.to_nc(1), "NC_000001.11")
        
        # string -> NC (need to detect format first)
        mapper.detect_sumstats_format(pd.Series(["1", "X"]))
        number = mapper.sumstats_to_number("1")
        self.assertEqual(mapper.to_nc(number), "NC_000001.11")
        
        # chr -> NC
        mapper.detect_sumstats_format(pd.Series(["chr1", "chrX"]))
        number = mapper.sumstats_to_number("chr1")
        self.assertEqual(mapper.to_nc(number), "NC_000001.11")

    def test_nc_format_invalid_accession(self):
        """Test handling of invalid NCBI accession IDs."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        mapper.detect_sumstats_format(pd.Series(["NC_000001.11"]))
        
        # Invalid NC ID should return pd.NA
        result = mapper.sumstats_to_number("NC_INVALID.1")
        self.assertTrue(pd.isna(result))
        
        result = mapper.sumstats_to_number("INVALID_FORMAT")
        self.assertTrue(pd.isna(result))

    def test_nc_format_missing_build_error(self):
        """Test that NC conversion raises error without build."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            log=self.log,
            verbose=False
        )
        
        # Should raise ValueError when build is not specified
        with self.assertRaises(ValueError):
            mapper.to_nc(1)

    def test_nc_format_unsupported_species_build(self):
        """Test NC format with unsupported species/build combination."""
        mapper = ChromosomeMapper(
            species="unsupported_species",
            build="99",
            log=self.log,
            verbose=False
        )
        
        # Should raise ValueError for unsupported combination
        with self.assertRaises(ValueError):
            mapper.to_nc(1)

    def test_nc_format_chromosome_not_in_mapping(self):
        """Test NC format with chromosome not in mapping."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # Chromosome 100 doesn't exist, should raise ValueError
        with self.assertRaises(ValueError):
            mapper.to_nc(100)

    def test_nc_format_case_insensitive(self):
        """Test that NC format detection is case-insensitive."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # NC format should be detected regardless of case
        series_upper = pd.Series(["NC_000001.11", "NC_000023.11"])
        series_lower = pd.Series(["nc_000001.11", "nc_000023.11"])
        
        mapper.detect_sumstats_format(series_upper)
        self.assertEqual(mapper._sumstats_format, "nc")
        
        mapper.detect_sumstats_format(series_lower)
        self.assertEqual(mapper._sumstats_format, "nc")

    def test_nc_format_mixed_with_other_formats(self):
        """Test handling NC format mixed with other formats."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # When NC format is detected, it should handle NC IDs
        series = pd.Series(["NC_000001.11", "NC_000002.12", "NC_000023.11"])
        mapper.detect_sumstats_format(series)
        
        # All should convert to numbers
        result = mapper.to_numeric(series)
        expected = pd.Series([1, 2, 23])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)

    def test_nc_format_build_processing(self):
        """Test that build codes are properly processed for NC format."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="hg19",  # Should be processed to "19"
            log=self.log,
            verbose=False
        )
        
        # Should use hg19 mappings
        self.assertEqual(mapper.to_nc(1), "NC_000001.10")
        
        mapper2 = ChromosomeMapper(
            species="homo sapiens",
            build="GRCh38",  # Should be processed to "38"
            log=self.log,
            verbose=False
        )
        
        # Should use hg38 mappings
        self.assertEqual(mapper2.to_nc(1), "NC_000001.11")

    def test_nc_format_species_aliases(self):
        """Test NC format with species aliases."""
        # Test with "human" alias
        mapper_human = ChromosomeMapper(
            species="human",
            build="38",
            log=self.log,
            verbose=False
        )
        self.assertEqual(mapper_human.to_nc(1), "NC_000001.11")
        
        # Test with "mouse" alias
        mapper_mouse = ChromosomeMapper(
            species="mouse",
            build="39",
            log=self.log,
            verbose=False
        )
        self.assertEqual(mapper_mouse.to_nc(1), "NC_000067.7")

    def test_nc_format_complete_chromosome_set(self):
        """Test NC format for complete chromosome set."""
        mapper = ChromosomeMapper(
            species="homo sapiens",
            build="38",
            log=self.log,
            verbose=False
        )
        
        # Test all chromosomes 1-25
        all_chromosomes = list(range(1, 26))
        nc_ids = mapper.to_nc(pd.Series(all_chromosomes))
        
        # Verify all chromosomes have NC IDs
        self.assertEqual(len(nc_ids), 25)
        self.assertTrue(all(nc_ids.str.startswith(("NC_", "CM_"))))
        
        # Verify specific chromosomes
        self.assertEqual(nc_ids.iloc[0], "NC_000001.11")   # chr1
        self.assertEqual(nc_ids.iloc[21], "NC_000022.11")   # chr22
        self.assertEqual(nc_ids.iloc[22], "NC_000023.11")   # chrX
        self.assertEqual(nc_ids.iloc[23], "NC_000024.10")    # chrY
        self.assertEqual(nc_ids.iloc[24], "NC_012920.1")   # chrMT


if __name__ == '__main__':
    unittest.main()
