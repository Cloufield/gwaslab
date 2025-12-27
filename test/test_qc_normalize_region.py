"""
Test file for _normalize_region function.

To run tests:
    python test/test_qc_normalize_region.py -v
    python -m unittest test.test_qc_normalize_region -v
    pytest test/test_qc_normalize_region.py -v
"""

import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
from gwaslab.qc.qc_normalize_args import _normalize_region, _parse_flanking
from gwaslab.info.g_Log import Log


def make_test_sumstats():
    """Create a test sumstats DataFrame for testing snpid:flanking format."""
    rows = [
        {"CHR": 1, "POS": 1000, "EA": "A", "NEA": "G", "SNPID": "1:1000:A:G", "rsID": "rs1000"},
        {"CHR": 1, "POS": 2000, "EA": "C", "NEA": "T", "SNPID": "1:2000:C:T", "rsID": "rs2000"},
        {"CHR": 2, "POS": 1500, "EA": "A", "NEA": "T", "SNPID": "2:1500:A:T", "rsID": "rs1500"},
        {"CHR": 2, "POS": 2500, "EA": "G", "NEA": "C", "SNPID": "2:2500:G:C", "rsID": "rs2500"},
    ]
    return pd.DataFrame(rows)


class TestParseFlanking(unittest.TestCase):
    """Test the _parse_flanking helper function."""
    
    def test_parse_base_pairs(self):
        """Test parsing base pairs (no kb suffix)."""
        self.assertEqual(_parse_flanking("500"), 500)
        self.assertEqual(_parse_flanking("1000"), 1000)
        self.assertEqual(_parse_flanking(" 250 "), 250)
    
    def test_parse_kilobases(self):
        """Test parsing kilobases with kb suffix."""
        self.assertEqual(_parse_flanking("500kb"), 500000)
        self.assertEqual(_parse_flanking("1kb"), 1000)
        self.assertEqual(_parse_flanking("0.5kb"), 500)
        self.assertEqual(_parse_flanking("500KB"), 500000)  # uppercase
        self.assertEqual(_parse_flanking("500 kb"), 500000)  # with space
        self.assertEqual(_parse_flanking(" 500 kb "), 500000)  # with spaces


class TestNormalizeRegionChrStartEnd(unittest.TestCase):
    """Test chr:start-end format."""
    
    def test_basic_chr_start_end(self):
        """Test basic chr:start-end format."""
        result = _normalize_region("chr1:1000-2000", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_chr_without_prefix(self):
        """Test chromosome without 'chr' prefix."""
        result = _normalize_region("1:1000-2000", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_chr_uppercase(self):
        """Test uppercase chromosome."""
        result = _normalize_region("CHR1:1000-2000", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_swapped_start_end(self):
        """Test that start > end gets swapped."""
        result = _normalize_region("chr1:2000-1000", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_large_numbers(self):
        """Test with large position numbers."""
        result = _normalize_region("chr1:12345678-23456789", verbose=False)
        self.assertEqual(result, (1, 12345678, 23456789))
    
    def test_float_positions(self):
        """Test that float positions are converted to int."""
        result = _normalize_region("chr1:1000.5-2000.7", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))


class TestNormalizeRegionChrPosFlanking(unittest.TestCase):
    """Test chr:pos:flanking format."""
    
    def test_basic_chr_pos_flanking(self):
        """Test basic chr:pos:flanking format."""
        result = _normalize_region("chr1:1500:500", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_chr_pos_flanking_with_kb(self):
        """Test chr:pos:flanking with kb suffix."""
        result = _normalize_region("chr1:1500:500kb", verbose=False)
        # 1500 - 500kb = 1500 - 500000 = -498500, 1500 + 500kb = 501500
        self.assertEqual(result, (1, -498500, 501500))
    
    def test_chr_pos_flanking_uppercase_kb(self):
        """Test chr:pos:flanking with uppercase KB."""
        result = _normalize_region("chr1:1500:500KB", verbose=False)
        # 1500 - 500kb = 1500 - 500000 = -498500, 1500 + 500kb = 501500
        self.assertEqual(result, (1, -498500, 501500))
    
    def test_chr_pos_flanking_with_space_kb(self):
        """Test chr:pos:flanking with space before kb."""
        result = _normalize_region("chr1:1500:500 kb", verbose=False)
        # 1500 - 500kb = 1500 - 500000 = -498500, 1500 + 500kb = 501500
        self.assertEqual(result, (1, -498500, 501500))
    
    def test_chr_pos_flanking_with_reasonable_kb(self):
        """Test chr:pos:flanking with reasonable kb that doesn't cause negative positions."""
        result = _normalize_region("chr1:1000000:500kb", verbose=False)
        # 1000000 - 500kb = 1000000 - 500000 = 500000, 1000000 + 500kb = 1500000
        self.assertEqual(result, (1, 500000, 1500000))
    
    def test_chr_pos_flanking_decimal_kb(self):
        """Test chr:pos:flanking with decimal kb."""
        result = _normalize_region("chr1:1500:0.5kb", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))  # 1500 - 500bp to 1500 + 500bp
    
    def test_chr_without_prefix_flanking(self):
        """Test chromosome without 'chr' prefix with flanking."""
        result = _normalize_region("1:1500:500", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_different_chromosomes(self):
        """Test different chromosomes."""
        result = _normalize_region("chr2:5000:1000", verbose=False)
        self.assertEqual(result, (2, 4000, 6000))
        result = _normalize_region("chr3:10000:2000", verbose=False)
        self.assertEqual(result, (3, 8000, 12000))


class TestNormalizeRegionSnpidFlanking(unittest.TestCase):
    """Test snpid:flanking format."""
    
    def setUp(self):
        """Set up test sumstats."""
        self.sumstats = make_test_sumstats()
    
    def test_rsid_flanking(self):
        """Test rsID:flanking format."""
        result = _normalize_region("rs1000:500", sumstats=self.sumstats, verbose=False)
        self.assertEqual(result, (1, 500, 1500))  # POS 1000, flanking 500
    
    def test_rsid_flanking_with_kb(self):
        """Test rsID:flanking with kb suffix."""
        result = _normalize_region("rs1000:500kb", sumstats=self.sumstats, verbose=False)
        # POS 1000, flanking 500kb: 1000 - 500000 = -499000, 1000 + 500000 = 501000
        self.assertEqual(result, (1, -499000, 501000))
    
    def test_rsid_flanking_with_reasonable_kb(self):
        """Test rsID:flanking with reasonable kb that doesn't cause negative positions."""
        # Use a position that's large enough
        sumstats = pd.DataFrame({
            "CHR": [1],
            "POS": [1000000],
            "rsID": ["rs1000"]
        })
        result = _normalize_region("rs1000:500kb", sumstats=sumstats, verbose=False)
        # POS 1000000, flanking 500kb: 1000000 - 500000 = 500000, 1000000 + 500000 = 1500000
        self.assertEqual(result, (1, 500000, 1500000))
    
    def test_snpid_string_flanking(self):
        """Test SNPID string:flanking format."""
        result = _normalize_region("1:1000:A:G:500", sumstats=self.sumstats, verbose=False)
        self.assertEqual(result, (1, 500, 1500))
    
    def test_snpid_string_flanking_with_kb(self):
        """Test SNPID string:flanking with kb suffix."""
        result = _normalize_region("1:1000:A:G:500kb", sumstats=self.sumstats, verbose=False)
        # POS 1000, flanking 500kb: 1000 - 500000 = -499000, 1000 + 500000 = 501000
        self.assertEqual(result, (1, -499000, 501000))
    
    def test_snpid_coordinate_format(self):
        """Test coordinate-based SNPID format (chr:pos:ea:nea:flanking)."""
        result = _normalize_region("1:2000:C:T:1000", sumstats=self.sumstats, verbose=False)
        self.assertEqual(result, (1, 1000, 3000))
    
    def test_snpid_coordinate_format_with_kb(self):
        """Test coordinate-based SNPID format with kb."""
        result = _normalize_region("2:1500:A:T:1kb", sumstats=self.sumstats, verbose=False)
        self.assertEqual(result, (2, 500, 2500))  # POS 1500, flanking 1kb = 1000bp
    
    def test_snpid_allele_swapping(self):
        """Test that allele swapping is handled correctly."""
        # The sumstats has EA=A, NEA=G at POS 1000
        # Searching with swapped alleles should still find it
        result = _normalize_region("1:1000:G:A:500", sumstats=self.sumstats, verbose=False)
        self.assertEqual(result, (1, 500, 1500))
    
    def test_snpid_not_found(self):
        """Test error when SNP is not found."""
        with self.assertRaises(ValueError) as context:
            _normalize_region("rs9999:500", sumstats=self.sumstats, verbose=False)
        self.assertIn("not found", str(context.exception))
    
    def test_snpid_flanking_no_sumstats(self):
        """Test error when sumstats is not provided."""
        with self.assertRaises(ValueError) as context:
            _normalize_region("rs1000:500", verbose=False)
        self.assertIn("requires sumstats", str(context.exception))


class TestNormalizeRegionTupleList(unittest.TestCase):
    """Test tuple/list input format."""
    
    def test_tuple_input(self):
        """Test tuple input."""
        result = _normalize_region((1, 1000, 2000), verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_list_input(self):
        """Test list input."""
        result = _normalize_region([2, 1500, 2500], verbose=False)
        self.assertEqual(result, (2, 1500, 2500))
    
    def test_tuple_with_chr_string(self):
        """Test tuple with chromosome as string."""
        result = _normalize_region(("chr1", 1000, 2000), verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_tuple_swapped_start_end(self):
        """Test tuple with swapped start and end."""
        result = _normalize_region((1, 2000, 1000), verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_tuple_float_positions(self):
        """Test tuple with float positions."""
        result = _normalize_region((1, 1000.5, 2000.7), verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_tuple_invalid_length(self):
        """Test error with invalid tuple length."""
        with self.assertRaises(ValueError) as context:
            _normalize_region((1, 1000), verbose=False)
        self.assertIn("must be a tuple/list of (chr, start, end)", str(context.exception))


class TestNormalizeRegionEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def test_none_input(self):
        """Test None input."""
        result = _normalize_region(None, verbose=False)
        self.assertIsNone(result)
    
    def test_invalid_format(self):
        """Test invalid format."""
        with self.assertRaises(ValueError) as context:
            _normalize_region("invalid_format", verbose=False)
        self.assertIn("must be in one of these formats", str(context.exception))
    
    def test_invalid_chromosome(self):
        """Test invalid chromosome."""
        # chr99 might be converted to int 99, which is technically valid
        # Let's test with a truly invalid chromosome string
        with self.assertRaises(ValueError) as context:
            _normalize_region("chrINVALID:1000-2000", verbose=False)
        self.assertIn("not recognized", str(context.exception))
    
    def test_empty_string(self):
        """Test empty string."""
        with self.assertRaises(ValueError):
            _normalize_region("", verbose=False)
    
    def test_whitespace_handling(self):
        """Test that whitespace is handled correctly."""
        result = _normalize_region("  chr1  :  1000  -  2000  ", verbose=False)
        self.assertEqual(result, (1, 1000, 2000))
    
    def test_custom_chr_dict(self):
        """Test with custom chromosome dictionary."""
        custom_dict = {"1": 1, "2": 2, "X": 23}
        result = _normalize_region("chr1:1000-2000", chr_dict=custom_dict, verbose=False)
        self.assertEqual(result, (1, 1000, 2000))


class TestNormalizeRegionWithSumstats(unittest.TestCase):
    """Test normalize_region with various sumstats configurations."""
    
    def test_sumstats_with_rsid_only(self):
        """Test with sumstats that only has rsID column."""
        sumstats = pd.DataFrame({
            "CHR": [1],
            "POS": [1000],
            "rsID": ["rs1000"]
        })
        result = _normalize_region("rs1000:500", sumstats=sumstats, verbose=False)
        self.assertEqual(result, (1, 500, 1500))
    
    def test_sumstats_with_snpid_only(self):
        """Test with sumstats that only has SNPID column."""
        sumstats = pd.DataFrame({
            "CHR": [1],
            "POS": [1000],
            "SNPID": ["1:1000:A:G"]
        })
        result = _normalize_region("1:1000:A:G:500", sumstats=sumstats, verbose=False)
        self.assertEqual(result, (1, 500, 1500))
    
    def test_sumstats_without_allele_columns(self):
        """Test with sumstats that doesn't have EA/NEA columns."""
        sumstats = pd.DataFrame({
            "CHR": [1],
            "POS": [1000],
            "rsID": ["rs1000"]
        })
        # Should still work for coordinate-based search without alleles
        result = _normalize_region("1:1000:500", sumstats=sumstats, verbose=False)
        # This should use chr:pos:flanking format, not snpid:flanking
        self.assertEqual(result, (1, 500, 1500))
        
        # Test coordinate-based snpid format without alleles - should search by chr:pos only
        sumstats2 = pd.DataFrame({
            "CHR": [1],
            "POS": [1000]
        })
        result2 = _normalize_region("1:1000:A:G:500", sumstats=sumstats2, verbose=False)
        # Should find by chr:pos only (no allele matching)
        self.assertEqual(result2, (1, 500, 1500))
    
    def test_custom_column_names(self):
        """Test with custom column names."""
        sumstats = pd.DataFrame({
            "chromosome": [1],
            "position": [1000],
            "rsID": ["rs1000"]
        })
        result = _normalize_region(
            "rs1000:500",
            sumstats=sumstats,
            chrom_col="chromosome",
            pos_col="position",
            verbose=False
        )
        self.assertEqual(result, (1, 500, 1500))


if __name__ == "__main__":
    unittest.main()

