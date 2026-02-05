"""
Test suite for FASTA reading functions with various chromosome naming formats.

Tests cover:
- Numeric chromosome names: "1", "2", "22", "23", "24", "25"
- String format: "1", "X", "Y", "MT"
- Chr-prefixed lowercase: "chr1", "chrX", "chrY", "chrMT"
- Chr-prefixed mixed case: "Chr1", "CHR1"
- Mixed formats in same file
"""

import os
import sys
import gzip
import unittest
import tempfile
import shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.io.io_fasta import (
    parse_fasta,
    load_fasta_auto,
    load_fasta_filtered,
    load_and_build_fasta_records,
    write_fasta,
    FastaRecord,
)
from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log
import pandas as pd


class TestReadFastaChromosomeFormats(unittest.TestCase):
    """Test suite for FASTA reading with various chromosome naming formats."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp(prefix="test_fasta_chr_formats_")
        self.log = Log()
        self.log.verbose = False

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _create_fasta_file(self, filename, records, gzipped=True):
        """
        Create a FASTA file with specified records.
        
        Parameters
        ----------
        filename : str
            Name of the file to create
        records : list of tuple
            List of (chromosome_name, sequence) tuples
        gzipped : bool
            Whether to gzip the file
        
        Returns
        -------
        str
            Path to created file
        """
        fasta_path = os.path.join(self.tmpdir, filename)
        
        lines = []
        for chrom_name, sequence in records:
            lines.append(f">{chrom_name}")
            # Wrap sequence at 60 characters
            for i in range(0, len(sequence), 60):
                lines.append(sequence[i:i+60])
        
        content = "\n".join(lines) + "\n" if lines else ""
        
        if gzipped:
            with gzip.open(fasta_path, "wt") as f:
                f.write(content)
        else:
            with open(fasta_path, "w") as f:
                f.write(content)
        
        return fasta_path

    def _make_sequence(self, length=100, seed=42):
        """Generate a random DNA sequence."""
        import random
        random.seed(seed)
        bases = "ACGT"
        return "".join(random.choice(bases) for _ in range(length))

    # ========================================================================
    # Tests for parse_fasta with numeric chromosome format
    # ========================================================================

    def test_parse_fasta_numeric_chromosome_1(self):
        """Test parsing FASTA with numeric chromosome '1'."""
        seq = self._make_sequence(100, seed=1)
        records = [("1", seq)]
        fasta_path = self._create_fasta_file("test_num1.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("1", result)
        self.assertEqual(result["1"], seq)

    def test_parse_fasta_numeric_chromosome_22(self):
        """Test parsing FASTA with numeric chromosome '22'."""
        seq = self._make_sequence(100, seed=22)
        records = [("22", seq)]
        fasta_path = self._create_fasta_file("test_num22.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("22", result)
        self.assertEqual(result["22"], seq)

    # ========================================================================
    # Tests for parse_fasta with string chromosome format (X, Y, MT)
    # ========================================================================

    def test_parse_fasta_string_chromosome_X(self):
        """Test parsing FASTA with string chromosome 'X'."""
        seq = self._make_sequence(100, seed=23)
        records = [("X", seq)]
        fasta_path = self._create_fasta_file("test_strX.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("X", result)
        self.assertEqual(result["X"], seq)

    def test_parse_fasta_string_chromosome_Y(self):
        """Test parsing FASTA with string chromosome 'Y'."""
        seq = self._make_sequence(100, seed=24)
        records = [("Y", seq)]
        fasta_path = self._create_fasta_file("test_strY.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("Y", result)
        self.assertEqual(result["Y"], seq)

    def test_parse_fasta_string_chromosome_MT(self):
        """Test parsing FASTA with string chromosome 'MT'."""
        seq = self._make_sequence(100, seed=25)
        records = [("MT", seq)]
        fasta_path = self._create_fasta_file("test_strMT.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("MT", result)
        self.assertEqual(result["MT"], seq)

    # ========================================================================
    # Tests for parse_fasta with chr-prefixed format (lowercase)
    # ========================================================================

    def test_parse_fasta_chr_lowercase_chr1(self):
        """Test parsing FASTA with chr-prefixed 'chr1'."""
        seq = self._make_sequence(100, seed=1)
        records = [("chr1", seq)]
        fasta_path = self._create_fasta_file("test_chr1.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("chr1", result)
        self.assertEqual(result["chr1"], seq)

    def test_parse_fasta_chr_lowercase_chrX(self):
        """Test parsing FASTA with chr-prefixed 'chrX'."""
        seq = self._make_sequence(100, seed=23)
        records = [("chrX", seq)]
        fasta_path = self._create_fasta_file("test_chrX.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("chrX", result)
        self.assertEqual(result["chrX"], seq)

    def test_parse_fasta_chr_lowercase_chrMT(self):
        """Test parsing FASTA with chr-prefixed 'chrMT'."""
        seq = self._make_sequence(100, seed=25)
        records = [("chrMT", seq)]
        fasta_path = self._create_fasta_file("test_chrMT.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("chrMT", result)
        self.assertEqual(result["chrMT"], seq)

    # ========================================================================
    # Tests for parse_fasta with chr-prefixed format (mixed case)
    # ========================================================================

    def test_parse_fasta_chr_mixed_case_Chr1(self):
        """Test parsing FASTA with mixed case 'Chr1'."""
        seq = self._make_sequence(100, seed=1)
        records = [("Chr1", seq)]
        fasta_path = self._create_fasta_file("test_Chr1.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("Chr1", result)
        self.assertEqual(result["Chr1"], seq)

    def test_parse_fasta_chr_mixed_case_CHR1(self):
        """Test parsing FASTA with uppercase 'CHR1'."""
        seq = self._make_sequence(100, seed=1)
        records = [("CHR1", seq)]
        fasta_path = self._create_fasta_file("test_CHR1.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("CHR1", result)
        self.assertEqual(result["CHR1"], seq)

    def test_parse_fasta_chr_mixed_case_ChrX(self):
        """Test parsing FASTA with mixed case 'ChrX'."""
        seq = self._make_sequence(100, seed=23)
        records = [("ChrX", seq)]
        fasta_path = self._create_fasta_file("test_ChrX.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("ChrX", result)
        self.assertEqual(result["ChrX"], seq)

    def test_parse_fasta_chr_mixed_case_CHRX(self):
        """Test parsing FASTA with uppercase 'CHRX'."""
        seq = self._make_sequence(100, seed=23)
        records = [("CHRX", seq)]
        fasta_path = self._create_fasta_file("test_CHRX.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("CHRX", result)
        self.assertEqual(result["CHRX"], seq)

    # ========================================================================
    # Tests for load_fasta_filtered with chromosome name normalization
    # ========================================================================

    def test_load_fasta_filtered_numeric_format(self):
        """Test load_fasta_filtered with numeric chromosome names."""
        seq1 = self._make_sequence(100, seed=1)
        seq2 = self._make_sequence(100, seed=2)
        records = [("1", seq1), ("2", seq2)]
        fasta_path = self._create_fasta_file("test_filtered_num.fa.gz", records)
        
        chromlist_set = {1, 2}
        chroms_in_sumstats_set = {1, 2}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2]))
        
        result = load_fasta_filtered(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            log=self.log,
            verbose=False
        )
        
        self.assertEqual(len(result), 2)
        self.assertIn(1, result)
        self.assertIn(2, result)

    def test_load_fasta_filtered_chr_prefixed_format(self):
        """Test load_fasta_filtered with chr-prefixed chromosome names."""
        seq1 = self._make_sequence(100, seed=1)
        seqX = self._make_sequence(100, seed=23)
        records = [("chr1", seq1), ("chrX", seqX)]
        fasta_path = self._create_fasta_file("test_filtered_chr.fa.gz", records)
        
        chromlist_set = {1, 23}
        chroms_in_sumstats_set = {1, 23}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 23]))
        
        result = load_fasta_filtered(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            log=self.log,
            verbose=False
        )
        
        self.assertEqual(len(result), 2)
        # Normalized to numeric format
        self.assertIn(1, result)
        self.assertIn(23, result)

    def test_load_fasta_filtered_Chr_mixed_case(self):
        """Test load_fasta_filtered with mixed case chr-prefixed names (Chr1, CHR1)."""
        seq1 = self._make_sequence(100, seed=1)
        seq2 = self._make_sequence(100, seed=2)
        records = [("Chr1", seq1), ("CHR2", seq2)]
        fasta_path = self._create_fasta_file("test_filtered_Chr.fa.gz", records)
        
        chromlist_set = {1, 2}
        chroms_in_sumstats_set = {1, 2}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2]))
        
        result = load_fasta_filtered(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            log=self.log,
            verbose=False
        )
        
        self.assertEqual(len(result), 2)
        self.assertIn(1, result)
        self.assertIn(2, result)

    def test_load_fasta_filtered_sex_chromosomes(self):
        """Test load_fasta_filtered with sex chromosomes (X, Y, MT)."""
        seqX = self._make_sequence(100, seed=23)
        seqY = self._make_sequence(100, seed=24)
        seqMT = self._make_sequence(100, seed=25)
        records = [("X", seqX), ("Y", seqY), ("MT", seqMT)]
        fasta_path = self._create_fasta_file("test_filtered_sex.fa.gz", records)
        
        chromlist_set = {23, 24, 25}
        chroms_in_sumstats_set = {23, 24, 25}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([23, 24, 25]))
        
        result = load_fasta_filtered(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            log=self.log,
            verbose=False
        )
        
        self.assertEqual(len(result), 3)
        self.assertIn(23, result)  # X
        self.assertIn(24, result)  # Y
        self.assertIn(25, result)  # MT

    def test_load_fasta_filtered_chrX_chrY_chrMT(self):
        """Test load_fasta_filtered with chr-prefixed sex chromosomes."""
        seqX = self._make_sequence(100, seed=23)
        seqY = self._make_sequence(100, seed=24)
        seqMT = self._make_sequence(100, seed=25)
        records = [("chrX", seqX), ("chrY", seqY), ("chrMT", seqMT)]
        fasta_path = self._create_fasta_file("test_filtered_chrXYMT.fa.gz", records)
        
        chromlist_set = {23, 24, 25}
        chroms_in_sumstats_set = {23, 24, 25}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([23, 24, 25]))
        
        result = load_fasta_filtered(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            log=self.log,
            verbose=False
        )
        
        self.assertEqual(len(result), 3)
        self.assertIn(23, result)  # chrX -> 23
        self.assertIn(24, result)  # chrY -> 24
        self.assertIn(25, result)  # chrMT -> 25

    # ========================================================================
    # Tests for load_and_build_fasta_records with chromosome name normalization
    # ========================================================================

    def test_load_and_build_fasta_records_numeric(self):
        """Test load_and_build_fasta_records with numeric chromosome names."""
        seq1 = self._make_sequence(100, seed=1)
        seq2 = self._make_sequence(100, seed=2)
        records = [("1", seq1), ("2", seq2)]
        fasta_path = self._create_fasta_file("test_build_num.fa.gz", records)
        
        chromlist_set = {1, 2}
        chroms_in_sumstats_set = {1, 2}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2]))
        
        record, starting_positions, records_len = load_and_build_fasta_records(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            pos_as_dict=True,
            log=self.log,
            verbose=False
        )
        
        self.assertIn(1, starting_positions)
        self.assertIn(2, starting_positions)
        self.assertEqual(records_len[1], 100)
        self.assertEqual(records_len[2], 100)

    def test_load_and_build_fasta_records_chr_prefixed(self):
        """Test load_and_build_fasta_records with chr-prefixed names."""
        seq1 = self._make_sequence(100, seed=1)
        seqX = self._make_sequence(100, seed=23)
        records = [("chr1", seq1), ("chrX", seqX)]
        fasta_path = self._create_fasta_file("test_build_chr.fa.gz", records)
        
        chromlist_set = {1, 23}
        chroms_in_sumstats_set = {1, 23}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 23]))
        
        record, starting_positions, records_len = load_and_build_fasta_records(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            pos_as_dict=True,
            log=self.log,
            verbose=False
        )
        
        self.assertIn(1, starting_positions)
        self.assertIn(23, starting_positions)
        self.assertEqual(records_len[1], 100)
        self.assertEqual(records_len[23], 100)

    def test_load_and_build_fasta_records_mixed_case(self):
        """Test load_and_build_fasta_records with mixed case chr prefix."""
        seq1 = self._make_sequence(100, seed=1)
        seq2 = self._make_sequence(100, seed=2)
        records = [("Chr1", seq1), ("CHR2", seq2)]
        fasta_path = self._create_fasta_file("test_build_mixed.fa.gz", records)
        
        chromlist_set = {1, 2}
        chroms_in_sumstats_set = {1, 2}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2]))
        
        record, starting_positions, records_len = load_and_build_fasta_records(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            pos_as_dict=True,
            log=self.log,
            verbose=False
        )
        
        self.assertIn(1, starting_positions)
        self.assertIn(2, starting_positions)

    # ========================================================================
    # Tests for load_fasta_auto
    # ========================================================================

    def test_load_fasta_auto_as_seqrecord(self):
        """Test load_fasta_auto returning FastaRecord objects."""
        seq1 = self._make_sequence(100, seed=1)
        seq2 = self._make_sequence(100, seed=2)
        records = [("chr1", seq1), ("chr2", seq2)]
        fasta_path = self._create_fasta_file("test_auto.fa.gz", records)
        
        result = list(load_fasta_auto(fasta_path, as_seqrecord=True))
        
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], FastaRecord)
        self.assertEqual(result[0].id, "chr1")
        self.assertEqual(result[1].id, "chr2")

    def test_load_fasta_auto_as_tuple(self):
        """Test load_fasta_auto returning tuples."""
        seq1 = self._make_sequence(100, seed=1)
        seq2 = self._make_sequence(100, seed=2)
        records = [("Chr1", seq1), ("ChrX", seq2)]
        fasta_path = self._create_fasta_file("test_auto_tuple.fa.gz", records)
        
        result = list(load_fasta_auto(fasta_path, as_seqrecord=False))
        
        self.assertEqual(len(result), 2)
        self.assertIsInstance(result[0], tuple)
        self.assertEqual(result[0][0], "Chr1")
        self.assertEqual(result[1][0], "ChrX")

    # ========================================================================
    # Tests for mixed formats in same file
    # ========================================================================

    def test_parse_fasta_mixed_formats(self):
        """Test parsing FASTA with mixed chromosome naming formats."""
        records = [
            ("1", self._make_sequence(100, seed=1)),
            ("2", self._make_sequence(100, seed=2)),
            ("X", self._make_sequence(100, seed=23)),
            ("Y", self._make_sequence(100, seed=24)),
            ("MT", self._make_sequence(100, seed=25)),
        ]
        fasta_path = self._create_fasta_file("test_mixed.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 5)
        self.assertIn("1", result)
        self.assertIn("2", result)
        self.assertIn("X", result)
        self.assertIn("Y", result)
        self.assertIn("MT", result)

    def test_parse_fasta_all_chr_prefixed(self):
        """Test parsing FASTA with all chr-prefixed names."""
        records = [
            ("chr1", self._make_sequence(100, seed=1)),
            ("chr2", self._make_sequence(100, seed=2)),
            ("chrX", self._make_sequence(100, seed=23)),
            ("chrY", self._make_sequence(100, seed=24)),
            ("chrMT", self._make_sequence(100, seed=25)),
        ]
        fasta_path = self._create_fasta_file("test_all_chr.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 5)
        self.assertIn("chr1", result)
        self.assertIn("chr2", result)
        self.assertIn("chrX", result)
        self.assertIn("chrY", result)
        self.assertIn("chrMT", result)

    # ========================================================================
    # Tests for edge cases
    # ========================================================================

    def test_empty_fasta(self):
        """Test parsing an empty FASTA file."""
        fasta_path = self._create_fasta_file("test_empty.fa.gz", [])
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 0)

    def test_uncompressed_fasta(self):
        """Test parsing an uncompressed FASTA file."""
        seq = self._make_sequence(100, seed=1)
        records = [("chr1", seq)]
        fasta_path = self._create_fasta_file("test_uncompressed.fa", records, gzipped=False)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        self.assertIn("chr1", result)
        self.assertEqual(result["chr1"], seq)

    def test_fasta_with_description(self):
        """Test parsing FASTA with description after chromosome name."""
        # Create FASTA with description line: ">chr1 some description"
        fasta_path = os.path.join(self.tmpdir, "test_desc.fa.gz")
        seq = self._make_sequence(100, seed=1)
        with gzip.open(fasta_path, "wt") as f:
            f.write(">chr1 Homo sapiens chromosome 1\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 1)
        # The record name should be just "chr1" (first word)
        self.assertIn("chr1", result)

    def test_all_standard_chromosomes(self):
        """Test FASTA with all standard human chromosomes."""
        records = []
        # Autosomes 1-22
        for i in range(1, 23):
            records.append((str(i), self._make_sequence(50, seed=i)))
        # Sex chromosomes and mitochondrial
        records.append(("X", self._make_sequence(50, seed=23)))
        records.append(("Y", self._make_sequence(50, seed=24)))
        records.append(("MT", self._make_sequence(50, seed=25)))
        
        fasta_path = self._create_fasta_file("test_all_chroms.fa.gz", records)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 25)
        # Check autosomes
        for i in range(1, 23):
            self.assertIn(str(i), result)
        # Check sex chromosomes and mitochondrial
        self.assertIn("X", result)
        self.assertIn("Y", result)
        self.assertIn("MT", result)

    def test_write_and_read_roundtrip(self):
        """Test writing and reading FASTA file roundtrip."""
        original_records = {
            "chr1": self._make_sequence(100, seed=1),
            "chrX": self._make_sequence(100, seed=23),
            "chrMT": self._make_sequence(100, seed=25),
        }
        
        fasta_path = os.path.join(self.tmpdir, "test_roundtrip.fa.gz")
        write_fasta(original_records, fasta_path)
        
        result = parse_fasta(fasta_path, as_dict=True)
        
        self.assertEqual(len(result), 3)
        for name, seq in original_records.items():
            self.assertIn(name, result)
            self.assertEqual(result[name], seq)

    def test_load_fasta_filtered_filtering(self):
        """Test that load_fasta_filtered correctly filters chromosomes."""
        records = [
            ("chr1", self._make_sequence(100, seed=1)),
            ("chr2", self._make_sequence(100, seed=2)),
            ("chr3", self._make_sequence(100, seed=3)),
        ]
        fasta_path = self._create_fasta_file("test_filter.fa.gz", records)
        
        # Only request chromosome 1
        chromlist_set = {1}
        chroms_in_sumstats_set = {1}
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1]))
        
        result = load_fasta_filtered(
            fasta_path,
            chromlist_set,
            chroms_in_sumstats_set,
            mapper=mapper,
            log=self.log,
            verbose=False
        )
        
        # Should only have chromosome 1
        self.assertEqual(len(result), 1)
        self.assertIn(1, result)
        self.assertNotIn(2, result)
        self.assertNotIn(3, result)


if __name__ == "__main__":
    unittest.main()
