"""
Test suite for read_gtf function with various chromosome naming formats.

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

from gwaslab.io.io_gtf import read_gtf


class TestReadGTFChromosomeFormats(unittest.TestCase):
    """Test suite for read_gtf with various chromosome naming formats."""

    def setUp(self):
        """Set up test fixtures."""
        self.tmpdir = tempfile.mkdtemp(prefix="test_gtf_chr_formats_")

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _create_gtf_file(self, filename, records, gzipped=True):
        """
        Create a GTF file with specified records.
        
        Parameters
        ----------
        filename : str
            Name of the file to create
        records : list of dict
            List of GTF records. Each dict should have keys:
            seqname, source, feature, start, end, score, strand, frame, attribute
        gzipped : bool
            Whether to gzip the file
        
        Returns
        -------
        str
            Path to created file
        """
        gtf_path = os.path.join(self.tmpdir, filename)
        
        lines = []
        for rec in records:
            line = "\t".join([
                str(rec.get("seqname", "1")),
                str(rec.get("source", "test")),
                str(rec.get("feature", "gene")),
                str(rec.get("start", 1000)),
                str(rec.get("end", 2000)),
                str(rec.get("score", ".")),
                str(rec.get("strand", "+")),
                str(rec.get("frame", ".")),
                str(rec.get("attribute", 'gene_id "GENE1"; gene_name "TestGene";')),
            ])
            lines.append(line)
        
        content = "\n".join(lines)
        
        if gzipped:
            with gzip.open(gtf_path, "wt") as f:
                f.write(content)
        else:
            with open(gtf_path, "w") as f:
                f.write(content)
        
        return gtf_path

    # ========================================================================
    # Tests for numeric chromosome format
    # ========================================================================

    def test_numeric_chromosome_1(self):
        """Test GTF with numeric chromosome '1'."""
        records = [
            {"seqname": "1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
            {"seqname": "1", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE2";'},
        ]
        gtf_path = self._create_gtf_file("test_num1.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 2)
        self.assertTrue(all(df["seqname"] == "1"))

    def test_numeric_chromosome_22(self):
        """Test GTF with numeric chromosome '22'."""
        records = [
            {"seqname": "22", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
        ]
        gtf_path = self._create_gtf_file("test_num22.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "22")

    def test_numeric_chromosome_sex_chroms(self):
        """Test GTF with numeric sex chromosomes (23, 24, 25)."""
        records = [
            {"seqname": "23", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEX";'},
            {"seqname": "24", "start": 3000, "end": 4000, "attribute": 'gene_id "GENEY";'},
            {"seqname": "25", "start": 5000, "end": 6000, "attribute": 'gene_id "GENEMT";'},
        ]
        gtf_path = self._create_gtf_file("test_num_sex.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 3)
        # After conversion: 23->X, 24->Y, 25->MT
        seqnames = df["seqname"].tolist()
        self.assertIn("X", seqnames)
        self.assertIn("Y", seqnames)
        self.assertIn("MT", seqnames)

    # ========================================================================
    # Tests for string chromosome format (X, Y, MT)
    # ========================================================================

    def test_string_chromosome_X(self):
        """Test GTF with string chromosome 'X'."""
        records = [
            {"seqname": "X", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEX";'},
        ]
        gtf_path = self._create_gtf_file("test_strX.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "X")

    def test_string_chromosome_Y(self):
        """Test GTF with string chromosome 'Y'."""
        records = [
            {"seqname": "Y", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEY";'},
        ]
        gtf_path = self._create_gtf_file("test_strY.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "Y")

    def test_string_chromosome_MT(self):
        """Test GTF with string chromosome 'MT'."""
        records = [
            {"seqname": "MT", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEMT";'},
        ]
        gtf_path = self._create_gtf_file("test_strMT.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "MT")

    # ========================================================================
    # Tests for chr-prefixed format (lowercase)
    # ========================================================================

    def test_chr_lowercase_chr1(self):
        """Test GTF with chr-prefixed 'chr1'."""
        records = [
            {"seqname": "chr1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
            {"seqname": "chr1", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE2";'},
        ]
        gtf_path = self._create_gtf_file("test_chr1.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 2)
        # After normalization, should be "1"
        self.assertTrue(all(df["seqname"] == "1"))

    def test_chr_lowercase_chrX(self):
        """Test GTF with chr-prefixed 'chrX'."""
        records = [
            {"seqname": "chrX", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEX";'},
        ]
        gtf_path = self._create_gtf_file("test_chrX.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "X")

    def test_chr_lowercase_chrY(self):
        """Test GTF with chr-prefixed 'chrY'."""
        records = [
            {"seqname": "chrY", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEY";'},
        ]
        gtf_path = self._create_gtf_file("test_chrY.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "Y")

    def test_chr_lowercase_chrMT(self):
        """Test GTF with chr-prefixed 'chrMT'."""
        records = [
            {"seqname": "chrMT", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEMT";'},
        ]
        gtf_path = self._create_gtf_file("test_chrMT.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "MT")

    # ========================================================================
    # Tests for chr-prefixed format (mixed case: Chr, CHR)
    # ========================================================================

    def test_chr_mixed_case_Chr1(self):
        """Test GTF with mixed case 'Chr1'."""
        records = [
            {"seqname": "Chr1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
        ]
        gtf_path = self._create_gtf_file("test_Chr1.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        # After normalization, should be "1"
        self.assertEqual(df["seqname"].iloc[0], "1")

    def test_chr_mixed_case_CHR1(self):
        """Test GTF with uppercase 'CHR1'."""
        records = [
            {"seqname": "CHR1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
        ]
        gtf_path = self._create_gtf_file("test_CHR1.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        # After normalization, should be "1"
        self.assertEqual(df["seqname"].iloc[0], "1")

    def test_chr_mixed_case_ChrX(self):
        """Test GTF with mixed case 'ChrX'."""
        records = [
            {"seqname": "ChrX", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEX";'},
        ]
        gtf_path = self._create_gtf_file("test_ChrX.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "X")

    def test_chr_mixed_case_CHRX(self):
        """Test GTF with uppercase 'CHRX'."""
        records = [
            {"seqname": "CHRX", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEX";'},
        ]
        gtf_path = self._create_gtf_file("test_CHRX.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "X")

    # ========================================================================
    # Tests for mixed formats in same file
    # ========================================================================

    def test_mixed_formats_numeric_and_string(self):
        """Test GTF with mixed numeric and string chromosome names."""
        records = [
            {"seqname": "1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
            {"seqname": "2", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE2";'},
            {"seqname": "X", "start": 5000, "end": 6000, "attribute": 'gene_id "GENEX";'},
            {"seqname": "Y", "start": 7000, "end": 8000, "attribute": 'gene_id "GENEY";'},
            {"seqname": "MT", "start": 9000, "end": 10000, "attribute": 'gene_id "GENEMT";'},
        ]
        gtf_path = self._create_gtf_file("test_mixed_num_str.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 5)
        seqnames = df["seqname"].tolist()
        self.assertIn("1", seqnames)
        self.assertIn("2", seqnames)
        self.assertIn("X", seqnames)
        self.assertIn("Y", seqnames)
        self.assertIn("MT", seqnames)

    def test_mixed_formats_chr_prefixed(self):
        """Test GTF with mixed chr-prefixed chromosome names."""
        records = [
            {"seqname": "chr1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
            {"seqname": "chr2", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE2";'},
            {"seqname": "chrX", "start": 5000, "end": 6000, "attribute": 'gene_id "GENEX";'},
            {"seqname": "chrY", "start": 7000, "end": 8000, "attribute": 'gene_id "GENEY";'},
            {"seqname": "chrMT", "start": 9000, "end": 10000, "attribute": 'gene_id "GENEMT";'},
        ]
        gtf_path = self._create_gtf_file("test_mixed_chr.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 5)
        seqnames = df["seqname"].tolist()
        # After normalization, chr prefix should be removed
        self.assertIn("1", seqnames)
        self.assertIn("2", seqnames)
        self.assertIn("X", seqnames)
        self.assertIn("Y", seqnames)
        self.assertIn("MT", seqnames)

    # ========================================================================
    # Tests for chrom filtering with various formats
    # ========================================================================

    def test_chrom_filter_numeric(self):
        """Test filtering by chromosome with numeric input."""
        records = [
            {"seqname": "1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1A";'},
            {"seqname": "1", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE1B";'},
            {"seqname": "2", "start": 5000, "end": 6000, "attribute": 'gene_id "GENE2";'},
        ]
        gtf_path = self._create_gtf_file("test_filter_num.gtf.gz", records)
        
        df = read_gtf(gtf_path, chrom="1")
        
        self.assertEqual(len(df), 2)
        self.assertTrue(all(df["seqname"] == "1"))

    def test_chrom_filter_chr_prefixed(self):
        """Test filtering by chromosome with chr-prefixed input on chr-prefixed GTF."""
        records = [
            {"seqname": "chr1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1A";'},
            {"seqname": "chr1", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE1B";'},
            {"seqname": "chr2", "start": 5000, "end": 6000, "attribute": 'gene_id "GENE2";'},
        ]
        gtf_path = self._create_gtf_file("test_filter_chr.gtf.gz", records)
        
        # Filter using "chr1" or "1" should both work
        df1 = read_gtf(gtf_path, chrom="chr1")
        df2 = read_gtf(gtf_path, chrom="1")
        
        self.assertEqual(len(df1), 2)
        self.assertEqual(len(df2), 2)

    def test_chrom_filter_X(self):
        """Test filtering by X chromosome with various input formats."""
        records = [
            {"seqname": "chrX", "start": 1000, "end": 2000, "attribute": 'gene_id "GENEX";'},
            {"seqname": "chr1", "start": 3000, "end": 4000, "attribute": 'gene_id "GENE1";'},
        ]
        gtf_path = self._create_gtf_file("test_filter_X.gtf.gz", records)
        
        # Filter using "X", "chrX", or "23" should all work
        df_X = read_gtf(gtf_path, chrom="X")
        df_chrX = read_gtf(gtf_path, chrom="chrX")
        df_23 = read_gtf(gtf_path, chrom="23")
        
        self.assertEqual(len(df_X), 1)
        self.assertEqual(len(df_chrX), 1)
        self.assertEqual(len(df_23), 1)
        self.assertEqual(df_X["seqname"].iloc[0], "X")

    # ========================================================================
    # Tests for edge cases
    # ========================================================================

    def test_empty_gtf_with_header(self):
        """Test reading a GTF file with only header comments (no records)."""
        # Create GTF with only comment lines (valid but empty)
        gtf_path = os.path.join(self.tmpdir, "test_empty_header.gtf.gz")
        with gzip.open(gtf_path, "wt") as f:
            f.write("#!genome-build GRCh38\n")
            f.write("#!genome-version GRCh38\n")
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 0)

    def test_uncompressed_gtf(self):
        """Test reading an uncompressed GTF file."""
        records = [
            {"seqname": "1", "start": 1000, "end": 2000, "attribute": 'gene_id "GENE1";'},
        ]
        gtf_path = self._create_gtf_file("test_uncompressed.gtf", records, gzipped=False)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["seqname"].iloc[0], "1")

    def test_usecols_with_seqname(self):
        """Test using usecols to select specific columns including seqname."""
        records = [
            {"seqname": "chr1", "feature": "gene", "start": 1000, "end": 2000, 
             "attribute": 'gene_id "GENE1"; gene_name "TestGene";'},
        ]
        gtf_path = self._create_gtf_file("test_usecols.gtf.gz", records)
        
        df = read_gtf(gtf_path, usecols=["seqname", "start", "end", "gene_id"])
        
        self.assertIn("seqname", df.columns)
        self.assertIn("start", df.columns)
        self.assertIn("end", df.columns)
        self.assertIn("gene_id", df.columns)
        self.assertEqual(df["seqname"].iloc[0], "1")

    def test_feature_filter(self):
        """Test filtering by feature type."""
        records = [
            {"seqname": "chr1", "feature": "gene", "start": 1000, "end": 5000, 
             "attribute": 'gene_id "GENE1";'},
            {"seqname": "chr1", "feature": "transcript", "start": 1000, "end": 5000, 
             "attribute": 'gene_id "GENE1"; transcript_id "TX1";'},
            {"seqname": "chr1", "feature": "exon", "start": 1000, "end": 2000, 
             "attribute": 'gene_id "GENE1"; transcript_id "TX1"; exon_id "EX1";'},
        ]
        gtf_path = self._create_gtf_file("test_feature_filter.gtf.gz", records)
        
        df = read_gtf(gtf_path, features={"gene"})
        
        self.assertEqual(len(df), 1)
        self.assertEqual(df["feature"].iloc[0], "gene")
        self.assertEqual(df["seqname"].iloc[0], "1")

    def test_all_chromosomes(self):
        """Test GTF with all standard human chromosomes."""
        records = []
        # Autosomes 1-22
        for i in range(1, 23):
            records.append({
                "seqname": str(i), 
                "start": i * 1000, 
                "end": i * 1000 + 1000, 
                "attribute": f'gene_id "GENE{i}";'
            })
        # Sex chromosomes and mitochondrial
        records.append({"seqname": "X", "start": 23000, "end": 24000, "attribute": 'gene_id "GENEX";'})
        records.append({"seqname": "Y", "start": 24000, "end": 25000, "attribute": 'gene_id "GENEY";'})
        records.append({"seqname": "MT", "start": 1, "end": 1000, "attribute": 'gene_id "GENEMT";'})
        
        gtf_path = self._create_gtf_file("test_all_chroms.gtf.gz", records)
        
        df = read_gtf(gtf_path)
        
        self.assertEqual(len(df), 25)
        seqnames = set(df["seqname"].tolist())
        # Check autosomes
        for i in range(1, 23):
            self.assertIn(str(i), seqnames)
        # Check sex chromosomes and mitochondrial
        self.assertIn("X", seqnames)
        self.assertIn("Y", seqnames)
        self.assertIn("MT", seqnames)


if __name__ == "__main__":
    unittest.main()
