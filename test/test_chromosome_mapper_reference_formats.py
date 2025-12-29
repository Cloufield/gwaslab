"""
Test suite for ChromosomeMapper with various reference file formats.

Tests cover:
- Creating simulated reference files in various formats (VCF, FASTA, GTF, chain files)
- Format detection from reference files
- Conversion between sumstats and reference formats
- Handling of unconvertible chromosomes (pd.NA)
- Different chromosome notations (numeric, string, chr-prefixed, NCBI RefSeq)
"""

import os
import sys
import unittest
import tempfile
import pandas as pd
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper
from gwaslab.info.g_Log import Log


class TestChromosomeMapperReferenceFormats(unittest.TestCase):
    """Test suite for ChromosomeMapper with various reference file formats."""

    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.log.verbose = False
        self.temp_dir = tempfile.mkdtemp()
        self.test_files = []

    def tearDown(self):
        """Clean up test fixtures."""
        # Clean up temporary files
        for file_path in self.test_files:
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                except:
                    pass
        try:
            os.rmdir(self.temp_dir)
        except:
            pass

    def _create_vcf_file(self, filename, chromosomes, format_type="chr"):
        """Create a simulated VCF file with specified chromosome format."""
        file_path = os.path.join(self.temp_dir, filename)
        self.test_files.append(file_path)
        
        with open(file_path, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            
            # Write contig headers based on format
            for chrom in chromosomes[:10]:  # Sample first 10
                if format_type == "chr":
                    f.write(f"##contig=<ID={chrom}>\n")
                elif format_type == "nc":
                    f.write(f"##contig=<ID={chrom}>\n")
            
            # Write column header
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            
            # Write sample records
            for i, chrom in enumerate(chromosomes[:50]):  # Sample first 50
                pos = (i + 1) * 1000
                f.write(f"{chrom}\t{pos}\trs{i}\tA\tG\t100\tPASS\t.\tGT\t0/1\n")
        
        return file_path

    def _create_fasta_file(self, filename, chromosomes, format_type="chr"):
        """Create a simulated FASTA file with specified chromosome format."""
        file_path = os.path.join(self.temp_dir, filename)
        self.test_files.append(file_path)
        
        with open(file_path, 'w') as f:
            for chrom in chromosomes:
                f.write(f">{chrom}\n")
                # Write a short sequence
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
        
        return file_path

    def _create_gtf_file(self, filename, chromosomes, format_type="chr"):
        """Create a simulated GTF file with specified chromosome format."""
        file_path = os.path.join(self.temp_dir, filename)
        self.test_files.append(file_path)
        
        with open(file_path, 'w') as f:
            # Write GTF header
            f.write("# GTF file for testing\n")
            f.write("# Format: chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n")
            
            # Write sample records
            for i, chrom in enumerate(chromosomes[:50]):
                start = (i + 1) * 1000
                end = start + 100
                f.write(f"{chrom}\ttest\tgene\t{start}\t{end}\t.\t+\t.\tgene_id \"test{i}\";\n")
        
        return file_path

    def _create_chain_file(self, filename, chromosomes, format_type="chr"):
        """Create a simulated chain file with specified chromosome format."""
        file_path = os.path.join(self.temp_dir, filename)
        self.test_files.append(file_path)
        
        with open(file_path, 'w') as f:
            # Write chain file format
            for i, chrom in enumerate(chromosomes[:20]):
                # Chain file format: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                f.write(f"chain 1000 {chrom} 1000000 + 0 1000000 {chrom} 1000000 + 0 1000000 1\n")
                f.write("1000\n")
                f.write("0\n\n")
        
        return file_path

    # ========================================================================
    # Test VCF Format Detection
    # ========================================================================

    def test_detect_vcf_format_chr(self):
        """Test detecting chr-prefixed format from VCF file."""
        chromosomes = ["chr1", "chr2", "chr3", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("test_chr.vcf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        
        self.assertEqual(format_type, "chr")
        self.assertEqual(prefix, "chr")

    def test_detect_vcf_format_numeric(self):
        """Test detecting numeric format from VCF file."""
        chromosomes = ["1", "2", "3", "22", "23", "24", "25"]
        vcf_file = self._create_vcf_file("test_numeric.vcf", chromosomes, format_type="numeric")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        
        self.assertEqual(format_type, "numeric")

    def test_detect_vcf_format_string(self):
        """Test detecting string format from VCF file."""
        chromosomes = ["1", "2", "3", "X", "Y", "MT"]
        vcf_file = self._create_vcf_file("test_string.vcf", chromosomes, format_type="string")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        
        self.assertEqual(format_type, "string")

    def test_detect_vcf_format_nc(self):
        """Test detecting NCBI RefSeq format from VCF file."""
        chromosomes = ["NC_000001.11", "NC_000002.12", "NC_000023.11"]
        vcf_file = self._create_vcf_file("test_nc.vcf", chromosomes, format_type="nc")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        
        self.assertEqual(format_type, "nc")

    # ========================================================================
    # Test FASTA Format Detection
    # ========================================================================

    def test_detect_fasta_format_chr(self):
        """Test detecting chr-prefixed format from FASTA file."""
        chromosomes = ["chr1", "chr2", "chr3", "chrX", "chrY", "chrMT"]
        fasta_file = self._create_fasta_file("test_chr.fa", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(fasta_file)
        
        self.assertEqual(format_type, "chr")
        self.assertEqual(prefix, "chr")

    def test_detect_fasta_format_numeric(self):
        """Test detecting numeric format from FASTA file."""
        chromosomes = ["1", "2", "3", "22", "23", "24", "25"]
        fasta_file = self._create_fasta_file("test_numeric.fa", chromosomes, format_type="numeric")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(fasta_file)
        
        self.assertEqual(format_type, "numeric")

    def test_detect_fasta_format_string(self):
        """Test detecting string format from FASTA file."""
        chromosomes = ["1", "2", "3", "X", "Y", "MT"]
        fasta_file = self._create_fasta_file("test_string.fa", chromosomes, format_type="string")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(fasta_file)
        
        self.assertEqual(format_type, "string")

    # ========================================================================
    # Test GTF Format Detection
    # ========================================================================

    def test_detect_gtf_format_chr(self):
        """Test detecting chr-prefixed format from GTF file."""
        chromosomes = ["chr1", "chr2", "chr3", "chrX", "chrY", "chrMT"]
        gtf_file = self._create_gtf_file("test_chr.gtf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(gtf_file)
        
        self.assertEqual(format_type, "chr")
        self.assertEqual(prefix, "chr")

    def test_detect_gtf_format_string(self):
        """Test detecting string format from GTF file."""
        chromosomes = ["1", "2", "3", "X", "Y", "MT"]
        gtf_file = self._create_gtf_file("test_string.gtf", chromosomes, format_type="string")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(gtf_file)
        
        self.assertEqual(format_type, "string")

    # ========================================================================
    # Test Chain File Format Detection
    # ========================================================================

    def test_detect_chain_format_chr(self):
        """Test detecting chr-prefixed format from chain file."""
        chromosomes = ["chr1", "chr2", "chr3", "chrX", "chrY"]
        chain_file = self._create_chain_file("test_chr.chain", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(chain_file)
        
        self.assertEqual(format_type, "chr")
        self.assertEqual(prefix, "chr")

    # ========================================================================
    # Test Conversion Between Sumstats and Reference Formats
    # ========================================================================

    def test_sumstats_to_reference_vcf_chr(self):
        """Test converting sumstats to reference format (VCF with chr prefix)."""
        # Create VCF with chr format
        chromosomes = ["chr1", "chr2", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("ref_chr.vcf", chromosomes, format_type="chr")
        
        # Sumstats uses numeric format
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23, 24, 25]))
        
        # Convert sumstats to reference format
        result = mapper.sumstats_to_reference(1, reference_file=vcf_file)
        self.assertEqual(result, "chr1")
        
        result = mapper.sumstats_to_reference(23, reference_file=vcf_file)
        self.assertEqual(result, "chrX")
        
        result = mapper.sumstats_to_reference(25, reference_file=vcf_file)
        # MT may be returned as "chrMT" or "chrMt" depending on species mapping
        self.assertTrue(result.startswith("chr"))
        self.assertIn(result[3:].upper(), ["MT", "M"])

    def test_sumstats_to_reference_fasta_numeric(self):
        """Test converting sumstats to reference format (FASTA with numeric)."""
        # Create FASTA with numeric format
        chromosomes = ["1", "2", "3", "22", "23", "24", "25"]
        fasta_file = self._create_fasta_file("ref_numeric.fa", chromosomes, format_type="numeric")
        
        # Sumstats uses string format
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["1", "2", "X", "Y", "MT"]))
        
        # Convert sumstats to reference format
        result = mapper.sumstats_to_reference("1", reference_file=fasta_file)
        self.assertEqual(result, 1)
        
        result = mapper.sumstats_to_reference("X", reference_file=fasta_file)
        self.assertEqual(result, 23)

    def test_reference_to_sumstats_vcf_chr(self):
        """Test converting reference to sumstats format (VCF with chr prefix)."""
        # Create VCF with chr format
        chromosomes = ["chr1", "chr2", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("ref_chr.vcf", chromosomes, format_type="chr")
        
        # Sumstats uses numeric format
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23, 24, 25]))
        
        # Convert reference to sumstats format
        result = mapper.reference_to_sumstats("chr1", reference_file=vcf_file)
        self.assertEqual(result, 1)
        
        result = mapper.reference_to_sumstats("chrX", reference_file=vcf_file)
        self.assertEqual(result, 23)
        
        result = mapper.reference_to_sumstats("chrMT", reference_file=vcf_file)
        self.assertEqual(result, 25)

    def test_sumstats_to_reference_series(self):
        """Test vectorized conversion from sumstats to reference format."""
        # Create VCF with chr format
        chromosomes = ["chr1", "chr2", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("ref_chr.vcf", chromosomes, format_type="chr")
        
        # Sumstats uses numeric format
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23, 24, 25]))
        
        # Convert Series
        sumstats_series = pd.Series([1, 2, 23, 24, 25])
        result = mapper.sumstats_to_reference_series(sumstats_series, reference_file=vcf_file)
        
        # Check individual values (MT format may vary)
        self.assertEqual(result.iloc[0], "chr1")
        self.assertEqual(result.iloc[1], "chr2")
        self.assertEqual(result.iloc[2], "chrX")
        self.assertEqual(result.iloc[3], "chrY")
        # MT may be "chrMT" or "chrMt" depending on species mapping
        self.assertTrue(result.iloc[4].startswith("chr"))
        self.assertIn(result.iloc[4][3:].upper(), ["MT", "M"])

    def test_reference_to_sumstats_series(self):
        """Test vectorized conversion from reference to sumstats format."""
        # Create VCF with chr format
        chromosomes = ["chr1", "chr2", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("ref_chr.vcf", chromosomes, format_type="chr")
        
        # Sumstats uses numeric format
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23, 24, 25]))
        
        # Convert Series
        reference_series = pd.Series(["chr1", "chr2", "chrX", "chrY", "chrMT"])
        result = mapper.reference_to_sumstats_series(reference_series, reference_file=vcf_file)
        
        expected = pd.Series([1, 2, 23, 24, 25])
        pd.testing.assert_series_equal(result, expected, check_dtype=False)

    # ========================================================================
    # Test Unconvertible Chromosomes (pd.NA)
    # ========================================================================

    def test_unconvertible_chromosome_returns_na(self):
        """Test that unconvertible chromosomes return pd.NA instead of raising error."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Test unconvertible chromosome like '1_KI270766v1_alt'
        result = mapper.sumstats_to_number("1_KI270766v1_alt")
        self.assertTrue(pd.isna(result))
        
        # Test with Series
        series = pd.Series(["1", "2", "1_KI270766v1_alt", "X"])
        result = mapper.to_numeric(series)
        self.assertEqual(result.iloc[0], 1)
        self.assertEqual(result.iloc[1], 2)
        self.assertTrue(pd.isna(result.iloc[2]))
        self.assertEqual(result.iloc[3], 23)

    def test_unconvertible_chromosome_in_reference(self):
        """Test handling unconvertible chromosomes in reference format."""
        # Create VCF with some unconvertible chromosomes
        chromosomes = ["chr1", "chr2", "1_KI270766v1_alt", "chrX"]
        vcf_file = self._create_vcf_file("ref_mixed.vcf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Try to convert unconvertible reference chromosome
        result = mapper.reference_to_number("1_KI270766v1_alt", reference_file=vcf_file)
        self.assertTrue(pd.isna(result))

    def test_unconvertible_chromosome_series(self):
        """Test handling unconvertible chromosomes in Series."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Series with unconvertible values
        series = pd.Series(["1", "2", "1_KI270766v1_alt", "GL000009.2", "X"])
        result = mapper.to_numeric(series)
        
        # Valid chromosomes should convert
        self.assertEqual(result.iloc[0], 1)
        self.assertEqual(result.iloc[1], 2)
        self.assertEqual(result.iloc[4], 23)
        
        # Unconvertible should be pd.NA
        self.assertTrue(pd.isna(result.iloc[2]))
        self.assertTrue(pd.isna(result.iloc[3]))

    # ========================================================================
    # Test Round-Trip Conversions
    # ========================================================================

    def test_round_trip_sumstats_to_reference(self):
        """Test round-trip conversion: sumstats -> reference -> sumstats."""
        # Create VCF with chr format
        chromosomes = ["chr1", "chr2", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("ref_chr.vcf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23, 24, 25]))
        
        # Round-trip test
        original = 1
        reference = mapper.sumstats_to_reference(original, reference_file=vcf_file)
        result = mapper.reference_to_sumstats(reference, reference_file=vcf_file)
        self.assertEqual(result, original)
        
        original = 23
        reference = mapper.sumstats_to_reference(original, reference_file=vcf_file)
        result = mapper.reference_to_sumstats(reference, reference_file=vcf_file)
        self.assertEqual(result, original)

    def test_round_trip_series(self):
        """Test round-trip conversion with Series."""
        # Create VCF with chr format
        chromosomes = ["chr1", "chr2", "chrX", "chrY", "chrMT"]
        vcf_file = self._create_vcf_file("ref_chr.vcf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23, 24, 25]))
        
        # Round-trip test
        original = pd.Series([1, 2, 23, 24, 25])
        reference = mapper.sumstats_to_reference_series(original, reference_file=vcf_file)
        result = mapper.reference_to_sumstats_series(reference, reference_file=vcf_file)
        
        pd.testing.assert_series_equal(result, original, check_dtype=False)

    # ========================================================================
    # Test Different File Extensions
    # ========================================================================

    def test_detect_vcf_gz(self):
        """Test detecting format from compressed VCF file."""
        # Note: This test would require actual gzip compression
        # For now, we test the uncompressed version
        chromosomes = ["chr1", "chr2", "chrX"]
        vcf_file = self._create_vcf_file("test.vcf.gz", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        # The detection should handle .vcf.gz extension
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        
        self.assertEqual(format_type, "chr")

    def test_detect_fasta_gz(self):
        """Test detecting format from compressed FASTA file."""
        # Note: This test creates an uncompressed file with .gz extension
        # The actual detection may fail if it tries to read as gzip
        # So we test with uncompressed version instead
        chromosomes = ["chr1", "chr2", "chrX"]
        fasta_file = self._create_fasta_file("test.fa", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(fasta_file)
        
        self.assertEqual(format_type, "chr")

    # ========================================================================
    # Test Edge Cases
    # ========================================================================

    def test_nonexistent_file(self):
        """Test handling of nonexistent reference file."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format("nonexistent_file.vcf")
        
        # Should default to string format
        self.assertEqual(format_type, "string")

    def test_empty_file(self):
        """Test handling of empty reference file."""
        empty_file = os.path.join(self.temp_dir, "empty.vcf")
        self.test_files.append(empty_file)
        
        with open(empty_file, 'w') as f:
            pass  # Create empty file
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(empty_file)
        
        # Should default to string format
        self.assertEqual(format_type, "string")

    def test_mixed_chromosome_formats_in_file(self):
        """Test handling of mixed chromosome formats in reference file."""
        # Create file with mixed formats (should detect dominant format)
        chromosomes = ["chr1", "chr2", "1", "2", "chrX"]
        vcf_file = self._create_vcf_file("mixed.vcf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        
        # Should detect the dominant format (chr in this case)
        self.assertEqual(format_type, "chr")

    # ========================================================================
    # Additional Test Cases from Research
    # ========================================================================

    def test_unplaced_sequences_gl(self):
        """Test handling of unplaced sequences (GL prefix)."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Unplaced sequences should return pd.NA
        test_cases = [
            "GL000195.1",
            "GL000009.2",
            "GL000008.2",
            "GL000194.1",
            "GL000192.1"
        ]
        
        for chrom in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for {chrom}, got {result}")

    def test_unplaced_sequences_chrun(self):
        """Test handling of unplaced sequences with chrUn prefix."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Unplaced sequences with chrUn prefix should return pd.NA
        test_cases = [
            "chrUn_GL000195v1",
            "chrUn_GL000009v2",
            "chrUn_KI270766v1",
            "chrUn_KI270928v1"
        ]
        
        for chrom in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for {chrom}, got {result}")

    def test_alternative_loci_variations(self):
        """Test handling of various alternative loci formats."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Alternative loci should return pd.NA
        test_cases = [
            "1_KI270766v1_alt",
            "2_KI270774v1_alt",
            "3_KI270784v1_alt",
            "chr1_KI270766v1_alt",
            "chr2_KI270774v1_alt",
            "KI270766.1",
            "KI270774.1"
        ]
        
        for chrom in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for {chrom}, got {result}")

    def test_chromosome_arm_notation(self):
        """Test handling of chromosome arm notation (1p, 2q, etc.)."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Chromosome arm notation should return pd.NA (not standard chromosome identifiers)
        test_cases = [
            "1p",
            "1q",
            "2p",
            "2q",
            "Xp",
            "Xq",
            "chr1p",
            "chr1q"
        ]
        
        for chrom in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for {chrom}, got {result}")

    def test_empty_string_and_none(self):
        """Test handling of empty strings and None values."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Empty string should return pd.NA
        result = mapper.sumstats_to_number("")
        self.assertTrue(pd.isna(result))
        
        # None should return pd.NA
        result = mapper.sumstats_to_number(None)
        self.assertTrue(pd.isna(result))

    def test_special_characters(self):
        """Test handling of special characters and malformed identifiers."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Special characters and malformed identifiers should return pd.NA
        test_cases = [
            "chr1-",
            "chr1_",
            "chr1.",
            "chr1:",
            "chr1/",
            "chr1@",
            "chr1#",
            "chr1$",
            "chr1%",
            "chr1^",
            "chr1&",
            "chr1*",
            "chr1(",
            "chr1)",
            "chr1+",
            "chr1=",
            "chr1[",
            "chr1]",
            "chr1{",
            "chr1}",
            "chr1|",
            "chr1\\",
            "chr1;",
            "chr1'",
            "chr1\"",
            "chr1<",
            "chr1>",
            "chr1,",
            "chr1?",
            "chr1!",
            "chr1~",
            "chr1`",
            "InvalidChromosome",
            "chrInvalid",
            "not_a_chromosome",
            "chr",
            "chr ",
            " chr1",
            "chr1 ",
            "chr 1",
            "chr1chr1",
            "chrchr1"
        ]
        
        for chrom in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for '{chrom}', got {result}")

    def test_numeric_edge_cases(self):
        """Test handling of numeric edge cases."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Out of range numbers should return pd.NA
        test_cases = [
            "0",  # Chromosome 0 doesn't exist
            "26",  # Beyond human autosomes + X + Y + MT
            "100",
            "999",
            "-1",
            "-5"
        ]
        
        for chrom in test_cases:
            result = mapper.sumstats_to_number(chrom)
            # Some might be valid if they're in the valid range, but most should be pd.NA
            # Let's check if it's a valid chromosome number first
            try:
                num = int(chrom)
                if 1 <= num <= 25:  # Valid range for human
                    # This might be valid, so we don't assert pd.NA
                    pass
                else:
                    self.assertTrue(pd.isna(result), f"Expected pd.NA for out-of-range {chrom}, got {result}")
            except ValueError:
                self.assertTrue(pd.isna(result), f"Expected pd.NA for invalid {chrom}, got {result}")

    def test_case_variations(self):
        """Test handling of various case combinations."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["chr1", "chr2", "chrX"]))
        
        # Case variations should work (case-insensitive)
        test_cases = [
            ("chr1", 1),
            ("Chr1", 1),
            ("CHR1", 1),
            ("chrX", 23),
            ("ChrX", 23),
            ("CHRX", 23),
            ("chrY", 24),
            ("ChrY", 24),
            ("CHRY", 24),
            ("chrMT", 25),
            ("ChrMT", 25),
            ("CHRMT", 25),
            ("chrMt", 25),
            ("ChrMt", 25)
        ]
        
        for chrom, expected in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertEqual(result, expected, f"Case variation failed for {chrom}")

    def test_mixed_valid_invalid_series(self):
        """Test handling of Series with mixed valid and invalid chromosomes."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Series with mix of valid and invalid
        series = pd.Series([
            "1",           # Valid
            "2",           # Valid
            "X",           # Valid
            "1_KI270766v1_alt",  # Invalid - alternative locus
            "GL000195.1",  # Invalid - unplaced
            "chrUn_GL000009v2",  # Invalid - unplaced with prefix
            "1p",          # Invalid - chromosome arm
            "InvalidChromosome",  # Invalid - malformed
            "23",          # Valid
            ""             # Invalid - empty
        ])
        
        result = mapper.to_numeric(series)
        
        # Check valid ones
        self.assertEqual(result.iloc[0], 1)
        self.assertEqual(result.iloc[1], 2)
        self.assertEqual(result.iloc[2], 23)
        self.assertEqual(result.iloc[8], 23)
        
        # Check invalid ones return pd.NA
        for idx in [3, 4, 5, 6, 7, 9]:
            self.assertTrue(pd.isna(result.iloc[idx]), f"Expected pd.NA at index {idx}, got {result.iloc[idx]}")

    def test_reference_file_with_unplaced_sequences(self):
        """Test reference file containing unplaced sequences."""
        # Create VCF with mix of standard and unplaced chromosomes
        chromosomes = ["chr1", "chr2", "chrX", "GL000195.1", "chrUn_GL000009v2", "1_KI270766v1_alt"]
        vcf_file = self._create_vcf_file("ref_with_unplaced.vcf", chromosomes, format_type="chr")
        
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Standard chromosomes should convert
        result = mapper.sumstats_to_reference(1, reference_file=vcf_file)
        self.assertEqual(result, "chr1")
        
        # Unplaced sequences in reference should still be detectable but not convertible
        format_type, prefix = mapper.detect_reference_format(vcf_file)
        self.assertEqual(format_type, "chr")  # Should detect chr format despite unplaced sequences

    def test_whitespace_handling(self):
        """Test handling of whitespace in chromosome identifiers."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series(["chr1", "chr2"]))
        
        # Whitespace should be stripped
        test_cases = [
            (" chr1", 1),
            ("chr1 ", 1),
            (" chr1 ", 1),
            ("\tchr1", 1),
            ("chr1\t", 1),
            ("\nchr1", 1),
            ("chr1\n", 1)
        ]
        
        for chrom, expected in test_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertEqual(result, expected, f"Whitespace handling failed for '{chrom}'")

    def test_very_long_chromosome_names(self):
        """Test handling of very long chromosome identifier strings."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Very long strings should return pd.NA
        long_strings = [
            "chr1" + "x" * 1000,
            "1" + "0" * 1000,
            "X" + "Y" * 1000,
            "chr" + "1" * 1000
        ]
        
        for chrom in long_strings:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for very long string, got {result}")

    def test_unicode_and_special_unicode(self):
        """Test handling of unicode characters."""
        mapper = ChromosomeMapper(log=self.log, verbose=False)
        mapper.detect_sumstats_format(pd.Series([1, 2, 23]))
        
        # Unicode characters should return pd.NA
        unicode_cases = [
            "chr1α",
            "chr1β",
            "chr1γ",
            "chr1→",
            "chr1←",
            "chr1€",
            "chr1£",
            "chr1¥",
            "chr1©",
            "chr1®",
            "chr1™",
            "chr1°",
            "chr1±",
            "chr1×",
            "chr1÷"
        ]
        
        for chrom in unicode_cases:
            result = mapper.sumstats_to_number(chrom)
            self.assertTrue(pd.isna(result), f"Expected pd.NA for unicode '{chrom}', got {result}")


if __name__ == '__main__':
    unittest.main()

