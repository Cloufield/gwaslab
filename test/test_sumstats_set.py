import os
import sys
import unittest
import tempfile
import shutil

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.g_Sumstats import Sumstats
from gwaslab.g_SumstatsSet import SumstatsSet


class TestSumstatsSet(unittest.TestCase):
    def setUp(self):
        rows1 = [
            {"CHR": 1, "POS": 100, "EA": "A", "NEA": "G", "BETA": 0.10, "SE": 0.02, "P": 1e-6, "SNPID": "1:100:A:G"},
            {"CHR": 2, "POS": 200, "EA": "C", "NEA": "T", "BETA": -0.05, "SE": 0.03, "P": 1e-4, "SNPID": "2:200:C:T"},
        ]
        rows2 = [
            {"CHR": 1, "POS": 100, "EA": "A", "NEA": "G", "BETA": 0.08, "SE": 0.02, "P": 5e-8, "SNPID": "1:100:A:G"},
            {"CHR": 3, "POS": 300, "EA": "G", "NEA": "C", "BETA": 0.02, "SE": 0.04, "P": 0.2,  "SNPID": "3:300:G:C"},
        ]

        df1 = pd.DataFrame(rows1)
        df2 = pd.DataFrame(rows2)

        self.gl1 = Sumstats(sumstats=df1, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", beta="BETA", se="SE", p="P", verbose=False)
        self.gl2 = Sumstats(sumstats=df2, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", beta="BETA", se="SE", p="P", verbose=False)

        self.sumstats_dic = {"StudyA": self.gl1, "StudyB": self.gl2}
        self.variant_set = [[1, 100], "2:200:C:T", "1:100:A:G"]

    def test_init_and_extract(self):
        sset = SumstatsSet(sumstats_dic=self.sumstats_dic, variant_set=self.variant_set, verbose=False)
        self.assertIn("STUDY", sset.data.columns)
        self.assertTrue({"CHR", "POS", "EA", "NEA", "BETA", "SE", "P"}.issubset(set(sset.data.columns)))
        self.assertEqual(len(sset.data), 3)
        self.assertEqual(sset.meta["gwaslab"]["genome_build"], "99")
        self.assertEqual(sset._build, "99")

    def test_inherited_basic_check(self):
        sset = SumstatsSet(sumstats_dic=self.sumstats_dic, variant_set=self.variant_set, verbose=False)
        sset.basic_check(normalize=False, remove_dup=False, verbose=False)
        self.assertIn("STATUS", sset.data.columns)

    def test_inherited_filter_value(self):
        sset = SumstatsSet(sumstats_dic=self.sumstats_dic, variant_set=self.variant_set, verbose=False)
        filtered = sset.filter_value("CHR == 1", inplace=False)
        self.assertIsInstance(filtered, SumstatsSet)
        self.assertEqual(len(filtered.data), 2)

    def test_plot_effect_runs(self):
        sset = SumstatsSet(sumstats_dic=self.sumstats_dic, variant_set=self.variant_set, verbose=False)
        sset.plot_effect(verbose=False, save=None, y="STUDY",
                         fig_kwargs={"figsize": (4, 3)}, group=["CHR","POS","STUDY"], y_sort=["CHR","POS","STUDY"])


class TestSumstatsSetGlobLoading(unittest.TestCase):
    """Tests for SumstatsSet glob pattern loading feature."""

    def setUp(self):
        """Create temporary directory with test sumstats files."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test data
        self.data1 = pd.DataFrame([
            {"CHR": 1, "POS": 100, "EA": "A", "NEA": "G", "BETA": 0.10, "SE": 0.02, "P": 1e-6, "SNPID": "rs123"},
            {"CHR": 2, "POS": 200, "EA": "C", "NEA": "T", "BETA": -0.05, "SE": 0.03, "P": 1e-4, "SNPID": "rs456"},
        ])
        self.data2 = pd.DataFrame([
            {"CHR": 1, "POS": 100, "EA": "A", "NEA": "G", "BETA": 0.08, "SE": 0.02, "P": 5e-8, "SNPID": "rs123"},
            {"CHR": 3, "POS": 300, "EA": "G", "NEA": "C", "BETA": 0.02, "SE": 0.04, "P": 0.2, "SNPID": "rs789"},
        ])
        self.data3 = pd.DataFrame([
            {"CHR": 1, "POS": 100, "EA": "A", "NEA": "G", "BETA": 0.12, "SE": 0.025, "P": 1e-5, "SNPID": "rs123"},
            {"CHR": 4, "POS": 400, "EA": "T", "NEA": "A", "BETA": 0.03, "SE": 0.05, "P": 0.5, "SNPID": "rs999"},
        ])
        
        # Write test files
        self.file1 = os.path.join(self.temp_dir, "study_EUR.txt")
        self.file2 = os.path.join(self.temp_dir, "study_EAS.txt")
        self.file3 = os.path.join(self.temp_dir, "study_AFR.txt")
        
        self.data1.to_csv(self.file1, sep="\t", index=False)
        self.data2.to_csv(self.file2, sep="\t", index=False)
        self.data3.to_csv(self.file3, sep="\t", index=False)

    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir)

    def test_glob_pattern_loading(self):
        """Test loading multiple files via glob pattern."""
        pattern = os.path.join(self.temp_dir, "study_*.txt")
        variant_set = ["rs123"]
        
        sset = SumstatsSet(
            pattern,
            variant_set=variant_set,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Should have loaded 3 studies
        studies = sset.data["STUDY"].unique()
        self.assertEqual(len(studies), 3)
        
        # Study names should be derived from filenames
        self.assertIn("study_EUR", studies)
        self.assertIn("study_EAS", studies)
        self.assertIn("study_AFR", studies)
        
        # Should have extracted rs123 from all 3 studies
        self.assertEqual(len(sset.data), 3)

    def test_glob_pattern_with_question_mark(self):
        """Test glob pattern with ? wildcard."""
        # Create additional files with single character variation
        file_a = os.path.join(self.temp_dir, "trait_A.txt")
        file_b = os.path.join(self.temp_dir, "trait_B.txt")
        self.data1.to_csv(file_a, sep="\t", index=False)
        self.data2.to_csv(file_b, sep="\t", index=False)
        
        pattern = os.path.join(self.temp_dir, "trait_?.txt")
        variant_set = ["rs123"]
        
        sset = SumstatsSet(
            pattern,
            variant_set=variant_set,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        studies = sset.data["STUDY"].unique()
        self.assertEqual(len(studies), 2)
        self.assertIn("trait_A", studies)
        self.assertIn("trait_B", studies)

    def test_glob_pattern_no_match_raises_error(self):
        """Test that FileNotFoundError is raised when no files match."""
        pattern = os.path.join(self.temp_dir, "nonexistent_*.txt")
        
        with self.assertRaises(FileNotFoundError) as context:
            SumstatsSet(
                pattern,
                variant_set=["rs123"],
                verbose=False
            )
        
        self.assertIn("No files match pattern", str(context.exception))

    def test_glob_single_file_match(self):
        """Test glob pattern that matches only one file."""
        pattern = os.path.join(self.temp_dir, "study_EUR.txt")
        variant_set = ["rs123", "rs456"]
        
        sset = SumstatsSet(
            pattern,
            variant_set=variant_set,
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        studies = sset.data["STUDY"].unique()
        self.assertEqual(len(studies), 1)
        self.assertIn("study_EUR", studies)
        self.assertEqual(len(sset.data), 2)

    def test_glob_pattern_without_variant_set(self):
        """Test loading all data when variant_set is None."""
        pattern = os.path.join(self.temp_dir, "study_*.txt")
        
        sset = SumstatsSet(
            pattern,
            variant_set=None,  # No variant filtering
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Should have loaded 3 studies
        studies = sset.data["STUDY"].unique()
        self.assertEqual(len(studies), 3)
        
        # Should have all variants from all studies (2 + 2 + 2 = 6)
        self.assertEqual(len(sset.data), 6)
        
        # Check STUDY column exists and has correct values
        self.assertIn("STUDY", sset.data.columns)
        self.assertIn("study_EUR", studies)
        self.assertIn("study_EAS", studies)
        self.assertIn("study_AFR", studies)

    def test_dict_input_without_variant_set(self):
        """Test loading all data from dict when variant_set is None."""
        rows1 = [
            {"CHR": 1, "POS": 100, "EA": "A", "NEA": "G", "BETA": 0.10, "SE": 0.02, "P": 1e-6, "SNPID": "rs1"},
        ]
        rows2 = [
            {"CHR": 2, "POS": 200, "EA": "C", "NEA": "T", "BETA": 0.05, "SE": 0.03, "P": 1e-4, "SNPID": "rs2"},
        ]
        
        gl1 = Sumstats(pd.DataFrame(rows1), chrom="CHR", pos="POS", ea="EA", nea="NEA", 
                       snpid="SNPID", beta="BETA", se="SE", p="P", verbose=False)
        gl2 = Sumstats(pd.DataFrame(rows2), chrom="CHR", pos="POS", ea="EA", nea="NEA",
                       snpid="SNPID", beta="BETA", se="SE", p="P", verbose=False)
        
        sset = SumstatsSet(
            {"study1": gl1, "study2": gl2},
            variant_set=None,
            verbose=False
        )
        
        # Should have all data from both studies
        self.assertEqual(len(sset.data), 2)
        self.assertIn("STUDY", sset.data.columns)
        studies = sset.data["STUDY"].unique()
        self.assertIn("study1", studies)
        self.assertIn("study2", studies)


class TestDeriveStudyName(unittest.TestCase):
    """Tests for the _derive_study_name static method."""

    def test_simple_txt_extension(self):
        """Test removing .txt extension."""
        result = SumstatsSet._derive_study_name("./data/study_EUR.txt")
        self.assertEqual(result, "study_EUR")

    def test_gz_extension(self):
        """Test removing .gz extension."""
        result = SumstatsSet._derive_study_name("./data/study_EUR.txt.gz")
        self.assertEqual(result, "study_EUR")

    def test_sumstats_gz_extension(self):
        """Test removing .sumstats.gz extension."""
        result = SumstatsSet._derive_study_name("/path/to/trait_1.sumstats.gz")
        self.assertEqual(result, "trait_1")

    def test_tsv_extension(self):
        """Test removing .tsv extension."""
        result = SumstatsSet._derive_study_name("analysis.tsv")
        self.assertEqual(result, "analysis")

    def test_csv_extension(self):
        """Test removing .csv extension."""
        result = SumstatsSet._derive_study_name("results.csv")
        self.assertEqual(result, "results")

    def test_bgz_extension(self):
        """Test removing .bgz extension."""
        result = SumstatsSet._derive_study_name("data.txt.bgz")
        self.assertEqual(result, "data")

    def test_no_extension(self):
        """Test file without recognized extension."""
        result = SumstatsSet._derive_study_name("./study_name")
        self.assertEqual(result, "study_name")

    def test_complex_path(self):
        """Test with complex directory path."""
        result = SumstatsSet._derive_study_name("/home/user/gwas/results/EUR_meta.sumstats.gz")
        self.assertEqual(result, "EUR_meta")


if __name__ == "__main__":
    unittest.main()
