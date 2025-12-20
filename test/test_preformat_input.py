import os
import sys
import unittest
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
import polars as pl

from gwaslab.io.io_preformat_input import _preformat
from gwaslab.io.io_preformat_input_polars import preformatp


class TestPreformatInput(unittest.TestCase):
    def setUp(self):
        rng = pd.Series(range(10))
        rows = []
        for i in rng:
            chr_ = 1
            pos_ = 1000 + i
            snpid = f"{chr_}:{pos_}_A_G"
            rows.append({"CHR": chr_, "POS": pos_, "EA": "A", "NEA": "G", "EAF": 0.2, "RAF": 0.3, "SNPID": snpid, "EXTRA1": i})
        self.df = pd.DataFrame(rows)

    def test_dataframe_preserves_columns(self):
        out = _preformat(sumstats=self.df, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", verbose=False)
        self.assertIn("EAF", out.columns)
        self.assertIn("RAF", out.columns)
        self.assertIn("EXTRA1", out.columns)

    def test_path_trims_unmapped_columns(self):
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        try:
            self.df.to_csv(path, sep="\t", index=False)
            out = _preformat(sumstats=path, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", eaf="EAF", verbose=False)
            self.assertIn("EAF", out.columns)
            self.assertIn("SNPID", out.columns)
            self.assertNotIn("EXTRA1", out.columns)
        finally:
            if os.path.exists(path):
                os.remove(path)


class TestPreformatInputPolars(unittest.TestCase):
    def setUp(self):
        self.raw_dir = os.path.join(os.path.dirname(__file__), "raw")
        self.dirty_sumstats_path = os.path.join(self.raw_dir, "dirty_sumstats.tsv")

    def test_load_dirty_sumstats_from_path(self):
        """Test loading dirty_sumstats.tsv using preformatp"""
        result = preformatp(
            sumstats=self.dirty_sumstats_path,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            n="N",
            ncase="N_CASE",
            ncontrol="N_CONTROL",
            other=["NOTE"],
            verbose=False
        )
        
        # Check that result is a polars DataFrame
        self.assertIsInstance(result, pl.DataFrame)
        
        # Check that it has data
        self.assertGreater(result.height, 0)
        
        # Check that expected columns are present
        expected_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P", "N", "N_CASE", "N_CONTROL", "NOTE"]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Column {col} not found in result")
        
        # Check that NOTE column is preserved
        self.assertIn("NOTE", result.columns)

    def test_load_dirty_sumstats_with_include(self):
        """Test loading with include parameter"""
        result = preformatp(
            sumstats=self.dirty_sumstats_path,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            include=["SNPID", "CHR", "POS", "EA", "NEA", "P"],
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertGreater(result.height, 0)
        
        # Should only have included columns
        included_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "P"]
        for col in included_cols:
            self.assertIn(col, result.columns)
        
        # Should not have excluded columns like BETA
        self.assertNotIn("BETA", result.columns)

    def test_load_dirty_sumstats_with_exclude(self):
        """Test loading with exclude parameter"""
        result = preformatp(
            sumstats=self.dirty_sumstats_path,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            exclude=["BETA", "SE"],
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertGreater(result.height, 0)
        
        # Should have basic columns
        self.assertIn("SNPID", result.columns)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        
        # Should not have excluded columns
        self.assertNotIn("BETA", result.columns)
        self.assertNotIn("SE", result.columns)

    def test_load_dirty_sumstats_from_pandas_dataframe(self):
        """Test loading from pandas DataFrame"""
        # First load as pandas
        df_pd = pd.read_csv(self.dirty_sumstats_path, sep="\t")
        
        result = preformatp(
            sumstats=df_pd,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertEqual(result.height, df_pd.shape[0])
        self.assertIn("SNPID", result.columns)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)

    def test_load_dirty_sumstats_creates_snpid_if_missing(self):
        """Test that SNPID is created if both rsID and SNPID are missing"""
        # Create a test file without SNPID
        import tempfile
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        try:
            # Create a minimal test file
            test_data = """CHR\tPOS\tEA\tNEA\tP
1\t100\tA\tG\t0.001
1\t200\tT\tC\t0.002"""
            with open(path, 'w') as f:
                f.write(test_data)
            
            result = preformatp(
                sumstats=path,
                chrom="CHR",
                pos="POS",
                ea="EA",
                nea="NEA",
                p="P",
                verbose=False
            )
            
            self.assertIsInstance(result, pl.DataFrame)
            self.assertIn("SNPID", result.columns)
            self.assertGreater(result.height, 0)
        finally:
            if os.path.exists(path):
                os.remove(path)


if __name__ == "__main__":
    unittest.main()

