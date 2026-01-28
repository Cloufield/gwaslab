"""
Test suite for get_novel with known loci (DataFrame or local file).
Uses known=DataFrame to avoid GWAS Catalog API calls. Covers:
- known as DataFrame, only_novel, output_known, if_get_lead
- ValueError when neither known nor efo provided
- Novel vs known by distance (windowsizekb_for_novel)
- show_child_traits keyword (no-op when using known)
"""

import os
import sys
import tempfile
import unittest
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.g_Sumstats import Sumstats


def _make_sumstats_two_leads():
    """Sumstats with two lead-like variants on chr1: 1e6 and 5e6."""
    df = pd.DataFrame({
        "CHR": [1, 1],
        "POS": [1000000, 5000000],
        "P": [1e-10, 1e-9],
        "MLOG10P": [10.0, 9.0],
        "SNPID": ["1:1000000_A_G", "1:5000000_A_G"],
        "EA": ["A", "A"],
        "NEA": ["G", "G"],
    })
    ss = Sumstats(
        sumstats=df,
        chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
        snpid="SNPID", ea="EA", nea="NEA",
        verbose=False,
    )
    ss.build = "38"
    return ss


def _make_known_at(chr, pos):
    """Known loci DataFrame with CHR, POS."""
    return pd.DataFrame({"CHR": [chr], "POS": [pos]})


class TestGetNovelWithKnownDataFrame(unittest.TestCase):
    """Test get_novel using known=DataFrame (no GWAS Catalog)."""

    def test_get_novel_returns_dataframe_with_novel_column(self):
        """get_novel(known=df) returns DataFrame with NOVEL column."""
        ss = _make_sumstats_two_leads()
        known = _make_known_at(1, 1000500)  # near first lead
        result = ss.get_novel(known=known, build="38", verbose=False)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("NOVEL", result.columns)
        self.assertEqual(len(result), 2)

    def test_get_novel_novel_vs_known_by_distance(self):
        """Variant within windowsizekb_for_novel of known is not novel; far one is novel."""
        ss = _make_sumstats_two_leads()
        # known at 1:1000500 -> 500 bp from 1:1000000, ~4 Mb from 1:5000000
        known = _make_known_at(1, 1000500)
        result = ss.get_novel(
            known=known,
            build="38",
            windowsizekb_for_novel=1000,  # 1000 kb
            verbose=False,
        )
        self.assertIsNotNone(result)
        row_1m = result[result["POS"] == 1000000]
        row_5m = result[result["POS"] == 5000000]
        self.assertEqual(len(row_1m), 1)
        self.assertEqual(len(row_5m), 1)
        self.assertFalse(row_1m.iloc[0]["NOVEL"], "variant at 1Mb should be 'known' (within 1000kb)")
        self.assertTrue(row_5m.iloc[0]["NOVEL"], "variant at 5Mb should be novel")

    def test_get_novel_only_novel_true(self):
        """only_novel=True returns only rows with NOVEL=True."""
        ss = _make_sumstats_two_leads()
        known = _make_known_at(1, 1000500)
        result = ss.get_novel(
            known=known,
            build="38",
            only_novel=True,
            windowsizekb_for_novel=1000,
            verbose=False,
        )
        self.assertIsNotNone(result)
        self.assertIn("NOVEL", result.columns)
        self.assertTrue(result["NOVEL"].all())
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["POS"], 5000000)

    def test_get_novel_output_known_true(self):
        """output_known=True returns (variants_df, known_df) tuple."""
        ss = _make_sumstats_two_leads()
        known = _make_known_at(1, 1000500)
        result = ss.get_novel(
            known=known,
            build="38",
            output_known=True,
            verbose=False,
        )
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        variants_df, known_df = result
        self.assertIsInstance(variants_df, pd.DataFrame)
        self.assertIsInstance(known_df, pd.DataFrame)
        self.assertIn("NOVEL", variants_df.columns)
        self.assertIn("CHR", known_df.columns)
        self.assertIn("POS", known_df.columns)

    def test_get_novel_if_get_lead_false(self):
        """if_get_lead=False uses input sumstats as leads (no lead extraction)."""
        ss = _make_sumstats_two_leads()
        known = _make_known_at(1, 1000500)
        result = ss.get_novel(
            known=known,
            build="38",
            if_get_lead=False,
            verbose=False,
        )
        self.assertIsNotNone(result)
        self.assertIn("NOVEL", result.columns)
        self.assertEqual(len(result), 2)

    def test_get_novel_requires_known_or_efo(self):
        """Raises ValueError when neither known nor efo is provided."""
        ss = _make_sumstats_two_leads()
        with self.assertRaises(ValueError) as ctx:
            ss.get_novel(build="38", verbose=False)
        self.assertIn("known", str(ctx.exception).lower())
        self.assertIn("efo", str(ctx.exception).lower())

    def test_get_novel_show_child_traits_kwarg_accepted(self):
        """show_child_traits can be passed when using known (no efo); should not error."""
        ss = _make_sumstats_two_leads()
        known = _make_known_at(1, 1000500)
        result = ss.get_novel(
            known=known,
            build="38",
            show_child_traits=False,
            verbose=False,
        )
        self.assertIsNotNone(result)
        self.assertIn("NOVEL", result.columns)

    def test_get_novel_has_distance_and_location_columns(self):
        """Result includes DISTANCE_TO_KNOWN and LOCATION_OF_KNOWN when applicable."""
        ss = _make_sumstats_two_leads()
        known = _make_known_at(1, 1000500)
        result = ss.get_novel(known=known, build="38", verbose=False)
        self.assertIn("DISTANCE_TO_KNOWN", result.columns)
        self.assertIn("LOCATION_OF_KNOWN", result.columns)

    def test_get_novel_known_from_path(self):
        """get_novel(known=path) reads CHR/POS from file."""
        ss = _make_sumstats_two_leads()
        known_df = _make_known_at(1, 1000500)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            known_df.to_csv(f.name, sep="\t", index=False)
            path = f.name
        try:
            result = ss.get_novel(known=path, build="38", verbose=False)
            self.assertIsNotNone(result)
            self.assertIn("NOVEL", result.columns)
            self.assertEqual(len(result), 2)
        finally:
            os.unlink(path)


class TestGetNovelEdgeCases(unittest.TestCase):
    """Edge cases: single variant, all novel, all known."""

    def test_single_lead_far_from_known_is_novel(self):
        """Single lead far from known -> one row, NOVEL=True."""
        df = pd.DataFrame({
            "CHR": [1], "POS": [10000000],
            "P": [1e-10], "MLOG10P": [10.0],
            "SNPID": ["1:10000000_A_G"], "EA": ["A"], "NEA": ["G"],
        })
        ss = Sumstats(
            sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
            snpid="SNPID", ea="EA", nea="NEA", verbose=False,
        )
        ss.build = "38"
        known = _make_known_at(1, 1000000)
        result = ss.get_novel(known=known, build="38", verbose=False)
        self.assertEqual(len(result), 1)
        self.assertTrue(result.iloc[0]["NOVEL"])

    def test_lead_same_pos_as_known_not_novel(self):
        """Lead at same position as known -> not novel."""
        df = pd.DataFrame({
            "CHR": [1], "POS": [1000000],
            "P": [1e-10], "MLOG10P": [10.0],
            "SNPID": ["1:1000000_A_G"], "EA": ["A"], "NEA": ["G"],
        })
        ss = Sumstats(
            sumstats=df, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P",
            snpid="SNPID", ea="EA", nea="NEA", verbose=False,
        )
        ss.build = "38"
        known = _make_known_at(1, 1000000)
        result = ss.get_novel(known=known, build="38", verbose=False)
        self.assertEqual(len(result), 1)
        self.assertFalse(result.iloc[0]["NOVEL"])


if __name__ == "__main__":
    unittest.main()
