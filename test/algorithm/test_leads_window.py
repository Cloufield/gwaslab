"""Tests for lead window clustering."""

import os
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.algorithm.leads.window import cluster_ids_from_positions
from gwaslab.g_Sumstats import Sumstats
from gwaslab.util.util_in_get_sig import _collect_leads_generic


class TestClusterIdsFromPositions(unittest.TestCase):
    def test_single_chromosome_clusters_by_window(self):
        chrom = np.array([1, 1, 1, 1])
        pos = np.array([1000, 1500, 6000, 6500])
        cluster_id = cluster_ids_from_positions(chrom, pos, window_bp=500)
        np.testing.assert_array_equal(cluster_id, [0, 0, 1, 1])

    def test_new_chromosome_starts_new_cluster(self):
        chrom = np.array([1, 1, 2, 2])
        pos = np.array([1000, 2000, 1000, 2000])
        cluster_id = cluster_ids_from_positions(chrom, pos, window_bp=500_000)
        np.testing.assert_array_equal(cluster_id, [0, 0, 1, 1])

    def test_nan_positions_raise(self):
        with self.assertRaises(ValueError):
            cluster_ids_from_positions(
                np.array([1, 1]),
                np.array([1000.0, np.nan]),
                window_bp=5000,
            )


class TestCollectLeadsGeneric(unittest.TestCase):
    def test_skips_nan_pos_without_error(self):
        df = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": pd.array([1_000_000, pd.NA, 2_000_000], dtype="Int64"),
                "P": [1e-10, 1e-11, 1e-9],
                "__SCALEDP": [1e-10, 1e-11, 1e-9],
                "__ID": [0, 1, 2],
            }
        )
        leads = _collect_leads_generic(
            df, "CHR", "POS", 500, "__SCALEDP", "__ID", maximize=False
        )
        self.assertEqual(leads, [0, 2])

    def test_get_lead_with_nan_pos_in_dirty_like_data(self):
        df = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": [1_000_000, None, 2_000_000],
                "P": [1e-10, 1e-12, 1e-9],
                "MLOG10P": [10.0, 12.0, 9.0],
                "SNPID": ["a", "b", "c"],
                "EA": ["A"] * 3,
                "NEA": ["G"] * 3,
            }
        )
        ss = Sumstats(
            sumstats=df,
            chrom="CHR",
            pos="POS",
            p="P",
            mlog10p="MLOG10P",
            snpid="SNPID",
            ea="EA",
            nea="NEA",
            verbose=False,
        )
        result = ss.get_lead(sig_level=5e-8, windowsizekb=500, verbose=False)
        self.assertIsNotNone(result)
        self.assertFalse(result["POS"].isna().any())
        self.assertEqual(len(result), 2)


if __name__ == "__main__":
    unittest.main()
