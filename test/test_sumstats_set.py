import os
import sys
import unittest

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


if __name__ == "__main__":
    unittest.main()
