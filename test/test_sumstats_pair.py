import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.g_Sumstats import Sumstats
from gwaslab.g_SumstatsPair import SumstatsPair


def make_sumstats(df, study_name):
    return Sumstats(
        sumstats=df,
        fmt=None,
        tab_fmt="tsv",
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        p="P",
        n=1000,
        study=study_name,
        verbose=False,
    )


class TestSumstatsPair(unittest.TestCase):
    def setUp(self):
        df1 = pd.DataFrame({
            "CHR": [1, 1, 2],
            "POS": [100, 200, 300],
            "EA": ["A", "G", "T"],
            "NEA": ["G", "C", "C"],
            "P": [0.05, 0.001, 0.2],
            "BETA": [0.1, -0.2, 0.3],
            "SE": [0.05, 0.1, 0.2],
        })
        df2 = pd.DataFrame({
            "CHR": [1, 1, 2],
            "POS": [100, 200, 300],
            "EA": ["A", "C", "T"],
            "NEA": ["G", "G", "C"],
            "P": [0.04, 0.01, 0.3],
            "BETA": [0.12, 0.25, -0.35],
            "SE": [0.06, 0.11, 0.25],
        })

        self.s1 = make_sumstats(df1, study_name="StudyA")
        self.s2 = make_sumstats(df2, study_name="StudyB")

    def test_init_and_merge(self):
        sp = SumstatsPair(self.s1, self.s2, verbose=False)
        self.assertEqual(sp.meta["gwaslab"]["group_name"], "StudyA_StudyB")
        self.assertIn("EA", sp.data.columns)
        self.assertIn("NEA", sp.data.columns)
        self.assertTrue(any(c.endswith("_2") for c in sp.data.columns))
        self.assertIn("P_1", sp.data.columns)
        self.assertIn("P_2", sp.data.columns)
        self.assertEqual(len(sp.data), 3)
        self.assertIn("N_1", sp.data.columns)
        self.assertIn("N_2", sp.data.columns)
        self.assertEqual(sp.ns, (1000, 1000))

    def test_filter_value_copy_and_inplace(self):
        sp = SumstatsPair(self.s1, self.s2, verbose=False)
        sp2 = sp.filter_value("P_1 < 0.05", inplace=False)
        self.assertIsInstance(sp2, SumstatsPair)
        self.assertLess(len(sp2.data), len(sp.data))

        sp.filter_value("P_1 < 0.05", inplace=True)
        self.assertLess(len(sp.data), 3)

    def test_compare_af_and_plot_miami_entry(self):
        sp = SumstatsPair(self.s1, self.s2, verbose=False)
        self.assertIn("P_1", sp.data.columns)
        self.assertIn("P_2", sp.data.columns)
        # ensure entry points exist; do not render files during tests
        self.assertTrue(callable(sp.compare_af))
        self.assertTrue(callable(sp.plot_miami))
        # call plot_miami; wrapper does not return fig/log, just ensure no exception
        sp.plot_miami(save=False, verbose=False)

    def test_stacked_mqq(self):
        sp = SumstatsPair(self.s1, self.s2, verbose=False)
        # stacked_mqq should call the plotting function without errors
        sp.stacked_mqq(titles=["StudyA", "StudyB"], vcfs=[None], save=False, verbose=False)


if __name__ == "__main__":
    unittest.main()
