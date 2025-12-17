import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
import numpy as np
from gwaslab.viz.viz_plot_compare_effect import compare_effect
from gwaslab.g_Sumstats import Sumstats


def make_sumstats_beta(n=120, seed=1357):
    rng = pd.Series(range(n)).sample(frac=1.0, random_state=seed)
    rows = []
    for i in range(n):
        chr_ = (i % 22) + 1
        pos_ = 1_000_000 + i * 100
        snpid = f"{chr_}:{pos_}_A_G"
        p = max(min((i + 1) / (n + 2), 0.999999), 1e-300)
        beta = ((i % 7) - 3) * 0.02
        se = 0.05
        rows.append({
            "SNPID": snpid,
            "CHR": chr_,
            "POS": pos_,
            "EA": "A",
            "NEA": "G",
            "P": p,
            "BETA": beta,
            "SE": se,
            "EAF": 0.3 + (i % 10) * 0.01,
        })
    return pd.DataFrame(rows)


class TestCompareEffect(unittest.TestCase):
    def setUp(self):
        df1 = make_sumstats_beta(n=150, seed=10)
        df2 = df1.copy()
        df2["BETA"] = df2["BETA"] * 1.05
        df2["SE"] = df2["SE"]
        self.gl1 = Sumstats(sumstats=df1, chrom="CHR", pos="POS", p="P", snpid="SNPID", ea="EA", nea="NEA", beta="BETA", se="SE", eaf="EAF", verbose=False)
        self.gl2 = Sumstats(sumstats=df2, chrom="CHR", pos="POS", p="P", snpid="SNPID", ea="EA", nea="NEA", beta="BETA", se="SE", eaf="EAF", verbose=False)

    def test_compare_effect_beta_basic(self):
        label = ["Study1", "Study2", "Both", "Not-significant"]
        result, fig, log = compare_effect(
            path1=self.gl1,
            path2=self.gl2,
            mode="beta",
            label=label,
            sig_level=1,
            verbose=False,
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)
        self.assertIn("EFFECT_1", result.columns)
        self.assertIn("EFFECT_2_aligned", result.columns)
        self.assertIn("indicator", result.columns)

    def test_compare_effect_scaled_mlog10p(self):
        df1 = self.gl1.data.copy()
        df2 = self.gl2.data.copy()
        df1["MLOG10P"] = -np.log10(df1["P"].clip(lower=1e-300))
        df2["MLOG10P"] = -np.log10(df2["P"].clip(lower=1e-300))

        label = ["Study1", "Study2", "Both", "Not-significant"]
        gl1s = Sumstats(sumstats=df1, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P", snpid="SNPID", ea="EA", nea="NEA", beta="BETA", se="SE", eaf="EAF", verbose=False)
        gl2s = Sumstats(sumstats=df2, chrom="CHR", pos="POS", p="P", mlog10p="MLOG10P", snpid="SNPID", ea="EA", nea="NEA", beta="BETA", se="SE", eaf="EAF", verbose=False)
        result, fig, log = compare_effect(
            path1=gl1s,
            path2=gl2s,
            mode="beta",
            label=label,
            sig_level=1,
            verbose=False,
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("MLOG10P_1", result.columns)
        self.assertIn("MLOG10P_2", result.columns)

    def test_compare_effect_with_maf_filter(self):
        label = ["Study1", "Study2", "Both", "Not-significant"]
        result, fig, log = compare_effect(
            path1=self.gl1,
            path2=self.gl2,
            mode="beta",
            label=label,
            eaf=["EAF", "EAF"],
            maf_level=0.05,
            sig_level=1,
            verbose=False,
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreaterEqual(len(result), 0)


if __name__ == "__main__":
    unittest.main()
