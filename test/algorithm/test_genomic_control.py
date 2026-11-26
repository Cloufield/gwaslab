import unittest

import numpy as np
import pandas as pd
import scipy.stats as stats

import gwaslab as gl
from gwaslab.algorithm.core.genomic_control import (
    lambda_gc_from_chisq,
    lambda_gc_from_mlog10p,
    lambda_gc_from_p,
    lambda_gc_from_z,
)


class TestGenomicControlQuantiles(unittest.TestCase):
    def setUp(self):
        # Deterministic probability grid aligned with np.quantile indices.
        u = np.linspace(1e-6, 1.0 - 1e-6, 10001)
        self.inflation = 1.25
        self.chisq = self.inflation * stats.chi2.ppf(u, df=1)
        self.z = np.sqrt(self.chisq)
        self.p = stats.chi2.sf(self.chisq, df=1)
        self.mlog10p = -stats.chi2.logsf(self.chisq, df=1) / np.log(10.0)

    def test_lambda_gc_matches_inflation_across_quantiles(self):
        for q in (0.25, 0.5, 0.9):
            with self.subTest(quantile=q):
                self.assertAlmostEqual(
                    lambda_gc_from_chisq(self.chisq, quantile=q),
                    self.inflation,
                    places=3,
                )
                self.assertAlmostEqual(
                    lambda_gc_from_z(self.z, quantile=q),
                    self.inflation,
                    places=3,
                )
                self.assertAlmostEqual(
                    lambda_gc_from_p(self.p, quantile=q),
                    self.inflation,
                    places=3,
                )
                self.assertAlmostEqual(
                    lambda_gc_from_mlog10p(self.mlog10p, quantile=q),
                    self.inflation,
                    places=3,
                )

    def test_sumstats_get_gc_non_median_level(self):
        df = pd.DataFrame(
            {
                "CHR": 1,
                "POS": np.arange(1, len(self.p) + 1),
                "P": self.p,
                "Z": self.z,
                "CHISQ": self.chisq,
                "MLOG10P": self.mlog10p,
            }
        )
        ss = gl.Sumstats(
            df,
            chrom="CHR",
            pos="POS",
            p="P",
            z="Z",
            chisq="CHISQ",
            mlog10p="MLOG10P",
            verbose=False,
        )
        for mode in ("P", "MLOG10P", "Z", "CHISQ"):
            with self.subTest(mode=mode):
                lam = ss.get_gc(mode=mode, level=0.9, verbose=False)
                self.assertAlmostEqual(lam, self.inflation, places=3)

    def test_empty_and_nan_inputs_return_nan(self):
        arr = np.array([np.nan, np.inf])
        self.assertTrue(np.isnan(lambda_gc_from_chisq(arr, quantile=0.9)))
        self.assertTrue(np.isnan(lambda_gc_from_z(arr, quantile=0.9)))
        self.assertTrue(np.isnan(lambda_gc_from_p(arr, quantile=0.9)))
        self.assertTrue(np.isnan(lambda_gc_from_mlog10p(arr, quantile=0.9)))


if __name__ == "__main__":
    unittest.main()
