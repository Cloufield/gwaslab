import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import random
import pandas as pd

from gwaslab.g_Sumstats import Sumstats


def make_sumstats(n=100):
    rng = random.Random(2025)
    rows = []
    for i in range(n):
        chr_ = rng.randint(1, 22)
        pos_ = rng.randint(1_000, 2_000_000)
        pval = max(min(rng.random(), 0.999999), 1e-300)
        snpid = f"{chr_}:{pos_}:A:G"
        rows.append({"CHR": chr_, "POS": pos_, "P": pval, "SNPID": snpid, "EA": "A", "NEA": "G"})
    return pd.DataFrame(rows)


class TestSumstatsChaining(unittest.TestCase):
    def setUp(self):
        self.df = make_sumstats()
        self.gl = Sumstats(sumstats=self.df, chrom="CHR", pos="POS", p="P", ea="EA", nea="NEA", snpid="SNPID", verbose=False)

    def test_fix_methods_chaining(self):
        out = self.gl.fix_id().fix_chr().fix_pos().fix_allele().check_sanity().check_data_consistency().remove_dup().sort_coordinate().sort_column()
        self.assertIsInstance(out, Sumstats)
        self.assertGreater(len(out.data), 0)

    def test_basic_check_chaining(self):
        out = self.gl.basic_check(remove_dup=True, verbose=False).sort_column()
        self.assertIsInstance(out, Sumstats)
        status = out.check_sumstats_qc_status()
        self.assertTrue(status.get("basic_check", {}).get("performed", False))

    def test_harmonize_chaining_ref_free(self):
        out = self.gl.harmonize(basic_check=True, verbose=False).sort_column()
        self.assertIsInstance(out, Sumstats)
        self.assertTrue(out.meta.get("is_harmonised", False))


if __name__ == "__main__":
    unittest.main()
