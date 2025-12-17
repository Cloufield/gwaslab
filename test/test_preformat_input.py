import os
import sys
import unittest
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.io.io_preformat_input import preformat


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
        out = preformat(sumstats=self.df, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", verbose=False)
        self.assertIn("EAF", out.columns)
        self.assertIn("RAF", out.columns)
        self.assertIn("EXTRA1", out.columns)

    def test_path_trims_unmapped_columns(self):
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        try:
            self.df.to_csv(path, sep="\t", index=False)
            out = preformat(sumstats=path, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", eaf="EAF", verbose=False)
            self.assertIn("EAF", out.columns)
            self.assertIn("SNPID", out.columns)
            self.assertNotIn("EXTRA1", out.columns)
        finally:
            if os.path.exists(path):
                os.remove(path)


if __name__ == "__main__":
    unittest.main()

