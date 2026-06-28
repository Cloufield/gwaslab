import os
import sys
import tempfile
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import gwaslab as gl
from gwaslab.bd.bd_sumstats_formats import detect_sumstats_format
from gwaslab.io.io_preformat_input import _preformat


PLINK2_DATA = """#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP
1\t1000\trs1\tG\tA\tA\tADD\t1000\t0.1\t0.01\t10.0\t0.001
1\t2000\trs2\tC\tT\tT\tADD\t1000\t0.15\t0.015\t10.0\t0.0005"""

PLINK2_WITH_META = """## PLINK2 GLM output
#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP
1\t1000\trs1\tG\tA\tA\tADD\t1000\t0.1\t0.01\t10.0\t0.001"""

VCF_DATA = """##fileformat=VCFv4.3
##FORMAT=<ID=ES,Number=1,Type=Float,Description="Effect Size">
##FORMAT=<ID=SE,Number=1,Type=Float,Description="Standard Error">
##FORMAT=<ID=LP,Number=1,Type=Float,Description="-log10 p-value">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t1000\trs1\tG\tA\t.\t.\t.\tES:SE:LP\t0.1:0.01:3.0
1\t2000\trs2\tC\tT\t.\t.\t.\tES:SE:LP\t0.15:0.015:3.3"""


class TestSumstatsFormatDetection(unittest.TestCase):
    def test_plink2_not_detected_as_vcf(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".glm.linear", delete=False) as f:
            f.write(PLINK2_DATA)
            path = f.name
        try:
            result = detect_sumstats_format(path)
            self.assertNotEqual(result["best_format"], "vcf")
            self.assertIn(result["best_format"], ("plink2", "plink2_linear"))
        finally:
            os.unlink(path)

    def test_vcf_still_detected(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(VCF_DATA)
            path = f.name
        try:
            result = detect_sumstats_format(path)
            self.assertEqual(result["best_format"], "vcf")
        finally:
            os.unlink(path)

    def test_sumstats_load_plink2_from_path(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".glm.linear", delete=False) as f:
            f.write(PLINK2_DATA)
            path = f.name
        try:
            ss = gl.Sumstats(path, fmt="plink2", verbose=False)
            self.assertGreater(len(ss.data), 0)
            self.assertIn("CHR", ss.data.columns)
            self.assertIn("POS", ss.data.columns)
            self.assertIn("BETA", ss.data.columns)
            self.assertEqual(len(ss.data), 2)
        finally:
            os.unlink(path)

    def test_plink2_with_leading_meta_lines(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".glm.linear", delete=False) as f:
            f.write(PLINK2_WITH_META)
            path = f.name
        try:
            result = _preformat(sumstats=path, fmt="plink2", verbose=False)
            self.assertEqual(len(result), 1)
            self.assertIn("CHR", result.columns)
            self.assertIn("BETA", result.columns)
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
