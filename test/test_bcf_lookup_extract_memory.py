"""Tests for streaming BCF lookup extraction (memory-safe path)."""
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "test/ref/1kg_eas_hg19.chr7_126253550_128253550.vcf.gz"


@unittest.skipUnless(
    REF.exists() and shutil.which("bcftools"),
    "needs test ref and bcftools",
)
class TestBcfLookupExtractStreaming(unittest.TestCase):
    def setUp(self):
        # Only chr7 positions inside the local test VCF window
        self.targets = pd.DataFrame(
            {"CHR": [7, 7, 7], "POS": [126253551, 126253567, 126253570]}
        )

    def test_extract_writes_lookup_without_concat(self):
        from gwaslab.hm.hm_assign_rsid import _extract_lookup_table_from_vcf_bcf

        with tempfile.NamedTemporaryFile(suffix=".lookup.tsv", delete=False) as tmp:
            out_path = tmp.name

        try:
            path, _ = _extract_lookup_table_from_vcf_bcf(
                vcf_path=str(REF),
                sumstats=self.targets,
                assign_cols=["rsID"],
                out_lookup=out_path,
                extract_threads=2,
                verbose=False,
                log_run_plan=False,
            )
            self.assertEqual(path, out_path)
            self.assertTrue(os.path.getsize(out_path) > 0)
            with open(out_path) as f:
                header = f.readline()
            self.assertIn("CHR", header)
            self.assertIn("rsID", header)
            df = pd.read_csv(out_path, sep="\t")
            self.assertGreater(len(df), 0)
        finally:
            if os.path.exists(out_path):
                os.remove(out_path)

    def test_extract_writes_parquet_lookup(self):
        from gwaslab.hm.hm_assign_rsid import _extract_lookup_table_from_vcf_bcf

        with tempfile.NamedTemporaryFile(suffix=".lookup.parquet", delete=False) as tmp:
            out_path = tmp.name

        try:
            path, _ = _extract_lookup_table_from_vcf_bcf(
                vcf_path=str(REF),
                sumstats=self.targets,
                assign_cols=["rsID"],
                out_lookup=out_path,
                extract_threads=2,
                verbose=False,
                log_run_plan=False,
            )
            self.assertEqual(path, out_path)
            self.assertTrue(os.path.getsize(out_path) > 0)
            df = pd.read_parquet(out_path)
            self.assertIn("CHR", df.columns)
            self.assertTrue("rsID" in df.columns or "ID" in df.columns)
            self.assertGreater(len(df), 0)
        finally:
            if os.path.exists(out_path):
                os.remove(out_path)
            build_path = f"{out_path}.build"
            if os.path.exists(build_path):
                os.remove(build_path)


if __name__ == "__main__":
    unittest.main()
