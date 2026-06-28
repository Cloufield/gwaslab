"""Parquet lookup extract + assign parity with TSV."""
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
class TestLookupParquetAssignParity(unittest.TestCase):
    def setUp(self):
        # Variants inside the chr7 test VCF window (126253550–128253550)
        self.targets = pd.DataFrame(
            {
                "CHR": [7, 7, 7],
                "POS": [126253551, 126253567, 126253570],
            }
        )
        self.sumstats = pd.DataFrame(
            {
                "CHR": [7, 7, 7],
                "POS": [126253551, 126253567, 126253570],
                "EA": ["T", "T", "A"],
                "NEA": ["G", "C", "G"],
            }
        )

    def test_tsv_and_parquet_assign_same_rsid(self):
        from gwaslab.hm.hm_assign_rsid import (
            _assign_from_lookup,
            _extract_lookup_table_from_vcf_bcf,
        )

        tsv_path = tempfile.NamedTemporaryFile(
            suffix=".lookup.tsv", delete=False
        ).name
        parquet_path = tempfile.NamedTemporaryFile(
            suffix=".lookup.parquet", delete=False
        ).name
        try:
            _extract_lookup_table_from_vcf_bcf(
                vcf_path=str(REF),
                sumstats=self.targets,
                assign_cols=["rsID"],
                out_lookup=tsv_path,
                extract_threads=1,
                verbose=False,
                log_run_plan=False,
            )
            _extract_lookup_table_from_vcf_bcf(
                vcf_path=str(REF),
                sumstats=self.targets,
                assign_cols=["rsID"],
                out_lookup=parquet_path,
                extract_threads=1,
                verbose=False,
                log_run_plan=False,
            )

            out_tsv = self.sumstats.copy()
            out_parquet = self.sumstats.copy()

            _assign_from_lookup(
                out_tsv,
                lookup_table=tsv_path,
                assign_cols=("rsID",),
                verbose=False,
                log_run_plan=False,
            )
            _assign_from_lookup(
                out_parquet,
                lookup_table=parquet_path,
                assign_cols=("rsID",),
                verbose=False,
                log_run_plan=False,
            )

            self.assertIn("rsID", out_tsv.columns)
            self.assertIn("rsID", out_parquet.columns)
            self.assertTrue(out_tsv["rsID"].equals(out_parquet["rsID"]))
            self.assertEqual(
                out_tsv["ALLELE_FLIPPED"].sum(),
                out_parquet["ALLELE_FLIPPED"].sum(),
            )
            self.assertTrue(out_tsv["rsID"].notna().all())
        finally:
            for p in (tsv_path, parquet_path, f"{parquet_path}.build"):
                if os.path.exists(p):
                    os.remove(p)


if __name__ == "__main__":
    unittest.main()
