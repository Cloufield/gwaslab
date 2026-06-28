"""Full-scale reference benchmark (skipped unless GWASLAB_RUN_FULL_BENCHMARK=1)."""
import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class TestReferenceRunFullBenchmark(unittest.TestCase):
    @unittest.skipUnless(
        os.environ.get("GWASLAB_RUN_FULL_BENCHMARK") == "1",
        "Set GWASLAB_RUN_FULL_BENCHMARK=1 to run",
    )
    def test_full_bbj_harmonize_manifest_exists(self):
        manifest = ROOT / "benchmarks/reference_run_calibration/manifests/yunye_local.json"
        self.assertTrue(manifest.is_file())
        import json
        with open(manifest) as f:
            cfg = json.load(f)
        self.assertTrue(Path(cfg["sumstats"]["path"]).exists())
        self.assertTrue(Path(cfg["refs"]["ref_seq"]).exists())
        self.assertTrue(Path(cfg["refs"]["ref_infer"]).exists())
        self.assertTrue(Path(cfg["refs"]["ref_rsid_vcf"]).exists())


if __name__ == "__main__":
    unittest.main()
