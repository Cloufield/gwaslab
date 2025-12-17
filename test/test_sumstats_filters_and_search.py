import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.g_Sumstats import Sumstats


def make_df():
    rows = [
        {"CHR": 1, "POS": 1000, "P": 5e-8, "EA": "A", "NEA": "G", "SNPID": "1:1000_G_A", "rsID": "rs1000", "EAF": 0.25, "INFO": 0.95},
        {"CHR": 1, "POS": 2000, "P": 0.05,  "EA": "AT", "NEA": "A", "SNPID": "1:2000_A_AT", "rsID": "rs2000", "EAF": 0.10, "INFO": 0.80},
        {"CHR": 2, "POS": 1500, "P": 1e-3, "EA": "C", "NEA": "T", "SNPID": "2:1500_T_C", "rsID": "rs1500", "EAF": 0.45, "INFO": 0.99},
        {"CHR": 2, "POS": 2500, "P": 0.2,  "EA": "G", "NEA": "GT", "SNPID": "2:2500_GT_G", "rsID": "rs2500", "EAF": 0.55, "INFO": 0.60},
    ]
    return pd.DataFrame(rows)


class TestSumstatsFiltersAndSearch(unittest.TestCase):
    def setUp(self):
        self.df = make_df()
        self.gl = Sumstats(sumstats=self.df, chrom="CHR", pos="POS", p="P", ea="EA", nea="NEA", snpid="SNPID", rsid="rsID", verbose=False)
        self.gl.set_build("19", verbose=False)

    def test_filter_region_in_out(self):
        bed_path = os.path.join(os.path.dirname(__file__), "_tmp_regions.bed")
        with open(bed_path, "w") as f:
            f.write("chr1\t900\t2100\n")
        kept = self.gl.filter_region_in(path=bed_path)
        self.assertTrue(all(kept.data["CHR"] == 1))
        self.assertTrue(kept.data["POS"].between(900, 2100).all())
        removed = self.gl.filter_region_out(path=bed_path)
        self.assertTrue(((removed.data["CHR"] != 1) | (~removed.data["POS"].between(900, 2100))).all())
        os.remove(bed_path)

    def test_filter_flanking_by_id(self):
        center = self.gl.data.iloc[0]["SNPID"]
        flank = self.gl.filter_flanking_by_id(snpid=center, windowsizekb=1)
        self.assertGreaterEqual(len(flank.data), 1)
        self.assertIn("CHR", flank.data.columns)
        self.assertIn("POS", flank.data.columns)

    def test_filter_value_query(self):
        kept = self.gl.filter_value(expr="P < 0.01")
        self.assertTrue((kept.data["P"] < 0.01).all())
        kept2 = self.gl.filter_in(gt={"INFO": 0.9})
        self.assertTrue((kept2.data["INFO"] > 0.9).all())
        removed = self.gl.filter_out(lt={"EAF": 0.2})
        self.assertTrue((removed.data["EAF"] >= 0.2).all())

    def test_filter_snp_and_indel(self):
        only_snps = self.gl.filter_snp()
        self.assertTrue(((only_snps.data["EA"].str.len() == 1) & (only_snps.data["NEA"].str.len() == 1)).all())
        only_indels = self.gl.filter_indel()
        self.assertTrue(((only_indels.data["EA"].str.len() != 1) | (only_indels.data["NEA"].str.len() != 1)).all())

    def test_filter_palindromic(self):
        out = self.gl.filter_palindromic(mode="out")
        self.assertGreaterEqual(len(self.gl.data), len(out.data))

    def test_search_variants_by_ids_and_coords(self):
        found = self.gl.search(snplist=["rs1000", "1:2000", [2, 1500], "2:2500:G:GT"])  # rsID, CHR:POS, [CHR,POS], CHR:POS:EA:NEA
        ids = set(found.data.get("rsID", pd.Series(dtype=str)).tolist()) | set(found.data.get("SNPID", pd.Series(dtype=str)).tolist())
        self.assertTrue("rs1000" in ids)
        self.assertTrue(any((found.data["CHR"] == 1) & (found.data["POS"] == 2000)))
        self.assertTrue(any((found.data["CHR"] == 2) & (found.data["POS"] == 1500)))
        self.assertTrue(any((found.data["CHR"] == 2) & (found.data["POS"] == 2500)))


if __name__ == "__main__":
    unittest.main()
