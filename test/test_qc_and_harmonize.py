import os
import sys
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
from gwaslab.g_Sumstats import Sumstats


RAW_DIR = os.path.join(os.path.dirname(__file__), "raw")


class TestBasicCheck(unittest.TestCase):
    def test_basic_check_on_dirty_sumstats(self):
        path = os.path.join(RAW_DIR, "dirty_sumstats.tsv")
        gl = Sumstats(sumstats=path, fmt=None, tab_fmt="tsv", snpid="SNPID", chrom="CHR", pos="POS", ea="EA", nea="NEA", p="P", eaf="EAF", verbose=False)
        gl.basic_check(remove_dup=True, verbose=False)
        status = gl.check_sumstats_qc_status()
        self.assertIn("basic_check", status)
        self.assertTrue(status["basic_check"].get("performed", False))
        self.assertGreater(len(gl.data), 0)


class TestRemoveDup(unittest.TestCase):
    def test_remove_dup_on_duplicate(self):
        path = os.path.join(RAW_DIR, "duplicate.tsv")
        gl = Sumstats(
            sumstats=path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="ALT",
            nea="REF",
            p="P",
            snpid="SNP",
            rsid="rsID",
            neaf="Frq",
            verbose=False,
        )
        before = len(gl.data)
        gl.remove_dup(mode="dsdrdc", keep="first", verbose=False)
        after = len(gl.data)
        self.assertLess(after, before)


class TestHarmonization(unittest.TestCase):
    def test_harmonize_on_to_harmonize(self):
        path = os.path.join(RAW_DIR, "to_harmonize.tsv")
        gl = Sumstats(sumstats=path, tab_fmt="tsv", chrom="CHR", pos="POS", ea="EA", nea="NEA", p="P", snpid="SNPID", verbose=False)
        # run basic_check within harmonize (ref-free path)
        gl.harmonize(basic_check=True, verbose=False)
        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertGreater(len(gl.data), 0)


class TestQCOutputMatches(unittest.TestCase):
    def test_clean_output_matches_correct(self):
        raw_path = os.path.join(os.path.dirname(__file__), "raw", "dirty_sumstats.tsv")
        correct_path = os.path.join(os.path.dirname(__file__), "correct", "clean_sumstats.tsv")
        out_dir = os.path.join(os.path.dirname(__file__), "output")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "clean_sumstats.tsv")

        gl = Sumstats(sumstats=raw_path, fmt="gwaslab", other=["NOTE"], verbose=False)
        gl.basic_check(remove=True, remove_dup=True, verbose=False)
        gl.data.to_csv(out_path, sep="\t", index=None)

        df_out = pd.read_csv(out_path, sep="\t")
        df_corr = pd.read_csv(correct_path, sep="\t")

        common_cols = [c for c in df_corr.columns if c in df_out.columns]
        df_out = df_out[common_cols]
        df_corr = df_corr[common_cols]

        sort_cols = [c for c in ["CHR", "POS", "SNPID"] if c in common_cols]
        if sort_cols:
            df_out = df_out.sort_values(by=sort_cols).reset_index(drop=True)
            df_corr = df_corr.sort_values(by=sort_cols).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_out, df_corr, check_dtype=False)


class TestHarmonizeWorkflowWithReferences(unittest.TestCase):
    def test_harmonize_bbj_chr7_with_ref_seq_rsid_and_infer(self):
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_seq_path = os.path.join(ref_dir, "chr7.fasta.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=True
        )

        gl.harmonize(
            basic_check=True,
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            threads=4,
            remove=False,
            verbose=True
        )

        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertIn("rsID", gl.data.columns)
        self.assertIn("EA", gl.data.columns)
        self.assertIn("NEA", gl.data.columns)
        self.assertGreater(len(gl.data), 0)


class TestFixIDWorkflow(unittest.TestCase):
    def test_fix_id_produces_expected_output(self):
        raw_path = os.path.join(os.path.dirname(__file__), "raw", "fixid.tsv")
        correct_path = os.path.join(os.path.dirname(__file__), "correct", "fixed_id.tsv")
        out_dir = os.path.join(os.path.dirname(__file__), "output")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "fixed_id.tsv")

        gl = Sumstats(sumstats=raw_path, fmt="gwaslab", verbose=False)
        gl.fix_id(fixchrpos=True, fixeanea=True, fixprefix=True, verbose=False)
        gl.data.to_csv(out_path, sep="\t", index=None)

        df_out = pd.read_csv(out_path, sep="\t")
        df_corr = pd.read_csv(correct_path, sep="\t")

        common_cols = [c for c in df_corr.columns if c in df_out.columns]
        df_out = df_out[common_cols]
        df_corr = df_corr[common_cols]

        sort_cols = [c for c in ["CHR", "POS", "SNPID"] if c in common_cols]
        if sort_cols:
            df_out = df_out.sort_values(by=sort_cols).reset_index(drop=True)
            df_corr = df_corr.sort_values(by=sort_cols).reset_index(drop=True)

        pd.testing.assert_frame_equal(df_out, df_corr, check_dtype=False)

    def test_harmonize_bbj_chr7_sweep_mode(self):
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_seq_path = os.path.join(ref_dir, "chr7.fasta.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")

        gl = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=True
        )

        gl.harmonize(
            basic_check=True,
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            threads=4,
            remove=False,
            verbose=True,
            sweep_mode=True,
        )

        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertIn("rsID", gl.data.columns)
        self.assertIn("EA", gl.data.columns)
        self.assertIn("NEA", gl.data.columns)
        self.assertGreater(len(gl.data), 0)


if __name__ == "__main__":
    unittest.main()
