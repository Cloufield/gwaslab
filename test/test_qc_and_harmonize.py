"""
Test file for QC and harmonization functions.

To see logs during test execution:
    python test/test_qc_and_harmonize.py -v
    python -m unittest test.test_qc_and_harmonize -v -s
    pytest test/test_qc_and_harmonize.py -v -s

The -s flag (or --capture=no in pytest) shows stdout/stderr, which includes gwaslab logs.
"""

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
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)
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
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)
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
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        gl.harmonize(
            basic_check=True,
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            threads=6,
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
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        gl.harmonize(
            basic_check=True,
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True,
            sweep_mode=True,
        )

        self.assertTrue(gl.meta.get("is_harmonised", False))
        self.assertIn("rsID", gl.data.columns)
        self.assertIn("EA", gl.data.columns)
        self.assertIn("NEA", gl.data.columns)
        self.assertGreater(len(gl.data), 0)


class TestCheckAF(unittest.TestCase):
    def test_check_af_with_reference_vcf(self):
        """Test check_af function to calculate DAF (difference in allele frequency)"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
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
            verbose=False
        )
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Run check_af
        gl.check_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that DAF column was created
        self.assertIn("DAF", gl.data.columns)
        
        # Check that DAF values are numeric (can be NaN for variants not found in VCF)
        daf_values = gl.data["DAF"].dropna()
        if len(daf_values) > 0:
            # DAF should be between -1 and 1 (difference between two frequencies)
            self.assertTrue(all(daf_values >= -1) and all(daf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)


class TestInferAF(unittest.TestCase):
    def test_infer_af_with_reference_vcf(self):
        """Test infer_af function to infer EAF from reference VCF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
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
            verbose=False
        )
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Run infer_af
        gl.infer_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that EAF column still exists
        self.assertIn("EAF", gl.data.columns)
        
        # Check that EAF values are numeric and between 0 and 1
        eaf_values = gl.data["EAF"].dropna()
        if len(eaf_values) > 0:
            self.assertTrue(all(eaf_values >= 0) and all(eaf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)
        
        # Check that at least some EAF values were inferred (if variants match the VCF region)
        # Note: infer_af only processes variants with good status codes and within the VCF region
        inferred_eaf_count = gl.data["EAF"].notna().sum()
        self.assertGreater(inferred_eaf_count, 0)


class TestConsistencyBetweenMethods(unittest.TestCase):
    """Test consistency between old and new methods (assign_rsid vs assign_rsid2, infer_strand vs infer_strand2)"""
    
    def test_assign_rsid_consistency(self):
        """Test that assign_rsid and assign_rsid2 produce consistent results"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl1.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl1.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl1.data = gl1.data[chr7_mask].reset_index(drop=True)
        # Use only every 5th variant for faster testing (1/5 of variants)
        gl1.data = gl1.data.iloc[::5].reset_index(drop=True)
        print(f"[test_assign_rsid_consistency] Using {len(gl1.data)} variants (every 5th variant, 1/5 of total)")
        gl1.assign_rsid(ref_rsid_vcf=ref_rsid_vcf_path, threads=6, verbose=False)
        result1 = gl1.data.copy()
        rsid1_count = result1["rsID"].notna().sum()
        
        # Test with assign_rsid2 (new method) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl2.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl2.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl2.data = gl2.data[chr7_mask].reset_index(drop=True)
        # Use only every 5th variant for faster testing (same as gl1, 1/5 of variants)
        gl2.data = gl2.data.iloc[::5].reset_index(drop=True)
        gl2.assign_rsid2(vcf_path=ref_rsid_vcf_path, threads=6, verbose=False)
        result2 = gl2.data.copy()
        rsid2_count = result2["rsID"].notna().sum()
        
        # Compare rsID assignments
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Check that both methods assigned rsIDs
        self.assertGreater(rsid1_count, 0, "assign_rsid should assign some rsIDs")
        self.assertGreater(rsid2_count, 0, "assign_rsid2 should assign some rsIDs")
        
        # Compare rsID assignments where both have values
        # Merge on CHR:POS to compare rsIDs
        merged = result1[["CHR", "POS", "rsID"]].merge(
            result2[["CHR", "POS", "rsID"]],
            on=["CHR", "POS"],
            suffixes=("_old", "_new"),
            how="inner"
        )
        
        # Check consistency: where both have rsIDs, they should match
        both_have_rsid = merged["rsID_old"].notna() & merged["rsID_new"].notna()
        if both_have_rsid.sum() > 0:
            matching = (merged.loc[both_have_rsid, "rsID_old"] == merged.loc[both_have_rsid, "rsID_new"]).sum()
            match_rate = matching / both_have_rsid.sum()
            print(f"\n[test_assign_rsid_consistency] rsID match rate: {match_rate:.2%} ({matching}/{both_have_rsid.sum()} variants)")
            # Allow for some differences due to different processing, but should be mostly consistent
            self.assertGreater(match_rate, 0.8, 
                             f"rsID assignments should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_assign_rsid_consistency] No variants with rsIDs in both methods to compare")
        
        # Check that the number of assigned rsIDs is similar (within 20%)
        if rsid1_count > 0 and rsid2_count > 0:
            ratio = min(rsid1_count, rsid2_count) / max(rsid1_count, rsid2_count)
            print(f"[test_assign_rsid_consistency] rsID count ratio: {ratio:.2%} (old: {rsid1_count}, new: {rsid2_count})")
            self.assertGreater(ratio, 0.99, 
                             f"Number of assigned rsIDs should be similar (ratio: {ratio:.2%})")
        else:
            print(f"[test_assign_rsid_consistency] rsID counts - old: {rsid1_count}, new: {rsid2_count}")
    
    def test_infer_strand_consistency(self):
        """Test that infer_strand (via harmonize) and infer_strand2 produce consistent results"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_seq_path = os.path.join(ref_dir, "chr7.fasta.gz")
        ref_rsid_vcf_path = os.path.join(ref_dir, "b157_2564.vcf.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl1.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl1.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl1.data = gl1.data[chr7_mask].reset_index(drop=True)
        # Use only every 5th variant for faster testing (1/5 of variants)
        gl1.data = gl1.data.iloc[::5].reset_index(drop=True)
        print(f"[test_infer_strand_consistency] Using {len(gl1.data)} variants (every 5th variant, 1/5 of total)")
        
        # Test with infer_strand via harmonize (old method, non-sweep mode)
        gl1.harmonize(
            basic_check=False,  # Already done above
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            ref_maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True,
            sweep_mode=False,  # Use old method
            infer_strand_kwargs={"daf_tolerance": 0.20}  # Pass daf_tolerance through kwargs
        )
        result1 = gl1.data.copy()
        status1 = result1["STATUS"].copy()
        
        # Test with infer_strand2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check first to normalize data
        gl2.basic_check(verbose=False)
        
        # File name suggests it's already chr7 only, so use all data
        # If needed, filter to chromosome 7 - handle both string and int CHR formats
        if "chr7" not in sumstats_path.lower():
            chr_col = gl2.data["CHR"]
            chr7_mask = (
                (chr_col == 7) | 
                (chr_col == "7") | 
                (chr_col.astype(str).str.strip() == "7")
            )
            gl2.data = gl2.data[chr7_mask].reset_index(drop=True)
        # Use only every 5th variant for faster testing (same as gl1, 1/5 of variants)
        gl2.data = gl2.data.iloc[::5].reset_index(drop=True)
        
        gl2.harmonize(
            basic_check=False,  # Already done above
            ref_seq=ref_seq_path,
            ref_rsid_vcf=ref_rsid_vcf_path,
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            maf_threshold=0.40,
            ref_maf_threshold=0.40,
            threads=6,
            remove=False,
            verbose=True,
            sweep_mode=True,  # Use new method
            infer_strand_kwargs={"daf_tolerance": 0.20}  # Pass daf_tolerance through kwargs
        )
        result2 = gl2.data.copy()
        status2 = result2["STATUS"].copy()
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare STATUS codes - the 7th digit should be consistent for strand inference
        # Extract 7th digit from STATUS (strand inference status)
        def get_strand_status(status_series):
            """Extract 7th digit from STATUS codes"""
            return status_series.astype(str).str[-1] if len(status_series) > 0 else pd.Series(dtype=str)
        
        strand_status1 = get_strand_status(status1)
        strand_status2 = get_strand_status(status2)
        
        # Compare strand status codes where both have valid status
        valid_mask = (strand_status1 != "nan") & (strand_status2 != "nan")
        if valid_mask.sum() > 0:
            matching = (strand_status1[valid_mask] == strand_status2[valid_mask]).sum()
            match_rate = matching / valid_mask.sum()
            print(f"\n[test_infer_strand_consistency] Strand STATUS match rate: {match_rate:.2%} ({matching}/{valid_mask.sum()} variants)")
            # Show distribution of status codes for debugging
            status1_dist = strand_status1[valid_mask].value_counts().sort_index()
            status2_dist = strand_status2[valid_mask].value_counts().sort_index()
            print(f"[test_infer_strand_consistency] Old method STATUS distribution: {dict(status1_dist)}")
            print(f"[test_infer_strand_consistency] New method STATUS distribution: {dict(status2_dist)}")
            # Allow for some differences, but should be mostly consistent
            self.assertGreater(match_rate, 0.99, 
                             f"Strand inference STATUS codes should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_infer_strand_consistency] No variants with valid STATUS in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_infer_strand_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


class TestCheckAF2(unittest.TestCase):
    def test_check_af2_with_reference_vcf(self):
        """Test check_af2 function (sweep mode) to calculate DAF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
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
            verbose=False
        )
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Run check_af2 (sweep mode)
        gl.check_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that DAF column was created
        self.assertIn("DAF", gl.data.columns)
        
        # Check that DAF values are numeric (can be NaN for variants not found in VCF)
        daf_values = gl.data["DAF"].dropna()
        if len(daf_values) > 0:
            # DAF should be between -1 and 1 (difference between two frequencies)
            self.assertTrue(all(daf_values >= -1) and all(daf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)


class TestInferAF2(unittest.TestCase):
    def test_infer_af2_with_reference_vcf(self):
        """Test infer_af2 function (sweep mode) to infer EAF from reference VCF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
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
            verbose=False
        )
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Store original EAF count
        original_eaf_count = gl.data["EAF"].notna().sum()

        # Run infer_af2 (sweep mode)
        gl.infer_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that EAF column still exists
        self.assertIn("EAF", gl.data.columns)
        
        # Check that EAF values are numeric and between 0 and 1
        eaf_values = gl.data["EAF"].dropna()
        if len(eaf_values) > 0:
            self.assertTrue(all(eaf_values >= 0) and all(eaf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)
        
        # Check that at least some EAF values were inferred (if variants match the VCF region)
        inferred_eaf_count = gl.data["EAF"].notna().sum()
        self.assertGreater(inferred_eaf_count, 0)


class TestInferEAFFromMAF2(unittest.TestCase):
    def test_infer_eaf_from_maf2_with_reference_vcf(self):
        """Test infer_eaf_from_maf2 function (sweep mode) to infer EAF from MAF"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
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
            verbose=False
        )
        # Use only 1/5 of variants for faster testing
        gl.data = gl.data.iloc[::5].reset_index(drop=True)

        # First harmonize to ensure proper status codes
        gl.harmonize(basic_check=True, verbose=False)

        # Check if MAF column exists, if not create it from EAF
        if "MAF" not in gl.data.columns or gl.data["MAF"].isna().all():
            # Create MAF from EAF if available
            if "EAF" in gl.data.columns:
                gl.data["MAF"] = gl.data["EAF"].apply(lambda x: min(x, 1-x) if pd.notna(x) else pd.NA)

        # Store original EAF count
        original_eaf_count = gl.data["EAF"].notna().sum()

        # Run infer_eaf_from_maf2 (sweep mode)
        gl.infer_eaf_from_maf2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )

        # Check that EAF column still exists
        self.assertIn("EAF", gl.data.columns)
        
        # Check that EAF values are numeric and between 0 and 1
        eaf_values = gl.data["EAF"].dropna()
        if len(eaf_values) > 0:
            self.assertTrue(all(eaf_values >= 0) and all(eaf_values <= 1))
        
        # Check that we have some data
        self.assertGreater(len(gl.data), 0)
        
        # Check that at least some EAF values were inferred (if variants match the VCF region and have MAF)
        inferred_eaf_count = gl.data["EAF"].notna().sum()
        self.assertGreater(inferred_eaf_count, 0)


class TestConsistencyCheckAF(unittest.TestCase):
    """Test consistency between check_af (ver1) and check_af2 (ver2)"""
    
    def test_check_af_consistency(self):
        """Test that check_af and check_af2 produce consistent DAF values"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl1.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl1.data = gl1.data[
            (gl1.data["CHR"] == 7) & 
            (gl1.data["POS"] >= 126253550) & 
            (gl1.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 5th variant for faster testing (1/5 of variants)
        gl1.data = gl1.data.iloc[::5].reset_index(drop=True)
        print(f"[test_check_af_consistency] Using {len(gl1.data)} variants (every 5th variant, 1/5 of total in chr7:126253550-128253550)")
        
        # Test with check_af (old method, normal mode)
        gl1.check_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result1 = gl1.data.copy()
        daf1 = result1["DAF"].copy()
        
        # Test with check_af2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl2.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl2.data = gl2.data[
            (gl2.data["CHR"] == 7) & 
            (gl2.data["POS"] >= 126253550) & 
            (gl2.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 5th variant for faster testing (same as gl1, 1/5 of variants)
        gl2.data = gl2.data.iloc[::5].reset_index(drop=True)
        
        gl2.check_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result2 = gl2.data.copy()
        daf2 = result2["DAF"].copy()
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare DAF values where both have values
        both_have_daf = daf1.notna() & daf2.notna()
        if both_have_daf.sum() > 0:
            daf1_values = daf1[both_have_daf]
            daf2_values = daf2[both_have_daf]
            
            # Calculate differences
            daf_diff = (daf1_values - daf2_values).abs()
            max_diff = daf_diff.max()
            mean_diff = daf_diff.mean()
            
            # Check that differences are small (within tolerance)
            # Allow for small numerical differences due to different processing methods
            tolerance = 0.01  # 1% tolerance
            within_tolerance = (daf_diff <= tolerance).sum()
            match_rate = within_tolerance / both_have_daf.sum()
            
            print(f"\n[test_check_af_consistency] DAF match rate (within {tolerance}): {match_rate:.2%} ({within_tolerance}/{both_have_daf.sum()} variants)")
            print(f"[test_check_af_consistency] Max DAF difference: {max_diff:.6f}, Mean DAF difference: {mean_diff:.6f}")
            
            # Should be mostly consistent
            self.assertGreater(match_rate, 0.95, 
                             f"DAF values should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_check_af_consistency] No variants with DAF in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_check_af_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


class TestConsistencyInferAF(unittest.TestCase):
    """Test consistency between infer_af (ver1) and infer_af2 (ver2)"""
    
    def test_infer_af_consistency(self):
        """Test that infer_af and infer_af2 produce consistent EAF values"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl1.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl1.data = gl1.data[
            (gl1.data["CHR"] == 7) & 
            (gl1.data["POS"] >= 126253550) & 
            (gl1.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 5th variant for faster testing (1/5 of variants)
        gl1.data = gl1.data.iloc[::5].reset_index(drop=True)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl1.data["EAF"] = pd.NA
        
        print(f"[test_infer_af_consistency] Using {len(gl1.data)} variants (every 5th variant, 1/5 of total in chr7:126253550-128253550)")
        
        # Test with infer_af (old method, normal mode)
        gl1.infer_af(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result1 = gl1.data.copy()
        eaf1 = result1["EAF"].copy()
        
        # Test with infer_af2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl2.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl2.data = gl2.data[
            (gl2.data["CHR"] == 7) & 
            (gl2.data["POS"] >= 126253550) & 
            (gl2.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 5th variant for faster testing (same as gl1, 1/5 of variants)
        gl2.data = gl2.data.iloc[::5].reset_index(drop=True)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl2.data["EAF"] = pd.NA
        
        gl2.infer_af2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        
        result2 = gl2.data.copy()
        eaf2 = result2["EAF"].copy()
        
        # ALLELE_FLIPPED is an internal temporary column and should be dropped in output
        if "ALLELE_FLIPPED" in result2.columns:
            print(f"[test_infer_af_consistency] WARNING: ALLELE_FLIPPED should be dropped but found in output")
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare EAF values where both have values
        both_have_eaf = eaf1.notna() & eaf2.notna()
        if both_have_eaf.sum() > 0:
            eaf1_values = eaf1[both_have_eaf]
            eaf2_values = eaf2[both_have_eaf]
            
            # Calculate differences
            eaf_diff = (eaf1_values - eaf2_values).abs()
            max_diff = eaf_diff.max()
            mean_diff = eaf_diff.mean()
            
            # Check if differences are close to 1 (suggesting a flip issue)
            close_to_one = (eaf_diff > 0.99) & (eaf_diff < 1.01)
            print(f"\n[test_infer_af_consistency] Differences close to 1 (flip issue?): {close_to_one.sum()} / {both_have_eaf.sum()} ({100*close_to_one.sum()/both_have_eaf.sum():.1f}%)")
            
            # ALLELE_FLIPPED is an internal temporary column and should be dropped in output
            if "ALLELE_FLIPPED" in result2.columns:
                print(f"[test_infer_af_consistency] WARNING: ALLELE_FLIPPED should be dropped but found in output")
            
            # Print out all variants that are different
            mismatched = eaf_diff > 0.01
            if mismatched.sum() > 0:
                print(f"\n[test_infer_af_consistency] Printing all {mismatched.sum()} variants with differences > 0.01:")
                print(f"{'CHR':<6} {'POS':<12} {'EA':<6} {'NEA':<6} {'EAF_old':<10} {'EAF_new':<10} {'Diff':<10}")
                print("-" * 70)
                
                # Get the actual indices from result2 for mismatched variants
                # Note: both_have_eaf is a boolean Series aligned with result2.index
                # mismatched is a boolean Series aligned with eaf_diff (which is also aligned with result2.index[both_have_eaf])
                mismatched_bool = pd.Series(False, index=result2.index)
                mismatched_bool.loc[result2.index[both_have_eaf][mismatched]] = True
                mismatched_indices = result2.index[mismatched_bool]
                
                mismatched_df = result2.loc[mismatched_indices, ["CHR", "POS", "EA", "NEA"]].copy()
                mismatched_df["EAF_old"] = eaf1_values.loc[mismatched_indices].values
                mismatched_df["EAF_new"] = eaf2_values.loc[mismatched_indices].values
                mismatched_df["Diff"] = eaf_diff[mismatched].values
                
                # ALLELE_FLIPPED is an internal temporary column and should be dropped
                # So we don't include it in the output comparison
                mismatched_df["FLIPPED"] = "N/A (internal column)"
                
                # Sort by difference (largest first)
                mismatched_df = mismatched_df.sort_values("Diff", ascending=False)
                
                # Print all mismatched variants
                for idx, row in mismatched_df.iterrows():
                    print(f"{int(row['CHR']):<6} {int(row['POS']):<12} {row['EA']:<6} {row['NEA']:<6} "
                          f"{row['EAF_old']:<10.6f} {row['EAF_new']:<10.6f} {row['Diff']:<10.6f}")
            else:
                print(f"\n[test_infer_af_consistency] No variants with differences > 0.01")
            
            # Check that differences are small (within tolerance)
            # Allow for small numerical differences due to different processing methods
            tolerance = 0.01  # 1% tolerance
            within_tolerance = (eaf_diff <= tolerance).sum()
            match_rate = within_tolerance / both_have_eaf.sum()
            
            print(f"\n[test_infer_af_consistency] EAF match rate (within {tolerance}): {match_rate:.2%} ({within_tolerance}/{both_have_eaf.sum()} variants)")
            print(f"[test_infer_af_consistency] Max EAF difference: {max_diff:.6f}, Mean EAF difference: {mean_diff:.6f}")
            
            # Should be mostly consistent - require 99% match rate
            self.assertGreater(match_rate, 0.99, 
                             f"EAF values should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_infer_af_consistency] No variants with EAF in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_infer_af_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


class TestConsistencyInferEAFFromMAF(unittest.TestCase):
    """Test consistency between infer_eaf_from_maf (ver1) and infer_eaf_from_maf2 (ver2)"""
    
    def test_infer_eaf_from_maf_consistency(self):
        """Test that infer_eaf_from_maf and infer_eaf_from_maf2 produce consistent EAF values"""
        ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        sumstats_path = os.path.join(ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        ref_infer_vcf_path = os.path.join(ref_dir, "1kg_eas_hg19.chr7_126253550_128253550.vcf.gz")
        
        # Load sumstats - load fresh each time
        gl1 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl1.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl1.data = gl1.data[
            (gl1.data["CHR"] == 7) & 
            (gl1.data["POS"] >= 126253550) & 
            (gl1.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 5th variant for faster testing (1/5 of variants)
        gl1.data = gl1.data.iloc[::5].reset_index(drop=True)
        
        # Check if MAF column exists, if not create it from EAF
        if "MAF" not in gl1.data.columns or gl1.data["MAF"].isna().all():
            # Create MAF from EAF if available (before dropping EAF)
            if "EAF" in gl1.data.columns:
                gl1.data["MAF"] = gl1.data["EAF"].apply(lambda x: min(x, 1-x) if pd.notna(x) else pd.NA)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl1.data["EAF"] = pd.NA
        
        print(f"[test_infer_eaf_from_maf_consistency] Using {len(gl1.data)} variants (every 5th variant, 1/5 of total in chr7:126253550-128253550)")
        
        # Test with infer_eaf_from_maf (old method, normal mode)
        gl1.infer_eaf_from_maf(
            ref_infer=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result1 = gl1.data.copy()
        eaf1 = result1["EAF"].copy()
        
        # Test with infer_eaf_from_maf2 (new method, sweep mode) - load fresh again
        gl2 = Sumstats(
            sumstats=sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            eaf="EAF",
            verbose=False
        )
        # Run basic_check and harmonize first to normalize data
        gl2.harmonize(basic_check=True, verbose=False)
        
        # Filter to only variants in chr7:126253550-128253550 (matching reference VCF region)
        gl2.data = gl2.data[
            (gl2.data["CHR"] == 7) & 
            (gl2.data["POS"] >= 126253550) & 
            (gl2.data["POS"] <= 128253550)
        ].reset_index(drop=True)
        
        # Use only every 5th variant for faster testing (same as gl1, 1/5 of variants)
        gl2.data = gl2.data.iloc[::5].reset_index(drop=True)
        
        # Check if MAF column exists, if not create it from EAF
        if "MAF" not in gl2.data.columns or gl2.data["MAF"].isna().all():
            # Create MAF from EAF if available (before dropping EAF)
            if "EAF" in gl2.data.columns:
                gl2.data["MAF"] = gl2.data["EAF"].apply(lambda x: min(x, 1-x) if pd.notna(x) else pd.NA)
        
        # Drop original EAF to ensure both methods actually infer EAF
        gl2.data["EAF"] = pd.NA
        
        gl2.infer_eaf_from_maf2(
            vcf_path=ref_infer_vcf_path,
            ref_alt_freq="AF",
            verbose=False
        )
        result2 = gl2.data.copy()
        eaf2 = result2["EAF"].copy()
        
        # Sort both dataframes for comparison
        result1 = result1.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        result2 = result2.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
        
        # Compare EAF values where both have values
        both_have_eaf = eaf1.notna() & eaf2.notna()
        if both_have_eaf.sum() > 0:
            eaf1_values = eaf1[both_have_eaf]
            eaf2_values = eaf2[both_have_eaf]
            
            # Calculate differences
            eaf_diff = (eaf1_values - eaf2_values).abs()
            max_diff = eaf_diff.max()
            mean_diff = eaf_diff.mean()
            
            # Check that differences are small (within tolerance)
            # Allow for small numerical differences due to different processing methods
            tolerance = 0.01  # 1% tolerance
            within_tolerance = (eaf_diff <= tolerance).sum()
            match_rate = within_tolerance / both_have_eaf.sum()
            
            print(f"\n[test_infer_eaf_from_maf_consistency] EAF match rate (within {tolerance}): {match_rate:.2%} ({within_tolerance}/{both_have_eaf.sum()} variants)")
            print(f"[test_infer_eaf_from_maf_consistency] Max EAF difference: {max_diff:.6f}, Mean EAF difference: {mean_diff:.6f}")
            
            # Should be mostly consistent - require 99% match rate
            self.assertGreater(match_rate, 0.99, 
                             f"EAF values should be consistent (match rate: {match_rate:.2%})")
        else:
            print(f"\n[test_infer_eaf_from_maf_consistency] No variants with EAF in both methods to compare")
        
        # Check that both methods processed similar numbers of variants
        print(f"[test_infer_eaf_from_maf_consistency] Variant counts - old: {len(result1)}, new: {len(result2)}")
        self.assertEqual(len(result1), len(result2), "Both methods should process the same number of variants")


if __name__ == "__main__":
    # Run with verbosity and show stdout/stderr to see logs
    # Usage: python test/test_qc_and_harmonize.py -v
    # Or: python -m unittest test.test_qc_and_harmonize -v -s
    import sys
    # Add -s flag to show stdout (logs) during test execution
    if "-s" not in sys.argv and "--buffer" not in sys.argv:
        sys.argv.append("-s")
    unittest.main(verbosity=2)
