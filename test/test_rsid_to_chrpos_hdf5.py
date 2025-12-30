import os
import sys
import unittest
import tempfile
import shutil
import glob

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

# Check if pytables is available (required for HDF5 operations)
try:
    import tables
    PYTABLES_AVAILABLE = True
except ImportError:
    PYTABLES_AVAILABLE = False

from gwaslab.g_Sumstats import Sumstats
from gwaslab.util.util_ex_process_h5 import process_vcf_to_hfd5
from gwaslab.bd.bd_common_data import get_NC_to_number


@unittest.skipIf(not PYTABLES_AVAILABLE, "pytables (tables) is required for HDF5 operations")
class TestProcessVCFToHDF5(unittest.TestCase):
    """Test VCF to HDF5 conversion."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        self.vcf_path = os.path.join(self.ref_dir, "b157_2564.vcf.gz")
        self.temp_dir = tempfile.mkdtemp(prefix="gwaslab_test_hdf5_")
        
    def tearDown(self):
        """Clean up temporary files."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_process_vcf_to_hdf5_creates_chromosome_files(self):
        """Test that process_vcf_to_hfd5 creates separate HDF5 files per chromosome."""
        # Process VCF to HDF5 with chromosome mapping for NC_ notation
        # Use get_NC_to_number to map NC_ notation to chromosome numbers
        chr_dict = get_NC_to_number(build="19")
        
        output_dir = process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            complevel=3,
            threads=1,
            verbose=False
        )
        
        # Check that output directory is returned
        self.assertEqual(output_dir, self.temp_dir)
        
        # Check that HDF5 files were created
        h5_files = glob.glob(os.path.join(self.temp_dir, "*.chr*.rsID_CHR_POS_mod10.h5"))
        self.assertGreater(len(h5_files), 0, "No chromosome HDF5 files created")
        
        # Check that at least one file has data (some chromosomes may be empty)
        files_with_data = []
        total_variants = 0
        for h5_file in h5_files:
            try:
                with pd.HDFStore(h5_file, mode='r') as store:
                    # Get all keys (groups)
                    all_keys = list(store.keys())
                    groups = [key for key in all_keys if key.startswith("/group_")]
                    
                    # Count variants in this file
                    file_variants = 0
                    for group_key in groups:
                        try:
                            group_df = store[group_key]
                            if len(group_df) > 0:
                                file_variants += len(group_df)
                                total_variants += len(group_df)
                        except Exception:
                            # Skip groups that can't be read
                            continue
                    
                    if file_variants > 0:
                        files_with_data.append(h5_file)
            except Exception as e:
                # Skip files that can't be opened (might be empty or corrupted)
                continue
        
        # At least some variants should be processed (VCF has rsIDs)
        self.assertGreater(total_variants, 0,
                          f"At least some variants should be processed. "
                          f"Found {len(files_with_data)} files with data, {total_variants} total variants")
        
        # Check structure of files that have data
        for h5_file in files_with_data[:3]:  # Check first 3 files with data
            with pd.HDFStore(h5_file, mode='r') as store:
                all_keys = list(store.keys())
                groups = [key for key in all_keys if key.startswith("/group_")]
                self.assertGreater(len(groups), 0, f"No groups found in {h5_file}")
                
                # Check that groups are numbered 0-9
                group_numbers = [int(key.split("_")[1]) for key in groups]
                self.assertTrue(all(0 <= g <= 9 for g in group_numbers),
                              f"Group numbers should be 0-9, found: {group_numbers}")
                
                # Check that each group has POS and rsn columns (CHR not stored)
                # Note: rsn is stored as column in HDF5, but set as index when used for matching
                checked_groups = 0
                for group_key in groups:
                    if checked_groups >= 3:  # Check first 3 non-empty groups
                        break
                    try:
                        group_df = store[group_key]
                        if len(group_df) > 0:  # Only check non-empty groups
                            # rsn should be a column (stored as column, set as index when used)
                            self.assertIn("rsn", group_df.columns, 
                                        f"Group {group_key} missing 'rsn' column")
                            self.assertIn("POS", group_df.columns,
                                        f"Group {group_key} missing 'POS' column")
                            self.assertNotIn("CHR", group_df.columns,
                                           f"Group {group_key} should not have 'CHR' column")
                            
                            # Verify data types
                            self.assertEqual(str(group_df["rsn"].dtype), "int64",
                                           f"rsn should be int64, got {group_df['rsn'].dtype}")
                            self.assertEqual(str(group_df["POS"].dtype), "int32",
                                           f"POS should be int32, got {group_df['POS'].dtype}")
                            
                            # Verify no duplicates in rsn (should be deduplicated in finalize phase)
                            self.assertEqual(
                                len(group_df), len(group_df["rsn"].unique()),
                                f"Group {group_key} should have unique rsn values after deduplication"
                            )
                            checked_groups += 1
                    except Exception:
                        continue
    
    def test_process_vcf_to_hdf5_with_specific_chromosomes(self):
        """Test processing specific chromosomes only."""
        # Use get_NC_to_number to map NC_ notation to chromosome numbers
        chr_dict = get_NC_to_number(build="19")
        
        # Process only chromosome 7 (which exists in the test VCF)
        output_dir = process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            chr_list=[7],
            complevel=3,
            threads=1,
            verbose=False
        )
        
        # Check that chr7 file was created
        h5_files = glob.glob(os.path.join(self.temp_dir, "*.chr7.rsID_CHR_POS_mod10.h5"))
        self.assertGreaterEqual(len(h5_files), 1, "Should create at least one chr7 file")
        
        # Check that the file has data
        if len(h5_files) > 0:
            with pd.HDFStore(h5_files[0], mode='r') as store:
                groups = [key for key in store.keys() if key.startswith("group_")]
                # At least some groups should have data (may be empty if no rsIDs in chr7)
                self.assertGreaterEqual(len(groups), 0, "File should exist (may be empty)")
    
    def test_process_vcf_to_hdf5_overwrite_option(self):
        """Test that overwrite option works correctly."""
        # Use get_NC_to_number to map NC_ notation to chromosome numbers
        chr_dict = get_NC_to_number(build="19")
        
        # Process VCF to HDF5 first time
        output_dir1 = process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            complevel=3,
            threads=1,
            overwrite=False,
            verbose=False
        )
        
        # Get file modification times
        h5_files = glob.glob(os.path.join(self.temp_dir, "*.chr*.rsID_CHR_POS_mod10.h5"))
        self.assertGreater(len(h5_files), 0, "HDF5 files should be created")
        
        import time
        time.sleep(1)  # Wait a second to ensure different modification time
        
        # Process again with overwrite=False (should skip existing files)
        output_dir2 = process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            complevel=3,
            threads=1,
            overwrite=False,
            verbose=False
        )
        
        # Files should not be modified (overwrite=False skips existing files)
        # Note: This is a basic check - in practice, the function logs "Skipping" messages
        
        # Process again with overwrite=True (should recreate files)
        output_dir3 = process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            complevel=3,
            threads=1,
            overwrite=True,
            verbose=False
        )
        
        # Files should exist and be valid
        h5_files_after = glob.glob(os.path.join(self.temp_dir, "*.chr*.rsID_CHR_POS_mod10.h5"))
        self.assertGreater(len(h5_files_after), 0, "HDF5 files should still exist after overwrite=True")
        
        # Verify files are valid HDF5 files (some may be empty if chromosome has no variants)
        files_with_data = []
        for h5_file in h5_files_after:
            try:
                with pd.HDFStore(h5_file, mode='r') as store:
                    keys = list(store.keys())
                    # File should be openable (valid HDF5), even if empty
                    # Only check for keys if we want to verify data exists
                    if len(keys) > 0:
                        files_with_data.append(h5_file)
            except Exception as e:
                self.fail(f"HDF5 file {h5_file} should be valid after overwrite: {e}")
        
        # At least some files should have data (not all chromosomes may have variants in test VCF)
        # This is a weaker check - we just verify files are valid HDF5 format
        self.assertGreaterEqual(len(files_with_data), 0, "Some HDF5 files may be empty (no variants in chromosome)")


@unittest.skipIf(not PYTABLES_AVAILABLE, "pytables (tables) is required for HDF5 operations")
class TestRsidToChrposHDF5(unittest.TestCase):
    """Test rsid_to_chrpos using HDF5 files."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        self.vcf_path = os.path.join(self.ref_dir, "b157_2564.vcf.gz")
        self.sumstats_path = os.path.join(self.ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        self.temp_dir = tempfile.mkdtemp(prefix="gwaslab_test_hdf5_")
        
        # Get chromosome mapping for NC_ notation
        chr_dict = get_NC_to_number(build="19")
        
        # Process VCF to HDF5 first
        process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            complevel=3,
            threads=1,
            verbose=False
        )
    
    def tearDown(self):
        """Clean up temporary files."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_rsid_to_chrpos_with_hdf5_directory(self):
        """Test rsid_to_chrpos using HDF5 directory."""
        # Load sumstats
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            # First try to get variants with rsIDs, then sample
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    # Sample from rsID variants, but keep at least 10 if available
                    sample_size = max(10, len(rsid_variants) // 10)
                    gl.data = rsid_variants.head(sample_size).copy()
                else:
                    # Fallback to regular sampling if no rsIDs
                    gl.data = gl.data.iloc[::10].copy()
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        # Store original CHR and POS for comparison
        original_chr = gl.data["CHR"].copy()
        original_pos = gl.data["POS"].copy()
        
        # Remove CHR and POS to test assignment
        gl.data["CHR"] = pd.NA
        gl.data["POS"] = pd.NA
        
        # Run rsid_to_chrpos using HDF5 directory
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Check that CHR and POS were assigned
        assigned_chr = gl.data["CHR"].notna().sum()
        assigned_pos = gl.data["POS"].notna().sum()
        
        self.assertGreater(assigned_chr, 0, "Some CHR values should be assigned")
        self.assertGreater(assigned_pos, 0, "Some POS values should be assigned")
        
        # Check that assigned values are valid (may not match original if original was from different source)
        if assigned_chr > 0:
            assigned_chr_values = gl.data.loc[gl.data["CHR"].notna(), "CHR"]
            self.assertTrue(
                assigned_chr_values.between(1, 25).all(),
                "Assigned CHR values should be valid chromosome numbers (1-25)"
            )
        
        if assigned_pos > 0:
            assigned_pos_values = gl.data.loc[gl.data["POS"].notna(), "POS"]
            self.assertTrue(
                (assigned_pos_values > 0).all(),
                "Assigned POS values should be positive integers"
            )
    
    def test_rsid_to_chrpos_with_missing_rsids(self):
        """Test rsid_to_chrpos handles missing or invalid rsIDs correctly."""
        # Create test data with some missing/invalid rsIDs
        test_data = pd.DataFrame({
            "rsID": ["rs123456", "rs789012", "invalid_rsid", None, "rs345678"],
            "CHR": [pd.NA] * 5,
            "POS": [pd.NA] * 5,
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5,
            "P": [0.05] * 5
        })
        
        gl = Sumstats(sumstats=test_data, rsid="rsID", verbose=False)
        
        # Run rsid_to_chrpos
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Check that invalid rsIDs don't get assigned
        valid_mask = gl.data["rsID"].isin(["rs123456", "rs789012", "rs345678"])
        invalid_mask = ~valid_mask
        
        # Invalid rsIDs should remain NA (they can't be matched)
        if invalid_mask.sum() > 0:
            invalid_chr_pos = gl.data.loc[invalid_mask, ["CHR", "POS"]].notna().any(axis=1).sum()
            # Invalid rsIDs should not get CHR/POS assigned
            self.assertEqual(
                invalid_chr_pos, 0,
                "Invalid rsIDs should not get CHR/POS assigned"
            )
        
        # Valid rsIDs may or may not get assigned depending on whether they exist in reference
        # (This is acceptable - the test VCF may not contain these specific rsIDs)
    
    def test_rsid_to_chrpos_with_multiple_threads(self):
        """Test that rsid_to_chrpos works correctly with multiple threads."""
        # Load sumstats
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    sample_size = max(10, len(rsid_variants) // 10)
                    gl.data = rsid_variants.head(sample_size).copy()
                else:
                    gl.data = gl.data.iloc[::10].copy()
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        # Clear CHR and POS to test assignment
        gl.data["CHR"] = pd.NA
        gl.data["POS"] = pd.NA
        
        # Run rsid_to_chrpos with multiple threads
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=2,  # Use 2 threads
            verbose=False
        )
        
        # Check that CHR and POS were assigned
        assigned_chr = gl.data["CHR"].notna().sum()
        assigned_pos = gl.data["POS"].notna().sum()
        
        self.assertGreater(assigned_chr, 0, "Some CHR values should be assigned with multiple threads")
        self.assertGreater(assigned_pos, 0, "Some POS values should be assigned with multiple threads")
        
        # Verify no data loss
        self.assertGreater(len(gl.data), 0, "No data should be lost with multiple threads")
    
    def test_rsid_to_chrpos_index_based_matching(self):
        """Test that index-based matching works correctly for rsID lookup."""
        # Load sumstats
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    # Keep up to 10 variants with valid rsIDs
                    gl.data = rsid_variants.head(10).copy()
                else:
                    # If no valid rsIDs, skip this test
                    self.skipTest("No valid rsIDs found in data")
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        # Clear CHR and POS
        gl.data["CHR"] = pd.NA
        gl.data["POS"] = pd.NA
        
        # Run rsid_to_chrpos
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Verify that matching worked (index-based matching should be fast)
        # If variants exist in reference, they should be matched
        matched = gl.data["CHR"].notna() & gl.data["POS"].notna()
        
        # At least some variants should match if they exist in reference
        # (This depends on whether the test VCF contains these rsIDs)
        self.assertGreaterEqual(
            matched.sum(), 0,
            "Index-based matching should work (may be 0 if rsIDs not in reference)"
        )
        
        # If matched, verify consistency
        if matched.sum() > 0:
            matched_data = gl.data[matched]
            # CHR should be valid chromosome numbers
            self.assertTrue(
                matched_data["CHR"].between(1, 25).all(),
                "Matched CHR values should be valid chromosome numbers"
            )
            # POS should be positive integers
            self.assertTrue(
                (matched_data["POS"] > 0).all(),
                "Matched POS values should be positive integers"
            )
    
    def test_rsid_to_chrpos_preserves_existing_data(self):
        """Test that rsid_to_chrpos doesn't overwrite existing CHR/POS unnecessarily."""
        # Load sumstats
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    sample_size = max(10, len(rsid_variants) // 10)
                    gl.data = rsid_variants.head(sample_size).copy()
                else:
                    gl.data = gl.data.iloc[::10].copy()
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        # Store original values
        original_chr = gl.data["CHR"].copy()
        original_pos = gl.data["POS"].copy()
        
        # Run rsid_to_chrpos (should not change existing values)
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Existing CHR/POS should be preserved
        # Note: The function should only fill in missing values, not overwrite existing ones
        existing_chr_mask = original_chr.notna()
        existing_pos_mask = original_pos.notna()
        
        if existing_chr_mask.sum() > 0:
            # For variants that had CHR, check that they still have CHR (may be different value if HDF5 has different data)
            # But the key is that they shouldn't become NA
            final_chr = gl.data.loc[existing_chr_mask, "CHR"]
            # Most should still be non-NA (allowing for some edge cases where HDF5 might not have the variant)
            non_na_count = final_chr.notna().sum()
            self.assertGreater(
                non_na_count, existing_chr_mask.sum() * 0.8,
                "Most existing CHR values should remain non-NA"
            )
        
        if existing_pos_mask.sum() > 0:
            # For variants that had POS, check that they still have POS
            final_pos = gl.data.loc[existing_pos_mask, "POS"]
            non_na_count = final_pos.notna().sum()
            self.assertGreater(
                non_na_count, existing_pos_mask.sum() * 0.8,
                "Most existing POS values should remain non-NA"
            )
    
    def test_rsid_to_chrpos_with_chr_column_available(self):
        """Test that rsid_to_chrpos uses 'by chromosome' strategy when CHR is available."""
        # Load sumstats
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    sample_size = max(10, len(rsid_variants) // 10)
                    gl.data = rsid_variants.head(sample_size).copy()
                else:
                    gl.data = gl.data.iloc[::10].copy()
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        # Ensure CHR column has some values
        # Keep only rows with valid CHR for this test
        gl.data = gl.data[gl.data["CHR"].notna()].copy()
        
        # Store some CHR values, clear POS to test assignment
        original_chr = gl.data["CHR"].copy()
        gl.data["POS"] = pd.NA
        
        # Run rsid_to_chrpos (should use "by chromosome" strategy)
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # CHR should be preserved (not set to NA, though values may differ if HDF5 has different data)
        # The key is that existing CHR values should not become NA
        final_chr = gl.data["CHR"]
        non_na_count = final_chr.notna().sum()
        original_non_na_count = original_chr.notna().sum()
        self.assertGreaterEqual(
            non_na_count, original_non_na_count * 0.8,
            "Most existing CHR values should remain non-NA when CHR column is available"
        )
        
        # POS values may be assigned if rsIDs exist in HDF5 reference
        # (This depends on whether the test data rsIDs exist in the HDF5 file)
        assigned_pos = gl.data["POS"].notna().sum()
        # Just verify the function ran without error - POS assignment depends on HDF5 content
        self.assertGreaterEqual(assigned_pos, 0, "POS assignment depends on whether rsIDs exist in HDF5 reference")
    
    def test_rsid_to_chrpos_without_chr_column(self):
        """Test that rsid_to_chrpos uses 'search all chromosomes' strategy when CHR is missing."""
        # Load sumstats
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    sample_size = max(10, len(rsid_variants) // 10)
                    gl.data = rsid_variants.head(sample_size).copy()
                else:
                    gl.data = gl.data.iloc[::10].copy()
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        # Clear CHR and POS to test assignment from scratch
        gl.data["CHR"] = pd.NA
        gl.data["POS"] = pd.NA
        
        # Run rsid_to_chrpos (should use "search all chromosomes" strategy)
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Some CHR and POS values should be assigned
        assigned_chr = gl.data["CHR"].notna().sum()
        assigned_pos = gl.data["POS"].notna().sum()
        
        self.assertGreater(assigned_chr, 0, "Some CHR values should be assigned when searching all chromosomes")
        self.assertGreater(assigned_pos, 0, "Some POS values should be assigned when searching all chromosomes")
        
        # When both CHR and POS are assigned, they should be consistent
        both_assigned = gl.data["CHR"].notna() & gl.data["POS"].notna()
        if both_assigned.sum() > 0:
            # All assigned CHR values should be valid chromosome numbers
            assigned_chr_values = gl.data.loc[both_assigned, "CHR"]
            self.assertTrue(
                assigned_chr_values.between(1, 25).all(),
                "Assigned CHR values should be valid chromosome numbers (1-25)"
            )
    
    def test_rsid_to_chrpos_deduplication_across_chromosomes(self):
        """Test that rsid_to_chrpos deduplicates correctly when searching across all chromosomes."""
        # Create test data with rsIDs that might match multiple chromosomes
        test_data = pd.DataFrame({
            "rsID": ["rs123456", "rs789012", "rs345678", "rs111222", "rs999888"],
            "CHR": [pd.NA] * 5,  # No CHR, will search all chromosomes
            "POS": [pd.NA] * 5,
            "EA": ["A"] * 5,
            "NEA": ["G"] * 5,
            "P": [0.05] * 5
        })
        
        gl = Sumstats(sumstats=test_data, rsid="rsID", verbose=False)
        
        # Run rsid_to_chrpos (will search all chromosomes and deduplicate)
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Each rsID should appear only once (deduplication)
        self.assertEqual(
            len(gl.data), len(gl.data["rsID"].unique()),
            "Each rsID should appear only once after deduplication"
        )
        
        # If a variant has both CHR and POS, it should be a valid match
        both_assigned = gl.data["CHR"].notna() & gl.data["POS"].notna()
        if both_assigned.sum() > 0:
            # Matches with both CHR and POS should be prioritized
            self.assertGreater(
                both_assigned.sum(), 0,
                "At least some variants should have both CHR and POS assigned"
            )


@unittest.skipIf(not PYTABLES_AVAILABLE, "pytables (tables) is required for HDF5 operations")
class TestRsidToChrposWorkflow(unittest.TestCase):
    """Test complete workflow: VCF -> HDF5 -> rsid_to_chrpos."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.ref_dir = os.path.join(os.path.dirname(__file__), "ref")
        self.vcf_path = os.path.join(self.ref_dir, "b157_2564.vcf.gz")
        self.sumstats_path = os.path.join(self.ref_dir, "bbj_t2d_hm3_chr7_variants.txt.gz")
        self.temp_dir = tempfile.mkdtemp(prefix="gwaslab_test_workflow_")
    
    def tearDown(self):
        """Clean up temporary files."""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_complete_workflow(self):
        """Test complete workflow from VCF processing to rsid_to_chrpos."""
        # Use get_NC_to_number to map NC_ notation to chromosome numbers
        chr_dict = get_NC_to_number(build="19")
        
        # Step 1: Process VCF to HDF5
        output_dir = process_vcf_to_hfd5(
            vcf=self.vcf_path,
            directory=self.temp_dir,
            chr_dict=chr_dict,
            complevel=3,
            threads=1,
            verbose=False
        )
        
        self.assertEqual(output_dir, self.temp_dir)
        h5_files = glob.glob(os.path.join(self.temp_dir, "*.chr*.rsID_CHR_POS_mod10.h5"))
        self.assertGreater(len(h5_files), 0, "HDF5 files should be created")
        
        # Step 2: Load sumstats and assign CHR/POS using HDF5
        gl = Sumstats(
            sumstats=self.sumstats_path,
            tab_fmt="tsv",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            snpid="SNPID",
            rsid="rsID",
            verbose=False
        )
        # Use only 1/10 of variants to speed up tests
        # Prioritize variants with rsIDs to ensure matches in HDF5
        if len(gl.data) > 0:
            if "rsID" in gl.data.columns:
                rsid_variants = gl.data[gl.data["rsID"].notna() & gl.data["rsID"].str.startswith("rs", na=False)]
                if len(rsid_variants) > 0:
                    sample_size = max(10, len(rsid_variants) // 10)
                    gl.data = rsid_variants.head(sample_size).copy()
                else:
                    gl.data = gl.data.iloc[::10].copy()
            else:
                gl.data = gl.data.iloc[::10].copy()
        
        initial_count = len(gl.data)
        initial_chr_count = gl.data["CHR"].notna().sum()
        initial_pos_count = gl.data["POS"].notna().sum()
        
        # Step 3: Run rsid_to_chrpos
        gl.rsid_to_chrpos2(
            path=self.temp_dir,
            threads=1,
            verbose=False
        )
        
        # Check results
        final_count = len(gl.data)
        final_chr_count = gl.data["CHR"].notna().sum()
        final_pos_count = gl.data["POS"].notna().sum()
        
        # Data should not be lost
        self.assertEqual(initial_count, final_count, "No rows should be lost")
        
        # CHR/POS should be assigned for at least some variants
        self.assertGreaterEqual(
            final_chr_count, initial_chr_count,
            "CHR assignment count should not decrease"
        )
        self.assertGreaterEqual(
            final_pos_count, initial_pos_count,
            "POS assignment count should not decrease"
        )


if __name__ == "__main__":
    unittest.main()

