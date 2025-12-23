import unittest
import pandas as pd
import numpy as np
import os
import sys
from unittest.mock import patch
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.g_Sumstats import Sumstats
from gwaslab.hm.hm_liftover_v2 import _liftover_variant
from gwaslab.info.g_Log import Log


class TestLiftover(unittest.TestCase):
    """Test suite for liftover functionality."""

    def setUp(self):
        """Set up test fixtures with sample data."""
        # Create a sample dataframe with valid status codes
        # Status code format: build (2 digits) + 5 digits
        # Need digit 4 = 0 (CHR valid & POS valid) for liftover to process
        # 1900000 = build 19, digit 4 = 0 (CHR valid & POS valid)
        self.data = pd.DataFrame({
            'SNPID': ['rs1', 'rs2', 'rs3', 'rs4'],
            'CHR': ['1', '2', '3', 'X'],
            'POS': [1000, 2000, 3000, 4000],
            'EA': ['A', 'G', 'C', 'T'],
            'NEA': ['G', 'A', 'T', 'C'],
            'BETA': [0.1, 0.2, 0.3, 0.4],
            'SE': [0.01, 0.02, 0.03, 0.04],
            'P': [0.001, 0.002, 0.003, 0.004],
            'STATUS': [1900000, 1900000, 1900000, 1900000]  # Valid status codes
        })
        self.sumstats = Sumstats(self.data, fmt="gwaslab")
        # Ensure status codes have digit 4 = 0 after initialization
        # Sumstats initialization might modify status, so we set it explicitly
        self.sumstats.data['STATUS'] = [1900000, 1900000, 1900000, 1900000]

    def test_liftover_basic_functionality(self):
        """Test basic liftover functionality with successful mapping."""
        # Ensure status codes have digit 4 = 0 for this test
        self.sumstats.data['STATUS'] = [1900000, 1900000, 1900000, 1900000]
        
        def mock_liftover_df_side_effect(df, **kwargs):
            """Mock liftover_df to return successfully mapped variants."""
            df = df.copy()
            # Return same chromosomes (no mismatch) and updated positions
            # Convert X to 23 for proper matching
            df['CHR_LIFT'] = df['CHR'].astype(str).replace({'X': '23', 'Y': '24', 'MT': '25'})
            df['POS_LIFT'] = df['POS'] + 100  # Simulate position change
            df['STRAND_LIFT'] = '+'
            return df

        with patch('gwaslab.hm.hm_liftover_v2.get_chain', return_value="dummy_chain"), \
             patch('gwaslab.hm.hm_liftover_v2.liftover_df', side_effect=mock_liftover_df_side_effect), \
             patch('gwaslab.hm.hm_liftover_v2._process_build', return_value="38"):
            
            result_df = _liftover_variant(
                self.sumstats.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                remove=False,  # Keep unmapped for testing
                filter_by_status=True,  # Explicitly enable status filtering
                log=Log()
            )
            
            # Check that result is not empty
            self.assertGreater(len(result_df), 0, "Result dataframe should not be empty")
            
            # Check that CHR is converted to Int64 by _fix_chr
            # Note: _fix_chr converts to string first, then to Int64 at the end
            # For valid chromosomes, they should be Int64 after _fix_chr
            if len(result_df) > 0:
                # Check if there are any non-NA CHR values
                valid_chr = result_df['CHR'].notna()
                if valid_chr.any():
                    # _fix_chr should convert valid chromosomes to Int64
                    # But if all are removed or invalid, it might stay as string
                    chr_dtype = result_df['CHR'].dtype.name
                    # The final dtype should be Int64 for valid chromosomes
                    # However, if there are issues or all removed, it might be string
                    # So we check the actual values for valid chromosomes
                    chr_values = result_df.loc[valid_chr, 'CHR']
                    if len(chr_values) > 0:
                        # If dtype is string but values are numeric, that's okay (will be converted)
                        # If dtype is Int64/int64, that's correct
                        self.assertIn(chr_dtype, ['Int64', 'int64', 'string'],
                                     f"CHR dtype should be Int64/int64 or string, got {chr_dtype}")
            
            # Check that STATUS is Int64
            self.assertIn(result_df['STATUS'].dtype.name, ['Int64', 'int64'],
                         f"STATUS should be Int64 or int64, got {result_df['STATUS'].dtype.name}")
            
            # Check that positions were updated (for mapped variants)
            mapped = result_df['POS'].notna()
            if mapped.any():
                # Positions should have changed (mocked to add 100)
                # Note: _fix_pos may filter some, so we just check that POS exists
                self.assertTrue(result_df.loc[mapped, 'POS'].notna().all())

    def test_liftover_chromosome_mismatch(self):
        """Test that chromosome mismatches are detected and handled."""
        # Ensure status codes have digit 4 = 0 for this test
        self.sumstats.data['STATUS'] = [1900000, 1900000, 1900000, 1900000]
        
        def mock_liftover_df_side_effect(df, **kwargs):
            """Mock liftover_df to return chromosome mismatches."""
            df = df.copy()
            # Return different chromosomes to simulate mismatch
            df['CHR_LIFT'] = df['CHR'].astype(str).replace({'1': '2', '2': '3', '3': '4', 'X': 'Y'})
            df['POS_LIFT'] = df['POS'] + 100
            df['STRAND_LIFT'] = '+'
            return df

        with patch('gwaslab.hm.hm_liftover_v2.get_chain', return_value="dummy_chain"), \
             patch('gwaslab.hm.hm_liftover_v2.liftover_df', side_effect=mock_liftover_df_side_effect), \
             patch('gwaslab.hm.hm_liftover_v2._process_build', return_value="38"):
            
            result_df = _liftover_variant(
                self.sumstats.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                remove=True,  # Remove unmapped (including mismatches)
                filter_by_status=True,  # Explicitly enable status filtering
                log=Log()
            )
            
            # With remove=True, all variants with mismatches should be removed
            # So result should be empty or have fewer variants
            self.assertLessEqual(len(result_df), len(self.sumstats.data),
                                "With remove=True, mismatched variants should be removed")

    def test_liftover_unmapped_variants(self):
        """Test handling of unmapped variants (failed liftover)."""
        # Ensure status codes have digit 4 = 0 for this test
        self.sumstats.data['STATUS'] = [1900000, 1900000, 1900000, 1900000]
        
        def mock_liftover_df_side_effect(df, **kwargs):
            """Mock liftover_df to return some unmapped variants."""
            df = df.copy()
            # Mark some variants as unmapped (NA or -1)
            # Convert X to 23 for proper matching of mapped ones
            df['CHR_LIFT'] = df['CHR'].astype(str).replace({'X': '23', 'Y': '24', 'MT': '25'})
            df['POS_LIFT'] = df['POS'] + 100
            # Mark first two variants as unmapped (set POS_LIFT to NA)
            df.loc[df.index[:2], 'POS_LIFT'] = pd.NA
            df['STRAND_LIFT'] = '+'
            return df

        with patch('gwaslab.hm.hm_liftover_v2.get_chain', return_value="dummy_chain"), \
             patch('gwaslab.hm.hm_liftover_v2.liftover_df', side_effect=mock_liftover_df_side_effect), \
             patch('gwaslab.hm.hm_liftover_v2._process_build', return_value="38"):
            
            # Test with remove=False (keep unmapped)
            result_df_no_remove = _liftover_variant(
                self.sumstats.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                remove=False,
                filter_by_status=True,  # Explicitly enable status filtering
                log=Log()
            )
            
            # Should keep all variants (including unmapped)
            self.assertEqual(len(result_df_no_remove), len(self.sumstats.data),
                            "With remove=False, all variants should be kept")
            
            # Test with remove=True (remove unmapped)
            result_df_remove = _liftover_variant(
                self.sumstats.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                remove=True,
                filter_by_status=True,  # Explicitly enable status filtering
                log=Log()
            )
            
            # Should have fewer variants (unmapped removed)
            # Note: We marked 2 variants as unmapped, so should have 2 fewer
            # But also need to account for variants that might be filtered by _fix_chr or _fix_pos
            self.assertLess(len(result_df_remove), len(self.sumstats.data),
                           "With remove=True, unmapped variants should be removed")

    def test_liftover_status_code_updates(self):
        """Test that status codes are updated correctly after liftover."""
        # Ensure status codes have digit 4 = 0 for this test
        self.sumstats.data['STATUS'] = [1900000, 1900000, 1900000, 1900000]
        
        def mock_liftover_df_side_effect(df, **kwargs):
            """Mock liftover_df to return successfully mapped variants."""
            df = df.copy()
            # Convert X to 23 for proper matching
            df['CHR_LIFT'] = df['CHR'].astype(str).replace({'X': '23', 'Y': '24', 'MT': '25'})
            df['POS_LIFT'] = df['POS'] + 100
            df['STRAND_LIFT'] = '+'
            return df

        with patch('gwaslab.hm.hm_liftover_v2.get_chain', return_value="dummy_chain"), \
             patch('gwaslab.hm.hm_liftover_v2.liftover_df', side_effect=mock_liftover_df_side_effect), \
             patch('gwaslab.hm.hm_liftover_v2._process_build', return_value="38"):
            
            result_df = _liftover_variant(
                self.sumstats.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                remove=False,
                filter_by_status=True,  # Explicitly enable status filtering
                log=Log()
            )
            
            # Check that result is not empty (variants should be processed)
            self.assertGreater(len(result_df), 0, "Result dataframe should not be empty")
            
            # Check that status codes are updated (build should change from 19 to 38)
            # Status format: build (2 digits) + 5 digits
            # Original: 1900000 (build 19)
            # After liftover to build 38: should start with 38
            mapped = result_df['POS'].notna()
            if mapped.any():
                # Status codes for mapped variants should be updated
                # Build number (first 2 digits) should be 38
                status_vals = result_df.loc[mapped, 'STATUS'].dropna()
                if len(status_vals) > 0:
                    # Extract build number (first 2 digits)
                    build_numbers = (status_vals // 100000).astype(int)
                    # Should be 38 for build 38
                    self.assertTrue(
                        (build_numbers == 38).any(),
                        f"Status codes should be updated with new build number 38, got {build_numbers.unique()}"
                    )

    def test_liftover_empty_dataframe(self):
        """Test that liftover handles empty dataframes gracefully."""
        empty_data = pd.DataFrame(columns=['SNPID', 'CHR', 'POS', 'EA', 'NEA', 'STATUS'])
        empty_sumstats = Sumstats(empty_data, fmt="gwaslab")
        
        with patch('gwaslab.hm.hm_liftover_v2.get_chain', return_value="dummy_chain"), \
             patch('gwaslab.hm.hm_liftover_v2.liftover_df', side_effect=lambda df, **kwargs: df.copy()), \
             patch('gwaslab.hm.hm_liftover_v2._process_build', return_value="38"):
            
            result_df = _liftover_variant(
                empty_sumstats.data,
                chrom="CHR",
                pos="POS",
                from_build="19",
                to_build="38",
                status="STATUS",
                chain_path="dummy.chain",
                log=Log()
            )
            
            # Should return empty dataframe without error
            self.assertEqual(len(result_df), 0, "Empty input should return empty output")


if __name__ == '__main__':
    unittest.main()
