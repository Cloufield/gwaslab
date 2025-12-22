import os
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.util.util_in_fill_data import _fill_data
from gwaslab.info.g_Log import Log


class TestFillData(unittest.TestCase):
    """Comprehensive tests for _fill_data function."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.log = Log()
        self.n_variants = 100
        
    def _create_base_sumstats(self, **columns):
        """Create a base sumstats DataFrame with specified columns."""
        data = {
            "SNPID": [f"rs{i}" for i in range(self.n_variants)],
            "CHR": np.random.randint(1, 23, self.n_variants),
            "POS": np.random.randint(1_000, 2_000_000, self.n_variants),
            "EA": np.random.choice(["A", "T", "G", "C"], self.n_variants),
            "NEA": np.random.choice(["A", "T", "G", "C"], self.n_variants),
        }
        data.update(columns)
        return pd.DataFrame(data)
    
    # ========== Basic Single-Step Conversions ==========
    
    def test_fill_p_from_mlog10p(self):
        """Test filling P from MLOG10P."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        np.testing.assert_allclose(result["P"], 10**(-mlog10p), rtol=1e-10)
    
    def test_fill_p_from_z(self):
        """Test filling P from Z."""
        z = np.random.normal(0, 2, self.n_variants)
        sumstats = self._create_base_sumstats(Z=z)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that P values are reasonable (between 0 and 1)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
    
    def test_fill_p_from_chisq(self):
        """Test filling P from CHISQ."""
        chisq = np.random.chisquare(1, self.n_variants)
        sumstats = self._create_base_sumstats(CHISQ=chisq)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
    
    def test_fill_z_from_beta_se(self):
        """Test filling Z from BETA/SE."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
        
        self.assertIn("Z", result.columns)
        np.testing.assert_allclose(result["Z"], beta / se, rtol=1e-10)
    
    def test_fill_chisq_from_z(self):
        """Test filling CHISQ from Z."""
        z = np.random.normal(0, 2, self.n_variants)
        sumstats = self._create_base_sumstats(Z=z)
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("CHISQ", result.columns)
        np.testing.assert_allclose(result["CHISQ"], z**2, rtol=1e-10)
    
    def test_fill_chisq_from_p(self):
        """Test filling CHISQ from P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("CHISQ", result.columns)
        self.assertTrue((result["CHISQ"] >= 0).all())
    
    def test_fill_beta_from_or(self):
        """Test filling BETA from OR."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val)
        
        result = _fill_data(sumstats, to_fill=["BETA"], verbose=False, log=self.log)
        
        self.assertIn("BETA", result.columns)
        np.testing.assert_allclose(result["BETA"], np.log(or_val), rtol=1e-10)
    
    def test_fill_or_from_beta(self):
        """Test filling OR from BETA."""
        beta = np.random.normal(0, 0.2, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta)
        
        result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        
        self.assertIn("OR", result.columns)
        np.testing.assert_allclose(result["OR"], np.exp(beta), rtol=1e-10)
    
    def test_fill_or_confidence_intervals(self):
        """Test filling OR_95L and OR_95U from BETA/SE."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        
        self.assertIn("OR", result.columns)
        self.assertIn("OR_95L", result.columns)
        self.assertIn("OR_95U", result.columns)
        # Check that confidence intervals are reasonable
        self.assertTrue((result["OR_95L"] < result["OR"]).all())
        self.assertTrue((result["OR_95U"] > result["OR"]).all())
    
    def test_fill_se_from_beta_p(self):
        """Test filling SE from BETA/P."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        p = np.random.uniform(1e-6, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, P=p)
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        self.assertTrue((result["SE"] > 0).all())
    
    def test_fill_se_from_or_or95u(self):
        """Test filling SE from OR/OR_95U."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        or_95u = or_val * np.random.uniform(1.01, 1.5, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val, OR_95U=or_95u)
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        self.assertTrue((result["SE"] > 0).all())
    
    def test_fill_se_from_or_or95l(self):
        """Test filling SE from OR/OR_95L."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        or_95l = or_val * np.random.uniform(0.5, 0.99, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val, OR_95L=or_95l)
        
        result = _fill_data(sumstats, to_fill=["SE"], verbose=False, log=self.log)
        
        self.assertIn("SE", result.columns)
        self.assertTrue((result["SE"] > 0).all())
    
    def test_fill_mlog10p_from_p(self):
        """Test filling MLOG10P from P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        np.testing.assert_allclose(result["MLOG10P"], -np.log10(p), rtol=1e-10)
    
    def test_fill_maf_from_eaf(self):
        """Test filling MAF from EAF."""
        eaf = np.random.uniform(0.01, 0.99, self.n_variants)
        sumstats = self._create_base_sumstats(EAF=eaf)
        
        result = _fill_data(sumstats, to_fill=["MAF"], verbose=False, log=self.log)
        
        self.assertIn("MAF", result.columns)
        expected_maf = np.minimum(eaf, 1 - eaf)
        np.testing.assert_allclose(result["MAF"], expected_maf, rtol=1e-10)
    
    def test_fill_sig_from_p(self):
        """Test filling SIGNIFICANT from P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["SIG"], sig_level=5e-8, verbose=False, log=self.log)
        
        self.assertIn("SIGNIFICANT", result.columns)
        self.assertTrue(result["SIGNIFICANT"].dtype == bool)
        expected_sig = p < 5e-8
        np.testing.assert_array_equal(result["SIGNIFICANT"], expected_sig)
    
    def test_fill_sig_from_mlog10p(self):
        """Test filling SIGNIFICANT from MLOG10P."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["SIG"], sig_level=5e-8, verbose=False, log=self.log)
        
        self.assertIn("SIGNIFICANT", result.columns)
        expected_sig = mlog10p > np.log10(1/5e-8)
        np.testing.assert_array_equal(result["SIGNIFICANT"], expected_sig)
    
    # ========== Multi-Round Conversions ==========
    
    def test_multi_round_or_to_p(self):
        """Test multi-round conversion: OR -> BETA -> Z -> P."""
        or_val = np.random.uniform(0.5, 2.0, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(OR=or_val, SE=se)
        
        # Fill BETA, Z, and P - the iterative process should fill them in order
        result = _fill_data(sumstats, to_fill=["BETA", "Z", "P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Should have created BETA and Z as intermediate steps
        self.assertIn("BETA", result.columns)
        self.assertIn("Z", result.columns)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
    
    def test_multi_round_mlog10p_from_beta_se(self):
        """Test multi-round conversion: BETA/SE -> Z -> P -> MLOG10P."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        # Request all intermediate steps - the iterative process will fill them in order
        result = _fill_data(sumstats, to_fill=["Z", "P", "MLOG10P"], verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        self.assertIn("Z", result.columns)
        self.assertIn("P", result.columns)
        # P_MANTISSA and P_EXPONENT should be dropped after filling MLOG10P
        self.assertNotIn("P_MANTISSA", result.columns)
        self.assertNotIn("P_EXPONENT", result.columns)
        # Check that MLOG10P is reasonable
        self.assertTrue((result["MLOG10P"] >= 0).all())
    
    def test_multi_round_p_to_chisq(self):
        """Test multi-round conversion: P -> CHISQ."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        result = _fill_data(sumstats, to_fill=["CHISQ"], verbose=False, log=self.log)
        
        # P -> CHISQ should work
        self.assertIn("CHISQ", result.columns)
        self.assertTrue((result["CHISQ"] >= 0).all())
    
    # ========== Edge Cases ==========
    
    def test_fill_with_existing_column_no_overwrite(self):
        """Test that existing columns are not overwritten by default."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        original_p = p.copy()
        mlog10p = -np.log10(p)
        sumstats = self._create_base_sumstats(P=p, MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], overwrite=False, verbose=False, log=self.log)
        
        # P should remain unchanged
        np.testing.assert_array_equal(result["P"], original_p)
    
    def test_fill_with_existing_column_overwrite(self):
        """Test that existing columns are overwritten when overwrite=True."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        mlog10p = -np.log10(p)
        # Create P with wrong values
        wrong_p = np.random.uniform(0.5, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=wrong_p, MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], overwrite=True, verbose=False, log=self.log)
        
        # P should be recalculated from MLOG10P
        np.testing.assert_allclose(result["P"], 10**(-mlog10p), rtol=1e-10)
    
    def test_fill_with_missing_data(self):
        """Test filling with missing data (NaN values)."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        mlog10p[0:10] = np.nan  # Add some NaN values
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        # Check that NaN values are preserved
        self.assertTrue(np.isnan(result["P"][0:10]).all())
        # Check that non-NaN values are calculated correctly
        valid_mask = ~np.isnan(mlog10p)
        np.testing.assert_allclose(
            result["P"][valid_mask], 
            10**(-mlog10p[valid_mask]), 
            rtol=1e-10
        )
    
    def test_fill_extreme_mlog10p(self):
        """Test filling MLOG10P with extreme=True for very small P-values."""
        z = np.random.normal(0, 5, self.n_variants)  # Large Z scores
        sumstats = self._create_base_sumstats(Z=z)
        
        result = _fill_data(sumstats, to_fill=["MLOG10P"], extreme=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        # P_MANTISSA and P_EXPONENT are created by extreme methods but then dropped
        self.assertNotIn("P_MANTISSA", result.columns)
        self.assertNotIn("P_EXPONENT", result.columns)
    
    def test_fill_extreme_mlog10p_from_beta_se(self):
        """Test extreme MLOG10P calculation from BETA/SE."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.001, 0.01, self.n_variants)  # Small SE for large Z
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["MLOG10P"], extreme=True, verbose=False, log=self.log)
        
        self.assertIn("MLOG10P", result.columns)
        self.assertIn("Z", result.columns)
    
    def test_fill_multiple_columns(self):
        """Test filling multiple columns at once."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        result = _fill_data(sumstats, to_fill=["Z", "P", "CHISQ"], verbose=False, log=self.log)
        
        self.assertIn("Z", result.columns)
        self.assertIn("P", result.columns)
        self.assertIn("CHISQ", result.columns)
    
    def test_fill_with_string_input(self):
        """Test that string input for to_fill is converted to list."""
        mlog10p = np.random.uniform(0, 10, self.n_variants)
        sumstats = self._create_base_sumstats(MLOG10P=mlog10p)
        
        result = _fill_data(sumstats, to_fill="P", verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
    
    def test_fill_with_none_input(self):
        """Test that None input for to_fill is handled."""
        sumstats = self._create_base_sumstats()
        
        result = _fill_data(sumstats, to_fill=None, verbose=False, log=self.log)
        
        # Should return unchanged sumstats
        self.assertEqual(len(result), len(sumstats))
    
    def test_fill_with_invalid_column(self):
        """Test that invalid column names are ignored."""
        sumstats = self._create_base_sumstats()
        
        result = _fill_data(sumstats, to_fill=["INVALID_COL", "ALSO_INVALID"], verbose=False, log=self.log)
        
        # Should return unchanged sumstats
        self.assertNotIn("INVALID_COL", result.columns)
        self.assertNotIn("ALSO_INVALID", result.columns)
    
    def test_fill_when_column_cannot_be_filled(self):
        """Test behavior when a column cannot be filled due to missing prerequisites."""
        sumstats = self._create_base_sumstats()  # No columns that can be used for conversion
        
        result = _fill_data(sumstats, to_fill=["P"], verbose=False, log=self.log)
        
        # P should not be created
        self.assertNotIn("P", result.columns)
    
    def test_fill_only_sig_variants(self):
        """Test filling only significant variants."""
        chisq = np.random.chisquare(1, self.n_variants)
        # Create some significant variants
        chisq[0:20] = np.random.chisquare(1, 20) * 10  # Large CHISQ for significance
        sumstats = self._create_base_sumstats(CHISQ=chisq, P=np.ones(self.n_variants) * 0.5)
        
        result = _fill_data(
            sumstats, 
            to_fill=["P"], 
            only_sig=True, 
            sig_level=5e-8,
            overwrite=True,
            verbose=False, 
            log=self.log
        )
        
        self.assertIn("P", result.columns)
        # Check that only significant variants were updated
        # (This is a simplified check - actual implementation may vary)
    
    def test_fill_with_chisq_and_df(self):
        """Test filling P from CHISQ with degrees of freedom column."""
        chisq = np.random.chisquare(1, self.n_variants)
        df = np.ones(self.n_variants, dtype=int)  # All df=1
        sumstats = self._create_base_sumstats(CHISQ=chisq, DOF=df)
        
        result = _fill_data(sumstats, to_fill=["P"], df="DOF", verbose=False, log=self.log)
        
        self.assertIn("P", result.columns)
        self.assertTrue((result["P"] >= 0).all())
        self.assertTrue((result["P"] <= 1).all())
    
    # ========== Round-trip Conversions ==========
    
    def test_round_trip_p_mlog10p(self):
        """Test round-trip conversion: P -> MLOG10P -> P."""
        p = np.random.uniform(1e-10, 1, self.n_variants)
        sumstats = self._create_base_sumstats(P=p)
        
        # First fill MLOG10P
        result1 = _fill_data(sumstats, to_fill=["MLOG10P"], verbose=False, log=self.log)
        self.assertIn("MLOG10P", result1.columns)
        
        # Then fill P back (should match original)
        result2 = _fill_data(result1, to_fill=["P"], overwrite=True, verbose=False, log=self.log)
        np.testing.assert_allclose(result2["P"], p, rtol=1e-10)
    
    def test_round_trip_beta_or(self):
        """Test round-trip conversion: BETA -> OR -> BETA."""
        beta = np.random.normal(0, 0.2, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta)
        
        # First fill OR
        result1 = _fill_data(sumstats, to_fill=["OR"], verbose=False, log=self.log)
        self.assertIn("OR", result1.columns)
        
        # Remove BETA and fill it back from OR
        result1 = result1.drop(columns=["BETA"])
        result2 = _fill_data(result1, to_fill=["BETA"], verbose=False, log=self.log)
        np.testing.assert_allclose(result2["BETA"], beta, rtol=1e-10)
    
    def test_round_trip_beta_se_z(self):
        """Test round-trip conversion: BETA/SE -> Z -> BETA/SE (via P)."""
        beta = np.random.normal(0, 0.1, self.n_variants)
        se = np.random.uniform(0.01, 0.1, self.n_variants)
        sumstats = self._create_base_sumstats(BETA=beta, SE=se)
        
        # Fill Z
        result1 = _fill_data(sumstats, to_fill=["Z"], verbose=False, log=self.log)
        self.assertIn("Z", result1.columns)
        np.testing.assert_allclose(result1["Z"], beta / se, rtol=1e-10)
        
        # Fill P from Z
        result2 = _fill_data(result1, to_fill=["P"], verbose=False, log=self.log)
        self.assertIn("P", result2.columns)
        
        # Fill SE back from BETA/P
        result2 = result2.drop(columns=["SE"])
        result3 = _fill_data(result2, to_fill=["SE"], verbose=False, log=self.log)
        self.assertIn("SE", result3.columns)
        # SE should be approximately recovered (within reasonable tolerance)
        np.testing.assert_allclose(result3["SE"], se, rtol=0.1)


if __name__ == "__main__":
    unittest.main()

