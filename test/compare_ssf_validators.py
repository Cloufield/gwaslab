#!/usr/bin/env python3
"""
Standalone script to compare new SSF validator with gwas-ssf CLI.
Run with: python test/compare_ssf_validators.py
"""
import os
import sys
import subprocess
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path

# Add src to path
ROOT = Path(__file__).parent.parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from gwaslab.extension.gwas_sumstats_tools.validate_ssf import validate_ssf_file
from gwaslab.info.g_Log import Log

# Path to gwas-ssf CLI
GWAS_SSF_CLI = "/home/yunye/anaconda3/envs/gwas-ssf-py311/bin/gwas-ssf"


def run_gwas_ssf_validate(filename: Path) -> tuple:
    """Run gwas-ssf validate and return (is_valid, stdout, stderr)."""
    try:
        result = subprocess.run(
            [GWAS_SSF_CLI, "validate", str(filename)],
            capture_output=True,
            text=True,
            timeout=60
        )
        is_valid = result.returncode == 0
        return is_valid, result.stdout, result.stderr
    except FileNotFoundError:
        return None, "", f"gwas-ssf CLI not found at {GWAS_SSF_CLI}"
    except Exception as e:
        return None, "", f"Error: {e}"


def create_test_file(data: pd.DataFrame, base_path: Path, compress=False) -> Path:
    """Create test SSF file."""
    if compress:
        output = base_path.with_suffix('.tsv.gz')
        data.to_csv(output, sep='\t', index=False, compression='gzip')
    else:
        output = base_path.with_suffix('.tsv')
        data.to_csv(output, sep='\t', index=False)
    return output


def test_case(name: str, data: pd.DataFrame, expected_new: bool = None, expected_orig: bool = None):
    """Run a single test case and compare results."""
    print(f"\n{'='*80}")
    print(f"Test: {name}")
    print('='*80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test"
        ssf_file = create_test_file(data, test_file)
        
        # New validator
        log = Log()
        new_valid, new_msg, new_errors = validate_ssf_file(
            ssf_file, log=log, verbose=False, minimum_rows=100000
        )
        
        # Original validator
        orig_valid, orig_stdout, orig_stderr = run_gwas_ssf_validate(ssf_file)
        
        # Print results
        print(f"New validator:    {'✓ PASS' if new_valid else '✗ FAIL'}")
        if not new_valid:
            print(f"  Message: {new_msg}")
            if new_errors:
                print(f"  Errors: {new_errors[:3]}")  # First 3 errors
        
        if orig_valid is not None:
            print(f"Original validator: {'✓ PASS' if orig_valid else '✗ FAIL'}")
            if orig_stderr and not orig_valid:
                print(f"  Stderr: {orig_stderr[:200]}")
        else:
            print(f"Original validator: ⚠ SKIPPED (CLI not available)")
        
        # Compare
        if orig_valid is not None:
            if new_valid == orig_valid:
                print(f"✓ Results MATCH")
            else:
                print(f"✗ Results DIFFER!")
                print(f"  New: {new_valid}, Original: {orig_valid}")
        
        # Check expectations (only assert if original validator is available)
        if expected_new is not None:
            assert new_valid == expected_new, f"Expected new validator to return {expected_new}, got {new_valid}"
        if expected_orig is not None and orig_valid is not None:
            # Only assert if original validator ran successfully
            if "Traceback" not in orig_stderr and "Error" not in orig_stderr[:50]:
                assert orig_valid == expected_orig, f"Expected original validator to return {expected_orig}, got {orig_valid}"


def main():
    """Run all comparison tests."""
    print("SSF Validator Comparison")
    print("Comparing new pandas-based validator with gwas-ssf CLI")
    
    # Test 1: Valid file
    n = 150000
    valid_data = pd.DataFrame({
        'chromosome': np.random.randint(1, 23, n),
        'base_pair_location': np.random.randint(1, 250000000, n),
        'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n),
        'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n),
        'beta': np.random.normal(0, 0.1, n),
        'standard_error': np.abs(np.random.normal(0, 0.05, n)),
        'effect_allele_frequency': np.random.uniform(0, 1, n),
        'p_value': np.random.uniform(1e-10, 1, n),
    })
    test_case("Valid SSF file", valid_data, expected_new=True, expected_orig=True)
    
    # Test 2: Wrong column order
    wrong_order_data = pd.DataFrame({
        'base_pair_location': np.random.randint(1, 250000000, n),  # Wrong position
        'chromosome': np.random.randint(1, 23, n),
        'effect_allele': np.random.choice(['A', 'T', 'G', 'C'], n),
        'other_allele': np.random.choice(['A', 'T', 'G', 'C'], n),
        'beta': np.random.normal(0, 0.1, n),
        'standard_error': np.abs(np.random.normal(0, 0.05, n)),
        'effect_allele_frequency': np.random.uniform(0, 1, n),
        'p_value': np.random.uniform(1e-10, 1, n),
    })
    test_case("Wrong column order", wrong_order_data, expected_new=False, expected_orig=False)
    
    # Test 3: Missing column
    missing_col_data = valid_data.drop(columns=['other_allele'])
    test_case("Missing required column", missing_col_data, expected_new=False, expected_orig=False)
    
    # Test 4: Insufficient rows
    small_data = valid_data.head(50000)
    test_case("Insufficient rows (<100k)", small_data, expected_new=False, expected_orig=False)
    
    # Test 5: Invalid chromosome
    invalid_chr_data = valid_data.copy()
    invalid_chr_data.loc[0, 'chromosome'] = 26
    invalid_chr_data.loc[1, 'chromosome'] = 0
    test_case("Invalid chromosome values", invalid_chr_data, expected_new=False, expected_orig=False)
    
    # Test 6: Invalid p-value
    invalid_pval_data = valid_data.copy()
    invalid_pval_data.loc[0, 'p_value'] = 1.5
    invalid_pval_data.loc[1, 'p_value'] = -0.1
    test_case("Invalid p-value (outside [0,1])", invalid_pval_data, expected_new=False, expected_orig=False)
    
    # Test 7: Very small p-values (scientific notation)
    small_pval_data = valid_data.copy()
    small_pval_data.loc[0, 'p_value'] = '1.23e-300'
    small_pval_data.loc[1, 'p_value'] = '5.67e-200'
    test_case("Very small p-values (scientific notation)", small_pval_data, expected_new=True, expected_orig=True)
    
    # Test 8: Missing autosomes
    missing_autosomes_data = valid_data.copy()
    missing_autosomes_data['chromosome'] = np.random.randint(1, 11, n)  # Only 1-10
    test_case("Missing autosomes (only 1-10)", missing_autosomes_data, expected_new=False, expected_orig=False)
    
    # Test 9: Invalid allele
    invalid_allele_data = valid_data.copy()
    invalid_allele_data.loc[0, 'effect_allele'] = 'X'  # Invalid nucleotide
    test_case("Invalid allele (non-ACTG)", invalid_allele_data, expected_new=False, expected_orig=False)
    
    # Test 10: Negative odds_ratio
    or_data = valid_data.copy()
    or_data = or_data.drop(columns=['beta'])
    or_data['odds_ratio'] = np.random.uniform(0.1, 10, n)
    or_data.loc[0, 'odds_ratio'] = -0.5  # Invalid
    test_case("Negative odds_ratio", or_data, expected_new=False, expected_orig=False)
    
    print(f"\n{'='*80}")
    print("Comparison complete!")
    print('='*80)


if __name__ == "__main__":
    main()

