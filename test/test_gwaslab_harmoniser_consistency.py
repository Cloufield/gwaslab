#!/usr/bin/env python3
"""
Test consistency between GWASLab and genetics-sumstat-harmoniser.

This script:
1. Creates toy sumstats with various scenarios
2. Creates a toy VCF reference file
3. Runs GWASLab harmonization
4. Runs genetics-sumstat-harmoniser
5. Compares the results (STATUS codes vs hm_code)
"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import sys
from pathlib import Path
from typing import Optional, List, Tuple

# Ensure we use the latest source code, not the installed package
# Add the source directory to the beginning of sys.path
_test_dir = Path(__file__).parent.resolve()
_project_root = _test_dir.parent.resolve()
_source_dir = _project_root / "src"
if _source_dir.exists() and str(_source_dir) not in sys.path:
    sys.path.insert(0, str(_source_dir))
    # Also clear any cached gwaslab modules to force reload
    modules_to_remove = [m for m in list(sys.modules.keys()) if m.startswith('gwaslab')]
    for m in modules_to_remove:
        del sys.modules[m]

# Import mapping function
def extract_status_digit(status: int, digit: int) -> int:
    """Extract a specific digit from a 7-digit status code."""
    power = 10 ** (7 - digit)
    return (status // power) % 10

def map_gwaslab_to_hm_code(status: int, ea: Optional[str] = None, nea: Optional[str] = None, 
                          vcf_ref: Optional[str] = None, vcf_alt: Optional[str] = None) -> Optional[int]:
    """
    Map GWASLab 7-digit status code to hm_code.
    
    For palindromic reverse strand variants (digit_7=2), the mapping depends on whether
    alleles were flipped during harmonization. If EA/allele information is provided,
    it will be used to determine the correct mapping.
    """
    digit_4 = extract_status_digit(status, 4)
    digit_6 = extract_status_digit(status, 6)
    digit_7 = extract_status_digit(status, 7)
    
    if digit_4 == 9:
        return 14
    if digit_6 == 9 or digit_7 == 9:
        return 14
    if digit_6 == 8:
        return 15
    
    if digit_7 != 0:
        if digit_7 == 8 and digit_6 != 8:
            return 17
        if digit_7 == 7:
            return 18
        if digit_7 == 7 or digit_7 == 8:
            return 9
        
        if digit_7 == 1:
            if digit_6 == 0:
                return 1
            elif digit_6 == 1:
                return 2
        if digit_7 == 2:
            # For reverse strand palindromic, check if alleles were flipped
            # The harmoniser checks if EA matches REF (flipped) or ALT (not flipped)
            if ea is not None and vcf_ref is not None and vcf_alt is not None:
                # If EA matches REF, alleles were flipped → hm_code=4
                if ea == vcf_ref:
                    return 4
                # If EA matches ALT, alleles not flipped → hm_code=3
                elif ea == vcf_alt:
                    return 3
            # Fallback: use digit_6 (but this may be incorrect if flip occurred)
            if digit_6 == 0:
                return 3  # Default assumption: not flipped
            elif digit_6 == 1:
                return 4  # Explicitly flipped
            elif digit_6 == 2:
                return 4  # Reverse complementary fixed (also flipped)
        if digit_7 == 5:
            if digit_6 == 0:
                return 5
            elif digit_6 == 3:
                return 6
            elif digit_6 == 4:
                return 7
            elif digit_6 == 5:
                return 8
        if digit_7 == 3:
            if digit_6 == 0:
                return 10
            elif digit_6 == 1:
                return 11
            elif digit_6 == 2 or digit_6 == 4:
                return 12
        if digit_7 == 4:
            return 11
        if digit_7 == 6:
            if digit_6 == 0 or digit_6 == 1 or digit_6 == 3:
                return 11
            elif digit_6 == 5:
                return 13
    else:
        if digit_6 == 0:
            return 10
        elif digit_6 == 1 or digit_6 == 3:
            return 11
        elif digit_6 == 2 or digit_6 == 4:
            return 12
        elif digit_6 == 5:
            return 13
    
    return None

def create_toy_vcf(vcf_file: str):
    """Create a toy VCF file with various test variants."""
    print(f"Creating toy VCF: {vcf_file}")
    
    # Create test variants covering different scenarios
    variants = [
        # Format: (CHROM, POS, ID, REF, ALT, AF_NFE)
        # Non-palindromic, forward strand, correct
        ("1", 1000, "rs1", "A", "G", 0.3),
        # Non-palindromic, forward strand, flipped
        ("1", 2000, "rs2", "C", "T", 0.4),
        # Non-palindromic, reverse strand, correct
        ("1", 3000, "rs3", "G", "A", 0.2),
        # Non-palindromic, reverse strand, flipped
        ("1", 4000, "rs4", "T", "C", 0.35),
        # Palindromic A/T, forward strand
        ("1", 5000, "rs5", "A", "T", 0.25),  # Low MAF for inference
        # Palindromic G/C, forward strand
        ("1", 6000, "rs6", "G", "C", 0.28),  # Low MAF for inference
        # Palindromic A/T, reverse strand
        ("1", 7000, "rs7", "A", "T", 0.75),  # High AF (reverse strand)
        # Palindromic G/C, reverse strand
        ("1", 8000, "rs8", "G", "C", 0.72),  # High AF (reverse strand)
        # Palindromic with high MAF (should be dropped)
        ("1", 9000, "rs9", "A", "T", 0.45),  # High MAF
        # Indel
        ("1", 10000, "rs10", "A", "AT", 0.3),
    ]
    
    with open(vcf_file, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.3\n")
        f.write("##contig=<ID=1,length=10000000>\n")
        f.write("##INFO=<ID=AF_NFE,Number=A,Type=Float,Description=\"Allele frequency in NFE\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Write variants
        for chrom, pos, rsid, ref, alt, af in variants:
            f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\tAF_NFE={af}\n")
    
    # Index with tabix
    subprocess.run(['bgzip', '-f', vcf_file], check=True)
    subprocess.run(['tabix', '-p', 'vcf', vcf_file + '.gz'], check=True)
    
    print(f"Created VCF with {len(variants)} variants")
    return vcf_file + '.gz'

def create_toy_sumstats(sumstats_file: str):
    """Create toy sumstats with various scenarios."""
    print(f"Creating toy sumstats: {sumstats_file}")
    
    # Create test cases covering different harmonization scenarios
    # Define variants with FASTA reference alleles
    # FASTA REF at each position: 1000=A, 2000=C, 3000=G, 4000=T, 5000=A, 6000=G, 7000=A, 8000=G, 9000=A, 10000=A, 20000=A
    # Format: (SNPID, CHR, POS, EA, NEA, EAF, BETA, SE, P)
    # IMPORTANT: NEA must always match FASTA REF for proper alignment
    # When swapping EA/NEA to align NEA with REF, EAF must be flipped: EAF_new = 1 - EAF_old
    sumstats = [
        # Case 1: Non-palindromic, forward strand, correct alleles (NEA matches REF)
        ("1:1000_A_G", 1, 1000, "G", "A", 0.3, 0.1, 0.01, 1e-5),  # NEA=A matches REF=A ✓
        # Case 2: Non-palindromic, forward strand, flipped alleles (swapped so NEA matches REF)
        ("1:2000_C_T", 1, 2000, "T", "C", 0.4, 0.1, 0.01, 1e-5),  # NEA=C matches REF=C ✓ (EAF: 0.6→0.4)
        # Case 3: Non-palindromic, reverse strand, correct alleles (NEA matches REF)
        ("1:3000_G_A", 1, 3000, "A", "G", 0.8, 0.1, 0.01, 1e-5),  # NEA=G matches REF=G ✓
        # Case 4: Non-palindromic, reverse strand, flipped alleles (swapped so NEA matches REF)
        ("1:4000_T_C", 1, 4000, "C", "T", 0.35, 0.1, 0.01, 1e-5),  # NEA=T matches REF=T ✓ (EAF: 0.65→0.35)
        # Case 5: Palindromic A/T, forward strand (low MAF, NEA matches REF)
        ("1:5000_A_T", 1, 5000, "T", "A", 0.25, 0.1, 0.01, 1e-5),  # NEA=A matches REF=A ✓
        # Case 6: Palindromic G/C, forward strand (low MAF, NEA matches REF)
        ("1:6000_G_C", 1, 6000, "C", "G", 0.28, 0.1, 0.01, 1e-5),  # NEA=G matches REF=G ✓
        # Case 7: Palindromic A/T, reverse strand (high EAF, swapped so NEA matches REF)
        ("1:7000_A_T", 1, 7000, "T", "A", 0.25, 0.1, 0.01, 1e-5),  # NEA=A matches REF=A ✓ (EAF: 0.75→0.25)
        # Case 8: Palindromic G/C, reverse strand (high EAF, swapped so NEA matches REF)
        ("1:8000_G_C", 1, 8000, "C", "G", 0.28, 0.1, 0.01, 1e-5),  # NEA=G matches REF=G ✓ (EAF: 0.72→0.28)
        # Case 9: Palindromic with high MAF (swapped so NEA matches REF)
        ("1:9000_A_T", 1, 9000, "T", "A", 0.55, 0.1, 0.01, 1e-5),  # NEA=A matches REF=A ✓ (EAF: 0.45→0.55)
        # Case 10: Indel (NEA matches REF)
        ("1:10000_A_AT", 1, 10000, "AT", "A", 0.3, 0.1, 0.01, 1e-5),  # NEA=A matches REF=A ✓
        # Case 11: Variant not in VCF (NEA matches REF)
        ("1:20000_A_G", 1, 20000, "G", "A", 0.3, 0.1, 0.01, 1e-5),  # NEA=A matches REF=A ✓
    ]
    
    df = pd.DataFrame(sumstats, columns=[
        'SNPID', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P'
    ])
    
    df.to_csv(sumstats_file, sep='\t', index=False)
    print(f"Created sumstats with {len(df)} variants")
    return df

def run_gwaslab_harmonize(sumstats_file: str, vcf_file: str, fasta_file: str, output_file: str):
    """Run GWASLab harmonization."""
    print(f"\nRunning GWASLab harmonization...")
    
    try:
        import gwaslab as gl
        
        # Load sumstats
        mysumstats = gl.Sumstats(sumstats_file, 
                                 snpid="SNPID",
                                 chrom="CHR",
                                 pos="POS",
                                 ea="EA",
                                 nea="NEA",
                                 eaf="EAF",
                                 beta="BETA",
                                 se="SE",
                                 p="P")
        
        # Harmonize
        mysumstats.harmonize(
            ref_rsid_vcf=vcf_file,
            ref_seq=fasta_file,
            ref_infer=vcf_file,
            ref_alt_freq="AF_NFE",
            maf_threshold=0.42,
            ref_maf_threshold=0.42,
            remove=False,
            verbose=True
        )
        
        # Save results
        mysumstats.data.to_csv(output_file, sep='\t', index=False)
        print(f"GWASLab harmonization completed. Results saved to {output_file}")
        return mysumstats.data
        
    except Exception as e:
        print(f"Error running GWASLab: {e}")
        import traceback
        traceback.print_exc()
        return None

def run_harmoniser(sumstats_file: str, vcf_file: str, output_file: str, stats_file: str):
    """Run genetics-sumstat-harmoniser."""
    print(f"\nRunning genetics-sumstat-harmoniser...")
    
    harmoniser_dir = Path('/home/yunye/work/github/genetics-sumstat-harmoniser')
    harmoniser_path = harmoniser_dir / 'sumstat_harmoniser' / 'main.py'
    
    cmd = [
        'python', str(harmoniser_path),
        '--sumstats', sumstats_file,
        '--vcf', vcf_file,
        '--hm_sumstats', output_file,
        '--hm_statfile', stats_file,
        '--chrom_col', 'CHR',
        '--pos_col', 'POS',
        '--effAl_col', 'EA',
        '--otherAl_col', 'NEA',
        '--eaf_col', 'EAF',
        '--beta_col', 'BETA',
        '--palin_mode', 'infer',
        '--infer_maf_threshold', '0.42',
        '--af_vcf_field', 'AF_NFE',
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
            cwd=str(harmoniser_dir)
        )
        print("Harmoniser completed successfully")
        
        if os.path.exists(output_file):
            df = pd.read_csv(output_file, sep='\t')
            print(f"Loaded {len(df)} variants from harmoniser output")
            if len(df) > 0:
                print(f"Columns: {df.columns.tolist()}")
                print(f"hm_code distribution: {df['hm_code'].value_counts().to_dict()}")
            return df
        else:
            print(f"Error: Output file not found: {output_file}")
            return None
            
    except subprocess.CalledProcessError as e:
        print(f"Error running harmoniser: {e}")
        print(f"Stderr: {e.stderr}")
        return None

def compare_results(df_gwaslab: pd.DataFrame, df_harmoniser: pd.DataFrame):
    """Compare GWASLab STATUS codes with harmoniser hm_code."""
    print("\n" + "=" * 80)
    print("Comparing Results")
    print("=" * 80)
    
    # Create key for matching - use position-based matching (CHR:POS) as it's more reliable
    # Both tools should have CHR and POS, and this avoids SNPID format differences
    def create_position_key(chr_col, pos_col):
        """Create position-based key from CHR and POS columns."""
        if chr_col is None or pos_col is None:
            return pd.Series('', index=chr_col.index if chr_col is not None else pos_col.index)
        return chr_col.astype(str) + ':' + pos_col.astype(str)
    
    # Get CHR and POS from both dataframes
    gwaslab_chr = df_gwaslab.get('CHR', None)
    gwaslab_pos = df_gwaslab.get('POS', None)
    
    # Try multiple possible column names for harmoniser
    harmoniser_chr = None
    harmoniser_pos = None
    for chr_col in ['CHR', 'chrom', 'hm_chrom', 'chromosome']:
        if chr_col in df_harmoniser.columns:
            harmoniser_chr = df_harmoniser[chr_col]
            break
    for pos_col in ['POS', 'pos', 'hm_pos', 'position', 'bp', 'BP']:
        if pos_col in df_harmoniser.columns:
            harmoniser_pos = df_harmoniser[pos_col]
            break
    
    # Create position-based keys
    if gwaslab_chr is not None and gwaslab_pos is not None:
        df_gwaslab['key'] = create_position_key(gwaslab_chr, gwaslab_pos)
    else:
        print("Warning: GWASLab missing CHR or POS, using SNPID as key")
        df_gwaslab['key'] = df_gwaslab.get('SNPID', pd.Series('', index=df_gwaslab.index)).astype(str)
    
    if harmoniser_chr is not None and harmoniser_pos is not None:
        df_harmoniser['key'] = create_position_key(harmoniser_chr, harmoniser_pos)
    else:
        print("Warning: Harmoniser missing CHR or POS, using SNPID as key")
        if 'SNPID' in df_harmoniser.columns:
            df_harmoniser['key'] = df_harmoniser['SNPID'].astype(str)
        else:
            df_harmoniser['key'] = pd.Series('', index=df_harmoniser.index)
    
    # Map GWASLab STATUS to hm_code
    # For palindromic reverse strand variants (digit_7=2), check EA vs REF/ALT
    # The harmoniser checks: after reverse complement, if EA == REF → hm_code=4, else hm_code=3
    # For palindromic A/T or G/C, reverse complement is just swapping alleles
    def map_with_allele_check(row):
        status = row.get('STATUS', 0)
        digit_7 = extract_status_digit(status, 7)
        
        # For reverse strand palindromic (digit_7=2), use EA vs REF/ALT to determine correct mapping
        if digit_7 == 2:
            ea_gwaslab = row.get('EA', None)
            nea_gwaslab = row.get('NEA', None)
            # Get VCF REF/ALT from harmoniser output
            # After harmonization: NEA = VCF REF, EA = VCF ALT
            key = row.get('key', '')
            if key and key in df_harmoniser['key'].values:
                harm_row = df_harmoniser[df_harmoniser['key'] == key].iloc[0]
                vcf_ref = harm_row.get('hm_other_allele', harm_row.get('NEA', None))  # NEA after harmonization = VCF REF
                vcf_alt = harm_row.get('hm_effect_allele', harm_row.get('EA', None))  # EA after harmonization = VCF ALT
                if ea_gwaslab is not None and nea_gwaslab is not None and vcf_ref is not None:
                    # The harmoniser takes reverse complement (swap for palindromic), then checks if EA == REF
                    # After swap: EA becomes NEA, NEA becomes EA
                    # So if NEA (which becomes EA after swap) matches REF → hm_code=4
                    if nea_gwaslab == vcf_ref:
                        return 4
                    # If EA matches ALT (which becomes EA after swap) → hm_code=3
                    elif vcf_alt is not None and ea_gwaslab == vcf_alt:
                        return 3
        
        # For all other cases, use standard mapping
        return map_gwaslab_to_hm_code(status)
    
    df_gwaslab['mapped_hm_code'] = df_gwaslab.apply(map_with_allele_check, axis=1)
    
    # Debug: print keys to understand matching issue
    print(f"\nDebug - Matching keys (CHR:POS):")
    print(f"GWASLab keys (first 5): {df_gwaslab['key'].head().tolist()}")
    print(f"Harmoniser keys (first 5): {df_harmoniser['key'].head().tolist()}")
    print(f"GWASLab CHR:POS (first 5): {df_gwaslab[['CHR', 'POS']].head().values.tolist()}")
    if harmoniser_chr is not None and harmoniser_pos is not None:
        print(f"Harmoniser CHR:POS (first 5): {pd.DataFrame({'CHR': harmoniser_chr, 'POS': harmoniser_pos}).head().values.tolist()}")
    print(f"Unique GWASLab keys: {df_gwaslab['key'].nunique()} / {len(df_gwaslab)}")
    print(f"Unique Harmoniser keys: {df_harmoniser['key'].nunique()} / {len(df_harmoniser)}")
    print(f"Overlapping keys: {len(set(df_gwaslab['key']) & set(df_harmoniser['key']))}")
    
    # Merge - include all relevant columns from both sources
    # Get all columns from harmoniser (except key which we'll use for merging)
    harmoniser_cols = ['key', 'hm_code']
    # Add EAF column if available
    if 'hm_eaf' in df_harmoniser.columns:
        harmoniser_cols.append('hm_eaf')
    elif 'EAF' in df_harmoniser.columns:
        harmoniser_cols.append('EAF')
    # Add other harmoniser columns (hm_effect_allele, hm_other_allele, etc.)
    for col in df_harmoniser.columns:
        if col.startswith('hm_') and col not in harmoniser_cols:
            harmoniser_cols.append(col)
    # Add original columns that might be preserved
    for col in ['CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P']:
        if col in df_harmoniser.columns and col not in harmoniser_cols:
            harmoniser_cols.append(col)
    
    # Rename columns before merge to avoid conflicts
    df_gwaslab_merge = df_gwaslab.copy()
    # Rename GWASLab columns with _gwaslab suffix
    gwaslab_rename = {}
    for col in df_gwaslab_merge.columns:
        if col in ['EAF', 'BETA', 'SE', 'P', 'EA', 'NEA', 'CHR', 'POS', 'STATUS']:
            gwaslab_rename[col] = f'{col}_gwaslab'
    df_gwaslab_merge = df_gwaslab_merge.rename(columns=gwaslab_rename)
    
    df_harmoniser_merge = df_harmoniser[harmoniser_cols].copy()
    # Rename harmoniser columns with _harmoniser suffix
    harmoniser_rename = {}
    for col in df_harmoniser_merge.columns:
        if col not in ['key', 'hm_code']:
            if col.startswith('hm_'):
                harmoniser_rename[col] = col.replace('hm_', 'hm_')  # Keep hm_ prefix
            elif col in ['CHR', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P']:
                harmoniser_rename[col] = f'{col}_harmoniser'
    df_harmoniser_merge = df_harmoniser_merge.rename(columns=harmoniser_rename)
    
    merged = df_gwaslab_merge.merge(
        df_harmoniser_merge,
        on='key',
        how='outer'  # Use outer join to show all variants from both sides
    )
    
    # Create comparison dataframe with all relevant columns
    comparison_dict = {}
    
    # Basic identifiers
    if 'SNPID' in merged.columns:
        comparison_dict['SNPID'] = merged['SNPID']
    else:
        comparison_dict['SNPID'] = merged.index
    
    # CHR and POS - try both sources
    if 'CHR_gwaslab' in merged.columns:
        comparison_dict['CHR'] = merged['CHR_gwaslab']
    elif 'CHR_harmoniser' in merged.columns:
        comparison_dict['CHR'] = merged['CHR_harmoniser']
    else:
        comparison_dict['CHR'] = pd.Series('', index=merged.index)
    
    if 'POS_gwaslab' in merged.columns:
        comparison_dict['POS'] = merged['POS_gwaslab']
    elif 'POS_harmoniser' in merged.columns:
        comparison_dict['POS'] = merged['POS_harmoniser']
    else:
        comparison_dict['POS'] = pd.Series('', index=merged.index)
    
    # GWASLab columns (skip SE and P)
    for col in ['EA', 'NEA', 'EAF', 'BETA', 'STATUS']:
        col_name = f'{col}_gwaslab'
        if col_name in merged.columns:
            comparison_dict[col_name] = merged[col_name]
        else:
            comparison_dict[col_name] = pd.Series(None, index=merged.index, dtype=object)
    
    if 'mapped_hm_code' in merged.columns:
        comparison_dict['mapped_hm_code'] = merged['mapped_hm_code']
    else:
        comparison_dict['mapped_hm_code'] = pd.Series(None, index=merged.index, dtype=object)
    
    # Harmoniser columns - try multiple possible column names
    if 'hm_effect_allele' in merged.columns:
        comparison_dict['EA_harmoniser'] = merged['hm_effect_allele']
    elif 'EA_harmoniser' in merged.columns:
        comparison_dict['EA_harmoniser'] = merged['EA_harmoniser']
    elif 'EA' in merged.columns:
        comparison_dict['EA_harmoniser'] = merged['EA']
    else:
        comparison_dict['EA_harmoniser'] = pd.Series('', index=merged.index)
    
    if 'hm_other_allele' in merged.columns:
        comparison_dict['NEA_harmoniser'] = merged['hm_other_allele']
    elif 'NEA_harmoniser' in merged.columns:
        comparison_dict['NEA_harmoniser'] = merged['NEA_harmoniser']
    elif 'NEA' in merged.columns:
        comparison_dict['NEA_harmoniser'] = merged['NEA']
    else:
        comparison_dict['NEA_harmoniser'] = pd.Series('', index=merged.index)
    
    if 'hm_eaf' in merged.columns:
        comparison_dict['EAF_harmoniser'] = merged['hm_eaf']
    elif 'EAF_harmoniser' in merged.columns:
        comparison_dict['EAF_harmoniser'] = merged['EAF_harmoniser']
    elif 'EAF' in merged.columns:
        comparison_dict['EAF_harmoniser'] = merged['EAF']
    else:
        comparison_dict['EAF_harmoniser'] = pd.Series(None, index=merged.index, dtype=float)
    
    # Add other harmoniser stats if available (skip SE and P)
    for col in ['BETA']:
        col_name = f'{col}_harmoniser'
        hm_col = f'hm_{col.lower()}'
        if hm_col in merged.columns:
            comparison_dict[col_name] = merged[hm_col]
        elif col_name in merged.columns:
            comparison_dict[col_name] = merged[col_name]
        elif col in merged.columns:
            comparison_dict[col_name] = merged[col]
        else:
            comparison_dict[col_name] = pd.Series(None, index=merged.index, dtype=float)
    
    if 'hm_code' in merged.columns:
        comparison_dict['hm_code'] = merged['hm_code']
    else:
        comparison_dict['hm_code'] = pd.Series(None, index=merged.index, dtype=object)
    
    # Add result comparisons (EA, NEA, EAF, BETA)
    comparison_dict['EA_match'] = (
        comparison_dict.get('EA_gwaslab', pd.Series('', index=merged.index)) == 
        comparison_dict.get('EA_harmoniser', pd.Series('', index=merged.index))
    )
    comparison_dict['NEA_match'] = (
        comparison_dict.get('NEA_gwaslab', pd.Series('', index=merged.index)) == 
        comparison_dict.get('NEA_harmoniser', pd.Series('', index=merged.index))
    )
    
    # EAF comparison (with tolerance for floating point)
    eaf_gwaslab = comparison_dict.get('EAF_gwaslab', pd.Series(None, index=merged.index, dtype=float))
    eaf_harmoniser = comparison_dict.get('EAF_harmoniser', pd.Series(None, index=merged.index, dtype=float))
    if isinstance(eaf_gwaslab, pd.Series) and isinstance(eaf_harmoniser, pd.Series):
        eaf_match = (eaf_gwaslab.isna() & eaf_harmoniser.isna()) | (
            eaf_gwaslab.notna() & eaf_harmoniser.notna() & 
            (abs(eaf_gwaslab - eaf_harmoniser) < 0.001)
        )
    else:
        eaf_match = pd.Series(False, index=merged.index)
    comparison_dict['EAF_match'] = eaf_match
    
    # BETA comparison (with tolerance)
    beta_gwaslab = comparison_dict.get('BETA_gwaslab', pd.Series(None, index=merged.index, dtype=float))
    beta_harmoniser = comparison_dict.get('BETA_harmoniser', pd.Series(None, index=merged.index, dtype=float))
    if isinstance(beta_gwaslab, pd.Series) and isinstance(beta_harmoniser, pd.Series):
        beta_match = (beta_gwaslab.isna() & beta_harmoniser.isna()) | (
            beta_gwaslab.notna() & beta_harmoniser.notna() & 
            (abs(beta_gwaslab - beta_harmoniser) < 0.001)
        )
    else:
        beta_match = pd.Series(False, index=merged.index)
    comparison_dict['BETA_match'] = beta_match
    
    # Overall match (all key fields match)
    comparison_dict['results_match'] = (
        comparison_dict['EA_match'] & 
        comparison_dict['NEA_match'] & 
        comparison_dict['EAF_match'] &
        comparison_dict['BETA_match']
    )
    
    # Also keep hm_code comparison for reference
    mapped_hm = comparison_dict['mapped_hm_code']
    harmoniser_hm = comparison_dict['hm_code']
    if isinstance(mapped_hm, pd.Series) and isinstance(harmoniser_hm, pd.Series):
        comparison_dict['hm_code_match'] = (mapped_hm == harmoniser_hm) | (mapped_hm.isna() & harmoniser_hm.isna())
    else:
        comparison_dict['hm_code_match'] = mapped_hm == harmoniser_hm
    
    comparison = pd.DataFrame(comparison_dict)
    
    # Print results
    print("\n" + "=" * 120)
    print("Detailed Comparison - Full Results from Both Sources")
    print("=" * 120)
    
    # Set pandas display options for better formatting
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 15)
    
    print(comparison.to_string(index=False))
    
    # Statistics - compare actual results
    total = len(comparison)
    ea_matched = comparison['EA_match'].sum()
    nea_matched = comparison['NEA_match'].sum()
    eaf_matched = comparison['EAF_match'].sum()
    beta_matched = comparison['BETA_match'].sum()
    results_matched = comparison['results_match'].sum()
    hm_code_matched = comparison['hm_code_match'].sum()
    
    print(f"\n" + "=" * 120)
    print("Summary - Result Comparisons:")
    print(f"  Total variants: {total}")
    print(f"  EA matches: {ea_matched} ({ea_matched/total*100:.1f}%)")
    print(f"  NEA matches: {nea_matched} ({nea_matched/total*100:.1f}%)")
    print(f"  EAF matches: {eaf_matched} ({eaf_matched/total*100:.1f}%)")
    print(f"  BETA matches: {beta_matched} ({beta_matched/total*100:.1f}%)")
    print(f"  All results match (EA, NEA, EAF, BETA): {results_matched} ({results_matched/total*100:.1f}%)")
    print(f"  hm_code matches: {hm_code_matched} ({hm_code_matched/total*100:.1f}%)")
    print("=" * 120)
    
    # Show mismatches in results
    mismatches = comparison[~comparison['results_match']]
    if len(mismatches) > 0:
        print(f"\nMismatches (results differ - EA, NEA, EAF, or BETA):")
        print("=" * 120)
        # Show which fields don't match
        mismatch_details = mismatches.copy()
        mismatch_details['mismatch_fields'] = mismatch_details.apply(
            lambda row: ', '.join([
                'EA' if not row['EA_match'] else '',
                'NEA' if not row['NEA_match'] else '',
                'EAF' if not row['EAF_match'] else '',
                'BETA' if not row['BETA_match'] else ''
            ]).strip(', '),
            axis=1
        )
        print(mismatch_details[['SNPID', 'CHR', 'POS', 'EA_gwaslab', 'EA_harmoniser', 
                                'NEA_gwaslab', 'NEA_harmoniser', 'EAF_gwaslab', 'EAF_harmoniser',
                                'BETA_gwaslab', 'BETA_harmoniser', 'mismatch_fields']].to_string(index=False))
        print("=" * 120)
    
    return comparison

def main():
    """Main test function."""
    print("=" * 80)
    print("GWASLab vs Harmoniser Consistency Test")
    print("=" * 80)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create toy data
        sumstats_file = os.path.join(tmpdir, 'toy_sumstats.txt')
        vcf_file = os.path.join(tmpdir, 'toy.vcf')
        fasta_file = os.path.join(tmpdir, 'toy.fasta')
        
        # Create toy sumstats
        df_sumstats = create_toy_sumstats(sumstats_file)
        
        # Create toy VCF
        vcf_gz = create_toy_vcf(vcf_file)
        
        # Create FASTA with actual sequences at variant positions
        # We need sequences that match the VCF reference alleles
        with open(fasta_file, 'w') as f:
            f.write(">1\n")
            # Create sequence with actual bases at variant positions
            # Positions: 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000
            sequence = ['N'] * 25000
            # Set reference alleles at variant positions
            sequence[999] = 'A'   # pos 1000: ref A
            sequence[1999] = 'C'  # pos 2000: ref C
            sequence[2999] = 'G'  # pos 3000: ref G
            sequence[3999] = 'T'  # pos 4000: ref T
            sequence[4999] = 'A'  # pos 5000: ref A (palindromic)
            sequence[5999] = 'G'  # pos 6000: ref G (palindromic)
            sequence[6999] = 'A'  # pos 7000: ref A (palindromic)
            sequence[7999] = 'G'  # pos 8000: ref G (palindromic)
            sequence[8999] = 'A'  # pos 9000: ref A (palindromic)
            sequence[9999] = 'A'  # pos 10000: ref A (indel)
            sequence[19999] = 'A' # pos 20000: ref A (not in VCF)
            f.write(''.join(sequence) + '\n')
        
        # Run GWASLab
        gwaslab_output = os.path.join(tmpdir, 'gwaslab_output.txt')
        df_gwaslab = run_gwaslab_harmonize(sumstats_file, vcf_gz, fasta_file, gwaslab_output)
        
        # Run harmoniser
        harmoniser_output = os.path.join(tmpdir, 'harmoniser_output.txt')
        harmoniser_stats = os.path.join(tmpdir, 'harmoniser_stats.txt')
        df_harmoniser = run_harmoniser(sumstats_file, vcf_gz, harmoniser_output, harmoniser_stats)
        
        # Compare results
        if df_gwaslab is not None and df_harmoniser is not None:
            comparison = compare_results(df_gwaslab, df_harmoniser)
            
            # Save comparison
            output_file = 'test/output/gwaslab_harmoniser_toy_comparison.tsv'
            os.makedirs('test/output', exist_ok=True)
            comparison.to_csv(output_file, sep='\t', index=False)
            print(f"\nComparison saved to: {output_file}")
        else:
            print("\nError: Could not complete comparison (one or both tools failed)")
    
    print("\n" + "=" * 80)
    print("Test completed!")
    print("=" * 80)

if __name__ == "__main__":
    main()

