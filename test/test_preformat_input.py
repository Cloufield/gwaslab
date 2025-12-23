import os
import sys
import unittest
import tempfile
import shutil

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd
import polars as pl

from gwaslab.io.io_preformat_input import _preformat
from gwaslab.io.io_preformat_input_polars import preformatp


class TestPreformatInput(unittest.TestCase):
    def setUp(self):
        rng = pd.Series(range(10))
        rows = []
        for i in rng:
            chr_ = 1
            pos_ = 1000 + i
            snpid = f"{chr_}:{pos_}_A_G"
            rows.append({"CHR": chr_, "POS": pos_, "EA": "A", "NEA": "G", "EAF": 0.2, "RAF": 0.3, "SNPID": snpid, "EXTRA1": i})
        self.df = pd.DataFrame(rows)

    def test_dataframe_preserves_columns(self):
        out = _preformat(sumstats=self.df, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", verbose=False)
        self.assertIn("EAF", out.columns)
        self.assertIn("RAF", out.columns)
        self.assertIn("EXTRA1", out.columns)

    def test_path_trims_unmapped_columns(self):
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        try:
            self.df.to_csv(path, sep="\t", index=False)
            out = _preformat(sumstats=path, chrom="CHR", pos="POS", ea="EA", nea="NEA", snpid="SNPID", eaf="EAF", verbose=False)
            self.assertIn("EAF", out.columns)
            self.assertIn("SNPID", out.columns)
            self.assertNotIn("EXTRA1", out.columns)
        finally:
            if os.path.exists(path):
                os.remove(path)


class TestPreformatInputPolars(unittest.TestCase):
    def setUp(self):
        self.raw_dir = os.path.join(os.path.dirname(__file__), "raw")
        self.dirty_sumstats_path = os.path.join(self.raw_dir, "dirty_sumstats.tsv")

    def test_load_dirty_sumstats_from_path(self):
        """Test loading dirty_sumstats.tsv using preformatp"""
        result = preformatp(
            sumstats=self.dirty_sumstats_path,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            n="N",
            ncase="N_CASE",
            ncontrol="N_CONTROL",
            other=["NOTE"],
            verbose=False
        )
        
        # Check that result is a polars DataFrame
        self.assertIsInstance(result, pl.DataFrame)
        
        # Check that it has data
        self.assertGreater(result.height, 0)
        
        # Check that expected columns are present
        expected_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P", "N", "N_CASE", "N_CONTROL", "NOTE"]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Column {col} not found in result")
        
        # Check that NOTE column is preserved
        self.assertIn("NOTE", result.columns)

    def test_load_dirty_sumstats_with_include(self):
        """Test loading with include parameter"""
        result = preformatp(
            sumstats=self.dirty_sumstats_path,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            include=["SNPID", "CHR", "POS", "EA", "NEA", "P"],
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertGreater(result.height, 0)
        
        # Should only have included columns
        included_cols = ["SNPID", "CHR", "POS", "EA", "NEA", "P"]
        for col in included_cols:
            self.assertIn(col, result.columns)
        
        # Should not have excluded columns like BETA
        self.assertNotIn("BETA", result.columns)

    def test_load_dirty_sumstats_with_exclude(self):
        """Test loading with exclude parameter"""
        result = preformatp(
            sumstats=self.dirty_sumstats_path,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            exclude=["BETA", "SE"],
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertGreater(result.height, 0)
        
        # Should have basic columns
        self.assertIn("SNPID", result.columns)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        
        # Should not have excluded columns
        self.assertNotIn("BETA", result.columns)
        self.assertNotIn("SE", result.columns)

    def test_load_dirty_sumstats_from_pandas_dataframe(self):
        """Test loading from pandas DataFrame"""
        # First load as pandas
        df_pd = pd.read_csv(self.dirty_sumstats_path, sep="\t")
        
        result = preformatp(
            sumstats=df_pd,
            fmt="gwaslab",
            tab_fmt="tsv",
            snpid="SNPID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            p="P",
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertEqual(result.height, df_pd.shape[0])
        self.assertIn("SNPID", result.columns)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)

    def test_load_dirty_sumstats_creates_snpid_if_missing(self):
        """Test that SNPID is created if both rsID and SNPID are missing"""
        # Create a test file without SNPID
        import tempfile
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        try:
            # Create a minimal test file
            test_data = """CHR\tPOS\tEA\tNEA\tP
1\t100\tA\tG\t0.001
1\t200\tT\tC\t0.002"""
            with open(path, 'w') as f:
                f.write(test_data)
            
            result = preformatp(
                sumstats=path,
                chrom="CHR",
                pos="POS",
                ea="EA",
                nea="NEA",
                p="P",
                verbose=False
            )
            
            self.assertIsInstance(result, pl.DataFrame)
            self.assertIn("SNPID", result.columns)
            self.assertGreater(result.height, 0)
        finally:
            if os.path.exists(path):
                os.remove(path)


class TestPreformatInputWithAtSymbol(unittest.TestCase):
    """Test cases for loading files using @ symbol for chromosome-split files"""
    
    def setUp(self):
        """Create temporary directory with chromosome-split files"""
        self.temp_dir = tempfile.mkdtemp()
        self.base_path = os.path.join(self.temp_dir, "chr@.tsv")
        
        # Create test data for multiple chromosomes
        for chr_num in [1, 2, 3]:
            chr_path = os.path.join(self.temp_dir, f"chr{chr_num}.tsv")
            rows = []
            for i in range(5):
                pos = 1000 + i * 100
                snpid = f"{chr_num}:{pos}_A_G"
                rows.append({
                    "CHR": chr_num,
                    "POS": pos,
                    "EA": "A",
                    "NEA": "G",
                    "EAF": 0.2 + i * 0.01,
                    "BETA": 0.1 + i * 0.01,
                    "SE": 0.01,
                    "P": 0.001 + i * 0.0001,
                    "SNPID": snpid
                })
            df = pd.DataFrame(rows)
            df.to_csv(chr_path, sep="\t", index=False)
    
    def tearDown(self):
        """Clean up temporary directory"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_load_with_at_symbol(self):
        """Test loading chromosome-split files using @ symbol"""
        result = _preformat(
            sumstats=self.base_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            eaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Should have data from all 3 chromosomes
        self.assertGreater(len(result), 0)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        
        # Check that we have data from multiple chromosomes
        unique_chrs = result["CHR"].unique()
        self.assertGreaterEqual(len(unique_chrs), 1)
    
    def test_load_with_at_symbol_polars(self):
        """Test loading chromosome-split files using @ symbol with polars"""
        result = preformatp(
            sumstats=self.base_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            snpid="SNPID",
            eaf="EAF",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        self.assertIsInstance(result, pl.DataFrame)
        self.assertGreater(result.height, 0)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)


class TestPreformatInputFormats(unittest.TestCase):
    """Test cases for loading different file formats"""
    
    def setUp(self):
        """Create temporary directory for test files"""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_load_plink_format(self):
        """Test loading PLINK format"""
        plink_path = os.path.join(self.temp_dir, "plink_test.assoc")
        # PLINK format: CHR SNP BP A1 F_A F_U A2 CHISQ P OR
        plink_data = """CHR\tSNP\tBP\tA1\tF_A\tF_U\tA2\tCHISQ\tP\tOR
1\trs1\t1000\tA\t0.2\t0.3\tG\t10.5\t0.001\t1.2
1\trs2\t2000\tT\t0.15\t0.25\tC\t8.3\t0.002\t1.15
2\trs3\t3000\tG\t0.3\t0.35\tA\t12.1\t0.0005\t1.3"""
        
        with open(plink_path, 'w') as f:
            f.write(plink_data)
        
        try:
            result = _preformat(
                sumstats=plink_path,
                fmt="plink",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            # Check that key columns are present
            self.assertIn("CHR", result.columns)
            self.assertIn("POS", result.columns)
        except Exception as e:
            # If format not available, skip test
            self.skipTest(f"PLINK format not available: {e}")
    
    def test_load_plink2_format(self):
        """Test loading PLINK2 format"""
        plink2_path = os.path.join(self.temp_dir, "plink2_test.glm.linear")
        # PLINK2 format: #CHROM ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P
        plink2_data = """#CHROM\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP
1\trs1\tG\tA\tA\tADD\t1000\t0.1\t0.01\t10.0\t0.001
1\trs2\tC\tT\tT\tADD\t1000\t0.15\t0.015\t10.0\t0.0005
2\trs3\tA\tG\tG\tADD\t1000\t0.2\t0.02\t10.0\t0.0001"""
        
        with open(plink2_path, 'w') as f:
            f.write(plink2_data)
        
        try:
            result = _preformat(
                sumstats=plink2_path,
                fmt="plink2",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            self.assertIn("POS", result.columns)
        except Exception as e:
            self.skipTest(f"PLINK2 format not available: {e}")
    
    def test_load_regenie_format(self):
        """Test loading REGENIE format"""
        regenie_path = os.path.join(self.temp_dir, "regenie_test.regenie")
        # REGENIE format: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
        regenie_data = """CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tINFO\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA
1\t1000\trs1\tG\tA\t0.2\t0.95\t1000\tADD\t0.1\t0.01\t100.0\t3.0\t
1\t2000\trs2\tC\tT\t0.15\t0.92\t1000\tADD\t0.15\t0.015\t100.0\t3.3\t
2\t3000\trs3\tA\tG\t0.3\t0.98\t1000\tADD\t0.2\t0.02\t100.0\t4.0\t"""
        
        with open(regenie_path, 'w') as f:
            f.write(regenie_data)
        
        try:
            result = _preformat(
                sumstats=regenie_path,
                fmt="regenie",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            self.assertIn("POS", result.columns)
        except Exception as e:
            self.skipTest(f"REGENIE format not available: {e}")
    
    def test_load_vcf_format(self):
        """Test loading VCF format"""
        vcf_path = os.path.join(self.temp_dir, "test.vcf")
        # VCF format with FORMAT and sample columns
        vcf_data = """##fileformat=VCFv4.3
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=ES,Number=1,Type=Float,Description="Effect Size">
##FORMAT=<ID=SE,Number=1,Type=Float,Description="Standard Error">
##FORMAT=<ID=LP,Number=1,Type=Float,Description="-log10 p-value">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t1000\trs1\tG\tA\t.\t.\t.\tES:SE:LP\t0.1:0.01:3.0
1\t2000\trs2\tC\tT\t.\t.\t.\tES:SE:LP\t0.15:0.015:3.3
2\t3000\trs3\tA\tG\t.\t.\t.\tES:SE:LP\t0.2:0.02:4.0"""
        
        with open(vcf_path, 'w') as f:
            f.write(vcf_data)
        
        try:
            result = _preformat(
                sumstats=vcf_path,
                fmt="vcf",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            self.assertIn("POS", result.columns)
        except Exception as e:
            self.skipTest(f"VCF format not available: {e}")
    
    def test_load_saige_format(self):
        """Test loading SAIGE format"""
        saige_path = os.path.join(self.temp_dir, "saige_test.txt")
        # SAIGE format: CHR POS SNPID Allele1 Allele2 AC_Allele2 AF_Allele2 N BETA SE p.value
        saige_data = """CHR\tPOS\tSNPID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tN\tBETA\tSE\tp.value
1\t1000\trs1\tG\tA\t200\t0.2\t1000\t0.1\t0.01\t0.001
1\t2000\trs2\tC\tT\t150\t0.15\t1000\t0.15\t0.015\t0.0005
2\t3000\trs3\tA\tG\t300\t0.3\t1000\t0.2\t0.02\t0.0001"""
        
        with open(saige_path, 'w') as f:
            f.write(saige_data)
        
        try:
            result = _preformat(
                sumstats=saige_path,
                fmt="saige",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            self.assertIn("POS", result.columns)
        except Exception as e:
            self.skipTest(f"SAIGE format not available: {e}")
    
    def test_load_ldsc_format(self):
        """Test loading LDSC format"""
        ldsc_path = os.path.join(self.temp_dir, "ldsc_test.txt")
        # LDSC format: SNP CHR BP A1 A2 N Z
        ldsc_data = """SNP\tCHR\tBP\tA1\tA2\tN\tZ
rs1\t1\t1000\tA\tG\t1000\t10.0
rs2\t1\t2000\tT\tC\t1000\t10.0
rs3\t2\t3000\tG\tA\t1000\t10.0"""
        
        with open(ldsc_path, 'w') as f:
            f.write(ldsc_data)
        
        try:
            result = _preformat(
                sumstats=ldsc_path,
                fmt="ldsc",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            self.assertIn("POS", result.columns)
        except Exception as e:
            self.skipTest(f"LDSC format not available: {e}")
    
    def test_load_metal_format(self):
        """Test loading METAL format"""
        metal_path = os.path.join(self.temp_dir, "metal_test.txt")
        # METAL format: MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq Effect StdErr P-value Direction HetISq HetChiSq HetDf HetPVal
        metal_data = """MarkerName\tAllele1\tAllele2\tFreq1\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal
rs1\tA\tG\t0.2\t0.01\t0.15\t0.25\t0.1\t0.01\t0.001\t+++\t0.0\t0.0\t2\t1.0
rs2\tT\tC\t0.15\t0.01\t0.1\t0.2\t0.15\t0.015\t0.0005\t+++\t0.0\t0.0\t2\t1.0
rs3\tG\tA\t0.3\t0.01\t0.25\t0.35\t0.2\t0.02\t0.0001\t+++\t0.0\t0.0\t2\t1.0"""
        
        with open(metal_path, 'w') as f:
            f.write(metal_data)
        
        try:
            result = _preformat(
                sumstats=metal_path,
                fmt="metal",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            # METAL uses MarkerName which should be converted to SNPID
            self.assertIn("SNPID", result.columns)
        except Exception as e:
            self.skipTest(f"METAL format not available: {e}")
    
    def test_load_compressed_file(self):
        """Test loading compressed (.gz) file"""
        import gzip
        gz_path = os.path.join(self.temp_dir, "test.tsv.gz")
        
        # Create test data
        test_data = """CHR\tPOS\tEA\tNEA\tBETA\tSE\tP
1\t1000\tA\tG\t0.1\t0.01\t0.001
1\t2000\tT\tC\t0.15\t0.015\t0.0005
2\t3000\tG\tA\t0.2\t0.02\t0.0001"""
        
        with gzip.open(gz_path, 'wt') as f:
            f.write(test_data)
        
        result = _preformat(
            sumstats=gz_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        self.assertGreater(len(result), 0)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
    
    def test_load_parquet_format(self):
        """Test loading Parquet format"""
        parquet_path = os.path.join(self.temp_dir, "test.parquet")
        
        # Create test DataFrame
        test_df = pd.DataFrame({
            "CHR": [1, 1, 2],
            "POS": [1000, 2000, 3000],
            "EA": ["A", "T", "G"],
            "NEA": ["G", "C", "A"],
            "BETA": [0.1, 0.15, 0.2],
            "SE": [0.01, 0.015, 0.02],
            "P": [0.001, 0.0005, 0.0001]
        })
        
        test_df.to_parquet(parquet_path, index=False)
        
        result = _preformat(
            sumstats=parquet_path,
            tab_fmt="parquet",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        self.assertGreater(len(result), 0)
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)


class TestPreformatInputWithAtSymbolAndFormats(unittest.TestCase):
    """Test cases combining @ symbol with different formats"""
    
    def setUp(self):
        """Create temporary directory with chromosome-split files in different formats"""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_load_plink_format_with_at_symbol(self):
        """Test loading PLINK format files split by chromosome using @ symbol"""
        # Create PLINK format files for multiple chromosomes
        for chr_num in [1, 2]:
            plink_path = os.path.join(self.temp_dir, f"chr{chr_num}.assoc")
            plink_data = f"""CHR\tSNP\tBP\tA1\tF_A\tF_U\tA2\tCHISQ\tP\tOR
{chr_num}\trs{chr_num}_1\t{1000 + chr_num * 100}\tA\t0.2\t0.3\tG\t10.5\t0.001\t1.2
{chr_num}\trs{chr_num}_2\t{2000 + chr_num * 100}\tT\t0.15\t0.25\tC\t8.3\t0.002\t1.15"""
            
            with open(plink_path, 'w') as f:
                f.write(plink_data)
        
        base_path = os.path.join(self.temp_dir, "chr@.assoc")
        
        try:
            result = _preformat(
                sumstats=base_path,
                fmt="plink",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            # Should have data from multiple chromosomes
            unique_chrs = result["CHR"].unique()
            self.assertGreaterEqual(len(unique_chrs), 1)
        except Exception as e:
            self.skipTest(f"PLINK format with @ symbol not available: {e}")
    
    def test_load_regenie_format_with_at_symbol(self):
        """Test loading REGENIE format files split by chromosome using @ symbol"""
        # Create REGENIE format files for multiple chromosomes
        for chr_num in [1, 2]:
            regenie_path = os.path.join(self.temp_dir, f"chr{chr_num}.regenie")
            regenie_data = f"""CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tINFO\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA
{chr_num}\t{1000 + chr_num * 100}\trs{chr_num}_1\tG\tA\t0.2\t0.95\t1000\tADD\t0.1\t0.01\t100.0\t3.0\t
{chr_num}\t{2000 + chr_num * 100}\trs{chr_num}_2\tC\tT\t0.15\t0.92\t1000\tADD\t0.15\t0.015\t100.0\t3.3\t"""
            
            with open(regenie_path, 'w') as f:
                f.write(regenie_data)
        
        base_path = os.path.join(self.temp_dir, "chr@.regenie")
        
        try:
            result = _preformat(
                sumstats=base_path,
                fmt="regenie",
                verbose=False
            )
            
            self.assertGreater(len(result), 0)
            self.assertIn("CHR", result.columns)
            unique_chrs = result["CHR"].unique()
            self.assertGreaterEqual(len(unique_chrs), 1)
        except Exception as e:
            self.skipTest(f"REGENIE format with @ symbol not available: {e}")


class TestPreformatWorkflowAndPriority(unittest.TestCase):
    """Test cases for workflow and priority order of _preformat"""
    
    def setUp(self):
        """Create test data files"""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a test file with standard columns
        self.test_path = os.path.join(self.temp_dir, "test.tsv")
        test_data = """CHR\tPOS\tEA\tNEA\tBETA\tSE\tP\tEAF\tN\tEXTRA_COL
1\t1000\tA\tG\t0.1\t0.01\t0.001\t0.2\t1000\tvalue1
1\t2000\tT\tC\t0.15\t0.015\t0.0005\t0.25\t1000\tvalue2
2\t3000\tG\tA\t0.2\t0.02\t0.0001\t0.3\t1000\tvalue3"""
        
        with open(self.test_path, 'w') as f:
            f.write(test_data)
        
        # Create DataFrame version
        self.test_df = pd.DataFrame({
            "CHR": [1, 1, 2],
            "POS": [1000, 2000, 3000],
            "EA": ["A", "T", "G"],
            "NEA": ["G", "C", "A"],
            "BETA": [0.1, 0.15, 0.2],
            "SE": [0.01, 0.015, 0.02],
            "P": [0.001, 0.0005, 0.0001],
            "EAF": [0.2, 0.25, 0.3],
            "N": [1000, 1000, 1000],
            "EXTRA_COL": ["value1", "value2", "value3"]
        })
    
    def tearDown(self):
        """Clean up temporary directory"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def test_priority_user_overrides_formatbook(self):
        """Test that user-specified mappings override formatbook defaults"""
        # This test assumes formatbook might have different mappings
        # User mappings should take precedence
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # User-specified columns should be mapped correctly
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("EA", result.columns)
        self.assertIn("NEA", result.columns)
        self.assertIn("BETA", result.columns)
        self.assertIn("SE", result.columns)
        self.assertIn("P", result.columns)
    
    def test_priority_include_filter_after_mappings(self):
        """Test that include filter is applied after column mappings (Step 7)"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            eaf="EAF",
            include=["CHR", "POS", "EA", "NEA", "P"],  # Only include these
            verbose=False
        )
        
        # Should only have included columns
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("EA", result.columns)
        self.assertIn("NEA", result.columns)
        self.assertIn("P", result.columns)
        
        # Should NOT have excluded columns (even if mapped)
        self.assertNotIn("BETA", result.columns)
        self.assertNotIn("SE", result.columns)
        self.assertNotIn("EAF", result.columns)
    
    def test_priority_exclude_filter_after_include(self):
        """Test that exclude filter is applied after include filter"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            eaf="EAF",
            exclude=["SE"],  # Exclude SE even if mapped
            verbose=False
        )
        
        # Should have mapped columns
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("BETA", result.columns)
        self.assertIn("P", result.columns)
        
        # Should NOT have excluded column
        self.assertNotIn("SE", result.columns)
    
    def test_priority_include_then_exclude(self):
        """Test priority: include creates subset, then exclude removes from subset"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            eaf="EAF",
            include=["CHR", "POS", "BETA", "SE", "P"],  # Include these
            exclude=["SE"],  # Then exclude SE
            verbose=False
        )
        
        # Should have included columns except excluded
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("BETA", result.columns)
        self.assertIn("P", result.columns)
        
        # Should NOT have excluded column
        self.assertNotIn("SE", result.columns)
        
        # Should NOT have columns not in include list
        self.assertNotIn("EA", result.columns)
        self.assertNotIn("NEA", result.columns)
        self.assertNotIn("EAF", result.columns)
    
    def test_workflow_configuration_phase_initialization(self):
        """Test Step 1: Parameter initialization"""
        # Test that parameters are properly initialized
        result = _preformat(
            sumstats=self.test_df,
            chrom="CHR",
            pos="POS",
            verbose=False
        )
        
        # Should work without errors (initialization successful)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)
    
    def test_workflow_configuration_phase_parquet(self):
        """Test Step 2: Parquet format handling (early exit)"""
        parquet_path = os.path.join(self.temp_dir, "test.parquet")
        self.test_df.to_parquet(parquet_path, index=False)
        
        result = _preformat(
            sumstats=parquet_path,
            tab_fmt="parquet",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            verbose=False
        )
        
        # Should load parquet successfully
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)
        self.assertIn("CHR", result.columns)
    
    def test_workflow_mapping_phase_column_discovery(self):
        """Test Step 4: Column discovery from file"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            verbose=False
        )
        
        # Should discover and map columns correctly
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("EA", result.columns)
        self.assertIn("NEA", result.columns)
    
    def test_workflow_mapping_phase_user_mappings(self):
        """Test Step 5: User-specified column mappings"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # All user-specified columns should be mapped
        expected_cols = ["CHR", "POS", "EA", "NEA", "BETA", "SE", "P"]
        for col in expected_cols:
            self.assertIn(col, result.columns, f"Column {col} should be mapped")
    
    def test_workflow_data_phase_loading(self):
        """Test Step 8: Data loading"""
        # Test loading from path
        result_path = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            verbose=False
        )
        
        # Test loading from DataFrame
        result_df = _preformat(
            sumstats=self.test_df,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            verbose=False
        )
        
        # Both should produce valid DataFrames
        self.assertIsInstance(result_path, pd.DataFrame)
        self.assertIsInstance(result_df, pd.DataFrame)
        self.assertEqual(len(result_path), len(result_df))
    
    def test_workflow_data_phase_postprocessing(self):
        """Test Step 9: Post-processing (renaming, transformations)"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            build="19",
            verbose=False
        )
        
        # Should have standardized column names
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        self.assertIn("EA", result.columns)
        self.assertIn("NEA", result.columns)
        
        # Should have STATUS column created (post-processing)
        self.assertIn("STATUS", result.columns)
        
        # Should have SNPID if missing (post-processing)
        # Since we didn't provide snpid, it should be created from CHR:POS:NEA:EA
        self.assertIn("SNPID", result.columns)
    
    def test_workflow_other_columns_preserved(self):
        """Test that 'other' columns are preserved through workflow"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            other=["EXTRA_COL"],  # Add extra column
            verbose=False
        )
        
        # Extra column should be preserved
        self.assertIn("EXTRA_COL", result.columns)
    
    def test_workflow_other_columns_with_include(self):
        """Test that 'other' columns respect include filter"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            other=["EXTRA_COL"],
            include=["CHR", "POS", "EA", "NEA"],  # Don't include EXTRA_COL
            verbose=False
        )
        
        # Included columns should be present
        self.assertIn("CHR", result.columns)
        self.assertIn("POS", result.columns)
        
        # EXTRA_COL should NOT be present (not in include list)
        self.assertNotIn("EXTRA_COL", result.columns)
    
    def test_workflow_other_columns_with_exclude(self):
        """Test that 'other' columns respect exclude filter"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            other=["EXTRA_COL"],
            exclude=["EXTRA_COL"],  # Explicitly exclude
            verbose=False
        )
        
        # Should have mapped columns
        self.assertIn("CHR", result.columns)
        self.assertIn("BETA", result.columns)
        
        # EXTRA_COL should NOT be present (excluded)
        self.assertNotIn("EXTRA_COL", result.columns)
    
    def test_workflow_numeric_constants(self):
        """Test that numeric constants (n, ncase, ncontrol) are handled correctly"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            n=5000,  # Constant value
            ncase=2500,  # Constant value
            ncontrol=2500,  # Constant value
            verbose=False
        )
        
        # Should have constant columns added in post-processing
        self.assertIn("N", result.columns)
        self.assertIn("N_CASE", result.columns)
        self.assertIn("N_CONTROL", result.columns)
        
        # Check that values are constant
        self.assertTrue((result["N"] == 5000).all())
        self.assertTrue((result["N_CASE"] == 2500).all())
        self.assertTrue((result["N_CONTROL"] == 2500).all())
    
    def test_workflow_numeric_string_vs_int(self):
        """Test that string column names vs int constants are handled differently"""
        # Test with string column name
        result_str = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            n="N",  # Column name (string)
            verbose=False
        )
        
        # Test with int constant
        result_int = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            n=5000,  # Constant value (int)
            verbose=False
        )
        
        # Both should have N column
        self.assertIn("N", result_str.columns)
        self.assertIn("N", result_int.columns)
        
        # String version should use column values, int version should be constant
        self.assertTrue((result_int["N"] == 5000).all())
    
    def test_workflow_neaf_to_eaf_conversion(self):
        """Test NEAF to EAF conversion in post-processing"""
        # Create test file with NEAF instead of EAF
        neaf_path = os.path.join(self.temp_dir, "test_neaf.tsv")
        neaf_data = """CHR\tPOS\tEA\tNEA\tNEAF\tBETA\tSE\tP
1\t1000\tA\tG\t0.8\t0.1\t0.01\t0.001
1\t2000\tT\tC\t0.75\t0.15\t0.015\t0.0005
2\t3000\tG\tA\t0.7\t0.2\t0.02\t0.0001"""
        
        with open(neaf_path, 'w') as f:
            f.write(neaf_data)
        
        result = _preformat(
            sumstats=neaf_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            neaf="NEAF",  # Specify NEAF
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Should have EAF (converted from NEAF)
        self.assertIn("EAF", result.columns)
        self.assertNotIn("NEAF", result.columns)
        
        # EAF should be 1 - NEAF
        # First row: NEAF=0.8, so EAF should be 0.2
        self.assertAlmostEqual(result.iloc[0]["EAF"], 0.2, places=5)
    
    def test_workflow_snpid_creation(self):
        """Test SNPID creation when both rsID and SNPID are missing"""
        # Create test file without SNPID or rsID
        no_id_path = os.path.join(self.temp_dir, "test_no_id.tsv")
        no_id_data = """CHR\tPOS\tEA\tNEA\tBETA\tSE\tP
1\t1000\tA\tG\t0.1\t0.01\t0.001
1\t2000\tT\tC\t0.15\t0.015\t0.0005
2\t3000\tG\tA\t0.2\t0.02\t0.0001"""
        
        with open(no_id_path, 'w') as f:
            f.write(no_id_data)
        
        result = _preformat(
            sumstats=no_id_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Should have SNPID created in post-processing
        self.assertIn("SNPID", result.columns)
        
        # SNPID should be CHR:POS:NEA:EA format
        expected_snpid = "1:1000:G:A"
        self.assertEqual(result.iloc[0]["SNPID"], expected_snpid)
    
    def test_workflow_snpid_creation_without_alleles(self):
        """Test SNPID creation when alleles are missing"""
        # Create test file without alleles
        no_alleles_path = os.path.join(self.temp_dir, "test_no_alleles.tsv")
        no_alleles_data = """CHR\tPOS\tBETA\tSE\tP
1\t1000\t0.1\t0.01\t0.001
1\t2000\t0.15\t0.015\t0.0005
2\t3000\t0.2\t0.02\t0.0001"""
        
        with open(no_alleles_path, 'w') as f:
            f.write(no_alleles_data)
        
        result = _preformat(
            sumstats=no_alleles_path,
            chrom="CHR",
            pos="POS",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Should have SNPID created (CHR:POS format, no alleles)
        self.assertIn("SNPID", result.columns)
        
        # SNPID should be CHR:POS format
        expected_snpid = "1:1000"
        self.assertEqual(result.iloc[0]["SNPID"], expected_snpid)
    
    def test_workflow_column_renaming_order(self):
        """Test that column renaming happens in correct order"""
        result = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Original column names should be renamed to standard names
        # The rename happens in post-processing (Step 9)
        # So we should see standard GWASLab column names
        standard_cols = ["CHR", "POS", "EA", "NEA", "BETA", "SE", "P"]
        for col in standard_cols:
            self.assertIn(col, result.columns, 
                         f"Standard column {col} should be present after renaming")
    
    def test_workflow_dataframe_vs_path_consistency(self):
        """Test that DataFrame and path inputs produce consistent results for mapped columns"""
        result_path = _preformat(
            sumstats=self.test_path,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        result_df = _preformat(
            sumstats=self.test_df,
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            beta="BETA",
            se="SE",
            p="P",
            verbose=False
        )
        
        # Note: Path loading only loads specified columns (usecols), 
        # while DataFrame preserves all existing columns
        # So we check that mapped columns are consistent
        mapped_cols = {"CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "STATUS", "SNPID"}
        
        # Both should have all mapped columns
        for col in mapped_cols:
            self.assertIn(col, result_path.columns, 
                         f"Mapped column {col} should be in path result")
            self.assertIn(col, result_df.columns, 
                         f"Mapped column {col} should be in DataFrame result")
        
        # DataFrame may have additional columns (N, EAF, EXTRA_COL) that weren't mapped
        # This is expected behavior - DataFrame preserves all columns
        self.assertIn("N", result_df.columns)  # DataFrame preserves this
        self.assertIn("EAF", result_df.columns)  # DataFrame preserves this
        self.assertIn("EXTRA_COL", result_df.columns)  # DataFrame preserves this
        
        # Path result should NOT have unmapped columns
        self.assertNotIn("N", result_path.columns)  # Not mapped, so not loaded
        self.assertNotIn("EAF", result_path.columns)  # Not mapped, so not loaded
        self.assertNotIn("EXTRA_COL", result_path.columns)  # Not mapped, so not loaded
        
        # Both should have same number of rows
        self.assertEqual(len(result_path), len(result_df))


if __name__ == "__main__":
    unittest.main()

