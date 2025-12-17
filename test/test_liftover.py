
import unittest
import pandas as pd
import numpy as np
import os
import sys
from unittest.mock import MagicMock, patch
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.g_Sumstats import Sumstats
from gwaslab.qc.qc_fix_sumstats import parallelizeliftovervariant

# Mock pyliftover to avoid installing/downloading chain files
class MockConverter:
    def __init__(self, from_build, to_build):
        self.from_build = from_build
        self.to_build = to_build

    def convert_coordinate(self, chrom, pos):
        # Mock conversion: just add 100 to pos
        return [(chrom, pos + 100, 1)]

    def query(self, chrom, pos):
         return self.convert_coordinate(chrom, pos)

    def __getitem__(self, key):
         # Support converter[chrom][pos] syntax which might be used in some implementations
         # This is a simplification. The actual pyliftover object might behave differently.
         # Based on gwaslab code: results = converter[chrom][pos_0_based]
         return MockConverterChrom(key)

class MockConverterChrom:
    def __init__(self, chrom):
        self.chrom = chrom
    def __getitem__(self, pos):
        # Return a list of tuples (chrom, pos, strand)
        return [(self.chrom, pos + 101, "+")]

class MockPool:
    def __init__(self, processes=1):
        pass
    
    def map(self, func, iterable):
        # Synchronous map
        return [func(item) for item in iterable]
    
    def starmap(self, func, iterable):
        return [func(*item) for item in iterable]

    def close(self):
        pass
        
    def join(self):
        pass

class TestLiftover(unittest.TestCase):

    def setUp(self):
        # Create a sample dataframe
        self.data = pd.DataFrame({
            'SNPID': ['rs1', 'rs2', 'rs3'],
            'CHR': ['1', '2', '3'],
            'POS': [1000, 2000, 3000],
            'EA': ['A', 'G', 'C'],
            'NEA': ['G', 'A', 'T'],
            'BETA': [0.1, 0.2, 0.3],
            'SE': [0.01, 0.02, 0.03],
            'P': [0.001, 0.002, 0.003],
            'STATUS': ['OK', 'OK', 'OK']
        })
        self.sumstats = Sumstats(self.data, fmt="gwaslab")

    @patch('gwaslab.qc.qc_fix_sumstats.get_lifter')
    @patch('gwaslab.qc.qc_fix_sumstats.ChainFile') 
    def test_liftover_basic(self, mock_chainfile, mock_get_lifter):
        # Mock the converter
        mock_converter = MockConverter('hg19', 'hg38')
        mock_get_lifter.return_value = mock_converter
        mock_chainfile.return_value = mock_converter 

        pass

    def test_categorical_columns_fix(self):
        # This test specifically targets the fix we implemented.
        # We want to ensure that if CHR or STATUS are categorical, they are converted to object/str
        # and the liftover proceeds (or at least the assignment doesn't fail).
        
        # Set columns to categorical
        self.sumstats.data['CHR'] = self.sumstats.data['CHR'].astype('category')
        self.sumstats.data['STATUS'] = self.sumstats.data['STATUS'].astype('category')
        
        # Verify they are categorical
        self.assertTrue(pd.api.types.is_categorical_dtype(self.sumstats.data['CHR']))
        self.assertTrue(pd.api.types.is_categorical_dtype(self.sumstats.data['STATUS']))
        
        # Mock Pool to avoid multiprocessing pickling issues with mocks
        with patch('gwaslab.qc.qc_fix_sumstats.Pool', side_effect=MockPool), \
             patch('gwaslab.qc.qc_fix_sumstats.liftover_variant') as mock_worker:
            
            # The worker should return a dataframe with the same index but potentially new values
            def side_effect(df, **kwargs):
                df = df.copy()
                # Modify values to something that likely isn't in the categories
                df['CHR'] = df['CHR'].astype(str) + '_new'
                df['STATUS'] = df['STATUS'].astype(str) + '_new'
                return df[['CHR', 'POS', 'STATUS']]
            
            mock_worker.side_effect = side_effect
            
            # Setup data to match the pattern expected by liftover
            self.sumstats.data['STATUS'] = 'xxx0xxx'
            self.sumstats.data['CHR'] = self.sumstats.data['CHR'].astype('category')
            self.sumstats.data['STATUS'] = self.sumstats.data['STATUS'].astype('category')

            with patch('gwaslab.qc.qc_fix_sumstats.get_chain', return_value="dummy_chain"), \
                 patch('gwaslab.qc.qc_fix_sumstats.ChainFile', return_value=MagicMock()), \
                 patch('gwaslab.qc.qc_fix_sumstats._process_build', return_value="hg19"), \
                 patch('gwaslab.qc.qc_fix_sumstats.fixchr', side_effect=lambda x, **kwargs: x), \
                 patch('gwaslab.qc.qc_fix_sumstats.fixpos', side_effect=lambda x, **kwargs: x):
                 
                 from gwaslab.qc.qc_fix_sumstats import parallelizeliftovervariant
                 from gwaslab.g_Log import Log
                 
                 result_df = parallelizeliftovervariant(
                     self.sumstats.data,
                     n_cores=1,
                     chrom="CHR",
                     pos="POS",
                     from_build="19",
                     to_build="38",
                     status="STATUS",
                     chain="dummy.chain", 
                     log=Log()
                 )
                 
                 # Check if CHR and STATUS are now object type (string)
                 # Since fixchr is mocked, it won't revert the type
                 self.assertTrue(result_df['CHR'].dtype == 'object' or result_df['CHR'].dtype == 'string')
                 self.assertTrue(result_df['STATUS'].dtype == 'object' or result_df['STATUS'].dtype == 'string')
                 
                 # Check if values are updated (based on our mock worker)
                 self.assertTrue(result_df['CHR'].iloc[0].endswith('_new'))

if __name__ == '__main__':
    unittest.main()
