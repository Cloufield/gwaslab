
import unittest
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gwaslab as gl
import os

class TestHighlightPinpoint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Create a dummy sumstats object
        cls.df = pd.DataFrame({
            'SNPID': [f'rs{i}' for i in range(1, 101)],
            'CHR': [1] * 50 + [2] * 50,
            'POS': list(range(1, 101)),
            'EA': ['A'] * 100,
            'NEA': ['G'] * 100,
            'P': np.random.uniform(0, 1, 100),
            'EAF': np.random.uniform(0.1, 0.9, 100),
            'BETA': np.random.normal(0, 1, 100),
            'SE': np.random.uniform(0.1, 0.2, 100)
        })
        cls.sumstats = gl.Sumstats(cls.df, fmt="gwaslab")
        cls.sumstats.basic_check()

    def tearDown(self):
        plt.close('all')

    def test_highlight_single_string(self):
        print("\nTesting highlight with single string...")
        # Should highlight rs1
        self.sumstats.plot_mqq(highlight="rs1", save="test_highlight_string.png")
        self.assertTrue(os.path.exists("test_highlight_string.png"))
        os.remove("test_highlight_string.png")

    def test_highlight_list_strings(self):
        print("\nTesting highlight with list of strings...")
        # Should highlight rs1 and rs2 in same color
        self.sumstats.plot_mqq(highlight=["rs1", "rs2"], save="test_highlight_list.png")
        self.assertTrue(os.path.exists("test_highlight_list.png"))
        os.remove("test_highlight_list.png")

    def test_highlight_grouped_single_color(self):
        print("\nTesting highlight with groups (single color)...")
        # Should highlight rs1/rs2 (group 1) and rs3/rs4 (group 2) with default cycling color
        self.sumstats.plot_mqq(highlight=[["rs1", "rs2"], ["rs3", "rs4"]], save="test_highlight_groups.png")
        self.assertTrue(os.path.exists("test_highlight_groups.png"))
        os.remove("test_highlight_groups.png")

    def test_highlight_grouped_multi_color(self):
        print("\nTesting highlight with groups (multiple colors)...")
        # Should highlight rs1/rs2 (red) and rs3/rs4 (blue)
        self.sumstats.plot_mqq(
            highlight=[["rs1", "rs2"], ["rs3", "rs4"]], 
            highlight_color=["red", "blue"],
            save="test_highlight_groups_color.png"
        )
        self.assertTrue(os.path.exists("test_highlight_groups_color.png"))
        os.remove("test_highlight_groups_color.png")

    def test_pinpoint_single_string(self):
        print("\nTesting pinpoint with single string...")
        self.sumstats.plot_mqq(pinpoint="rs10", save="test_pinpoint_string.png")
        self.assertTrue(os.path.exists("test_pinpoint_string.png"))
        os.remove("test_pinpoint_string.png")

    def test_pinpoint_list_strings(self):
        print("\nTesting pinpoint with list of strings...")
        self.sumstats.plot_mqq(pinpoint=["rs10", "rs11"], save="test_pinpoint_list.png")
        self.assertTrue(os.path.exists("test_pinpoint_list.png"))
        os.remove("test_pinpoint_list.png")

    def test_pinpoint_grouped_multi_color(self):
        print("\nTesting pinpoint with groups (multiple colors)...")
        self.sumstats.plot_mqq(
            pinpoint=[["rs10", "rs11"], ["rs12", "rs13"]], 
            pinpoint_color=["green", "purple"],
            save="test_pinpoint_groups.png"
        )
        self.assertTrue(os.path.exists("test_pinpoint_groups.png"))
        os.remove("test_pinpoint_groups.png")

    def test_mixed_highlight_pinpoint(self):
        print("\nTesting mixed highlight and pinpoint...")
        self.sumstats.plot_mqq(
            highlight=["rs1", "rs2"],
            pinpoint=[["rs10"], ["rs20"]],
            pinpoint_color=["orange", "cyan"],
            save="test_mixed.png"
        )
        self.assertTrue(os.path.exists("test_mixed.png"))
        os.remove("test_mixed.png")

if __name__ == '__main__':
    unittest.main()
