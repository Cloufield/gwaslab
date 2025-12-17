import unittest
import json
import os
import sys
import os
sys.path.insert(0, os.path.abspath("src"))
from unittest.mock import patch, mock_open
from src.gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

class TestLoadConfigNoPlots(unittest.TestCase):
    def setUp(self):
        self.pm = VizParamsManager()

    def test_load_config_no_plots(self):
        """Test that no_plots from config file is correctly loaded and applied."""
        
        # Mock config data
        config_data = {
            "args": {
                "arg_should_be_banned": {
                    "plots": ["all"],
                    "default": "val",
                    "no_plots": ["plot_A"]
                }
            }
        }
        
        # Mock registry data (optional, but good for testing interaction)
        registry_data = {
            "registry": {
                "plot_A": {
                    "allowed": ["arg_should_be_banned"],
                    "defaults": {"arg_should_be_banned": "reg_val"}
                },
                "plot_B": {
                    "allowed": [],
                    "defaults": {}
                }
            }
        }

        # Mock open to return different content based on filename
        def side_effect(filename, *args, **kwargs):
            if "viz_aux_params.txt" in filename:
                return mock_open(read_data=json.dumps(config_data)).return_value
            elif "viz_aux_params_registry.txt" in filename:
                return mock_open(read_data=json.dumps(registry_data)).return_value
            return mock_open()()

        with patch("builtins.open", side_effect=side_effect):
            with patch("os.path.exists", return_value=True):
                load_viz_config(self.pm)
        
        print("Registry keys:", self.pm._registry.keys())
        if "plot_B" in self.pm._registry:
            print("plot_B allowed:", self.pm._registry["plot_B"].get("allowed"))

        # Verify arg_should_be_banned is NOT in plot_A allowed
        allowed_A = self.pm.allowed("plot_A")
        self.assertNotIn("arg_should_be_banned", allowed_A)
        
        # Verify arg_should_be_banned IS in plot_B allowed (because "all" in plots, and not in no_plots)
        allowed_B = self.pm.allowed("plot_B")
        self.assertIn("arg_should_be_banned", allowed_B)

if __name__ == '__main__':
    unittest.main()
