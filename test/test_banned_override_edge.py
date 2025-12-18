import unittest
import json
import os
import sys
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.viz.viz_aux_params import VizParamsManager, _apply_no_plots

class TestBannedOverrideEdgeCase(unittest.TestCase):
    def setUp(self):
        self.pm = VizParamsManager()
        
    def test_banned_overrides_registry_orphan(self):
        """
        Test that no_plots overrides registry even if the plot 
        is NOT touched by any other args in the config.
        """
        
        # 1. Simulate Registry Loading (Base Layer)
        # Register plot_Orphan with "arg_orphan" allowed
        self.pm._registry["plot_Orphan"] = {
            "allowed": {"arg_orphan"},
            "defaults": {"arg_orphan": "reg_val"}
        }
        
        # 2. Simulate Args Loading (Update Layer)
        # Register arg_orphan with no_plots=["plot_Orphan"]
        # And NO other args for plot_Orphan
        self.pm.register_arg("arg_orphan", plots=["all"], no_plots=["plot_Orphan"])
        
        # 3. Compile
        self.pm.compile_from_kwargs(merge=True)
        
        # 4. Apply no_plots exclusions
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed = self.pm.allowed("plot_Orphan")
        
        # arg_orphan should be GONE
        # If my hypothesis is correct, this will FAIL with current implementation
        self.assertNotIn("arg_orphan", allowed)

if __name__ == "__main__":
    unittest.main()
