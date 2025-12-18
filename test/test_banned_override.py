import unittest
import json
import os
import sys
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.viz.viz_aux_params import VizParamsManager, _apply_no_plots

class TestBannedOverride(unittest.TestCase):
    def setUp(self):
        self.pm = VizParamsManager()
        
    def test_banned_overrides_registry(self):
        """Test that no_plots in args overrides allowed in registry."""
        
        # 1. Simulate Registry Loading (Base Layer)
        # Register plot_A with "arg1" allowed
        self.pm._registry["plot_A"] = {
            "allowed": {"arg1", "arg2"},
            "defaults": {"arg1": "reg_val"}
        }
        
        # 2. Simulate Args Loading (Update Layer)
        # Register arg1 with no_plots=["plot_A"]
        # This should REMOVE arg1 from plot_A, even though it's in the registry
        self.pm.register_arg("arg1", plots=["plot_B"], no_plots=["plot_A"])
        
        # Register arg2 normally (should remain)
        self.pm.register_arg("arg2", plots=["plot_A"])
        
        # 3. Compile (Trigger the override logic)
        self.pm.compile_from_kwargs(merge=True)
        
        # 4. Apply no_plots exclusions
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed_A = self.pm.allowed("plot_A")
        
        # arg1 should be GONE
        self.assertNotIn("arg1", allowed_A)
        
        # arg2 should be PRESENT
        self.assertIn("arg2", allowed_A)
        
        # Defaults for arg1 should be kept (implementation keeps defaults for banned args)
        defaults_A = self.pm.defaults("plot_A")
        # arg1 default is kept even though it's banned, so it can be replaced with default value
        self.assertIn("arg1", defaults_A)

if __name__ == "__main__":
    unittest.main()
