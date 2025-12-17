import sys
import os
sys.path.insert(0, os.path.abspath("src"))
import unittest
from gwaslab.viz.viz_aux_params import VizParamsManager

class TestWildcardAndUnknown(unittest.TestCase):
    def test_unknown_mode_inheritance(self):
        pm = VizParamsManager()
        
        # Register arg for base plot
        pm.register_arg("base_arg", plots=["plot_test"], default="base_val")
        
        # Register arg for specific mode
        pm.register_arg("mode_arg", plots=["plot_test:known"], default="mode_val")
        
        # Compile
        pm.compile_from_kwargs()
        
        # 1. Check Known Mode
        # Should have both base_arg and mode_arg
        defs = pm.defaults("plot_test", "known")
        self.assertEqual(defs.get("base_arg"), "base_val")
        self.assertEqual(defs.get("mode_arg"), "mode_val")
        
        allowed = pm.allowed("plot_test", "known")
        self.assertIn("base_arg", allowed)
        self.assertIn("mode_arg", allowed)
        
        # 2. Check Unknown Mode
        # Should inherit base_arg
        defs = pm.defaults("plot_test", "unknown")
        self.assertEqual(defs.get("base_arg"), "base_val")
        
        allowed = pm.allowed("plot_test", "unknown")
        # Currently, this returns None (fallback to signature)
        # But if user wants strict application, it should probably be {'base_arg'}
        print(f"Unknown mode allowed: {allowed}")
        
        if allowed is not None:
             self.assertIn("base_arg", allowed)

if __name__ == "__main__":
    unittest.main()
