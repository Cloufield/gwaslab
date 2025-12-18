import unittest
import sys
import os
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.viz.viz_aux_params import VizParamsManager, _apply_no_plots

class TestNoPlots(unittest.TestCase):
    def setUp(self):
        self.pm = VizParamsManager()

    def test_no_plots_explicit(self):
        """Test explicit exclusion of a plot."""
        # Register arg for all plots, but exclude "plot_A"
        self.pm.register_arg(
            "arg1", 
            plots=["all"], 
            default="val1",
            no_plots=["plot_A"]
        )
        
        # Simulate knowing about plot_A and plot_B
        # We need to manually register them or trick compile_from_kwargs to see them
        # compile_from_kwargs discovers plots from the registry or arg definitions.
        # Since "all" is used, we need at least one other arg or registry entry to define the universe of plots.
        
        # Let's register a dummy arg for plot_A and plot_B to make them "known"
        self.pm.register_arg("dummy_A", plots=["plot_A"])
        self.pm.register_arg("dummy_B", plots=["plot_B"])
        
        self.pm.compile_from_kwargs()
        _apply_no_plots(self.pm)
        
        # Check plot_A (should NOT have arg1)
        allowed_A = self.pm.allowed("plot_A")
        self.assertNotIn("arg1", allowed_A if allowed_A else [])
        
        # Check plot_B (should have arg1)
        allowed_B = self.pm.allowed("plot_B")
        self.assertIn("arg1", allowed_B)

    def test_no_plots_wildcard(self):
        """Test exclusion with wildcard modes."""
        # Register arg for plot_C:mode1 and plot_C:mode2
        # Exclude plot_C:*
        
        self.pm.register_arg(
            "arg2",
            plots=["plot_C"], # implies plot_C:*
            default="val2",
            no_plots=["plot_C:mode1"]
        )
        
        # Register dummy to ensure modes are known
        self.pm.register_arg("dummy_C1", plots=["plot_C:mode1"])
        self.pm.register_arg("dummy_C2", plots=["plot_C:mode2"])
        
        self.pm.compile_from_kwargs()
        _apply_no_plots(self.pm)
        
        # Check plot_C:mode1 (should NOT have arg2)
        allowed_C1 = self.pm.allowed("plot_C", "mode1")
        self.assertNotIn("arg2", allowed_C1 if allowed_C1 else [])
        
        # Check plot_C:mode2 (should have arg2)
        allowed_C2 = self.pm.allowed("plot_C", "mode2")
        self.assertIn("arg2", allowed_C2)

    def test_no_plots_base_exclusion(self):
        """Test exclusion of base plot affecting unknown modes."""
        # Register arg for plot_D
        # Exclude plot_D (base)
        
        self.pm.register_arg(
            "arg3",
            plots=["plot_D"],
            default="val3",
            no_plots=["plot_D"] # Excludes base plot_D
        )
        
        self.pm.compile_from_kwargs()
        _apply_no_plots(self.pm)
        
        # Check plot_D (base) - should be excluded
        allowed_D = self.pm.allowed("plot_D")
        self.assertNotIn("arg3", allowed_D if allowed_D else [])
        
        # Check plot_D:unknown (inherits from base) - should NOT see arg3 because base doesn't have it
        allowed_D_unk = self.pm.allowed("plot_D", "unknown")
        # allowed() for unknown mode returns base allowed.
        self.assertNotIn("arg3", allowed_D_unk if allowed_D_unk else [])

if __name__ == "__main__":
    unittest.main()
