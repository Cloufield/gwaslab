import unittest
import os
import json
import sys
import os
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

class TestRegistryPriorityAndUnknownModes(unittest.TestCase):
    def setUp(self):
        self.pm = VizParamsManager()
        self.registry_file = "src/gwaslab/viz/viz_aux_params_registry.txt"
        self.args_file = "src/gwaslab/viz/viz_aux_params.txt"
        
        # Ensure directories exist
        os.makedirs(os.path.dirname(self.registry_file), exist_ok=True)
        os.makedirs(os.path.dirname(self.args_file), exist_ok=True)
        
        # Backup existing files if they exist
        self.backup_registry = None
        self.backup_args = None
        
        if os.path.exists(self.registry_file):
            with open(self.registry_file, 'r') as f:
                self.backup_registry = f.read()
                
        if os.path.exists(self.args_file):
            with open(self.args_file, 'r') as f:
                self.backup_args = f.read()

    def tearDown(self):
        # Restore files
        if self.backup_registry:
            with open(self.registry_file, 'w') as f:
                f.write(self.backup_registry)
        elif os.path.exists(self.registry_file):
            os.remove(self.registry_file)
            
        if self.backup_args:
            with open(self.args_file, 'w') as f:
                f.write(self.backup_args)
        elif os.path.exists(self.args_file):
            os.remove(self.args_file)

    def test_registry_priority_and_args_update(self):
        """
        Test that:
        1. Registry is loaded first (Base Layer)
        2. Args are loaded next (Update Layer)
        3. Args update the registry (overwrite defaults)
        """
        # 1. Create dummy registry
        registry_content = {
            "registry": {
                "plot_test": {
                    "allowed": ["reg_arg"],
                    "defaults": {"reg_arg": "reg_val", "shared_arg": "reg_shared"}
                }
            }
        }
        with open(self.registry_file, 'w') as f:
            json.dump(registry_content, f)
            
        # 2. Create dummy args
        args_content = {
            "args": {
                "args_arg": {
                    "plots": ["plot_test"],
                    "default": "args_val"
                },
                "shared_arg": {
                    "plots": ["plot_test"],
                    "default": "args_shared"  # Should overwrite reg_shared
                }
            }
        }
        with open(self.args_file, 'w') as f:
            json.dump(args_content, f)
            
        # 3. Load config
        load_viz_config(self.pm, path=self.args_file)
        
        # 4. Verify
        defaults = self.pm.defaults("plot_test")
        
        # Registry-only arg should be present
        self.assertEqual(defaults.get("reg_arg"), "reg_val")
        
        # Args-only arg should be present
        self.assertEqual(defaults.get("args_arg"), "args_val")
        
        # Shared arg should be overwritten by args (Update Layer)
        self.assertEqual(defaults.get("shared_arg"), "args_shared")
        
    def test_unknown_mode_inheritance(self):
        """
        Test that unknown modes inherit from base plot.
        User request: "plot_mqq" mean the arg will be applied to modes like "plot_mqq:m" and all modes including known and unknown ones!
        """
        # 1. Register arg for base plot
        self.pm.register_arg("base_arg", plots=["plot_unknown_test"], default="base_val")
        self.pm.compile_from_kwargs()
        
        # 2. Check Unknown Mode (should inherit from base)
        # "plot_unknown_test:unknown_mode" is NOT in registry, so it should fallback to "plot_unknown_test"
        
        # Check defaults
        defs = self.pm.defaults("plot_unknown_test", "unknown_mode")
        self.assertEqual(defs.get("base_arg"), "base_val")
        
        # Check allowed
        allowed = self.pm.allowed("plot_unknown_test", "unknown_mode")
        self.assertIsNotNone(allowed)
        self.assertIn("base_arg", allowed)

if __name__ == '__main__':
    unittest.main()
