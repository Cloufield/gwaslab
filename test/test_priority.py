import sys
import os
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.viz.viz_aux_params import VizParamsManager

def test_priority():
    pm = VizParamsManager()
    
    # 1. Register Arg with Global Default "Global"
    # plots=["plotA"] now implies wildcard "plotA:*" due to recent change
    pm.register_arg("arg1", plots=["plotA"], default="Global")
    
    # 2. Pre-seed registry to ensure wildcard works (mimicking load_viz_config)
    # We need both base and mode to exist
    pm._registry["plotA"] = {"allowed": set(), "defaults": {}}
    pm._registry["plotA:mode1"] = {"allowed": set(), "defaults": {}}
    
    # 3. Compile (populates registry from args)
    pm.compile_from_kwargs()
    
    # 4. Simulate loading Registry file with Plot-level override
    # This happens AFTER compile_from_kwargs in load_viz_config
    # We set a default for the BASE plot "plotA"
    pm._registry["plotA"]["defaults"]["arg1"] = "RegistryBase"
    
    # 5. Check defaults for Mode1
    # Expected: "RegistryBase" (inherited from plotA because mode1 has no specific override)
    # Current Bug: compile_from_kwargs puts "Global" into mode1, blocking inheritance
    
    d = pm.defaults("plotA", "mode1")
    val = d.get('arg1')
    print(f"Arg1 default for mode1: {val}")
    
    if val == "RegistryBase":
        print("PASS: RegistryBase override inherited.")
    elif val == "Global":
        print("FAIL: Global default blocked inheritance.")
    else:
        print(f"FAIL: Got {val}")

if __name__ == "__main__":
    test_priority()