import sys
import os
sys.path.insert(0, os.path.abspath("src"))
from gwaslab.viz.viz_aux_params import VizParamsManager

def test_expansion():
    import gwaslab
    print(f"Loaded gwaslab from: {gwaslab.__file__}")
    pm = VizParamsManager()
    
    # Simulate known pairs
    # In a real scenario, these would come from the registry or be inferred
    # For this test, we rely on how compile_from_kwargs handles expansion.
    # But wait, compile_from_kwargs relies on "known_plot_mode_pairs" which is 
    # calculated in `load_viz_config` but NOT stored in `pm` or passed to `compile_from_kwargs`.
    
    # Wait, looking at the code:
    # `compile_from_kwargs` Step 1 calculates `known_plot_mode_pairs` from `self._registry` and `self._arg_map`.
    
    # Let's set up a scenario:
    # 1. We have a plot "plotA" with mode "mode1" explicitly registered via some argument.
    # 2. We register a NEW argument "arg_universal" for "plotA" (implying all modes).
    # 3. We expect "arg_universal" to appear in "plotA:mode1".

    # Setup
    pm.register_arg("arg_specific", plots=["plotA:mode1"], default="spec")
    pm.register_arg("arg_universal", plots=["plotA"], default="univ")
    
    # Run compilation
    # pm.compile_from_kwargs() # This is the internal method, we should simulate what load_viz_config does
    
    # Manually pre-seed registry with known pairs because load_viz_config does that
    known_pairs = [("plotA", "mode1"), ("plotA", None)]
    for p, m in known_pairs:
        key = pm._km(p, m)
        if key not in pm._registry:
            pm._registry[key] = {"allowed": set(), "defaults": {}}
            
    pm.compile_from_kwargs()
    
    # Check results
    print("--- Testing Wildcard Expansion ---")
    
    # Check if plotA:mode1 got arg_universal
    allowed_mode1 = pm.allowed("plotA", "mode1")
    if allowed_mode1 and "arg_universal" in allowed_mode1:
        print(f"Test 1 (Universal arg in specific mode): PASS")
    else:
        print(f"Test 1 (Universal arg in specific mode): FAIL - allowed: {allowed_mode1}")

    # Check defaults
    defaults_mode1 = pm.defaults("plotA", "mode1")
    if defaults_mode1.get("arg_universal") == "univ":
        print(f"Test 2 (Universal default in specific mode): PASS")
    else:
        print(f"Test 2 (Universal default in specific mode): FAIL - defaults: {defaults_mode1}")

    # Check if base plot got it too
    allowed_base = pm.allowed("plotA", None)
    if allowed_base and "arg_universal" in allowed_base:
        print(f"Test 3 (Universal arg in base plot): PASS")
    else:
        print(f"Test 3 (Universal arg in base plot): FAIL - allowed: {allowed_base}")

if __name__ == "__main__":
    test_expansion()
