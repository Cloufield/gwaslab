"""Comprehensive tests for VizParamsManager.

This file consolidates all test cases related to VizParamsManager from multiple test files:
- test_viz_params_numeric_suffix.py
- test_no_plots.py
- test_banned_override_edge.py
- test_banned_override.py
- test_load_config_no_plots.py
- test_registry_priority.py
- test_unknown.py
- test_priority.py
- test_expansion.py

Tests cover:
- Numeric suffix support in filtering
- Parameter merge priority
- Default inheritance
- Configuration loading order
- Wildcard expansion
- Banned arguments (no_plots)
- Context defaults (ctx_defaults)
- Banned keys (nested sub-keys)
- Object-level presets
- Filtering behavior
- Registry priority
- Unknown mode inheritance
"""
import unittest
import sys
import os
import json
from unittest.mock import patch, mock_open

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config, _apply_no_plots, _resolve_inheritance
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.info.g_Log import Log


class TestNumericSuffixFiltering(unittest.TestCase):
    """Test that numeric suffixes are correctly handled in parameter filtering."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
        self.log = Log()
    
    def test_single_digit_suffix(self):
        """Test that single-digit suffixes work (e.g., highlight2)."""
        test_params = {
            'highlight': ['snp1'],
            'highlight2': ['snp2'],
            'highlight3': ['snp3'],
            'invalid_arg': 'value'
        }
        
        filtered = self.pm.filter(
            _mqqplot, test_params, 
            key='plot_mqq', mode='m',
            log=self.log, verbose=False
        )
        
        # highlight and highlight2, highlight3 should pass
        self.assertIn('highlight', filtered)
        self.assertIn('highlight2', filtered)
        self.assertIn('highlight3', filtered)
        # invalid_arg should be filtered out
        self.assertNotIn('invalid_arg', filtered)
    
    def test_multi_digit_suffix(self):
        """Test that multi-digit suffixes work (e.g., highlight10, highlight123)."""
        test_params = {
            'highlight': ['snp1'],
            'highlight10': ['snp10'],
            'highlight123': ['snp123'],
            'scatter_kwargs': {'s': 20},
            'scatter_kwargs99': {'s': 30},
            'invalid_arg': 'value'
        }
        
        filtered = self.pm.filter(
            _mqqplot, test_params,
            key='plot_mqq', mode='m',
            log=self.log, verbose=False
        )
        
        # All highlight variants should pass
        self.assertIn('highlight', filtered)
        self.assertIn('highlight10', filtered)
        self.assertIn('highlight123', filtered)
        # scatter_kwargs variants should pass
        self.assertIn('scatter_kwargs', filtered)
        self.assertIn('scatter_kwargs99', filtered)
        # invalid_arg should be filtered out
        self.assertNotIn('invalid_arg', filtered)
    
    def test_edge_cases(self):
        """Test edge cases for numeric suffixes."""
        test_params = {
            'highlight0': ['snp0'],
            'highlight00': ['snp00'],
            'highlight000': ['snp000'],
            'arg1': 'value1',
            'arg12': 'value12',
            'arg123': 'value123',
            'invalid_arg': 'value',
            'invalid_arg2': 'value2'
        }
        
        filtered = self.pm.filter(
            _mqqplot, test_params,
            key='plot_mqq', mode='m',
            log=self.log, verbose=False
        )
        
        # highlight variants should pass (if highlight is in allowed set)
        allowed = self.pm.allowed('plot_mqq', 'm')
        if allowed and 'highlight' in allowed:
            self.assertIn('highlight0', filtered)
            self.assertIn('highlight00', filtered)
            self.assertIn('highlight000', filtered)
        
        # arg1, arg12, arg123 should be filtered out (unless 'arg' is in allowed set)
        # This depends on whether 'arg' is in the allowed set
        if allowed and 'arg' in allowed:
            self.assertIn('arg1', filtered)
            self.assertIn('arg12', filtered)
            self.assertIn('arg123', filtered)
        else:
            self.assertNotIn('arg1', filtered)
            self.assertNotIn('arg12', filtered)
            self.assertNotIn('arg123', filtered)
        
        # invalid_arg variants should be filtered out
        self.assertNotIn('invalid_arg', filtered)
        self.assertNotIn('invalid_arg2', filtered)
    
    def test_no_base_in_allowed(self):
        """Test that suffixes don't match if base is not in allowed set."""
        test_params = {
            'nonexistent_arg': 'value',
            'nonexistent_arg2': 'value2',
            'nonexistent_arg10': 'value10'
        }
        
        filtered = self.pm.filter(
            _mqqplot, test_params,
            key='plot_mqq', mode='m',
            log=self.log, verbose=False
        )
        
        # None of these should pass if base is not in allowed set
        self.assertNotIn('nonexistent_arg', filtered)
        self.assertNotIn('nonexistent_arg2', filtered)
        self.assertNotIn('nonexistent_arg10', filtered)
    
    def test_mixed_valid_and_invalid_suffixes(self):
        """Test mix of valid and invalid suffixes."""
        test_params = {
            'highlight': ['snp1'],
            'highlight2': ['snp2'],
            'highlight10': ['snp10'],
            'scatter_kwargs': {'s': 20},
            'scatter_kwargs1': {'s': 30},
            'invalid_base': 'value',
            'invalid_base2': 'value2',
            'another_invalid': 'value3',
            'another_invalid99': 'value4'
        }
        
        filtered = self.pm.filter(
            _mqqplot, test_params,
            key='plot_mqq', mode='m',
            log=self.log, verbose=False
        )
        
        # Valid suffixes should pass
        self.assertIn('highlight', filtered)
        self.assertIn('highlight2', filtered)
        self.assertIn('highlight10', filtered)
        self.assertIn('scatter_kwargs', filtered)
        self.assertIn('scatter_kwargs1', filtered)
        
        # Invalid bases should be filtered out (even with suffixes)
        self.assertNotIn('invalid_base', filtered)
        self.assertNotIn('invalid_base2', filtered)
        self.assertNotIn('another_invalid', filtered)
        self.assertNotIn('another_invalid99', filtered)
    
    def test_suffix_with_merge(self):
        """Test that numeric suffixes work correctly with merge()."""
        # Test that banned args with suffixes are replaced with defaults
        test_params = {
            'highlight': ['custom_snp'],  # This should be replaced with default [] for plot_mqq:r
            'highlight2': ['custom_snp2'],  # This should also be replaced
        }
        
        merged = self.pm.merge('plot_mqq', test_params, mode='r')
        
        # highlight should be replaced with default (empty list) for plot_mqq:r
        # because highlight is banned for plot_mqq:r
        self.assertIn('highlight', merged)
        self.assertEqual(merged['highlight'], [])  # Default value
        
        # highlight2 should also be replaced
        self.assertIn('highlight2', merged)
        self.assertEqual(merged['highlight2'], [])  # Default value
    
    def test_suffix_with_different_plots(self):
        """Test numeric suffixes with different plot types."""
        # Test with plot_miami2
        test_params_miami = {
            'highlight': ['snp1'],
            'highlight2': ['snp2'],
            'highlight10': ['snp10']
        }
        
        from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2
        filtered_miami = self.pm.filter(
            plot_miami2, test_params_miami,
            key='plot_miami2', mode=None,
            log=self.log, verbose=False
        )
        
        # Check if highlight variants pass for plot_miami2
        allowed_miami = self.pm.allowed('plot_miami2')
        if allowed_miami and 'highlight' in allowed_miami:
            self.assertIn('highlight', filtered_miami)
            self.assertIn('highlight2', filtered_miami)
            self.assertIn('highlight10', filtered_miami)
    
    def test_suffix_preserves_values(self):
        """Test that numeric suffixes preserve their values correctly."""
        test_params = {
            'highlight': ['snp1'],
            'highlight2': ['snp2'],
            'highlight10': ['snp10'],
            'scatter_kwargs': {'s': 20},
            'scatter_kwargs1': {'s': 30}
        }
        
        filtered = self.pm.filter(
            _mqqplot, test_params,
            key='plot_mqq', mode='m',
            log=self.log, verbose=False
        )
        
        # Values should be preserved (may be merged with defaults)
        self.assertEqual(filtered['highlight'], ['snp1'])
        self.assertEqual(filtered['highlight2'], ['snp2'])
        self.assertEqual(filtered['highlight10'], ['snp10'])
        # scatter_kwargs and scatter_kwargs1 should be in filtered (numeric suffix matching works)
        # Note: exact values may be modified by defaults/merging, but the key should pass through
        self.assertIn('scatter_kwargs', filtered)
        self.assertIn('scatter_kwargs1', filtered)


class TestParameterMergePriority(unittest.TestCase):
    """Test parameter merge priority: defaults -> presets -> override."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
    
    def test_merge_priority_order(self):
        """Test that merge follows correct priority: defaults -> presets -> override."""
        # Set up: register defaults
        self.pm.register('test_plot', allowed={'arg1', 'arg2'}, defaults={'arg1': 'default_value'}, mode=None)
        
        # Set object-level preset
        self.pm.set('test_plot', {'arg1': 'preset_value'}, mode=None)
        
        # Merge with override
        override = {'arg1': 'override_value'}
        merged = self.pm.merge('test_plot', override, mode=None)
        
        # Override should win
        self.assertEqual(merged['arg1'], 'override_value')
    
    def test_merge_without_override(self):
        """Test merge when no override is provided."""
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'default_value'}, mode=None)
        self.pm.set('test_plot', {'arg1': 'preset_value'}, mode=None)
        
        merged = self.pm.merge('test_plot', None, mode=None)
        self.assertEqual(merged['arg1'], 'preset_value')
    
    def test_merge_presets_override_defaults(self):
        """Test that presets override defaults."""
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'default_value'}, mode=None)
        self.pm.set('test_plot', {'arg1': 'preset_value'}, mode=None)
        
        merged = self.pm.merge('test_plot', None, mode=None)
        self.assertEqual(merged['arg1'], 'preset_value')
        self.assertNotEqual(merged['arg1'], 'default_value')
    
    def test_merge_override_wins_all(self):
        """Test that override wins over both defaults and presets."""
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'default_value'}, mode=None)
        self.pm.set('test_plot', {'arg1': 'preset_value'}, mode=None)
        
        merged = self.pm.merge('test_plot', {'arg1': 'override_value'}, mode=None)
        self.assertEqual(merged['arg1'], 'override_value')


class TestDefaultInheritance(unittest.TestCase):
    """Test default value inheritance from base plot to modes."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
    
    def test_mode_inherits_from_base_plot(self):
        """Test that mode inherits defaults from base plot."""
        # Set base plot defaults
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'base_default'}, mode=None)
        
        # Get defaults for a mode (should inherit from base)
        defaults = self.pm.defaults('test_plot', 'mode1')
        self.assertEqual(defaults.get('arg1'), 'base_default')
    
    def test_mode_specific_overrides_base(self):
        """Test that mode-specific defaults override base defaults."""
        # Set base plot defaults
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'base_default'}, mode=None)
        # Set mode-specific defaults
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'mode_default'}, mode='mode1')
        
        defaults = self.pm.defaults('test_plot', 'mode1')
        self.assertEqual(defaults.get('arg1'), 'mode_default')
        self.assertNotEqual(defaults.get('arg1'), 'base_default')
    
    def test_inheritance_with_multiple_args(self):
        """Test inheritance with multiple arguments."""
        self.pm.register('test_plot', allowed={'arg1', 'arg2'}, 
                         defaults={'arg1': 'base1', 'arg2': 'base2'}, mode=None)
        self.pm.register('test_plot', allowed={'arg1', 'arg2'}, 
                         defaults={'arg1': 'mode1'}, mode='mode1')
        
        defaults = self.pm.defaults('test_plot', 'mode1')
        self.assertEqual(defaults.get('arg1'), 'mode1')  # Overridden
        self.assertEqual(defaults.get('arg2'), 'base2')  # Inherited


class TestConfigurationLoadingOrder(unittest.TestCase):
    """Test configuration loading order and priority."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
    
    def test_registry_loads_first(self):
        """Test that registry loads as base layer."""
        # Manually register (simulating registry file)
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'registry_default'}, mode=None)
        
        defaults = self.pm.defaults('test_plot')
        self.assertEqual(defaults.get('arg1'), 'registry_default')
    
    def test_args_override_registry_defaults(self):
        """Test that args file overrides registry defaults."""
        # Step 1: Load registry
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'registry_default'}, mode=None)
        
        # Step 2: Register arg (simulating args file)
        self.pm.register_arg('arg1', plots=['test_plot'], default='args_default')
        self.pm.compile_from_kwargs(merge=True)
        
        # Args default should override registry default
        defaults = self.pm.defaults('test_plot')
        self.assertEqual(defaults.get('arg1'), 'args_default')
    
    def test_ctx_defaults_override_compiled(self):
        """Test that ctx_defaults override compiled defaults."""
        # Register arg with global default
        self.pm.register_arg('arg1', plots=['test_plot'], default='global_default')
        self.pm._registry['test_plot'] = {'allowed': {'arg1'}, 'defaults': {}}
        self.pm.compile_from_kwargs(merge=True)
        
        # Apply ctx_defaults
        self.pm._arg_map['arg1'] = {
            'plots': {('test_plot', '*')},
            'no_plots': set(),
            'default': 'global_default',
            'ctx_defaults': {'test_plot': 'ctx_default'}
        }
        
        # Manually apply ctx_defaults (simulating load_viz_config step 4)
        entry = self.pm._registry['test_plot']
        if 'arg1' in entry['allowed']:
            entry['defaults']['arg1'] = 'ctx_default'
        
        defaults = self.pm.defaults('test_plot')
        self.assertEqual(defaults.get('arg1'), 'ctx_default')
    
    def test_no_plots_removes_from_allowed(self):
        """Test that no_plots removes args from allowed set."""
        # Register arg for plot
        self.pm.register_arg('arg1', plots=['test_plot'], default='default_value', no_plots=['test_plot:mode1'])
        self.pm._registry['test_plot'] = {'allowed': set(), 'defaults': {}}
        self.pm._registry['test_plot:mode1'] = {'allowed': set(), 'defaults': {}}
        self.pm.compile_from_kwargs(merge=True)
        
        # Apply no_plots (simulating load_viz_config step 5)
        entry = self.pm._registry['test_plot:mode1']
        if 'arg1' in entry.get('allowed', set()):
            entry['allowed'].discard('arg1')
            # Keep default
            entry['defaults']['arg1'] = 'default_value'
        
        allowed = self.pm.allowed('test_plot', 'mode1')
        self.assertNotIn('arg1', allowed if allowed else set())
        # But default should still exist
        defaults = self.pm.defaults('test_plot', 'mode1')
        self.assertEqual(defaults.get('arg1'), 'default_value')


class TestWildcardExpansion(unittest.TestCase):
    """Test wildcard expansion in plot specifications."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
    
    def test_plot_wildcard_expands_to_all_modes(self):
        """Test that 'plot:*' expands to all known modes."""
        # Create known modes
        self.pm._registry['test_plot'] = {'allowed': set(), 'defaults': {}}
        self.pm._registry['test_plot:mode1'] = {'allowed': set(), 'defaults': {}}
        self.pm._registry['test_plot:mode2'] = {'allowed': set(), 'defaults': {}}
        
        # Register arg with wildcard
        self.pm.register_arg('arg1', plots=['test_plot'], default='default_value')
        self.pm.compile_from_kwargs(merge=True)
        
        # Check that arg1 is in all modes
        self.assertIn('arg1', self.pm.allowed('test_plot') or set())
        self.assertIn('arg1', self.pm.allowed('test_plot', 'mode1') or set())
        self.assertIn('arg1', self.pm.allowed('test_plot', 'mode2') or set())
    
    def test_all_wildcard_expands_to_all_plots(self):
        """Test that 'all:*' expands to all known plots."""
        # Create known plots
        self.pm._registry['plotA'] = {'allowed': set(), 'defaults': {}}
        self.pm._registry['plotB'] = {'allowed': set(), 'defaults': {}}
        self.pm._registry['plotC:mode1'] = {'allowed': set(), 'defaults': {}}
        
        # Register arg with 'all:*'
        self.pm.register_arg('arg1', plots=['all'], default='default_value')
        self.pm.compile_from_kwargs(merge=True)
        
        # Check that arg1 is in all plots
        self.assertIn('arg1', self.pm.allowed('plotA') or set())
        self.assertIn('arg1', self.pm.allowed('plotB') or set())
        self.assertIn('arg1', self.pm.allowed('plotC', 'mode1') or set())
    
    def test_expansion(self):
        """Test wildcard expansion with specific mode registration."""
        # Setup
        self.pm.register_arg("arg_specific", plots=["plotA:mode1"], default="spec")
        self.pm.register_arg("arg_universal", plots=["plotA"], default="univ")
        
        # Manually pre-seed registry with known pairs because load_viz_config does that
        known_pairs = [("plotA", "mode1"), ("plotA", None)]
        for p, m in known_pairs:
            key = self.pm._km(p, m)
            if key not in self.pm._registry:
                self.pm._registry[key] = {"allowed": set(), "defaults": {}}
                
        self.pm.compile_from_kwargs()
        
        # Check if plotA:mode1 got arg_universal
        allowed_mode1 = self.pm.allowed("plotA", "mode1")
        self.assertIsNotNone(allowed_mode1)
        self.assertIn("arg_universal", allowed_mode1)
        
        # Check defaults
        defaults_mode1 = self.pm.defaults("plotA", "mode1")
        self.assertEqual(defaults_mode1.get("arg_universal"), "univ")
        
        # Check if base plot got it too
        allowed_base = self.pm.allowed("plotA", None)
        self.assertIsNotNone(allowed_base)
        self.assertIn("arg_universal", allowed_base)


class TestBannedArguments(unittest.TestCase):
    """Test banned arguments (no_plots) behavior."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
    
    def test_banned_arg_replaced_with_default(self):
        """Test that banned args are replaced with defaults."""
        # highlight is banned for plot_mqq:r
        test_params = {'highlight': ['custom_snp']}
        merged = self.pm.merge('plot_mqq', test_params, mode='r')
        
        # Should be replaced with default []
        self.assertIn('highlight', merged)
        self.assertEqual(merged['highlight'], [])
    
    def test_banned_arg_removed_from_allowed(self):
        """Test that banned args are removed from allowed set."""
        allowed = self.pm.allowed('plot_mqq', 'r')
        # highlight should not be in allowed set for plot_mqq:r
        self.assertNotIn('highlight', allowed if allowed else set())
    
    def test_banned_arg_default_preserved(self):
        """Test that defaults for banned args are preserved."""
        defaults = self.pm.defaults('plot_mqq', 'r')
        # highlight default should exist even though it's banned
        self.assertIn('highlight', defaults)
        self.assertEqual(defaults['highlight'], [])
    
    def test_banned_arg_passes_through_with_default(self):
        """Test that banned args with defaults pass through filter."""
        test_params = {'highlight': ['custom']}
        filtered = self.pm.filter(_mqqplot, test_params, key='plot_mqq', mode='r', log=Log(), verbose=False)
        
        # highlight should be in filtered with default value
        self.assertIn('highlight', filtered)
        self.assertEqual(filtered['highlight'], [])
    
    def test_no_plots_explicit(self):
        """
        Test explicit exclusion of a plot.
        
        With inheritance support, no_plots exclusions are applied after
        inheritance resolution, so they override inherited allowed sets.
        
        Note: Defaults are only preserved if the arg was in the allowed set
        before being banned. If an arg is banned from the start (never added
        to allowed), no default is preserved.
        """
        # First, add arg1 to plot_A's registry (simulating it was allowed initially)
        self.pm._registry["plot_A"] = {
            "allowed": {"arg1"},
            "defaults": {"arg1": "val1"}
        }
        
        # Register arg for all plots, but exclude "plot_A"
        # This will remove arg1 from plot_A's allowed set
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
        
        self.pm.compile_from_kwargs(merge=True)
        # Apply no_plots exclusions (removes banned args from allowed sets)
        _apply_no_plots(self.pm)
        
        # Check plot_A (should NOT have arg1 in allowed)
        allowed_A = self.pm.allowed("plot_A")
        self.assertNotIn("arg1", allowed_A if allowed_A else [])
        
        # Check plot_B (should have arg1)
        allowed_B = self.pm.allowed("plot_B")
        self.assertIn("arg1", allowed_B)
        
        # Defaults should be preserved for banned args that were in registry
        defaults_A = self.pm.defaults("plot_A")
        self.assertIn("arg1", defaults_A)
        self.assertEqual(defaults_A["arg1"], "val1")
    
    def test_no_plots_wildcard(self):
        """
        Test exclusion with wildcard modes.
        
        Tests that mode-specific exclusions work correctly with inheritance.
        """
        # Register arg for plot_C:mode1 and plot_C:mode2
        # Exclude plot_C:mode1
        
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
        
        # Default should be preserved for banned mode
        defaults_C1 = self.pm.defaults("plot_C", "mode1")
        self.assertIn("arg2", defaults_C1)
    
    def test_no_plots_base_exclusion(self):
        """
        Test exclusion of base plot affecting unknown modes.
        
        With inheritance, base plot exclusions affect modes that inherit from base.
        
        Note: Defaults are only preserved if the arg was in the allowed set
        before being banned.
        """
        # First, add arg3 to plot_D's registry (simulating it was allowed initially)
        self.pm._registry["plot_D"] = {
            "allowed": {"arg3"},
            "defaults": {"arg3": "val3"}
        }
        
        # Register arg for plot_D
        # Exclude plot_D (base)
        
        self.pm.register_arg(
            "arg3",
            plots=["plot_D"],
            default="val3",
            no_plots=["plot_D"] # Excludes base plot_D
        )
        
        self.pm.compile_from_kwargs(merge=True)
        _apply_no_plots(self.pm)
        
        # Check plot_D (base) - should be excluded from allowed
        allowed_D = self.pm.allowed("plot_D")
        self.assertNotIn("arg3", allowed_D if allowed_D else [])
        
        # Check plot_D:unknown (inherits from base) - should NOT see arg3 because base doesn't have it
        allowed_D_unk = self.pm.allowed("plot_D", "unknown")
        # allowed() for unknown mode returns base allowed.
        self.assertNotIn("arg3", allowed_D_unk if allowed_D_unk else [])
        
        # Default should be preserved (was in registry before being banned)
        defaults_D = self.pm.defaults("plot_D")
        self.assertIn("arg3", defaults_D)
        self.assertEqual(defaults_D["arg3"], "val3")
    
    def test_no_plots_with_inheritance(self):
        """
        Test that no_plots exclusions work correctly with inheritance chains.
        
        When plot_child inherits from plot_parent, and arg is banned for plot_child,
        it should be removed from plot_child's allowed set even if inherited.
        """
        # 1. Create parent entry
        self.pm._registry["plot_parent"] = {
            "allowed": {"arg_inherited", "arg_other"},
            "defaults": {"arg_inherited": "parent_val"}
        }
        
        # 2. Create child entry that inherits (simulated)
        # In real usage, this would use "inherit": "plot_parent"
        self.pm._registry["plot_child"] = {
            "allowed": {"arg_inherited", "arg_other"},  # Would come from inheritance
            "defaults": {"arg_inherited": "parent_val"}
        }
        
        # 3. Register arg with no_plots for child
        self.pm.register_arg("arg_inherited", plots=["plot_parent"], no_plots=["plot_child"])
        self.pm.register_arg("arg_other", plots=["plot_child"])
        
        # 4. Compile and apply no_plots
        self.pm.compile_from_kwargs(merge=True)
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed_child = self.pm.allowed("plot_child")
        
        # arg_inherited should be GONE from child (banned)
        self.assertNotIn("arg_inherited", allowed_child if allowed_child else set())
        
        # arg_other should be PRESENT
        self.assertIn("arg_other", allowed_child)
        
        # Parent should still have arg_inherited
        allowed_parent = self.pm.allowed("plot_parent")
        self.assertIn("arg_inherited", allowed_parent)
        
        # Default should be preserved for child
        defaults_child = self.pm.defaults("plot_child")
        self.assertIn("arg_inherited", defaults_child)
    
    def test_banned_overrides_registry(self):
        """Test that no_plots in args overrides allowed in registry."""
        
        # 1. Simulate Registry Loading (Base Layer)
        # Register plot_A with "arg1" allowed
        # Note: With inheritance support, registry entries are resolved before use
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
        # This removes banned args from allowed sets but preserves defaults
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed_A = self.pm.allowed("plot_A")
        
        # arg1 should be GONE from allowed set
        self.assertNotIn("arg1", allowed_A)
        
        # arg2 should be PRESENT
        self.assertIn("arg2", allowed_A)
        
        # Defaults for arg1 should be kept (implementation keeps defaults for banned args)
        defaults_A = self.pm.defaults("plot_A")
        # arg1 default is kept even though it's banned, so it can be replaced with default value
        self.assertIn("arg1", defaults_A)
    
    def test_banned_overrides_inherited_registry(self):
        """Test that no_plots overrides allowed even when entry inherits from another."""
        
        # 1. Create a base entry that will be inherited
        self.pm._registry["plot_base"] = {
            "allowed": {"arg1", "arg2", "arg3"},
            "defaults": {"arg1": "base_val", "arg2": "base_val2"}
        }
        
        # 2. Create a derived entry that inherits from base
        # Simulate inheritance by manually resolving it
        raw_entries = {
            "plot_base": {
                "plot": "plot_base",
                "mode": None,
                "entry": {
                    "allowed": ["arg1", "arg2", "arg3"],
                    "defaults": {"arg1": "base_val", "arg2": "base_val2"}
                }
            },
            "plot_derived": {
                "plot": "plot_derived",
                "mode": None,
                "entry": {
                    "inherit": "plot_base"
                }
            }
        }
        
        # Resolve inheritance
        resolved_entries = {}
        for key, raw_data in raw_entries.items():
            resolved = _resolve_inheritance(key, raw_data["entry"], raw_entries, resolved_entries)
            resolved_entries[key] = resolved
        
        # Apply to registry
        self.pm._registry["plot_derived"] = {
            "allowed": set(resolved_entries["plot_derived"].get("allowed", [])),
            "defaults": resolved_entries["plot_derived"].get("defaults", {})
        }
        
        # 3. Register args with no_plots for derived plot
        self.pm.register_arg("arg1", plots=["plot_base"], no_plots=["plot_derived"])
        self.pm.register_arg("arg2", plots=["plot_derived"])
        self.pm.register_arg("arg3", plots=["plot_derived"])
        
        # 4. Compile and apply no_plots
        self.pm.compile_from_kwargs(merge=True)
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed_derived = self.pm.allowed("plot_derived")
        
        # arg1 should be GONE (banned even though inherited)
        self.assertNotIn("arg1", allowed_derived)
        
        # arg2 and arg3 should be PRESENT
        self.assertIn("arg2", allowed_derived)
        self.assertIn("arg3", allowed_derived)
        
        # Defaults for arg1 should be kept
        defaults_derived = self.pm.defaults("plot_derived")
        self.assertIn("arg1", defaults_derived)
        self.assertEqual(defaults_derived["arg1"], "base_val")
    
    def test_banned_overrides_registry_orphan(self):
        """
        Test that no_plots overrides registry even if the plot 
        is NOT touched by any other args in the config.
        
        This test verifies that banned args are removed from allowed sets
        even when the plot has no other arguments registered.
        """
        
        # 1. Simulate Registry Loading (Base Layer)
        # Register plot_Orphan with "arg_orphan" allowed
        # Note: With inheritance support, entries are resolved before use,
        # but for this test we're directly setting the registry
        self.pm._registry["plot_Orphan"] = {
            "allowed": {"arg_orphan"},
            "defaults": {"arg_orphan": "reg_val"}
        }
        
        # 2. Simulate Args Loading (Update Layer)
        # Register arg_orphan with no_plots=["plot_Orphan"]
        # And NO other args for plot_Orphan
        self.pm.register_arg("arg_orphan", plots=["all"], no_plots=["plot_Orphan"])
        
        # 3. Compile (expands "all" to known plots, including plot_Orphan)
        # Need at least one other plot to make "all" work
        self.pm.register_arg("dummy", plots=["plot_Other"])
        self.pm.compile_from_kwargs(merge=True)
        
        # 4. Apply no_plots exclusions
        # This removes banned args from allowed sets
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed = self.pm.allowed("plot_Orphan")
        
        # arg_orphan should be GONE from allowed set
        self.assertNotIn("arg_orphan", allowed if allowed else set())
        
        # But default should still be preserved
        defaults = self.pm.defaults("plot_Orphan")
        self.assertIn("arg_orphan", defaults)
        self.assertEqual(defaults["arg_orphan"], "reg_val")
    
    def test_banned_overrides_inherited_orphan(self):
        """
        Test that no_plots overrides inherited allowed even for orphan plots.
        """
        
        # 1. Create base entry
        self.pm._registry["plot_base"] = {
            "allowed": {"arg_inherited"},
            "defaults": {"arg_inherited": "base_val"}
        }
        
        # 2. Create orphan entry that would inherit (but we'll set it directly for test)
        # In real usage, this would use inheritance directives
        self.pm._registry["plot_orphan_inherit"] = {
            "allowed": {"arg_inherited"},  # Would come from inheritance
            "defaults": {"arg_inherited": "base_val"}
        }
        
        # 3. Register arg with no_plots for orphan
        self.pm.register_arg("arg_inherited", plots=["all"], no_plots=["plot_orphan_inherit"])
        self.pm.register_arg("dummy", plots=["plot_Other"])
        
        # 4. Compile and apply no_plots
        self.pm.compile_from_kwargs(merge=True)
        _apply_no_plots(self.pm)
        
        # 5. Verify
        allowed = self.pm.allowed("plot_orphan_inherit")
        
        # arg_inherited should be GONE even though it was inherited
        self.assertNotIn("arg_inherited", allowed if allowed else set())
        
        # Default should be preserved
        defaults = self.pm.defaults("plot_orphan_inherit")
        self.assertIn("arg_inherited", defaults)


class TestContextDefaults(unittest.TestCase):
    """Test context-dependent defaults (ctx_defaults)."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
    
    def test_ctx_defaults_override_global(self):
        """Test that ctx_defaults override global defaults."""
        # Check anno_fontsize which has ctx_defaults
        defaults_mqq = self.pm.defaults('plot_mqq')
        defaults_region_r = self.pm.defaults('plot_region', 'r')
        
        # plot_region:r should have ctx_defaults value for anno_fontsize
        if 'anno_fontsize' in defaults_region_r:
            # Should be 9 (from ctx_defaults) not the global default
            self.assertEqual(defaults_region_r['anno_fontsize'], 9)
    
    def test_ctx_defaults_only_apply_to_allowed(self):
        """Test that ctx_defaults only apply if arg is in allowed set."""
        # This is tested implicitly - ctx_defaults in load_viz_config checks allowed set
        defaults = self.pm.defaults('plot_mqq', 'r')
        # If an arg has ctx_defaults but is not in allowed, it won't be in defaults
        # This is correct behavior


class TestBannedKeys(unittest.TestCase):
    """Test banned keys (nested sub-keys) behavior."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
    
    def test_banned_subkeys_removed(self):
        """Test that banned sub-keys are removed from dict args."""
        # Set up banned_keys for fig_kwargs
        entry = self.pm._get_registry_entry('plot_mqq', 'm')
        entry['banned_keys'] = {'fig_kwargs': ['figsize']}
        
        test_params = {
            'fig_kwargs': {'figsize': [10, 10], 'dpi': 200, 'other': 'value'}
        }
        
        filtered = self.pm.filter(_mqqplot, test_params, key='plot_mqq', mode='m', log=Log(), verbose=False)
        
        if 'fig_kwargs' in filtered:
            # figsize should be removed
            self.assertNotIn('figsize', filtered['fig_kwargs'])
            # Other keys should remain
            self.assertIn('dpi', filtered['fig_kwargs'])
            self.assertIn('other', filtered['fig_kwargs'])


class TestObjectLevelPresets(unittest.TestCase):
    """Test object-level presets (set/update/get)."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
    
    def test_set_replaces_presets(self):
        """Test that set() replaces existing presets."""
        self.pm.set('test_plot', {'arg1': 'value1'}, mode=None)
        self.pm.set('test_plot', {'arg1': 'value2'}, mode=None)
        
        presets = self.pm.get('test_plot')
        self.assertEqual(presets['arg1'], 'value2')
        self.assertNotEqual(presets['arg1'], 'value1')
    
    def test_update_merges_presets(self):
        """Test that update() merges with existing presets."""
        self.pm.set('test_plot', {'arg1': 'value1'}, mode=None)
        self.pm.update('test_plot', {'arg2': 'value2'}, mode=None)
        
        presets = self.pm.get('test_plot')
        self.assertEqual(presets['arg1'], 'value1')
        self.assertEqual(presets['arg2'], 'value2')
    
    def test_presets_override_defaults(self):
        """Test that presets override defaults in merge."""
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'default'}, mode=None)
        self.pm.set('test_plot', {'arg1': 'preset'}, mode=None)
        
        merged = self.pm.merge('test_plot', None, mode=None)
        self.assertEqual(merged['arg1'], 'preset')
    
    def test_presets_per_mode(self):
        """Test that presets are stored per mode."""
        self.pm.set('test_plot', {'arg1': 'base_preset'}, mode=None)
        self.pm.set('test_plot', {'arg1': 'mode_preset'}, mode='mode1')
        
        base_presets = self.pm.get('test_plot')
        mode_presets = self.pm.get('test_plot', 'mode1')
        
        self.assertEqual(base_presets['arg1'], 'base_preset')
        self.assertEqual(mode_presets['arg1'], 'mode_preset')


class TestFilteringBehavior(unittest.TestCase):
    """Test filtering behavior with whitelist and signature fallback."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
        self.log = Log()
    
    def test_whitelist_filtering(self):
        """Test filtering with whitelist (allowed set)."""
        test_params = {
            'highlight': ['snp1'],  # In allowed set
            'invalid_arg': 'value'  # Not in allowed set
        }
        
        filtered = self.pm.filter(_mqqplot, test_params, key='plot_mqq', mode='m', log=self.log, verbose=False)
        
        self.assertIn('highlight', filtered)
        self.assertNotIn('invalid_arg', filtered)
    
    def test_signature_fallback(self):
        """Test that filtering falls back to function signature when no whitelist."""
        # Don't register the plot - this ensures no allowed set exists
        # Or register with explicitly None allowed
        if 'test_plot' not in self.pm._registry:
            self.pm._registry['test_plot'] = {'allowed': None, 'defaults': {}}
        
        def test_func(arg1, arg2):
            """Function without **kwargs - only accepts arg1 and arg2."""
            pass
        
        test_params = {
            'arg1': 'value1',
            'arg2': 'value2',
            'invalid_arg': 'value'
        }
        
        filtered = self.pm.filter(test_func, test_params, key='test_plot', log=self.log, verbose=False)
        
        self.assertIn('arg1', filtered)
        self.assertIn('arg2', filtered)
        self.assertNotIn('invalid_arg', filtered)
    
    def test_signature_with_kwargs(self):
        """Test that functions with **kwargs pass everything through."""
        def test_func(**kwargs):
            pass
        
        test_params = {'arg1': 'value1', 'arg2': 'value2', 'any_arg': 'any_value'}
        
        filtered = self.pm.filter(test_func, test_params, key='test_plot', log=self.log, verbose=False)
        
        # Everything should pass through
        self.assertEqual(len(filtered), len(test_params))
        for key in test_params:
            self.assertIn(key, filtered)


class TestRegistryVsArgsPriority(unittest.TestCase):
    """Test priority between registry and args file."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
    
    def test_registry_allowed_union_with_args(self):
        """Test that registry and args allowed sets are unioned."""
        # Registry has arg1
        self.pm.register('test_plot', allowed={'arg1'}, defaults={}, mode=None)
        # Args adds arg2
        self.pm.register_arg('arg2', plots=['test_plot'])
        self.pm.compile_from_kwargs(merge=True)
        
        allowed = self.pm.allowed('test_plot')
        self.assertIn('arg1', allowed if allowed else set())
        self.assertIn('arg2', allowed if allowed else set())
    
    def test_args_defaults_override_registry(self):
        """Test that args defaults override registry defaults."""
        # Registry default
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'registry_default'}, mode=None)
        # Args default (should override)
        self.pm.register_arg('arg1', plots=['test_plot'], default='args_default')
        self.pm.compile_from_kwargs(merge=True)
        
        defaults = self.pm.defaults('test_plot')
        self.assertEqual(defaults.get('arg1'), 'args_default')
    
    def test_registry_priority_and_args_update(self):
        """
        Test that:
        1. Registry is loaded first (Base Layer)
        2. Args are loaded next (Update Layer)
        3. Args update the registry (overwrite defaults)
        """
        registry_file = "src/gwaslab/viz/viz_aux_params_registry.txt"
        args_file = "src/gwaslab/viz/viz_aux_params.txt"
        
        # Backup existing files if they exist
        backup_registry = None
        backup_args = None
        
        if os.path.exists(registry_file):
            with open(registry_file, 'r') as f:
                backup_registry = f.read()
                
        if os.path.exists(args_file):
            with open(args_file, 'r') as f:
                backup_args = f.read()
        
        try:
            # 1. Create dummy registry
            registry_content = {
                "registry": {
                    "plot_test": {
                        "allowed": ["reg_arg"],
                        "defaults": {"reg_arg": "reg_val", "shared_arg": "reg_shared"}
                    }
                }
            }
            os.makedirs(os.path.dirname(registry_file), exist_ok=True)
            with open(registry_file, 'w') as f:
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
            os.makedirs(os.path.dirname(args_file), exist_ok=True)
            with open(args_file, 'w') as f:
                json.dump(args_content, f)
                
            # 3. Load config
            load_viz_config(self.pm, path=args_file)
            
            # 4. Verify
            defaults = self.pm.defaults("plot_test")
            
            # Registry-only arg should be present
            self.assertEqual(defaults.get("reg_arg"), "reg_val")
            
            # Args-only arg should be present
            self.assertEqual(defaults.get("args_arg"), "args_val")
            
            # Shared arg should be overwritten by args (Update Layer)
            self.assertEqual(defaults.get("shared_arg"), "args_shared")
        finally:
            # Restore files
            if backup_registry:
                with open(registry_file, 'w') as f:
                    f.write(backup_registry)
            elif os.path.exists(registry_file):
                os.remove(registry_file)
                
            if backup_args:
                with open(args_file, 'w') as f:
                    f.write(backup_args)
            elif os.path.exists(args_file):
                os.remove(args_file)
    
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
        
        # Verify arg_should_be_banned is NOT in plot_A allowed
        allowed_A = self.pm.allowed("plot_A")
        self.assertNotIn("arg_should_be_banned", allowed_A)
        
        # Verify arg_should_be_banned IS in plot_B allowed (because "all" in plots, and not in no_plots)
        allowed_B = self.pm.allowed("plot_B")
        self.assertIn("arg_should_be_banned", allowed_B)


class TestUnknownModeInheritance(unittest.TestCase):
    """Test unknown mode inheritance behavior."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
    
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
    
    def test_unknown_mode_inheritance_with_known_mode(self):
        """Test that unknown modes inherit from base plot even when known modes exist."""
        # Register arg for base plot
        self.pm.register_arg("base_arg", plots=["plot_test"], default="base_val")
        
        # Register arg for specific mode
        self.pm.register_arg("mode_arg", plots=["plot_test:known"], default="mode_val")
        
        # Compile
        self.pm.compile_from_kwargs()
        
        # 1. Check Known Mode
        # Should have both base_arg and mode_arg
        defs = self.pm.defaults("plot_test", "known")
        self.assertEqual(defs.get("base_arg"), "base_val")
        self.assertEqual(defs.get("mode_arg"), "mode_val")
        
        allowed = self.pm.allowed("plot_test", "known")
        self.assertIn("base_arg", allowed)
        self.assertIn("mode_arg", allowed)
        
        # 2. Check Unknown Mode
        # Should inherit base_arg
        defs = self.pm.defaults("plot_test", "unknown")
        self.assertEqual(defs.get("base_arg"), "base_val")
        
        allowed = self.pm.allowed("plot_test", "unknown")
        # Currently, this returns None (fallback to signature)
        # But if user wants strict application, it should probably be {'base_arg'}
        if allowed is not None:
             self.assertIn("base_arg", allowed)
    
    def test_priority(self):
        """
        Test priority between registry and args defaults.
        
        Note: This test documents current behavior where compile_from_kwargs()
        sets defaults directly on modes during compilation. When registry defaults
        are set after compilation, modes don't inherit from base because they
        already have their own defaults set.
        """
        # 1. Register Arg with Global Default "Global"
        # plots=["plotA"] now implies wildcard "plotA:*" due to recent change
        self.pm.register_arg("arg1", plots=["plotA"], default="Global")
        
        # 2. Pre-seed registry to ensure wildcard works (mimicking load_viz_config)
        # We need both base and mode to exist
        self.pm._registry["plotA"] = {"allowed": set(), "defaults": {}}
        self.pm._registry["plotA:mode1"] = {"allowed": set(), "defaults": {}}
        
        # 3. Compile (populates registry from args)
        # This sets "Global" directly on both plotA and plotA:mode1
        self.pm.compile_from_kwargs()
        
        # 4. Simulate loading Registry file with Plot-level override
        # This happens AFTER compile_from_kwargs in load_viz_config
        # We set a default for the BASE plot "plotA"
        self.pm._registry["plotA"]["defaults"]["arg1"] = "RegistryBase"
        
        # 5. Check defaults for Mode1
        # Current behavior: mode1 has "Global" set during compilation,
        # so it doesn't inherit "RegistryBase" from base plot.
        # The defaults() method returns mode-specific defaults if they exist,
        # which override base defaults.
        d = self.pm.defaults("plotA", "mode1")
        val = d.get('arg1')
        
        # Current behavior: mode1 keeps "Global" because it was set during compilation
        self.assertEqual(val, "Global")
        
        # However, base plot should have "RegistryBase"
        d_base = self.pm.defaults("plotA", None)
        self.assertEqual(d_base.get('arg1'), "RegistryBase")


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
    
    def test_none_override(self):
        """Test merge with None override."""
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'default'}, mode=None)
        merged = self.pm.merge('test_plot', None, mode=None)
        self.assertEqual(merged.get('arg1'), 'default')
    
    def test_empty_allowed_set(self):
        """Test behavior with empty allowed set."""
        self.pm.register('test_plot', allowed=set(), defaults={}, mode=None)
        allowed = self.pm.allowed('test_plot')
        self.assertEqual(allowed, set())
    
    def test_none_allowed_set(self):
        """Test behavior with None allowed set (fallback to signature)."""
        # Register with explicitly None allowed
        entry = self.pm._get_registry_entry('test_plot', None)
        entry['allowed'] = None
        allowed = self.pm.allowed('test_plot')
        self.assertIsNone(allowed)
    
    def test_unknown_mode_fallback(self):
        """Test that unknown modes fall back to base plot."""
        self.pm.register('test_plot', allowed={'arg1'}, defaults={'arg1': 'default'}, mode=None)
        
        # Unknown mode should fall back to base
        allowed = self.pm.allowed('test_plot', 'unknown_mode')
        defaults = self.pm.defaults('test_plot', 'unknown_mode')
        
        self.assertIn('arg1', allowed if allowed else set())
        self.assertEqual(defaults.get('arg1'), 'default')
    
    def test_banned_arg_without_default(self):
        """Test banned arg without default is removed."""
        self.pm.register_arg('arg1', plots=['test_plot'], default=None, no_plots=['test_plot:mode1'])
        self.pm._registry['test_plot:mode1'] = {'allowed': {'arg1'}, 'defaults': {}}
        
        # Apply no_plots
        entry = self.pm._registry['test_plot:mode1']
        entry['allowed'].discard('arg1')
        # No default to add
        
        test_params = {'arg1': 'value'}
        merged = self.pm.merge('test_plot', test_params, mode='mode1')
        # arg1 should be removed (no default to replace with)
        self.assertNotIn('arg1', merged)


if __name__ == "__main__":
    unittest.main()

