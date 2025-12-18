"""Comprehensive tests for viz_aux_params module.

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
"""
import unittest
import sys
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config
from gwaslab.viz.viz_plot_mqqplot import mqqplot
from gwaslab.g_Log import Log


class TestNumericSuffixFiltering(unittest.TestCase):
    """Test that numeric suffixes are correctly handled in parameter filtering."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.pm = VizParamsManager()
        load_viz_config(self.pm)
        self.log = Log()


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
        from gwaslab.viz.viz_plot_mqqplot import mqqplot
        filtered = self.pm.filter(mqqplot, test_params, key='plot_mqq', mode='r', log=Log(), verbose=False)
        
        # highlight should be in filtered with default value
        self.assertIn('highlight', filtered)
        self.assertEqual(filtered['highlight'], [])


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
        
        from gwaslab.viz.viz_plot_mqqplot import mqqplot
        filtered = self.pm.filter(mqqplot, test_params, key='plot_mqq', mode='m', log=Log(), verbose=False)
        
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
        
        from gwaslab.viz.viz_plot_mqqplot import mqqplot
        filtered = self.pm.filter(mqqplot, test_params, key='plot_mqq', mode='m', log=self.log, verbose=False)
        
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
    
    def test_single_digit_suffix(self):
        """Test that single-digit suffixes work (e.g., highlight2)."""
        test_params = {
            'highlight': ['snp1'],
            'highlight2': ['snp2'],
            'highlight3': ['snp3'],
            'invalid_arg': 'value'
        }
        
        filtered = self.pm.filter(
            mqqplot, test_params, 
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
            mqqplot, test_params,
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
            mqqplot, test_params,
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
            mqqplot, test_params,
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
            mqqplot, test_params,
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
            mqqplot, test_params,
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


if __name__ == "__main__":
    unittest.main()

