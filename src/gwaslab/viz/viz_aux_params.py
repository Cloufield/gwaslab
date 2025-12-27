"""Visualization parameter management for gwaslab plots.

This module provides a centralized parameter management system for gwaslab visualization
functions. It manages per-plot and per-mode parameter surfaces, allowing fine-grained
control over which arguments are allowed, their default values, and how they are filtered
before being passed to plotting functions.

Key Features:
- Centralized registry of allowed parameters and defaults per plot:mode combination
- Parameter filtering with whitelist support and function signature fallback
- Banned argument handling: banned args are replaced with defaults rather than removed
- Automatic default inclusion: registered defaults are included even when not explicitly passed
- None value handling: None values preserve defaults when defaults exist
- Nested sub-key filtering for dict-type arguments (e.g., removing 'figsize' from fig_kwargs)
- Registry inheritance: registry entries can inherit from other entries
- Configuration via JSON files: `viz_aux_params.txt` (argument catalog) and
  `viz_aux_params_registry.txt` (explicit registry with inheritance support)
- Numeric suffix support: arguments like "highlight2" match base "highlight" in filtering
- Context-specific defaults and descriptions for mode-specific customizations

The system is loaded via `load_viz_config()` which processes both configuration files
and applies them in a specific order to build the final parameter registry.

================================================================================
USAGE RULES AND PRIORITIES
================================================================================

1. PARAMETER MERGE PRIORITY (from lowest to highest precedence)
   -----------------------------------------------------------
   When calling merge(key, override, mode):
   
   a) Registered defaults (from registry files)
      - Base plot defaults: `plot_mqq` → applies to all modes
      - Mode-specific defaults: `plot_mqq:r` → overrides base for mode 'r'
      - Inheritance: Mode inherits from base plot if not explicitly set
      - ALL registered defaults are included in the base, ensuring they're available
        even if not explicitly passed in override
   
   b) Object-level presets (via set()/update())
      - Override registered defaults
      - Stored per plot:mode combination
   
   c) Call-time kwargs (override parameter)
      - Highest precedence
      - Banned args are replaced with defaults (see Rule 4)
      - None values preserve defaults: if override contains None and a default exists,
        the default is kept rather than being overridden with None
      - If override contains None and no default exists, None is used
   
   Example:
     registry default: {'highlight': [], 'region_ld_colors': ['#E4E4E4', ...]}
     object preset: {'highlight': ['snp1']}
     call-time: {'highlight': ['snp2'], 'region_ld_colors': None}
     Result: {'highlight': ['snp2'], 'region_ld_colors': ['#E4E4E4', ...]}
            # call-time wins for highlight, None preserves default for region_ld_colors


2. CONFIGURATION LOADING ORDER (load_viz_config workflow)
   -------------------------------------------------------
   The system loads and applies configurations in this strict order:
   
   Step 1: Base Registry (viz_aux_params_registry.txt)
          - Loads explicit allowed keys, defaults, and banned_keys
          - Acts as foundation layer
   
   Step 2: Argument Definitions (viz_aux_params.txt)
          - Registers arguments with plots associations
          - Stores global defaults and no_plots exclusions
   
   Step 3: Compilation & Merge
          - Expands argument definitions into per-plot:mode registry entries
          - Merges with base registry (union for allowed, overwrite for defaults)
   
   Step 4: Context-Dependent Defaults (ctx_defaults)
          - Applies per-plot:mode default overrides
          - Only applies to arguments already in allowed set
          - Overrides both registry and compiled defaults
   
   Step 5: No-Plots Exclusions (no_plots)
          - Removes arguments from allowed sets for specific plots/modes
          - Keeps defaults for banned args (so they can be replaced with defaults)
          - Final override: bans apply even if arg was in registry
   
   Step 6: Banned Keys (banned_keys)
          - Applies nested sub-key bans (e.g., remove 'figsize' from fig_kwargs)
          - Can be global or context-specific (ctx_banned_keys)
   
   Priority within steps:
     - Later steps override earlier steps
     - ctx_defaults override compiled defaults
     - no_plots overrides allowed sets (but preserves defaults)


3. PLOT/MODE PARSING RULES
   ------------------------
   Plot specifications are parsed as follows:
   
   String formats:
     "plot_name"        → (plot_name, "*")     # Wildcard: all modes
     "plot_name:mode"   → (plot_name, mode)    # Specific mode
     "plot_name:"       → (plot_name, "")      # Empty mode (treated as None)
   
   List/Tuple formats:
     ["plot_name"]      → (plot_name, "*")     # Wildcard
     ["plot_name", "r"] → (plot_name, "r")     # Specific mode
     [("plot", "r")]    → (plot, "r")          # Already a tuple
   
   Special wildcards:
     "all:*"            → Applies to ALL known plots and modes
     "plot:*"           → Applies to plot with mode=None AND all specific modes
   
   Context keys (ctx_defaults, ctx_desc):
     "plot:mode"        → (plot, mode)
     "plot"             → (plot, None)         # Base plot only


4. BANNED ARGUMENTS (no_plots) BEHAVIOR
   --------------------------------------
   When an argument is banned for a plot:mode:
   
   a) It is REMOVED from the allowed set
      - Will not pass whitelist filtering
      - Cannot be explicitly passed by users
   
   b) Its DEFAULT VALUE is PRESERVED
      - Default remains in registry.defaults
      - Used to replace any attempt to pass the banned arg
   
   c) Replacement logic (in merge/filter):
      - If user passes banned arg → replaced with default value
      - If no default exists → arg is removed entirely
      - Priority: registry default > global default > remove
      - Numeric suffixes are handled: "highlight2" is banned if "highlight" is banned
        * "highlight2" → replaced with default for "highlight"
        * "highlight10" → replaced with default for "highlight"
   
   Example:
     highlight is banned for plot_mqq:r
     User passes: {'highlight': ['snp1'], 'highlight2': ['snp2']}
     Result: {'highlight': [], 'highlight2': []}  # Both replaced with default []
   
   d) Banned args still pass through filter:
      - After whitelist filtering, banned args with defaults are added back
      - Ensures function receives default value even if not in allowed set


5. WILDCARD EXPANSION RULES
   -------------------------
   When an argument has wildcard plot specification:
   
   a) "plot_name:*" expansion:
      - Applies to: (plot_name, None) + all (plot_name, specific_mode)
      - Example: "plot_mqq:*" → applies to plot_mqq, plot_mqq:r, plot_mqq:m, etc.
   
   b) "all:*" expansion:
      - Applies to ALL known (plot, mode) pairs in registry
      - Used for truly global arguments
   
   c) Expansion happens during compilation:
      - Wildcards expand to all known plot:mode pairs
      - no_plots exclusions are checked during expansion
      - Banned targets are skipped during expansion
   
   d) Mode inheritance:
      - Mode-specific entries inherit from base plot
      - Base plot defaults apply to all modes unless overridden
      - Mode-specific allowed sets do NOT inherit (must be explicit)


6. DEFAULT VALUE INHERITANCE
   --------------------------
   Default values follow this inheritance:
   
   a) Mode inherits from base plot:
      - plot_mqq:r inherits defaults from plot_mqq
      - Mode-specific defaults override base plot defaults
   
   b) Inheritance order (defaults() method):
      1. Base plot defaults (plot_mqq)
      2. Mode-specific defaults (plot_mqq:r) - overrides base
   
   c) Allowed sets do NOT inherit:
      - Each plot:mode must have explicit allowed set
      - Falls back to function signature if not in registry


7. FILTERING RULES
   ---------------
   When filtering parameters (filter() method):
   
   a) Banned args are replaced FIRST (before whitelist)
      - _apply_no_plots_filter() runs first
      - Replaces banned args with defaults
   
   b) Whitelist filtering (if allowed set exists):
      - Only args in allowed set pass through
      - Numeric suffixes supported: args ending with digits match base name
        * "highlight2" matches "highlight" if "highlight" is allowed
        * "highlight10" matches "highlight" (strips all trailing digits)
        * "scatter_kwargs1" matches "scatter_kwargs"
        * Works with any number of trailing digits: "arg123" → "arg"
      - Banned args with defaults are added back after filtering
      - Registered defaults for non-banned args are added back if missing
        (ensures defaults are available even if not in params)
   
   c) Signature fallback (if no allowed set):
      - Uses function signature to determine valid params
      - Functions with **kwargs pass everything through
      - Registered defaults for valid function parameters are added back if missing
        (ensures defaults are available even if not in params)
   
   d) Nested sub-key filtering:
      - banned_keys removes nested keys from dict-type args
      - Example: fig_kwargs['figsize'] removed if 'figsize' in banned_keys['fig_kwargs']
   
   e) Default preservation:
      - Registered defaults are preserved for valid (non-banned) parameters
      - Defaults are included even if the parameter wasn't in the original params dict
      - This ensures functions receive default values from the registry when parameters
        are not explicitly passed or are passed as None


8. REGISTRY KEY FORMAT
   --------------------
   Registry keys use format: "plot_name" or "plot_name:mode"
   
   - "plot_mqq"      → Base plot (mode=None)
   - "plot_mqq:r"    → Specific mode
   - Keys are case-sensitive
   - Mode can be any string (e.g., "r", "m", "qq", "mqq")


9. CONFIG FILE STRUCTURE
   ----------------------
   viz_aux_params.txt (argument catalog):
     {
       "args": {
         "arg_name": {
           "plots": ["plot1", "plot2:mode"],  # Where arg applies
           "no_plots": ["plot3:r"],           # Where arg is banned
           "default": value,                   # Global default
           "ctx_defaults": {                  # Context-specific defaults
             "plot:mode": value
           },
           "banned_keys": ["subkey"],          # Nested sub-keys to ban
           "ctx_banned_keys": {                # Context-specific banned sub-keys
             "plot:mode": ["subkey"]
           },
           "desc": "description",              # Argument description
           "ctx_desc": {                       # Context-specific descriptions
             "plot:mode": "description"
           },
           "add_to_docstring": false
         }
       }
     }
   
   viz_aux_params_registry.txt (explicit registry with inheritance):
     {
       "registry": {
         "plot:mode": {
           "allowed": ["arg1", "arg2"],         # Whitelist
           "defaults": {                       # Default values
             "arg1": value1
           },
           "banned_keys": {                    # Nested sub-key bans
             "arg_name": ["subkey1", "subkey2"]
           },
           "inherit": "parent_plot:mode",      # Inherit all (allowed, defaults, banned_keys)
           "allowed_inherit": "parent",        # Inherit only allowed list
           "defaults_inherit": "parent",        # Inherit only defaults
           "banned_keys_inherit": "parent"     # Inherit only banned_keys
         }
       }
     }
     
     Inheritance directives:
     - "inherit": Inherits allowed, defaults, and banned_keys from parent entry
       (child values override parent, lists are merged)
     - "allowed_inherit": Inherits only the allowed list (merged with child's allowed)
     - "defaults_inherit": Inherits only defaults (child overrides parent)
     - "banned_keys_inherit": Inherits only banned_keys (lists are merged)
     - Multiple inheritance directives can be used together
     - Inheritance is resolved recursively with cycle detection


10. BEST PRACTICES
    --------------
    - Use registry file for stable, explicit configurations
    - Use args file for dynamic, argument-centric definitions
    - Use wildcards ("plot:*") for arguments that apply to most modes
    - Use no_plots to explicitly ban args rather than omitting from plots
    - Keep defaults for banned args so functions receive valid values
    - Use ctx_defaults for mode-specific customizations
    - Use banned_keys for fine-grained control over nested dict args
"""
import inspect
import os
import json


def _parse_plot_mode(plot_spec):
    """Parse plot specification into (plot, mode) tuple.
    
    Args:
        plot_spec: str like "plot:mode" or "plot", or tuple/list [plot, mode]
    
    Returns:
        tuple: (plot_name, mode_name) where mode_name is "*" for wildcard or None
    """
    if isinstance(plot_spec, str):
        if ":" in plot_spec:
            plot_name, mode_name = plot_spec.split(":", 1)
            return (plot_name, mode_name)
        else:
            return (plot_spec, "*")  # Wildcard: apply to all modes
    elif isinstance(plot_spec, (list, tuple)):
        if len(plot_spec) == 2:
            return (plot_spec[0], plot_spec[1])
        elif len(plot_spec) == 1:
            return (plot_spec[0], "*")
    return None


def _parse_context_key(ctx_key):
    """Parse context key into (plot, mode) tuple.
    
    Args:
        ctx_key: str like "plot:mode" or "plot"
    
    Returns:
        tuple: (plot_name, mode_name or None)
    """
    if ":" in ctx_key:
        return ctx_key.split(":", 1)
    return (ctx_key, None)


def _is_banned(plot, mode, no_plots_set):
    """Check if (plot, mode) is banned by no_plots set."""
    return (plot, mode) in no_plots_set or (plot, "*") in no_plots_set


class VizParamsManager:
    """Parameter manager for visualization functions.
    
    Organizes and exposes relevant parameters for `viz_plot_*` functions.
    Supports per-plot and per-mode registrations of allowed keys and default values.
    """
    
    def __init__(self):
        self._store = {}  # Object-level presets: {key:mode -> {params}}
        self._registry = {}  # Registry: {key:mode -> {allowed: set, defaults: dict, banned_keys: dict}}
        self._arg_doc = {}  # Argument documentation
        self._arg_map = {}  # Argument definitions: {name -> {plots: set, no_plots: set, default, ctx_defaults}}
        self._arg_desc = {}  # Argument descriptions: {name -> desc}
        self._arg_ctx_desc = {}  # Context-specific descriptions: {name -> {key:mode -> desc}}
    
    def _km(self, key, mode):
        """Get registry key for plot:mode."""
        return key if mode is None else f"{key}:{mode}"
    
    def _get_registry_entry(self, key, mode):
        """Get or create registry entry for plot:mode."""
        km = self._km(key, mode)
        if km not in self._registry:
            self._registry[km] = {"allowed": set(), "defaults": {}, "banned_keys": {}}
        return self._registry[km]
    
    # ========== Core Registration Methods ==========
    
    def register(self, key, allowed=None, defaults=None, mode=None):
        """Register allowed keys and defaults for a plot and optional mode."""
        entry = self._get_registry_entry(key, mode)
        if allowed is not None:
            entry["allowed"] = set(allowed)
        if defaults is not None:
            entry["defaults"].update(defaults)
    
    def set(self, key, params, mode=None):
        """Set object-level presets for a plot/mode, replacing existing values."""
        km = self._km(key, mode)
        if params is None:
            self._store[km] = {}
        elif isinstance(params, dict):
            self._store[km] = dict(params)
        else:
            raise TypeError("params must be a dict or None")
    
    def update(self, key, params, mode=None):
        """Update object-level presets for a plot/mode."""
        km = self._km(key, mode)
        if params is not None and isinstance(params, dict):
            self._store.setdefault(km, {}).update(params)
        elif params is not None:
            raise TypeError("params must be a dict or None")
    
    def get(self, key, mode=None):
        """Get object-level presets for a plot/mode."""
        return dict(self._store.get(self._km(key, mode), {}))
    
    def defaults(self, key, mode=None):
        """Get registered default values for a plot/mode (inherits from plot-level)."""
        base = {}
        plot_key = self._km(key, None)
        if plot_key in self._registry:
            base.update(self._registry[plot_key].get("defaults", {}))
        if mode is not None:
            mode_key = self._km(key, mode)
            if mode_key in self._registry:
                base.update(self._registry[mode_key].get("defaults", {}))
        return base
    
    def allowed(self, key, mode=None):
        """Get the registered allowed key set for a plot/mode."""
        km = self._km(key, mode)
        if km in self._registry:
            return self._registry[km].get("allowed")
        if mode is not None:
            base_km = self._km(key, None)
            if base_km in self._registry:
                return self._registry[base_km].get("allowed")
        return None
    
    # ========== Merge and Filter Methods ==========
    
    def merge(self, key, override, mode=None):
        """Merge defaults, object presets, and call-time kwargs.
        
        Merge order: registered defaults -> object presets -> override kwargs.
        
        Behavior:
        - Banned args from no_plots are replaced with their default values (if defaults exist)
        - Banned args without defaults are removed from the result
        - If override contains None values and there are registered defaults, the defaults are preserved
        - If override contains None values and there are no defaults, None is used
        - All registered defaults for the key/mode are included in the base, ensuring they're available
          even if not explicitly passed in override
        
        Parameters
        ----------
        key : str
            Plot identifier (e.g., "plot_mqq", "plot_region")
        override : dict or None
            Call-time kwargs to merge. If None, returns only defaults and object presets.
        mode : str or None, optional
            Mode identifier (e.g., "r", "m", "qq")
        
        Returns
        -------
        dict
            Merged parameters with defaults, object presets, and override values.
        """
        # Check for deprecated *_args parameters
        if override is not None:
            for param_key in override:
                if param_key.endswith("_args"):
                    suggested_name = param_key.replace("_args", "_kwargs")
                    raise ValueError(
                        f"The parameter '{param_key}' is deprecated. "
                        f"Did you mean '{suggested_name}'? Please use '{suggested_name}' instead."
                    )
        
        base = self.defaults(key, mode)
        base.update(self.get(key, mode))
        if override is None:
            return base
        
        merged = dict(base)
        filtered_override = self._apply_no_plots_filter(override, key, mode)
        
        # Remove args that were in override but removed by _apply_no_plots_filter
        # These are banned args without defaults - they should not be in the final result
        for k in list(merged.keys()):
            if k in override and k not in filtered_override:
                # This arg was removed by _apply_no_plots_filter (banned with no default)
                # Remove it from merged to ensure it doesn't pass through
                merged.pop(k, None)
        
        # Update merged dict with filtered override values
        # Preserve defaults when override value is None (don't override defaults with None)
        for k, v in filtered_override.items():
            # If override value is None and there's a default in merged, keep the default
            # If override value is not None, use it (overrides default)
            # If override value is None and no default exists, use None
            if v is not None or k not in merged:
                merged[k] = v
        
        return merged
    
    def filter(self, func, params, key=None, mode=None, log=None, verbose=True):
        """Filter merged params by whitelist or function signature.
        
        Filters parameters to only include those allowed by the whitelist (if available)
        or those that match the function signature. Preserves registered defaults for
        valid parameters even if they weren't in the original params dict.
        
        Behavior:
        - Banned args are replaced with their default values (handled by _apply_no_plots_filter)
        - If whitelist exists: only allowed keys pass through
        - If no whitelist: uses function signature to determine valid parameters
        - Registered defaults for valid (non-banned) parameters are added back if missing
        - Banned args with defaults are handled by _add_banned_defaults and pass through
        
        Parameters
        ----------
        func : callable
            Target function to filter parameters for
        params : dict
            Merged parameters (typically from merge())
        key : str or None, optional
            Plot identifier for context
        mode : str or None, optional
            Mode identifier for context
        log : Log or None, optional
            Logger instance for verbose output
        verbose : bool, default True
            Whether to log filtered parameters
        
        Returns
        -------
        dict
            Filtered parameters that are valid for the function, with defaults preserved
        """
        # Replace banned args with defaults
        params = self._apply_no_plots_filter(params, key, mode)
        
        # Get defaults to preserve them after filtering
        # This ensures registered defaults are available even if not in params
        plot_defaults = self.defaults(key, mode)
        
        # Helper to check if an arg is banned (no_plots)
        def is_banned_arg(arg_name):
            if not self._arg_map or key is None:
                return False
            arg_def = self._arg_map.get(arg_name)
            if arg_def:
                p_key, m_key = key, mode
                return _is_banned(p_key, m_key, arg_def.get("no_plots", set()))
            return False
        
        allowed_keys = self.allowed(key, mode)
        if allowed_keys is not None:
            # Whitelist-based filtering: only allowed keys pass through
            filtered = self._filter_by_whitelist(params, allowed_keys, key, mode, log, verbose)
            # Add defaults for banned args so they pass through with default values
            # (banned args are not in allowed_keys but should still get defaults)
            self._add_banned_defaults(filtered, allowed_keys, key, mode)
            # Add back defaults for non-banned args that are in allowed set but missing
            # This ensures defaults are available even if not explicitly passed
            for arg_name, default_val in plot_defaults.items():
                if arg_name in allowed_keys and arg_name not in filtered and not is_banned_arg(arg_name):
                    filtered[arg_name] = default_val
            return filtered
        
        # Signature-based filtering: use function signature to determine valid parameters
        filtered = self._filter_by_signature(func, params, key, log, verbose)
        func_signature = inspect.signature(func)
        func_param_names = set(func_signature.parameters.keys())
        # Add back defaults for non-banned args that are valid function parameters but missing
        # This ensures defaults are available even if not explicitly passed
        for arg_name, default_val in plot_defaults.items():
            if arg_name in func_param_names and arg_name not in filtered and not is_banned_arg(arg_name):
                filtered[arg_name] = default_val
        return filtered
    
    def _apply_no_plots_filter(self, params, key, mode):
        """Replace banned arguments with their default values.
        
        Also handles numeric suffixes: if "highlight2" is passed and "highlight" is banned,
        "highlight2" will be replaced with the default for "highlight".
        """
        if not self._arg_map or key is None:
            return params
        
        filtered = dict(params)
        plot_defaults = self.defaults(key, mode)
        p_key, m_key = key, mode
        
        for arg_name in list(filtered.keys()):
            # Check exact match first
            arg_def = self._arg_map.get(arg_name)
            if arg_def and _is_banned(p_key, m_key, arg_def.get("no_plots", set())):
                # Replace with default value
                if arg_name in plot_defaults:
                    filtered[arg_name] = plot_defaults[arg_name]
                elif arg_def.get("default") is not None:
                    filtered[arg_name] = arg_def["default"]
                else:
                    filtered.pop(arg_name, None)
                continue
            
            # Check if arg_name has numeric suffix and base is banned
            if arg_name:
                base = arg_name.rstrip('0123456789')
                if base != arg_name:  # Has numeric suffix
                    base_def = self._arg_map.get(base)
                    if base_def and _is_banned(p_key, m_key, base_def.get("no_plots", set())):
                        # Replace with default value for base arg
                        if base in plot_defaults:
                            filtered[arg_name] = plot_defaults[base]
                        elif base_def.get("default") is not None:
                            filtered[arg_name] = base_def["default"]
                        else:
                            filtered.pop(arg_name, None)
        
        return filtered
    
    def _add_banned_defaults(self, filtered, allowed_keys, key, mode):
        """Add default values for banned args that have defaults."""
        if not self._arg_map:
            return
        
        plot_defaults = self.defaults(key, mode)
        p_key, m_key = key, mode
        
        for arg_name, default_val in plot_defaults.items():
            if arg_name not in allowed_keys and arg_name not in filtered:
                arg_def = self._arg_map.get(arg_name)
                if arg_def and _is_banned(p_key, m_key, arg_def.get("no_plots", set())):
                    filtered[arg_name] = default_val
    
    def _filter_by_whitelist(self, params, allowed_keys, key, mode, log, verbose):
        """Filter parameters using a whitelist and remove banned subkeys.
        
        Supports numeric suffixes: "highlight2" matches "highlight" if "highlight" is allowed.
        Strips all trailing digits (e.g., "highlight10" → "highlight").
        """
        target_name = key or "unknown_plot"
        
        def is_allowed_key(arg_name):
            if arg_name in allowed_keys:
                return True
            # Check if arg ends with digits and base name is allowed
            if arg_name:
                # Strip all trailing digits
                base = arg_name.rstrip('0123456789')
                if base != arg_name and base in allowed_keys:
                    return True
            return False
        
        filtered_params = {k: v for k, v in params.items() if is_allowed_key(k)}
        dropped_keys = [k for k in params.keys() if not is_allowed_key(k)]
        
        if log and dropped_keys:
            log.write(f"Filtered out args for `{target_name}`: {', '.join(dropped_keys)}", verbose=verbose)
        
        # Remove banned nested sub-keys
        entry = self._get_registry_entry(key, mode)
        banned_subkeys_map = entry.get("banned_keys", {})
        removed_subkey_paths = []
        
        for arg_name, banned_list in banned_subkeys_map.items():
            if arg_name in filtered_params and isinstance(filtered_params[arg_name], dict):
                for subkey in banned_list:
                    if subkey in filtered_params[arg_name]:
                        filtered_params[arg_name].pop(subkey, None)
                        removed_subkey_paths.append(f"{arg_name}.{subkey}")
        
        if log and removed_subkey_paths:
            log.write(f"Filtered out kwargs sub-keys for `{target_name}`: {', '.join(removed_subkey_paths)}", verbose=verbose)
        
        return filtered_params
    
    def _filter_by_signature(self, func, params, key, log, verbose):
        """Filter parameters based on the target function's signature."""
        func_signature = inspect.signature(func)
        
        # If function accepts **kwargs, pass through unchanged
        if any(p.kind == inspect.Parameter.VAR_KEYWORD for p in func_signature.parameters.values()):
            return params
        
        func_param_names = set(func_signature.parameters.keys())
        target_name = key or func.__name__
        
        filtered_params = {k: v for k, v in params.items() if k in func_param_names}
        dropped_keys = [k for k in params.keys() if k not in func_param_names]
        
        if log and dropped_keys:
            log.write(f"Filtered out args for `{target_name}`: {', '.join(dropped_keys)}", verbose=verbose)
        
        return filtered_params
    
    def public(self, key, mode=None):
        """Return exposed parameters (merged and filtered) for a plot/mode."""
        merged = self.merge(key, None, mode)
        allowed = self.allowed(key, mode)
        if allowed is None:
            return merged
        return {k: v for k, v in merged.items() if k in allowed}
    
    # ========== Argument Registration Methods ==========
    
    def register_arg(self, name, plots=None, default=None, no_plots=None):
        """Register a single argument with plots and global default."""
        entry = self._arg_map.setdefault(name, {"plots": set(), "no_plots": set(), "default": None})
        
        if plots:
            for p in plots:
                parsed = _parse_plot_mode(p)
                if parsed:
                    entry["plots"].add(parsed)
        
        if no_plots:
            for p in no_plots:
                parsed = _parse_plot_mode(p)
                if parsed:
                    entry["no_plots"].add(parsed)
        
        if default is not None:
            entry["default"] = default
        
        self._arg_map[name] = entry
    
    def attach_arg(self, name, key, mode=None):
        """Attach an existing argument to an additional plot/mode."""
        entry = self._arg_map.setdefault(name, {"plots": set(), "default": None})
        entry["plots"].add((key, mode))
    
    def args_for_plot(self, key, mode=None):
        """List argument names associated with a given plot/mode."""
        return {arg for arg, entry in self._arg_map.items() if (key, mode) in entry.get("plots", set())}
    
    # ========== Compilation Method ==========
    
    def compile_from_kwargs(self, merge=True):
        """Compile plot/mode allowed lists and defaults from registered args."""
        # Step 1: Collect known (plot, mode) pairs
        known_pairs = self._collect_known_plot_mode_pairs()
        
        # Step 2: Compile registry from args
        compiled_registry = self._compile_registry_from_args(known_pairs)
        
        # Step 3: Merge into main registry
        self._merge_compiled_registry(compiled_registry, merge)
    
    def _collect_known_plot_mode_pairs(self):
        """Collect all known (plot, mode) pairs from registry and arg_map."""
        known_pairs = set()
        
        # From existing registry
        for registry_key in self._registry.keys():
            if ":" in registry_key:
                plot_name, mode_name = registry_key.split(":", 1)
                known_pairs.add((plot_name, mode_name))
            else:
                known_pairs.add((registry_key, None))
        
        # From argument map definitions
        for entry in self._arg_map.values():
            for plot_name, mode_name in entry.get("plots", set()):
                if mode_name != "*":
                    known_pairs.add((plot_name, mode_name))
                elif plot_name != "all":
                    known_pairs.add((plot_name, None))
        
        return known_pairs
    
    def _compile_registry_from_args(self, known_pairs):
        """Compile registry entries from argument definitions."""
        compiled_registry = {}
        
        def add_to_registry(plot_name, mode_name, arg_name, default_value):
            km = self._km(plot_name, mode_name)
            if km not in compiled_registry:
                compiled_registry[km] = {"allowed": set(), "defaults": {}}
            compiled_registry[km]["allowed"].add(arg_name)
            if default_value is not None:
                compiled_registry[km]["defaults"][arg_name] = default_value
        
        for arg_name, arg_entry in self._arg_map.items():
            global_default = arg_entry.get("default")
            no_plots_set = arg_entry.get("no_plots", set())
            
            for plot_name, mode_name in arg_entry.get("plots", set()):
                if mode_name == "*":
                    # Wildcard: apply to all known modes for this plot
                    if plot_name == "all":
                        targets = list(known_pairs)
                    else:
                        targets = [(k, m) for (k, m) in known_pairs if k == plot_name]
                        if (plot_name, None) not in targets:
                            targets.append((plot_name, None))
                    
                    for p, m in targets:
                        if not _is_banned(p, m, no_plots_set):
                            add_to_registry(p, m, arg_name, global_default)
                else:
                    # Explicit mapping
                    if not _is_banned(plot_name, mode_name, no_plots_set):
                        add_to_registry(plot_name, mode_name, arg_name, global_default)
        
        return compiled_registry
    
    def _merge_compiled_registry(self, compiled_registry, merge):
        """Merge compiled registry into main registry."""
        for registry_key, compiled_data in compiled_registry.items():
            if merge and registry_key in self._registry:
                entry = self._registry[registry_key]
                if entry.get("allowed") is None:
                    entry["allowed"] = set()
                entry["allowed"].update(compiled_data["allowed"])
                entry["defaults"].update(compiled_data["defaults"])
            else:
                self._registry[registry_key] = {
                    "allowed": set(compiled_data["allowed"]),
                    "defaults": dict(compiled_data["defaults"])
                }
    
    # ========== Description Methods ==========
    
    def set_arg_desc(self, name, desc):
        """Set human-readable description for an argument."""
        self._arg_desc[name] = str(desc)
    
    def get_arg_desc(self, name):
        """Get description for an argument, or empty string if missing."""
        return self._arg_desc.get(name, "")
    
    def set_arg_desc_for(self, name, key, mode, desc):
        """Set plot/mode-specific description for an argument."""
        km = self._km(key, mode)
        self._arg_ctx_desc.setdefault(name, {})[km] = str(desc)
    
    def get_arg_desc_for(self, name, key, mode):
        """Get plot/mode-specific description for an argument with fallback."""
        km = self._km(key, mode)
        return self._arg_ctx_desc.get(name, {}).get(km) or self.get_arg_desc(name)
    
    def na_kwargs(self):
        """List argument names not attached to any plot/mode."""
        return sorted([n for n, e in self._arg_map.items() if len(e.get("plots", set())) == 0])


def load_viz_config(pm, path=None, merge=True):
    """Load visualization parameter settings from config files.
    
    Workflow:
    1. Load registry (base layer)
    2. Load args and register them
    3. Compile args into registry
    4. Apply ctx_defaults (context-dependent defaults)
    5. Apply no_plots (ban args, but keep defaults)
    6. Apply banned_keys (nested sub-keys)
    """
    # Step 1: Load Registry (Base Layer)
    _load_registry_file(pm)
    
    # Step 2: Load Args (Update Layer)
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "viz_aux_params.txt")
    
    with open(path, "r", encoding="utf-8") as f:
        config_data = json.loads(f.read())
    raw_args = config_data.get("args", {})
    
    # Register all args
    for arg_name, arg_entry in raw_args.items():
        pm.register_arg(
            arg_name,
            plots=arg_entry.get("plots", []),
            default=arg_entry.get("default"),
            no_plots=arg_entry.get("no_plots")
        )
        
        # Descriptions
        pm._arg_doc[arg_name] = bool(arg_entry.get("add_to_docstring", False))
        if "desc" in arg_entry:
            pm.set_arg_desc(arg_name, arg_entry["desc"])
        
        # Context descriptions
        for ctx_key, desc_text in arg_entry.get("ctx_desc", {}).items():
            p_name, m_name = _parse_context_key(ctx_key)
            pm.set_arg_desc_for(arg_name, p_name, m_name, desc_text)
        
        # Context defaults
        if "ctx_defaults" in arg_entry:
            pm._arg_map[arg_name]["ctx_defaults"] = arg_entry["ctx_defaults"]
    
    # Step 3: Compile and Update Registry
    pm.compile_from_kwargs(merge=merge)
    
    # Step 4: Apply ctx_defaults
    _apply_ctx_defaults(pm, raw_args)
    
    # Step 5: Apply no_plots (ban args but keep defaults)
    _apply_no_plots(pm)
    
    # Step 6: Apply banned_keys
    _apply_banned_keys(pm, raw_args)


def _load_registry_file(pm):
    """Load registry from viz_aux_params_registry.txt with inheritance support."""
    registry_file_path = os.path.join(os.path.dirname(__file__), "viz_aux_params_registry.txt")
    if not os.path.exists(registry_file_path):
        return
    
    try:
        with open(registry_file_path, "r", encoding="utf-8") as rf:
            data = json.loads(rf.read())
            registry_source = data.get("registry", {})
            
            # First pass: Load all entries without inheritance
            raw_entries = {}
            for key, entry in registry_source.items():
                if ":" in key:
                    plot_name, mode_name = key.split(":", 1)
                else:
                    plot_name, mode_name = key, None
                
                raw_entries[key] = {
                    "plot": plot_name,
                    "mode": mode_name,
                    "entry": dict(entry)  # Make a copy
                }
            
            # Second pass: Resolve inheritance
            resolved_entries = {}
            for key, raw_data in raw_entries.items():
                resolved = _resolve_inheritance(key, raw_data["entry"], raw_entries, resolved_entries)
                resolved_entries[key] = resolved
            
            # Third pass: Apply resolved entries to registry
            for key, resolved_entry in resolved_entries.items():
                raw_data = raw_entries[key]
                plot_name = raw_data["plot"]
                mode_name = raw_data["mode"]
                target = pm._get_registry_entry(plot_name, mode_name)
                
                if "allowed" in resolved_entry and resolved_entry["allowed"] is not None:
                    target["allowed"] = set(resolved_entry["allowed"])
                if "defaults" in resolved_entry:
                    target["defaults"].update(resolved_entry["defaults"])
                if "banned_keys" in resolved_entry:
                    target.setdefault("banned_keys", {}).update(resolved_entry["banned_keys"])
    except Exception:
        pass


def _resolve_inheritance(key, entry, raw_entries, resolved_entries, visited=None):
    """Resolve inheritance for a registry entry.
    
    Supports:
    - "inherit": inherits allowed, banned_keys, and defaults from another entry
    - "allowed_inherit": inherits only allowed list
    - "banned_keys_inherit": inherits only banned_keys
    - "defaults_inherit": inherits only defaults
    
    Args:
        key: Registry key (e.g., "plot_manhattan")
        entry: Entry dictionary with potential inheritance directives
        raw_entries: All raw entries from registry file
        resolved_entries: Cache of already resolved entries
        visited: Set of keys being resolved (for cycle detection)
    
    Returns:
        Resolved entry dictionary with inheritance applied
    """
    if visited is None:
        visited = set()
    
    if key in resolved_entries:
        return resolved_entries[key]
    
    if key in visited:
        # Cycle detected, return entry without inheritance
        return dict(entry)
    
    visited.add(key)
    resolved = dict(entry)  # Start with a copy
    
    # Handle general inheritance (inherits all)
    if "inherit" in entry:
        parent_key = entry["inherit"]
        if parent_key in raw_entries:
            parent_entry = _resolve_inheritance(parent_key, raw_entries[parent_key]["entry"], 
                                               raw_entries, resolved_entries, visited.copy())
            # Inherit allowed
            if "allowed" not in resolved or resolved["allowed"] is None:
                resolved["allowed"] = parent_entry.get("allowed", [])
            elif isinstance(resolved["allowed"], list):
                # Merge: combine parent and child allowed lists
                parent_allowed = set(parent_entry.get("allowed", []))
                child_allowed = set(resolved["allowed"])
                resolved["allowed"] = list(parent_allowed | child_allowed)
            
            # Inherit defaults
            if "defaults" not in resolved:
                resolved["defaults"] = {}
            parent_defaults = parent_entry.get("defaults", {})
            # Merge defaults (child overrides parent)
            resolved["defaults"] = {**parent_defaults, **resolved["defaults"]}
            
            # Inherit banned_keys
            if "banned_keys" not in resolved:
                resolved["banned_keys"] = {}
            parent_banned = parent_entry.get("banned_keys", {})
            # Merge banned_keys (child overrides parent)
            merged_banned = dict(parent_banned)
            for banned_key, banned_list in resolved.get("banned_keys", {}).items():
                if banned_key in merged_banned:
                    # Merge lists
                    merged_banned[banned_key] = list(set(merged_banned[banned_key]) | set(banned_list))
                else:
                    merged_banned[banned_key] = banned_list
            resolved["banned_keys"] = merged_banned
        
        # Remove inherit directive from resolved entry
        resolved.pop("inherit", None)
    
    # Handle specific inheritance directives
    if "allowed_inherit" in entry:
        parent_key = entry["allowed_inherit"]
        if parent_key in raw_entries:
            parent_entry = _resolve_inheritance(parent_key, raw_entries[parent_key]["entry"],
                                               raw_entries, resolved_entries, visited.copy())
            parent_allowed = parent_entry.get("allowed", [])
            if "allowed" not in resolved or resolved["allowed"] is None:
                resolved["allowed"] = parent_allowed
            elif isinstance(resolved["allowed"], list):
                # Merge: combine parent and child
                resolved["allowed"] = list(set(parent_allowed) | set(resolved["allowed"]))
        resolved.pop("allowed_inherit", None)
    
    if "banned_keys_inherit" in entry:
        parent_key = entry["banned_keys_inherit"]
        if parent_key in raw_entries:
            parent_entry = _resolve_inheritance(parent_key, raw_entries[parent_key]["entry"],
                                               raw_entries, resolved_entries, visited.copy())
            parent_banned = parent_entry.get("banned_keys", {})
            if "banned_keys" not in resolved:
                resolved["banned_keys"] = {}
            # Merge banned_keys
            merged_banned = dict(parent_banned)
            for banned_key, banned_list in resolved.get("banned_keys", {}).items():
                if banned_key in merged_banned:
                    merged_banned[banned_key] = list(set(merged_banned[banned_key]) | set(banned_list))
                else:
                    merged_banned[banned_key] = banned_list
            resolved["banned_keys"] = merged_banned
        resolved.pop("banned_keys_inherit", None)
    
    if "defaults_inherit" in entry:
        parent_key = entry["defaults_inherit"]
        if parent_key in raw_entries:
            parent_entry = _resolve_inheritance(parent_key, raw_entries[parent_key]["entry"],
                                               raw_entries, resolved_entries, visited.copy())
            parent_defaults = parent_entry.get("defaults", {})
            if "defaults" not in resolved:
                resolved["defaults"] = {}
            # Merge defaults (child overrides parent)
            resolved["defaults"] = {**parent_defaults, **resolved["defaults"]}
        resolved.pop("defaults_inherit", None)
    
    visited.remove(key)
    resolved_entries[key] = resolved
    return resolved


def _apply_ctx_defaults(pm, raw_args):
    """Apply context-specific default overrides."""
    for arg_name, arg_entry in raw_args.items():
        if "ctx_defaults" not in arg_entry:
            continue
        
        for ctx_key, ctx_value in arg_entry["ctx_defaults"].items():
            ctx_plot, ctx_mode = _parse_context_key(ctx_key)
            entry = pm._registry.get(pm._km(ctx_plot, ctx_mode))
            
            if entry and entry.get("allowed") and arg_name in entry["allowed"]:
                entry["defaults"][arg_name] = ctx_value


def _apply_no_plots(pm):
    """Apply no_plots exclusions (remove from allowed but keep defaults)."""
    for registry_key in list(pm._registry.keys()):
        p_key, m_key = _parse_context_key(registry_key) if ":" in registry_key else (registry_key, None)
        entry = pm._registry[registry_key]
        current_allowed = entry.get("allowed", set())
        
        if not current_allowed:
            continue
        
        args_to_remove = set()
        for arg_name in current_allowed:
            arg_def = pm._arg_map.get(arg_name)
            if arg_def and _is_banned(p_key, m_key, arg_def.get("no_plots", set())):
                args_to_remove.add(arg_name)
        
        if args_to_remove:
            entry["allowed"] -= args_to_remove
            # Ensure defaults exist for banned args
            for banned_arg in args_to_remove:
                if banned_arg not in entry["defaults"]:
                    arg_def = pm._arg_map.get(banned_arg)
                    if arg_def and arg_def.get("default") is not None:
                        entry["defaults"][banned_arg] = arg_def["default"]


def _apply_banned_keys(pm, raw_args):
    """Apply banned_keys rules for nested sub-keys."""
    # Collect all known plot/mode pairs
    known_pairs = set()
    for k in pm._registry.keys():
        if ":" in k:
            p, m = k.split(":", 1)
            known_pairs.add((p, m))
        else:
            known_pairs.add((k, None))
    
    for arg_name, arg_entry in raw_args.items():
        # Get target contexts for this arg
        normalized_plots = pm._arg_map.get(arg_name, {}).get("plots", set())
        target_contexts = []
        for p_name, m_name in normalized_plots:
            if m_name == "*":
                if p_name == "all":
                    target_contexts.extend(known_pairs)
                else:
                    target_contexts.extend([(k, m) for (k, m) in known_pairs if k == p_name])
            else:
                target_contexts.append((p_name, m_name))
        
        # Apply global banned_keys
        global_banned = arg_entry.get("banned_keys", [])
        if global_banned:
            for p_name, m_name in target_contexts:
                entry = pm._get_registry_entry(p_name, m_name)
                entry.setdefault("banned_keys", {}).setdefault(arg_name, [])
                entry["banned_keys"][arg_name] = list(set(entry["banned_keys"][arg_name]) | set(global_banned))
        
        # Apply context-specific banned_keys
        ctx_banned_map = arg_entry.get("ctx_banned_keys", {})
        if isinstance(ctx_banned_map, dict):
            for ctx_key, banned_list in ctx_banned_map.items():
                p_name, m_name = _parse_context_key(ctx_key)
                entry = pm._get_registry_entry(p_name, m_name)
                entry.setdefault("banned_keys", {}).setdefault(arg_name, [])
                entry["banned_keys"][arg_name] = list(set(entry["banned_keys"][arg_name]) | set(banned_list))
