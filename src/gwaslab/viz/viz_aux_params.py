"""Visualization parameter management for gwaslab plots.

Provides a centralized registry of per-plot and per-mode parameter surfaces:
- Allowed keys and defaults are stored under `key[:mode]`.
- Filtering preserves only whitelisted params (or those in the target function
  signature when no whitelist exists) and removes banned nested sub-keys inside
  dict-type args.
- Configuration is loaded from `viz_aux_params.txt` (argument catalog) and
  `viz_aux_params_registry.txt` (registry) via `load_viz_config`.
- Arg-centric APIs (`register_arg`, `compile_from_kwargs`) can generate registry
  entries from a single catalog; plot-centric APIs (`register`, `merge`, `filter`)
  manage calls to plotting functions.

Priority of parameters during a call:
1) Call-time kwargs (`override`)
3) Registered defaults (`register`/compiled from config)

Config files:
- `viz_aux_params.txt` defines the argument catalog:
  plots/modes associations (`plots`), global defaults (`default`), per-context
  overrides (`ctx_defaults`), per-context descriptions (`ctx_desc`), and nested
  `banned_keys` rules (global and context-specific). These entries are compiled
  into per-plot/mode `allowed` and `defaults` via `compile_from_kwargs`.
- `viz_aux_params_registry.txt` provides the explicit registry:
  final per-plot/mode `allowed`, `defaults`, and `banned_keys`. When present,
  it replaces compiled `allowed` and `banned_keys` and merges `defaults` (if
  `merge=True`), ensuring stable and precise parameter surfaces.
"""
import inspect
import os
import json


class VizParamsManager:
    """Parameter manager for visualization functions.

    use files to supply parameters

    Organizes and exposes relevant parameters for `viz_plot_*` functions.
    Supports per-plot and per-mode (e.g., "qq", "r", "b") registrations
    of allowed keys and default values, and merges them with object-level
    presets and call-time kwargs. Irrelevant parameters are filtered out
    before passing to the underlying plotting function.

    Priority
    1) Call-time kwargs (`override`) have highest precedence.
    2) Object-level presets set via `set()`/`update()` override defaults.
    3) Registered defaults provide base values when nothing else is specified.

    Config Precedence (from `viz_aux_params.txt`)
    1) `args` section loads argument-to-plot mappings and base defaults
       via `register_arg()` and `compile_from_kwargs()`.
    2) `args.ctx_defaults` applies per-plot/mode default overrides on top
       of compiled defaults.
    3) `registry` section is merged last and wins on conflicts:
       - `allowed` extends/overrides prior allowed keys
       - `defaults` override previous defaults
       - `banned_keys` replaces earlier maps when provided in `registry`.

    Args Entry Keys (from `viz_aux_params.txt`)
    plots : List of `[plot, mode]` pairs this arg applies to; `null` mode means None.
            Drives whitelist compilation (`compile_from_kwargs`) and global `banned_keys` propagation.
    desc : Human-readable argument description; used for documentation.
    default : Global default value attached when compiling into per-plot/mode registry.
    ctx_defaults : Per-context overrides mapping `plot:mode` → value; merged after compilation.
    ctx_desc : Per-context description mapping `plot:mode` → text; used by `get_arg_desc_for`.
    banned_keys : List of nested sub-keys to remove from this arg's dict across listed plots.
    ctx_banned_keys : Per-context nested sub-keys to remove, mapping `plot:mode` → list.
    add_to_docstring : Boolean flag; when true, include this arg in `build_docstring_for_plot` output.

    Notes
    - The complete catalog, including plot/mode-specific context descriptions,
      is loaded from `viz_aux_params.txt` via `load_viz_config()`.
    - Use `build_docstring_for_plot(key, mode)` to render per-plot/mode
      parameter documentation assembled from the config.

    Attributes
    ----------
    _store : dict
        Object-level presets keyed by `<key>:<mode>` or `<key>`.
    _registry : dict
        Allowed keys and defaults for plots keyed by `<key>:<mode>` or `<key>`.

    Usage
    -----
    - Register plot/mode: `register('plot_region', allowed={...}, defaults={...}, mode='r')`
    - Set object presets: `set('plot_region', {'title': 'Locus'}, mode='r')`
    - Merge and filter on call: `filter(func, merge('plot_region', kwargs, mode='r'), key='plot_region', mode='r')`
    """
    def __init__(self):
        self._store = {}
        self._registry = {}
        self._arg_doc = {}

    def _km(self, key, mode):
        if mode is None:
            return key
        return f"{key}:{mode}"

    def register(self, key, allowed=None, defaults=None, mode=None):
        """Register allowed keys and defaults for a plot and optional mode.

        Parameters
        ----------
        key : str
            Plot identifier (e.g., 'plot_mqq').
        allowed : set or list or None
            Whitelist of exposed parameter names for this plot/mode.
            If None, filtering falls back to function signature.
        defaults : dict or None
            Default parameter values for this plot/mode.
        mode : str or None
            Optional sub-mode (e.g., 'qq', 'r', 'b').
        """
        km = self._km(key, mode)
        a = set(allowed) if allowed is not None else None
        d = dict(defaults) if defaults is not None else {}
        self._registry[km] = {"allowed": a, "defaults": d}

    def set(self, key, params, mode=None):
        """Set object-level presets for a plot/mode, replacing existing values.

        Parameters
        ----------
        key : str
            Plot identifier.
        params : dict or None
            Preset values; None resets to empty.
        mode : str or None
            Optional sub-mode.
        """
        km = self._km(key, mode)
        if params is None:
            self._store[km] = {}
        elif isinstance(params, dict):
            self._store[km] = dict(params)
        else:
            raise TypeError("params must be a dict or None")

    def update(self, key, params, mode=None):
        """Update object-level presets for a plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        params : dict or None
            Values to update; None keeps existing.
        mode : str or None
            Optional sub-mode.
        """
        km = self._km(key, mode)
        base = self._store.get(km, {})
        if params is None:
            self._store[km] = base
        elif isinstance(params, dict):
            base.update(params)
            self._store[km] = base
        else:
            raise TypeError("params must be a dict or None")

    def get(self, key, mode=None):
        """Get object-level presets for a plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        dict
            A copy of stored presets (may be empty).
        """
        km = self._km(key, mode)
        return dict(self._store.get(km, {}))

    def defaults(self, key, mode=None):
        """Get registered default values for a plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        dict
            A copy of registered defaults (may be empty).
        """
        # 1. Start with empty base
        base = {}

        # 2. Update with plot-level defaults (if key is registered)
        plot_key = self._km(key, None)
        if plot_key in self._registry:
            plot_defaults = self._registry[plot_key].get("defaults", {})
            base.update(plot_defaults)

        # 3. Update with mode-specific defaults (if mode is provided)
        if mode is not None:
            # Check for specific mode
            mode_key = self._km(key, mode)
            if mode_key in self._registry:
                mode_defaults = self._registry[mode_key].get("defaults", {})
                base.update(mode_defaults)
            # If specific mode not found, check if base plot has defaults
            # (Inheritance for unknown modes)
            elif plot_key in self._registry:
                 # Re-apply base defaults (already done in step 2, but just to be sure)
                 pass
        
        return base
        
    def allowed(self, key, mode=None):
        """Get the registered allowed key set for a plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        set or None
            Allowed keys; None means use function signature filtering.
        """
        # If mode is unknown (not in registry), fallback to plot-level allowed
        # This is for "plot_mqq" matching "plot_mqq:unknown"
        
        km = self._km(key, mode)
        if km in self._registry:
            return self._registry[km].get("allowed")
            
        # If specific mode not found, check base plot
        if mode is not None:
            base_km = self._km(key, None)
            if base_km in self._registry:
                return self._registry[base_km].get("allowed")
                 
        return None

    def merge(self, key, override, mode=None):
        """Merge defaults, object presets, and call-time kwargs.

        Merge order: registered defaults -> object presets -> override kwargs.

        Parameters
        ----------
        key : str
            Plot identifier.
        override : dict or None
            Call-time kwargs to override.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        dict
            Merged parameter dictionary.
        """
        base = self.defaults(key, mode)
        base.update(self.get(key, mode))
        if override is None:
            return base
        merged = dict(base)
        merged.update(override)
        return merged


    def filter(self, func, params, key=None, mode=None, log=None, verbose=True):
        """Filter merged params by whitelist or function signature.

        Parameters
        ----------
        func : callable
            Target plotting function to inspect when whitelist is absent.
        params : dict
            Parameters to be filtered.
        key : str or None
            Plot identifier, used to look up whitelist.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        dict
            Parameters safe to pass to `func`.
        """
        # 0) Enforce no_plots exclusions at filter time (defensive)
        params = self._apply_no_plots_filter(params, key, mode)
        # 1) Look up whitelist for this plot/mode
        allowed_keys = self.allowed(key, mode)
        # 2) If whitelist exists, filter by allowed keys
        if allowed_keys is not None:
            return self._filter_by_whitelist(params, allowed_keys, key, mode, log, verbose)
        
        # 3) Fallback: filter by function signature
        return self._filter_by_signature(func, params, key, log, verbose)

    def _apply_no_plots_filter(self, params, key, mode):
        """Remove arguments banned by no_plots for the given plot/mode."""
        if not hasattr(self, "_arg_map"):
            return params
        if key is None:
            return params
        p_key, m_key = key, mode
        filtered = dict(params)
        for arg_name in list(filtered.keys()):
            arg_def = self._arg_map.get(arg_name)
            if not arg_def:
                continue
            no_plots = arg_def.get("no_plots", set())
            if (p_key, m_key) in no_plots or (p_key, "*") in no_plots:
                filtered.pop(arg_name, None)
        return filtered

    def _filter_by_whitelist(self, params, allowed_keys, key, mode, log, verbose):
        """Filter parameters using a whitelist and remove banned subkeys."""
        target_name = key or "unknown_plot"
        
        # Helper for numeric suffix handling (e.g. "arg1" matches "arg")
        def is_allowed_key(arg_name):
            if arg_name in allowed_keys:
                return True
            if arg_name[-1:].isdigit():
                base = arg_name[:-1]
                return base in allowed_keys
            return False

        # Apply whitelist
        filtered_params = {k: v for k, v in params.items() if is_allowed_key(k)}
        dropped_keys = [k for k in params.keys() if not is_allowed_key(k)]
        
        if log is not None and dropped_keys:
            log.write(f"Filtered out args for `{target_name}`: {', '.join(dropped_keys)}", verbose=verbose)

        # Remove banned nested sub-keys inside dict-type args
        registry_entry = self._registry.get(self._km(key, mode), {})
        banned_subkeys_map = registry_entry.get("banned_keys", {})
        removed_subkey_paths = []
        
        for arg_name, banned_list in banned_subkeys_map.items():
            if arg_name in filtered_params and isinstance(filtered_params[arg_name], dict):
                for subkey in banned_list:
                    if subkey in filtered_params[arg_name]:
                        filtered_params[arg_name].pop(subkey, None)
                        removed_subkey_paths.append(f"{arg_name}.{subkey}")
                        
        if log is not None and removed_subkey_paths:
            log.write(f"Filtered out kwargs sub-keys for `{target_name}`: {', '.join(removed_subkey_paths)}", verbose=verbose)
            
        return filtered_params

    def _filter_by_signature(self, func, params, key, log, verbose):
        """Filter parameters based on the target function's signature."""
        func_signature = inspect.signature(func)
        
        # If function accepts **kwargs, pass through unchanged
        has_var_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in func_signature.parameters.values())
        if has_var_kwargs:
            return params
            
        # Apply signature-based filtering
        func_param_names = set(func_signature.parameters.keys())
        target_name = key or func.__name__
        
        filtered_params = {k: v for k, v in params.items() if k in func_param_names}
        dropped_keys = [k for k in params.keys() if k not in func_param_names]
        
        if log is not None and dropped_keys:
            log.write(f"Filtered out args for `{target_name}`: {', '.join(dropped_keys)}", verbose=verbose)
            
        return filtered_params

    def public(self, key, mode=None):
        """Return exposed parameters (merged and filtered) for a plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        dict
            Public-facing parameters with defaults applied.
        """
        merged = self.merge(key, None, mode)
        allowed = self.allowed(key, mode)
        if allowed is None:
            return merged
        return {k: v for k, v in merged.items() if k in allowed}


    def register_arg(self, name, plots=None, default=None, no_plots=None):
        """Register a single argument with plots and global default.

        Parameters
        ----------
        name : str
            Argument name.
        plots : list, optional
            List of plot strings.
            Format: "plot_name" or "plot_name:mode".
            If a string is provided without mode (e.g., "plot_mqq"),
            it applies to "plot_mqq" with mode=None.
        default : Any, optional
            Default value to apply for plots that include this argument.
        no_plots : list, optional
            List of plot strings to exclude.
        """
        if not hasattr(self, "_arg_map"):
            self._arg_map = {}
        if not hasattr(self, "_arg_desc"):
            self._arg_desc = {}
        entry = self._arg_map.get(name, {"plots": set(), "no_plots": set(), "default": None})
        if "no_plots" not in entry:
            entry["no_plots"] = set()
            
        if plots:
            for p in plots:
                if isinstance(p, str):
                    if ":" in p:
                        # "plot:mode" -> (plot, mode)
                        plot_name, mode_name = p.split(":", 1)
                        entry["plots"].add((plot_name, mode_name))
                    else:
                        # "plot" -> (plot, *) (Wildcard: apply to all modes)
                        entry["plots"].add((p, "*"))
                elif isinstance(p, (list, tuple)):
                    if len(p) == 2:
                        entry["plots"].add((p[0], p[1]))
                    elif len(p) == 1:
                        # ["plot"] -> (plot, *)
                        entry["plots"].add((p[0], "*"))
        
        if no_plots:
            for p in no_plots:
                if isinstance(p, str):
                    if ":" in p:
                        plot_name, mode_name = p.split(":", 1)
                        entry["no_plots"].add((plot_name, mode_name))
                    else:
                        entry["no_plots"].add((p, "*"))
                elif isinstance(p, (list, tuple)):
                    if len(p) == 2:
                        entry["no_plots"].add((p[0], p[1]))
                    elif len(p) == 1:
                        entry["no_plots"].add((p[0], "*"))

        # If the argument is registered with default value, add it to global defaults
        if default is not None:
            entry["default"] = default
            
        self._arg_map[name] = entry

    def attach_arg(self, name, key, mode=None):
        """Attach an existing argument to an additional plot/mode.

        Parameters
        ----------
        name : str
            Argument name.
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.
        """
        if not hasattr(self, "_arg_map"):
            self._arg_map = {}
        entry = self._arg_map.get(name, {"plots": set(), "default": None})
        entry["plots"].add((key, mode))
        self._arg_map[name] = entry

    def args_for_plot(self, key, mode=None):
        """List argument names associated with a given plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        set
            Names of arguments associated with `(key, mode)`.
        """
        if not hasattr(self, "_arg_map"):
            self._arg_map = {}
        out = set()
        for arg, entry in self._arg_map.items():
            if (key, mode) in entry["plots"]:
                out.add(arg)
        return out

    def compile_from_kwargs(self, merge=True):
        """Compile plot/mode allowed lists and defaults from registered args.

        Parameters
        ----------
        merge : bool, default=True
            If True, merge into existing registry entries; otherwise replace.
        """
        if not hasattr(self, "_arg_map"):
            self._arg_map = {}
        
        # Dictionary to store the compiled registry before merging
        # Structure: { "plot:mode": { "allowed": set(), "defaults": {} } }
        compiled_registry = {}

        # Helper to add an argument to the compiled registry
        def _add_to_compiled_registry(plot_name, mode_name, argument_name, default_value):
            registry_key = self._km(plot_name, mode_name)
            
            if registry_key not in compiled_registry:
                compiled_registry[registry_key] = {"allowed": set(), "defaults": {}}
            
            compiled_registry[registry_key]["allowed"].add(argument_name)
            if default_value is not None:
                compiled_registry[registry_key]["defaults"][argument_name] = default_value

        # Step 1: Identify all known (plot, mode) pairs
        # We need this to resolve wildcards like "plot_name:*" or "all:*"
        known_plot_mode_pairs = set()

        # 1a. From existing registry
        for registry_key in self._registry.keys():
            if ":" in registry_key:
                plot_name, mode_name = registry_key.split(":", 1)
                known_plot_mode_pairs.add((plot_name, mode_name))
            else:
                known_plot_mode_pairs.add((registry_key, None))

        # 1b. From argument map definitions
        for argument_entry in self._arg_map.values():
            for plot_name, mode_name in argument_entry["plots"]:
                if mode_name != "*":
                    known_plot_mode_pairs.add((plot_name, mode_name))
                else:
                    # If an argument is registered for "plot_name" (wildcard), 
                    # we should at least acknowledge the base plot "plot_name" exists.
                    if plot_name != "all":
                        known_plot_mode_pairs.add((plot_name, None))

        # Step 2: Iterate through each argument and populate the compiled registry
        for argument_name, argument_entry in self._arg_map.items():
            global_default = argument_entry["default"]
            no_plots_set = argument_entry.get("no_plots", set())
            
            def is_banned(p, m):
                if (p, m) in no_plots_set: return True
                if (p, "*") in no_plots_set: return True
                return False
            
            # 2a. Determine target (plot, mode) pairs
            
            for plot_name, mode_name in argument_entry["plots"]:
                 if mode_name == "*":
                     # Wildcard handling: apply to all known modes for this plot
                     # plus the base plot itself (None mode)
                     
                     # 1. Base plot
                     if plot_name == "all":
                         # "all:*" -> apply to everything
                         targets = list(known_plot_mode_pairs)
                     else:
                         # "plot:*" -> apply to (plot, None) and all (plot, specific_mode)
                         targets = [(k, m) for (k, m) in known_plot_mode_pairs if k == plot_name]
                         # Ensure base is included if not already in known pairs (though it should be)
                         if (plot_name, None) not in targets:
                             targets.append((plot_name, None))
                     
                     for p, m in targets:
                         if not is_banned(p, m):
                             _add_to_compiled_registry(p, m, argument_name, global_default)
                 else:
                     # Explicit mapping
                     if not is_banned(plot_name, mode_name):
                         _add_to_compiled_registry(plot_name, mode_name, argument_name, global_default)

            # 2b. Apply context-specific default overrides (ctx_defaults)
            if "ctx_defaults" in argument_entry:
                for ctx_key, ctx_value in argument_entry["ctx_defaults"].items():
                    if ":" in ctx_key:
                        ctx_plot, ctx_mode = ctx_key.split(":", 1)
                    else:
                        ctx_plot, ctx_mode = ctx_key, None
                    
                    # Update default if the argument is already allowed for this context
                    # (Note: Logic assumes ctx_defaults only apply if arg is allowed)
                    registry_key = self._km(ctx_plot, ctx_mode)
                    if registry_key in compiled_registry:
                         if argument_name in compiled_registry[registry_key]["allowed"]:
                             compiled_registry[registry_key]["defaults"][argument_name] = ctx_value

        # Step 3: Merge compiled results into the main registry
        for registry_key, compiled_data in compiled_registry.items():
            new_allowed = compiled_data["allowed"]
            new_defaults = compiled_data["defaults"]
            
            if merge and registry_key in self._registry:
                existing_entry = self._registry[registry_key]
                
                # Update allowed set
                if existing_entry["allowed"] is None:
                    existing_entry["allowed"] = set()
                existing_entry["allowed"].update(new_allowed)
                
                # Update defaults - OVERWRITE existing defaults
                # User request: "use src/gwaslab/viz/viz_aux_params.txt to update it"
                for arg_name, default_val in new_defaults.items():
                    existing_entry["defaults"][arg_name] = default_val
            else:
                # Create new entry
                self._registry[registry_key] = {
                    "allowed": set(new_allowed), 
                    "defaults": dict(new_defaults)
                }
        
        # Step 4: Force-apply no_plots exclusions across ENTIRE registry
        # This handles cases where registry has an allowed key, but args say "no_plots" for it,
        # even if that arg wasn't being actively added to this plot in Step 2.
        
        # 4a. Build a list of all known plot/mode keys from registry + compiled
        all_keys = set(self._registry.keys())
        
        for registry_key in all_keys:
            if ":" in registry_key:
                p_key, m_key = registry_key.split(":", 1)
            else:
                p_key, m_key = registry_key, None
                
            entry = self._registry[registry_key]
            current_allowed = entry.get("allowed", set())
            if not current_allowed:
                continue
                
            # Check every argument currently allowed to see if it should be banned
            # We iterate over args instead of allowed keys because we need to check no_plots from arg definition
            
            # Optimization: Only check args that are actually in current_allowed
            args_to_remove = set()
            
            for arg_name in current_allowed:
                arg_def = self._arg_map.get(arg_name)
                if not arg_def:
                    continue
                
                no_plots = arg_def.get("no_plots", set())
                if not no_plots:
                    continue
                    
                is_banned = False
                # Check exact match
                if (p_key, m_key) in no_plots: 
                    is_banned = True
                # Check wildcard match (plot:*)
                elif (p_key, "*") in no_plots: 
                    is_banned = True
                
                if is_banned:
                    args_to_remove.add(arg_name)
            
            # Apply removals
            if args_to_remove:
                entry["allowed"] -= args_to_remove
                # Also remove from defaults if present
                for banned_arg in args_to_remove:
                    entry["defaults"].pop(banned_arg, None)

    def na_kwargs(self):
        """List argument names not attached to any plot/mode.

        Returns
        -------
        list
            Sorted list of orphan argument names for manual fixing.
        """
        if not hasattr(self, "_arg_map"):
            self._arg_map = {}
        return sorted([n for n, e in self._arg_map.items() if len(e["plots"]) == 0])

    def set_arg_desc(self, name, desc):
        """Set human-readable description for an argument.

        Parameters
        ----------
        name : str
            Argument name.
        desc : str
            Description text for documentation reuse.
        """
        if not hasattr(self, "_arg_desc"):
            self._arg_desc = {}
        self._arg_desc[name] = str(desc)

    def get_arg_desc(self, name):
        """Get description for an argument, or empty string if missing."""
        if not hasattr(self, "_arg_desc"):
            self._arg_desc = {}
        return self._arg_desc.get(name, "")

    def set_arg_desc_for(self, name, key, mode, desc):
        """Set plot/mode-specific description for an argument."""
        if not hasattr(self, "_arg_ctx_desc"):
            self._arg_ctx_desc = {}
        km = self._km(key, mode)
        per_arg = self._arg_ctx_desc.get(name, {})
        per_arg[km] = str(desc)
        self._arg_ctx_desc[name] = per_arg

    def get_arg_desc_for(self, name, key, mode):
        """Get plot/mode-specific description for an argument with fallback."""
        if not hasattr(self, "_arg_ctx_desc"):
            self._arg_ctx_desc = {}
        km = self._km(key, mode)
        d = self._arg_ctx_desc.get(name, {}).get(km)
        if d:
            return d
        return self.get_arg_desc(name)


def load_viz_config(pm, path=None, merge=True):
    """Load visualization parameter settings from config files.

    This function initializes the VizParamsManager with configuration from two sources,
    establishing a layered parameter system where argument definitions can update and
    expand upon a base registry.

    Workflow & Priorities
    ---------------------
    1. Base Registry Layer (`viz_aux_params_registry.txt`)
       - Loads the foundational configuration: explicit allowed keys, base defaults,
         and banned sub-keys for plots/modes.
       - Acts as the initial state of the parameter registry.

    2. Argument Definition Layer (`viz_aux_params.txt`)
       - Loads argument definitions, including:
         - Plot associations (wildcards supported)
         - Global defaults and context-specific defaults (`ctx_defaults`)
         - Plot exclusions (`no_plots`)
         - Descriptions for documentation
       - These definitions are compiled into the registry via `compile_from_kwargs`.

    3. Compilation & Merge (Update Step)
       - The arguments from Step 2 are expanded into per-plot/mode settings.
       - These new settings are merged into the Base Registry from Step 1.
       - **Priority Rule**: Defaults defined in the Argument Layer (Step 2) OVERWRITE
         defaults from the Base Registry (Step 1). This allows `viz_aux_params.txt`
         to dynamically update and refine the static registry.
       - Allowed keys are unioned (Registry + Args).
       - **Override Rule**: If an argument is listed in `no_plots` (Args Layer), it
         is actively removed from the allowed set, overriding any presence in the
         Registry Layer.

    4. Banned Keys Application
       - Applies `banned_keys` rules from the Argument Layer to the registry,
         ensuring nested dictionary sub-keys are removed where specified.

    Parameters
    ----------
    pm : VizParamsManager
        The parameter manager instance to populate.
    path : str, optional
        Path to the argument definition file (default: `viz_aux_params.txt`).
    merge : bool, default=True
        Whether to merge compiled arguments into the existing registry.
        If False, the registry is replaced by the compiled arguments (ignoring Step 1).
    """
    
    # --- Step 1: Load Registry (Base Layer) ---
    registry_file_path = os.path.join(os.path.dirname(__file__), "viz_aux_params_registry.txt")
    if os.path.exists(registry_file_path):
        try:
            with open(registry_file_path, "r", encoding="utf-8") as rf:
                data = json.loads(rf.read())
                registry_source = data.get("registry", {})
                
                for key, entry in registry_source.items():
                    if key not in pm._registry:
                        pm._registry[key] = {"allowed": set(), "defaults": {}}
                    
                    target = pm._registry[key]
                    
                    if "allowed" in entry and entry["allowed"] is not None:
                        target["allowed"] = set(entry["allowed"])
                    
                    if "defaults" in entry:
                        target["defaults"].update(entry["defaults"])
                        
                    if "banned_keys" in entry:
                        if "banned_keys" not in target:
                            target["banned_keys"] = {}
                        target["banned_keys"].update(entry["banned_keys"])
        except Exception:
            pass

    # --- Step 2: Load Args (Update Layer) ---
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "viz_aux_params.txt")
        
    with open(path, "r", encoding="utf-8") as f:
        config_data = json.loads(f.read())
    raw_args = config_data.get("args", {})
    
    # Register args
    for arg_name, arg_entry in raw_args.items():
        raw_plots = arg_entry.get("plots", [])
        default_value = arg_entry.get("default", None)

        # Pass through raw plot specs to let register_arg handle parsing of
        # strings like "plot:mode", singletons, and wildcards.
        pm.register_arg(arg_name, plots=raw_plots, default=default_value, no_plots=arg_entry.get("no_plots"))
        
        # Descriptions
        pm._arg_doc[arg_name] = bool(arg_entry.get("add_to_docstring", False))
        if "desc" in arg_entry:
            pm.set_arg_desc(arg_name, arg_entry["desc"])
        
        # Context descriptions
        ctx_desc = arg_entry.get("ctx_desc", {})
        for ctx_key, desc_text in ctx_desc.items():
            if ":" in ctx_key:
                p_name, m_name = ctx_key.split(":", 1)
            else:
                p_name, m_name = ctx_key, None
            pm.set_arg_desc_for(arg_name, p_name, m_name, desc_text)
        
        # Context defaults
        if "ctx_defaults" in arg_entry:
            pm._arg_map[arg_name]["ctx_defaults"] = arg_entry["ctx_defaults"]

    # --- Step 3: Compile and Update Registry ---
    # This applies args to the registry, overwriting defaults where they overlap
    pm.compile_from_kwargs(merge=merge)
    
    # --- Step 4: Apply banned_keys from Args ---
    # Reconstruct known pairs from registry (now populated by both sources)
    known_plot_mode_pairs = set()
    for k in pm._registry.keys():
        if ":" in k:
            p, m = k.split(":", 1)
            known_plot_mode_pairs.add((p, m))
        else:
            known_plot_mode_pairs.add((k, None))

    for arg_name, arg_entry in raw_args.items():
        global_banned = arg_entry.get("banned_keys", [])
        ctx_banned_map = arg_entry.get("ctx_banned_keys", {})
        
        # Targets for this arg
        normalized_plots = pm._arg_map[arg_name]["plots"]
        target_contexts = []
        for p_name, m_name in normalized_plots:
            if m_name == "*":
                if p_name == "all":
                    target_contexts.extend(list(known_plot_mode_pairs))
                else:
                    target_contexts.extend([(k, m) for (k, m) in known_plot_mode_pairs if k == p_name])
            else:
                target_contexts.append((p_name, m_name))
        
        # Apply global banned
        if global_banned:
            for p_name, m_name in target_contexts:
                reg_key = pm._km(p_name, m_name)
                if reg_key not in pm._registry:
                     pm._registry[reg_key] = {"allowed": set(), "defaults": {}}
                if "banned_keys" not in pm._registry[reg_key]:
                    pm._registry[reg_key]["banned_keys"] = {}
                target_banned = pm._registry[reg_key]["banned_keys"]
                target_banned.setdefault(arg_name, [])
                target_banned[arg_name] = list(set(target_banned[arg_name]) | set(global_banned))
                
        # Apply context banned
        ctx_banned_items = []
        if isinstance(ctx_banned_map, dict):
            ctx_banned_items = list(ctx_banned_map.items())
        elif isinstance(ctx_banned_map, list):
            for item in ctx_banned_map:
                if isinstance(item, dict):
                    ctx_banned_items.extend(item.items())
                elif isinstance(item, (list, tuple)) and len(item) == 2:
                    ctx_banned_items.append((item[0], item[1]))

        for ctx_key, banned_list in ctx_banned_items:
            if ":" in ctx_key:
                p_name, m_name = ctx_key.split(":", 1)
            else:
                p_name, m_name = ctx_key, None
            reg_key = pm._km(p_name, m_name)
            if reg_key not in pm._registry:
                pm._registry[reg_key] = {"allowed": set(), "defaults": {}}
            if "banned_keys" not in pm._registry[reg_key]:
                pm._registry[reg_key]["banned_keys"] = {}
            target_banned = pm._registry[reg_key]["banned_keys"]
            target_banned.setdefault(arg_name, [])
            target_banned[arg_name] = list(set(target_banned[arg_name]) | set(banned_list))
            if "banned_keys" not in pm._registry[reg_key]:
                pm._registry[reg_key]["banned_keys"] = {}
            target_banned = pm._registry[reg_key]["banned_keys"]
            target_banned.setdefault(arg_name, [])
            target_banned[arg_name] = list(set(target_banned[arg_name]) | set(banned_list))
