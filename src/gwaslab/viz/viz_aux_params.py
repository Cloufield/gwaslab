"""Visualization parameter management for gwaslab plots.

Provides a centralized registry of per-plot and per-mode parameter surfaces:
- Allowed keys and defaults are stored under `key[:mode]`.
- Object-level presets overlay registered defaults; call-time kwargs override both.
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
2) Object-level presets (`set`/`update`)
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
        self._global_defaults = {}

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
        km = self._km(key, mode)
        reg = self._registry.get(km)
        base = {}
        if reg is not None:
            base = dict(reg.get("defaults", {}))
        for arg, val in getattr(self, "_global_defaults", {}).items():
            if arg not in base:
                base[arg] = val
        return base

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
        km = self._km(key, mode)
        reg = self._registry.get(km)
        if reg is None:
            return None
        return reg.get("allowed")

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
        # 1) Look up whitelist for this plot/mode; None means signature-based filtering
        allowed_keys = self.allowed(key, mode)
        if allowed_keys is not None:
            target_name = key or func.__name__
            # 2) Define predicate: direct allow or numeric-suffix allow (e.g., "arg1", "arg2")
            def is_allowed_key(arg_name):
                if arg_name in allowed_keys:
                    return True
                # numeric suffix handling
                if arg_name[-1:].isdigit():
                    base = arg_name[:-1]
                    return base in allowed_keys
                return False
            # 3) Apply whitelist to params; record dropped keys for logging
            filtered_params = {k: v for k, v in params.items() if is_allowed_key(k)}
            dropped_keys = [k for k in params.keys() if not is_allowed_key(k)]
            if log is not None and dropped_keys:
                log.write(f"Filtered out args for `{target_name}`: {', '.join(dropped_keys)}", verbose=verbose)
            # 4) Remove banned nested sub-keys inside dict-type args
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
        # 5) No whitelist: inspect function signature to determine allowed params
        func_signature = inspect.signature(func)
        # If function accepts **kwargs, pass through unchanged to allow downstream forwarding
        has_var_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in func_signature.parameters.values())
        if has_var_kwargs:
            return params
        # 6) Apply signature-based filtering; record dropped keys for logging
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


    def register_arg(self, name, plots=None, default=None):
        """Register a single argument and its plot/mode associations.

        Parameters
        ----------
        name : str
            Argument name to register.
        plots : list or set of (key, mode) or str, optional
            Plot/mode pairs this argument applies to. If an element is a str,
            it is treated as `(str, None)`.
        default : Any, optional
            Default value to apply for plots that include this argument.
        """
        if not hasattr(self, "_arg_map"):
            self._arg_map = {}
        if not hasattr(self, "_arg_desc"):
            self._arg_desc = {}
        entry = self._arg_map.get(name, {"plots": set(), "default": None})
        if plots:
            for p in plots:
                if isinstance(p, tuple):
                    entry["plots"].add((p[0], p[1]))
                else:
                    entry["plots"].add((p, None))
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
        compiled = {}
        for arg, entry in self._arg_map.items():
            for key, mode in entry["plots"]:
                km = self._km(key, mode)
                if km not in compiled:
                    compiled[km] = {"allowed": set(), "defaults": {}}
                compiled[km]["allowed"].add(arg)
                if entry["default"] is not None:
                    compiled[km]["defaults"][arg] = entry["default"]
        for km, obj in compiled.items():
            allowed = obj["allowed"]
            defaults = obj["defaults"]
            if merge and km in self._registry:
                reg = self._registry[km]
                if reg["allowed"] is None:
                    reg["allowed"] = set()
                reg["allowed"].update(allowed)
                reg["defaults"].update(defaults)
            else:
                self._registry[km] = {"allowed": set(allowed), "defaults": dict(defaults)}

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

    def build_docstring_for_plot(self, key, mode=None):
        """Construct docstring section from registered args for a plot/mode.

        Parameters
        ----------
        key : str
            Plot identifier.
        mode : str or None
            Optional sub-mode.

        Returns
        -------
        str
            A docstring-like parameter block with names and descriptions.
        """
        names = sorted([n for n in self.args_for_plot(key, mode) if self._arg_doc.get(n, True)])
        lines = []
        for n in names:
            d = self.get_arg_desc_for(n, key, mode)
            if d:
                lines.append(f"{n} : {d}")
            else:
                lines.append(f"{n} : ")
        return "\n".join(lines)

def init_viz_registry(pm):
    """Register allowed keys and defaults for major visualization functions.

    Parameters
    ----------
    pm : VizParamsManager
        The parameter manager to populate.

    Notes
    -----
    - Keys follow `Sumstats` wrapper names (e.g., 'plot_region', 'plot_qq').
    - Modes separate distinct parameter surfaces (e.g., 'r', 'qq', 'b').
    - Defaults dictionaries are empty by design; users set presets per-object.
    """
    return

def init_viz_arg_entries(pm):
    """Seed argument-to-plot mappings to bootstrap registry construction.

    This allows a workflow where arguments are curated in a single place, then
    translated into per-plot/mode allowed lists via `compile_from_kwargs()`.
    """
    return

def load_viz_config(pm, path=None, merge=True):
    """Load visualization parameter settings from config files.

    Primary args map is loaded from `viz_aux_params.txt`.
    Plot registry (allowed/defaults/banned_keys) is loaded from
    `viz_aux_params_registry.txt` if present, otherwise falls back
    to the `registry` embedded in `viz_aux_params.txt`.
    """
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "viz_aux_params.txt")
    with open(path, "r", encoding="utf-8") as f:
        data = json.loads(f.read())
    # Discover known plot:mode pairs from external registry or embedded registry
    reg_path = os.path.join(os.path.dirname(__file__), "viz_aux_params_registry.txt")
    known_plot_pairs = []
    try:
        if os.path.exists(reg_path):
            with open(reg_path, "r", encoding="utf-8") as rf:
                reg_container = json.loads(rf.read())
                reg_src = reg_container.get("registry", {})
        else:
            reg_src = data.get("registry", {})
    except Exception:
        reg_src = data.get("registry", {})
    for km in reg_src.keys():
        if ":" in km:
            key, mode = km.split(":", 1)
        else:
            key, mode = km, None
        known_plot_pairs.append((key, mode))
    # Fallback: derive known plots from args if registry is empty
    if not known_plot_pairs:
        args_probe = data.get("args", {})
        seen = set()
        for _, entry in args_probe.items():
            for p in entry.get("plots", []):
                if isinstance(p, str) and p != "all":
                    seen.add((p, None))
                elif isinstance(p, (list, tuple)):
                    if len(p) == 1 and p[0] != "all":
                        seen.add((p[0], None))
                    elif len(p) >= 2:
                        seen.add((p[0], p[1]))
        known_plot_pairs = list(seen)

    args = data.get("args", {})
    for name, entry in args.items():
        plots = entry.get("plots", [])
        default = entry.get("default", None)
        expanded_plots = []
        for p in plots:
            if isinstance(p, str):
                if p == "all":
                    expanded_plots.extend(known_plot_pairs)
                    if default is not None:
                        pm._global_defaults[name] = default
                else:
                    expanded_plots.extend([(k, m) for (k, m) in known_plot_pairs if k == p])
            elif isinstance(p, (list, tuple)):
                if len(p) == 1:
                    if p[0] == "all":
                        expanded_plots.extend(known_plot_pairs)
                        if default is not None:
                            pm._global_defaults[name] = default
                    else:
                        expanded_plots.extend([(k, m) for (k, m) in known_plot_pairs if k == p[0]])
                elif len(p) >= 2:
                    key, mode = p[0], p[1]
                    if mode is None:
                        expanded_plots.extend([(k, m) for (k, m) in known_plot_pairs if k == key])
                    else:
                        expanded_plots.append((key, mode))
            # ignore other cases
        pm.register_arg(name, plots=expanded_plots, default=default)
        pm._arg_doc[name] = bool(entry.get("add_to_docstring", False))
        desc = entry.get("desc")
        if desc is not None:
            pm.set_arg_desc(name, desc)
        ctx = entry.get("ctx_desc", {})
        for km, d in ctx.items():
            if ":" in km:
                key, mode = km.split(":", 1)
            else:
                key, mode = km, None
            pm.set_arg_desc_for(name, key, mode, d)
    # compile allowed/defaults from args
    pm.compile_from_kwargs(merge=merge)
    # apply per-plot/mode defaults for args (override compiled defaults)
    for name, entry in args.items():
        ctx_defaults = entry.get("ctx_defaults", {})
        for km, val in ctx_defaults.items():
            if ":" in km:
                key, mode = km.split(":", 1)
            else:
                key, mode = km, None
            kmm = pm._km(key, mode)
            if kmm not in pm._registry:
                pm._registry[kmm] = {"allowed": set(), "defaults": {}}
            pm._registry[kmm]["defaults"][name] = val
    # attach banned_keys from args (global and context-specific) into registry
    # global banned_keys apply to all plots/modes listed in the arg's plots
    # ctx_banned_keys uses "plot:mode" keys to target specific contexts
    for name, entry in args.items():
        g_banned = entry.get("banned_keys", [])
        ctx_banned = entry.get("ctx_banned_keys", {})
        plots = entry.get("plots", [])
        expanded_plots = []
        for p in plots:
            if isinstance(p, str):
                if p == "all":
                    expanded_plots.extend(known_plot_pairs)
                else:
                    expanded_plots.extend([(k, m) for (k, m) in known_plot_pairs if k == p])
            elif isinstance(p, (list, tuple)):
                if len(p) == 1:
                    if p[0] == "all":
                        expanded_plots.extend(known_plot_pairs)
                    else:
                        expanded_plots.extend([(k, m) for (k, m) in known_plot_pairs if k == p[0]])
                elif len(p) >= 2:
                    key, mode = p[0], p[1]
                    if mode is None:
                        expanded_plots.extend([(k, m) for (k, m) in known_plot_pairs if k == key])
                    else:
                        expanded_plots.append((key, mode))
        # global banned keys apply to all listed plots/modes for this arg
        if g_banned:
            for key, mode in expanded_plots:
                kmm = pm._km(key, mode)
                if kmm not in pm._registry:
                    pm._registry[kmm] = {"allowed": set(), "defaults": {}}
                if "banned_keys" not in pm._registry[kmm]:
                    pm._registry[kmm]["banned_keys"] = {}
                pm._registry[kmm]["banned_keys"].setdefault(name, [])
                pm._registry[kmm]["banned_keys"][name] = list(set(pm._registry[kmm]["banned_keys"][name]) | set(g_banned))
        # context-specific banned keys
        if isinstance(ctx_banned, dict):
            iterable = ctx_banned.items()
        elif isinstance(ctx_banned, list):
            tmp = []
            for item in ctx_banned:
                if isinstance(item, dict):
                    tmp.extend(item.items())
                elif isinstance(item, (list, tuple)) and len(item) == 2:
                    tmp.append((item[0], item[1]))
            iterable = tmp
        else:
            iterable = []
        for km, val in iterable:
            if ":" in km:
                key, mode = km.split(":", 1)
            else:
                key, mode = km, None
            kmm = pm._km(key, mode)
            if kmm not in pm._registry:
                pm._registry[kmm] = {"allowed": set(), "defaults": {}}
            if "banned_keys" not in pm._registry[kmm]:
                pm._registry[kmm]["banned_keys"] = {}
            pm._registry[kmm]["banned_keys"].setdefault(name, [])
            pm._registry[kmm]["banned_keys"][name] = list(set(pm._registry[kmm]["banned_keys"][name]) | set(val))
    # Prefer external registry file when available
    reg_path = os.path.join(os.path.dirname(__file__), "viz_aux_params_registry.txt")
    reg = {}
    try:
        if os.path.exists(reg_path):
            with open(reg_path, "r", encoding="utf-8") as rf:
                reg_container = json.loads(rf.read())
                reg = reg_container.get("registry", {})
        else:
            reg = data.get("registry", {})
    except Exception:
        reg = data.get("registry", {})
    for km, obj in reg.items():
        allowed = obj.get("allowed")
        defaults = obj.get("defaults", {})
        banned_map = obj.get("banned_keys", {})
        if ":" in km:
            key, mode = km.split(":", 1)
        else:
            key, mode = km, None
        kmm = pm._km(key, mode)
        if merge and kmm in pm._registry:
            existing = pm._registry[kmm]
            # Replace allowed set with external registry definition to avoid unintended unions
            if allowed is not None:
                existing["allowed"] = set(allowed)
            # Merge defaults from external registry
            existing["defaults"].update(dict(defaults))
            # Replace banned_keys with external registry mapping
            if banned_map:
                existing["banned_keys"] = dict(banned_map)
        else:
            pm.register(key, allowed=set(allowed) if allowed is not None else None, defaults=defaults, mode=mode)
            # attach banned_keys to registry entry
            pm._registry[kmm]["banned_keys"] = dict(banned_map)

def _append_args_catalog_to_docstring():
    try:
        path = os.path.join(os.path.dirname(__file__), "viz_aux_params.txt")
        with open(path, "r", encoding="utf-8") as f:
            data = json.loads(f.read())
        args = data.get("args", {})
        lines = []
        lines.append("\nArgs (Full Catalog)")
        for name, entry in args.items():
            desc = entry.get("desc", "")
            lines.append(f"{name} : {desc}")
        VizParamsManager.__doc__ = (VizParamsManager.__doc__ or "") + "\n" + "\n".join(lines)
    except Exception:
        pass

_append_args_catalog_to_docstring()

def _sort_viz_params_file(path=None):
    try:
        if path is None:
            path = os.path.join(os.path.dirname(__file__), "viz_aux_params.txt")
        with open(path, "r", encoding="utf-8") as f:
            data = json.loads(f.read())
        args = data.get("args", {})
        registry = data.get("registry", {})
        sorted_args = {k: args[k] for k in sorted(args.keys())}
        sorted_registry = {k: registry[k] for k in sorted(registry.keys())}
        new_data = {"args": sorted_args, "registry": sorted_registry}

        def _dump_inline(obj):
            if isinstance(obj, dict):
                return json.dumps(obj, ensure_ascii=False, separators=(",", ": "))
            if isinstance(obj, list):
                return "[" + ", ".join(_dump_inline(x) for x in obj) + "]"
            return json.dumps(obj, ensure_ascii=False)

        def _dump_dict(obj, indent=2, level=0):
            sp = " " * (indent * level)
            lines = []
            for k in obj:
                val = obj[k]
                if isinstance(val, dict):
                    if len(val) == 0:
                        v = "{}"
                    else:
                        v = _dump_dict(val, indent, level + 1)
                elif isinstance(val, list):
                    v = _dump_inline(val)
                else:
                    v = json.dumps(val, ensure_ascii=False)
                lines.append((k, v))
            body = ",\n".join("{}\"{}\": {}".format(" " * (indent * (level + 1)), k, v) for k, v in lines)
            return "{\n" + body + "\n" + sp + "}"

        text = _dump_dict(new_data, indent=2, level=0)
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
    except Exception:
        pass

_sort_viz_params_file()
