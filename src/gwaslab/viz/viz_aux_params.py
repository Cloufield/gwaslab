import inspect
import os
import json


class VizParamsManager:
    """Parameter manager for visualization functions.

    Organizes and exposes relevant parameters for `viz_plot_*` functions.
    Supports per-plot and per-mode (e.g., "qq", "r", "b") registrations
    of allowed keys and default values, and merges them with object-level
    presets and call-time kwargs. Irrelevant parameters are filtered out
    before passing to the underlying plotting function.

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
        km = self._km(key, mode)
        reg = self._registry.get(km)
        if reg is None:
            return {}
        return dict(reg.get("defaults", {}))

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
        allowed = self.allowed(key, mode)
        if allowed is not None:
            name = key or func.__name__
            filtered = {k: v for k, v in params.items() if k in allowed}
            dropped = [k for k in params.keys() if k not in allowed]
            if log is not None and dropped:
                log.write(f"Filtered out args for `{name}`: {', '.join(dropped)}", verbose=verbose)
            reg = self._registry.get(self._km(key, mode), {})
            banned_map = reg.get("banned_keys", {})
            removed_subkeys = []
            for arg_name, banned in banned_map.items():
                if arg_name in filtered and isinstance(filtered[arg_name], dict):
                    for bk in banned:
                        if bk in filtered[arg_name]:
                            filtered[arg_name].pop(bk, None)
                            removed_subkeys.append(f"{arg_name}.{bk}")
            if log is not None and removed_subkeys:
                log.write(f"Filtered out kwargs sub-keys for `{name}`: {', '.join(removed_subkeys)}", verbose=verbose)
            return filtered
        sig = inspect.signature(func)
        # If function accepts **kwargs, pass through without filtering to allow downstream forwarding
        has_var_kw = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())
        if has_var_kw:
            return params
        names = set(sig.parameters.keys())
        name = key or func.__name__
        filtered = {k: v for k, v in params.items() if k in names}
        dropped = [k for k in params.keys() if k not in names]
        if log is not None and dropped:
            log.write(f"Filtered out args for `{name}`: {', '.join(dropped)}", verbose=verbose)
        return filtered

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

    def compile_from_args(self, merge=True):
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

    def na_args(self):
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
    translated into per-plot/mode allowed lists via `compile_from_args()`.
    """
    return

def load_viz_config(pm, path=None, merge=True):
    """Load visualization parameter settings from a config file."""
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "viz_aux_params.txt")
    with open(path, "r", encoding="utf-8") as f:
        data = json.loads(f.read())
    args = data.get("args", {})
    for name, entry in args.items():
        plots = entry.get("plots", [])
        default = entry.get("default", None)
        pm.register_arg(name, plots=[(p[0], p[1]) for p in plots], default=default)
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
    pm.compile_from_args(merge=merge)
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
        # global banned keys apply to all listed plots/modes for this arg
        if g_banned:
            for p in plots:
                key, mode = p[0], p[1]
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
    reg = data.get("registry", {})
    for km, obj in reg.items():
        allowed = obj.get("allowed")
        defaults = obj.get("defaults", {})
        banned_map = obj.get("banned_keys", {})
        if ":" in km:
            key, mode = km.split(":", 1)
        else:
            key, mode = km, None
        if merge and pm._km(key, mode) in pm._registry:
            existing = pm._registry[pm._km(key, mode)]
            if allowed is not None:
                if existing["allowed"] is None:
                    existing["allowed"] = set()
                existing["allowed"].update(set(allowed))
            existing["defaults"].update(dict(defaults))
            if banned_map:
                existing["banned_keys"] = dict(banned_map)
        else:
            pm.register(key, allowed=set(allowed) if allowed is not None else None, defaults=defaults, mode=mode)
            # attach banned_keys to registry entry
            pm._registry[pm._km(key, mode)]["banned_keys"] = dict(banned_map)
"""
模块描述
========

本模块提供可复用的可视化参数管理器，用于统一管理 `viz_plot_*` 中大量重叠/跨模式的参数，
并在调用具体绘图函数时按绘图和模式只暴露相关参数。核心目标是：

- 将参数定义和暴露范围集中化，避免调用处传入无效或冗余的 kwargs；
- 支持两种工作流：
  1) 基于绘图/模式的注册（白名单 + 默认值），适合直接按绘图组织；
  2) 基于参数条目的整理（为每个参数附加“潜在绘图/模式”），适合先收敛参数再自动生成每个绘图的参数清单；
- 合并顺序明确：注册默认值 → 对象级预设 → 调用时 kwargs；
- 过滤策略明确：优先按白名单过滤；若未注册白名单，回退到函数签名过滤。

关键概念
--------
- `key`：绘图标识，通常与 `Sumstats` 的包装器同名，如 `plot_region`、`plot_qq` 等。
- `mode`：同一绘图下的子模式，如 `r`（区域）、`qq`、`b`（密度/二分类）、`q`（定量）等。
- `_registry`：按 `key:mode` 记录“允许参数集合（allowed）”与“默认值（defaults）”。
- `_store`：对象级预设，按 `key:mode` 存储，覆盖 `_registry` 的默认值。
- `_arg_map`：参数条目字典，记录参数名 → 适用的 (key, mode) 集合与默认值，用于从“参数为中心”自动编译到 `_registry`。

两种工作流
----------
1) 绘图/模式注册工作流：
   - 使用 `register(key, allowed, defaults, mode)` 直接声明白名单与默认值；
   - 调用时通过 `merge(...) + filter(...)` 合并与裁剪参数；
   - 适合已明确各绘图的参数面时快速集成。

2) 参数条目整理工作流：
   - 使用 `register_arg(name, plots=[(key, mode), ...], default=...)` 为参数建立“潜在绘图/模式”映射；
   - 通过 `compile_from_args()` 自动生成每个绘图/模式的 allowed 与 defaults；
   - 使用 `na_args()` 列出未关联任何绘图的参数，手动修复后再次编译；
   - 适合先在一个地方统一列出参数，再让系统构建每个绘图的白名单。

集成与示例
----------
- 在 `Sumstats.__init__` 中创建并初始化：
  - `self.viz_params = VizParamsManager()`；
  - 调用 `init_viz_registry(self.viz_params)`（绘图/模式注册流），或调用
    `init_viz_arg_entries(self.viz_params); self.viz_params.compile_from_args()`（参数条目流）。
- 在 `Sumstats.plot_*` 包装器中：
  - `params = self.viz_params.merge('plot_region', kwargs, mode='r')`；
  - `params = self.viz_params.filter(_plot_fn, params, key='plot_region', mode='r')`；
  - 然后将 `params` 传给底层绘图函数；必要时补齐 `build` 保持兼容。

查询与暴露
----------
- 使用 `public(key, mode)` 获取“对外暴露”的参数（已合并并按白名单裁剪），便于生成 UI/文档。
- 使用 `args_for_plot(key, mode)` 查看某绘图/模式当前挂接的参数名集合。

注意事项
--------
- 若某绘图/模式未注册白名单，系统会回退到函数签名过滤，这样仍能安全调用但建议补充注册以保证稳定性；
- 默认值设计为“可为空”，鼓励在对象级设定偏好（`set/update`），以减少调用处的冗长参数；
- `compile_from_args(merge=True)` 会并入现有注册表，若希望覆盖，请把 `merge=False`（可在需要时扩展）。
"""
