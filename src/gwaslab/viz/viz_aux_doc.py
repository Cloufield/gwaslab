"""Build Sumstats plot docstrings from VizParamsManager registry.
"""

from __future__ import annotations

import re
from typing import Any, Callable, FrozenSet, Optional, Sequence

from gwaslab.info.g_numpy_doc import format_parameter, format_section, infer_type_annotation
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config

_DOC_PM = VizParamsManager()
load_viz_config(_DOC_PM)

_SUMSTATS_HIDE: FrozenSet[str] = frozenset({"insumstats", "log", "object"})
_DOC_ALIAS_HIDE: FrozenSet[str] = frozenset({
    "fig_args",
    "figargs",
    "scatter_args",
    "scatterargs",
    "save_args",
    "saveargs",
    "qq_scatter_args",
    "anno_args",
    "anno_args_single",
    "region_anno_bbox_args",
    "region_title_args",
    "highlight_anno_args",
})
_ARG_DESC_ALIASES = {
    "fig_args": "fig_kwargs",
    "figargs": "fig_kwargs",
    "scatter_args": "scatter_kwargs",
    "scatterargs": "scatter_kwargs",
    "save_args": "save_kwargs",
    "saveargs": "save_kwargs",
    "qq_scatter_args": "qq_scatter_kwargs",
    "anno_args": "anno_kwargs",
    "anno_args_single": "anno_kwargs_single",
    "region_anno_bbox_args": "region_anno_bbox_kwargs",
    "region_title_args": "region_title_kwargs",
    "highlight_anno_args": "highlight_anno_kwargs",
}
_SECTION_HEADING = re.compile(
    r"\n\s*Parameters\n\s*-+\n|\nParameters\n-+\n",
    re.MULTILINE,
)
_RETURNS_HEADING = re.compile(
    r"\n\s*Returns\n\s*-+\n|\nReturns\n-+\n",
    re.MULTILINE,
)
_NOTES_HEADING = re.compile(
    r"\n\s*Notes\n\s*-+\n|\nNotes\n-+\n",
    re.MULTILINE,
)
_NEXT_SECTION_AFTER_RETURNS = re.compile(
    r"\n\s*(Examples|Raises|Warns|References|Attributes|Notes|See Also)\n\s*-+\n"
    r"|\n(Examples|Raises|Warns|References|Attributes|Notes|See Also)\n-+\n",
    re.MULTILINE | re.IGNORECASE,
)
_PRIORITY_PARAMS: Sequence[str] = (
    "mode",
    "region",
    "build",
    "sig_level",
    "anno_sig_level",
    "skip",
    "cut",
    "marker_size",
    "vcf_path",
    "ld_block",
    "ld_link",
)

_REGISTRY_TYPE_OVERRIDES: dict[str, str] = {
    "mode": "str, optional",
    "region": "tuple or str, optional",
    "build": "str, optional",
    "vcf_path": "str, optional",
    "fig_kwargs": "dict, optional",
    "scatter_kwargs": "dict, optional",
    "save_kwargs": "dict, optional",
    "qq_scatter_kwargs": "dict, optional",
    "anno_kwargs": "dict, optional",
    "title": "str, optional",
    "anno": "DataFrame or list, optional",
    "sig_level": "float, optional",
    "anno_sig_level": "float, optional",
    "marker_size": "list, optional",
    "windowsizekb": "int, optional",
    "skip": "int, optional",
    "cut": "int, optional",
    "verbose": "bool, optional",
    "save": "bool, optional",
    "ld_block": "bool, optional",
    "ld_link": "bool, optional",
}


def _registry_type_annotation(name: str, default_val: Any, entry: dict) -> Optional[str]:
    if entry.get("type"):
        return str(entry["type"])
    if name in _REGISTRY_TYPE_OVERRIDES:
        override = _REGISTRY_TYPE_OVERRIDES[name]
        if default_val is not None and "optional" in override:
            return infer_type_annotation(default_val, has_default=True)
        return override
    if name.endswith("_kwargs") or name.endswith("_args"):
        return infer_type_annotation(default_val, has_default=default_val is not None)
    if name.endswith("_path") or name.endswith("_file") or name.endswith("_dir"):
        return "str, optional" if default_val is None else infer_type_annotation(default_val)
    if name.endswith("_level") or name.endswith("_thr"):
        return infer_type_annotation(default_val, has_default=default_val is not None)
    if name.endswith("_color") or name.endswith("_colors"):
        return "str or list, optional"
    if name.endswith("_size") or name.endswith("_fontsize"):
        return infer_type_annotation(default_val, has_default=default_val is not None)
    return None


def _resolve_param_desc(pm: VizParamsManager, name: str, key: str, mode: Optional[str]) -> str:
    desc = pm.get_arg_desc_for(name, key, mode).strip()
    if desc:
        return desc
    canonical = _ARG_DESC_ALIASES.get(name)
    if canonical:
        base = pm.get_arg_desc_for(canonical, key, mode).strip()
        if base:
            return f"Deprecated alias for ``{canonical}``. {base}"
    return name


def get_doc_hidden_params() -> FrozenSet[str]:
    """Return parameter names omitted from generated Sumstats plot docs.
"""
    return _SUMSTATS_HIDE | _DOC_ALIAS_HIDE


def get_doc_params_manager() -> VizParamsManager:
    """Return the module-level VizParamsManager used for doc generation.
"""
    return _DOC_PM


def _sort_param_names(names: Sequence[str]) -> list[str]:
    priority = {name: index for index, name in enumerate(_PRIORITY_PARAMS)}

    def sort_key(name: str) -> tuple[int, str]:
        return (priority.get(name, len(_PRIORITY_PARAMS)), name)

    return sorted(names, key=sort_key)


def _extract_impl_sections(doc: str) -> tuple[str, str, str]:
    if not doc:
        return "", "", ""

    param_match = _SECTION_HEADING.search(doc)
    returns_match = _RETURNS_HEADING.search(doc)
    notes_match = _NOTES_HEADING.search(doc)

    summary_end = param_match.start() if param_match else len(doc)
    summary = doc[:summary_end].strip()

    returns = ""
    if returns_match:
        search_from = returns_match.end()
        next_section = _NEXT_SECTION_AFTER_RETURNS.search(doc, search_from)
        returns_end = next_section.start() if next_section else len(doc)
        returns = doc[search_from:returns_end].strip()

    notes = ""
    if notes_match:
        notes = doc[notes_match.end():].strip()

    return summary, returns, notes


def _banned_subkeys(pm: VizParamsManager, key: str, mode: Optional[str], arg_name: str) -> list[str]:
    entry = pm._get_registry_entry(key, mode)
    banned = entry.get("banned_keys", {}).get(arg_name, [])
    return list(banned) if banned else []


def _sanitize_viz_returns(returns: str) -> str:
    """Normalize impl Returns for public plot docs (drop internal log, fix layout).
"""
    if not returns:
        return ""
    out: list[str] = []
    for line in returns.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if re.match(r"^log\s*:", stripped, re.IGNORECASE):
            continue
        name_type = re.match(r"^(\w+)\s*:\s*(.+)$", stripped)
        if name_type:
            name, type_part = name_type.group(1), name_type.group(2).strip()
            if name.lower() == "log":
                continue
            out.append(type_part)
            continue
        if out and not stripped.startswith("Optional") and "optional" not in stripped.lower():
            out.append(f"    {stripped}")
        else:
            out.append(stripped)
    return "\n".join(out).strip()


def build_viz_docstring(
    pm: VizParamsManager,
    key: str,
    mode: Optional[str],
    impl_func: Callable[..., Any],
    *,
    wrapper_doc: str = "",
    hide_params: Optional[FrozenSet[str]] = None,
    include_internal: bool = False,
    positional_params: Optional[Sequence[tuple[str, str, str]]] = None,
) -> str:
    """Build a numpy-style docstring from registry allowed keys and catalog descriptions.

Parameters
----------
positional_params : sequence of (name, type, description), optional
    Required positional arguments prepended before registry kwargs.
"""
    hide = hide_params if hide_params is not None else get_doc_hidden_params()
    allowed = pm.allowed(key, mode)
    if not allowed:
        allowed = set()

    param_names = [
        name
        for name in _sort_param_names(sorted(allowed))
        if (include_internal or not name.startswith("_"))
        and name not in hide
    ]

    defaults = pm.defaults(key, mode)
    lines: list[str] = []

    wrapper = (wrapper_doc or "").strip()
    summary, returns, notes = _extract_impl_sections(impl_func.__doc__ or "")

    if wrapper:
        lines.append(wrapper)
    elif summary:
        lines.append(summary)

    if param_names or positional_params:
        lines.append("")
        lines.append("Parameters")
        lines.append("----------")
        if positional_params:
            for name, type_ann, desc in positional_params:
                lines.append(f"{name} : {type_ann}")
                lines.append(f"    {desc.strip()}")
        for name in param_names:
            desc = _resolve_param_desc(pm, name, key, mode)
            default_val = defaults.get(name, pm._arg_map.get(name, {}).get("default"))
            entry = pm._arg_map.get(name, {})
            type_ann = _registry_type_annotation(name, default_val, entry)

            banned = _banned_subkeys(pm, key, mode, name)
            if banned:
                desc = f"{desc} Do not set in ``{name}``: {', '.join(banned)}."

            lines.extend(format_parameter(name, desc, default_val, type_annotation=type_ann))

    if returns:
        sanitized = _sanitize_viz_returns(returns)
        if sanitized:
            lines.extend(format_section("Returns", sanitized.splitlines()))

    if notes:
        lines.extend(format_section("Notes", notes.splitlines()))

    return "\n".join(lines).strip() + "\n"


def build_gl_viz_docstring(
    key: str,
    impl_func: Callable[..., Any],
    mode: Optional[str] = None,
    *,
    summary: str = "",
    positional_params: Optional[Sequence[tuple[str, str, str]]] = None,
) -> str:
    """Build a registry-aligned docstring for gl.* plot wrappers.
"""
    return build_viz_docstring(
        _DOC_PM,
        key,
        mode,
        impl_func,
        wrapper_doc=summary,
        positional_params=positional_params,
    )
