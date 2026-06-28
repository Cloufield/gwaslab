"""Shared NumPy (numpydoc) docstring helpers for GWASLab.
"""

from __future__ import annotations

import inspect
import re
from typing import Any, Callable, FrozenSet, Optional, Sequence

try:
    from docstring_parser import DocstringParam, parse as parse_docstring
except ImportError:  # pragma: no cover
    DocstringParam = Any  # type: ignore[misc, assignment]
    parse_docstring = None  # type: ignore[assignment]

_NUMPY_SECTIONS = frozenset({
    "Parameters",
    "Returns",
    "Yields",
    "Raises",
    "Warns",
    "Notes",
    "References",
    "Examples",
    "Attributes",
})

_SECTION_HEADING = re.compile(
    r"^(\s*)(Parameters|Returns|Yields|Raises|Warns|Notes|References|Examples|Attributes|Args):?\s*$",
    re.MULTILINE,
)
_GOOGLE_ARGS = re.compile(r"^(\s*)Args:\s*$", re.MULTILINE)
_GOOGLE_RETURNS = re.compile(r"^(\s*)Returns:\s*$", re.MULTILINE)
_GOOGLE_RAISES = re.compile(r"^(\s*)Raises:\s*$", re.MULTILINE)
_GOOGLE_YIELDS = re.compile(r"^(\s*)Yields:\s*$", re.MULTILINE)
_SECTION_UNDERLINE = re.compile(
    r"^(\s*)(Parameters|Returns|Yields|Raises|Warns|Notes|References|Examples|Attributes)\n\s*-+\s*$",
    re.MULTILINE,
)
_PARAM_LINE = re.compile(r"^(\s*)(\w+)\s*:\s*(.+)$", re.MULTILINE)
_GOOGLE_PARAM_LINE = re.compile(r"^(\s*)(\w+)\s*(?:\(([^)]+)\))?\s*:\s*(.+)$", re.MULTILINE)
_PROSE_DEFAULT = re.compile(r"\.\s*Default is [^.]+(?:\.|$)\s*", re.IGNORECASE)
_METH_NOTE = re.compile(r"When called via :meth:", re.IGNORECASE)
_DEFAULT_DROP_PARAMS = frozenset({
    "sumstats_obj",
    "log",
    "insumstats",
    "insumstats_or_dataframe",
    "meta",
})


def _section_underline(title: str) -> str:
    return "-" * len(title)


def infer_type_annotation(default_val: Any, *, has_default: bool = True) -> str:
    """Infer a numpydoc type annotation from a default value.
"""
    if has_default and default_val is not None:
        base: str
        if isinstance(default_val, bool):
            base = "bool"
        elif isinstance(default_val, int) and not isinstance(default_val, bool):
            base = "int"
        elif isinstance(default_val, float):
            base = "float"
        elif isinstance(default_val, str):
            base = "str"
        elif isinstance(default_val, dict):
            base = "dict"
        elif isinstance(default_val, (list, tuple)):
            base = "list"
        else:
            base = "Any"
        return f"{base}, default {repr(default_val)}"
    if default_val is None and has_default:
        return "Any, optional"
    return "Any, optional"


def format_parameter(
    name: str,
    desc: str,
    default: Any = ...,
    *,
    type_annotation: Optional[str] = None,
) -> list[str]:
    """Return two lines: ``name : type`` and an indented description.
"""
    if type_annotation is None:
        has_default = default is not ...
        default_val = None if default is ... else default
        type_annotation = infer_type_annotation(default_val, has_default=has_default)

    text = desc.strip()
    if text and not text.endswith("."):
        text += "."

    return [f"{name} : {type_annotation}", f"    {text}"] if text else [f"{name} : {type_annotation}"]


def _fix_section_header_gaps(text: str) -> str:
    """Remove blank lines between a section title and its underline.
"""
    pattern = re.compile(
        r"^([ \t]*)(Parameters|Returns|Yields|Raises|Warns|Notes|References|Examples|Attributes)[ \t]*\n(?:[ \t]*\n)+([ \t]*-+)[ \t]*$",
        re.MULTILINE,
    )

    def repl(match: re.Match[str]) -> str:
        indent, title = match.group(1), match.group(2)
        return f"{indent}{title}\n{_section_underline(title)}"

    return pattern.sub(repl, text)


def _format_parameter_body(body: str) -> list[str]:
    """Normalize Parameters section lines to column-0 names and 4-space descriptions.
"""
    lines: list[str] = []
    for raw in body.splitlines():
        stripped = raw.strip()
        if not stripped:
            continue
        param_match = re.match(r"(\w+)\s*:\s*(.+)", stripped)
        if param_match:
            lines.append(f"{param_match.group(1)} : {param_match.group(2)}")
        else:
            lines.append(f"    {stripped}")
    return lines


def _normalize_parameter_sections(text: str) -> str:
    pattern = re.compile(
        r"^([ \t]*Parameters\n[ \t]*-+\n)(.*?)(?=^[ \t]*(Returns|Yields|Raises|Warns|Notes|References|Examples|Attributes)\n[ \t]*-+\n|\Z)",
        re.MULTILINE | re.DOTALL,
    )

    def repl(match: re.Match[str]) -> str:
        header = match.group(1)
        body = match.group(2)
        formatted = _format_parameter_body(body)
        if not formatted:
            return header.rstrip() + "\n"
        return header + "\n".join(formatted) + "\n"

    return pattern.sub(repl, text)


def format_section(title: str, body_lines: Sequence[str]) -> list[str]:
    """Format a numpydoc section with correct underline length.
"""
    if title not in _NUMPY_SECTIONS and title != "Args":
        title = title.rstrip(":")

    lines: list[str] = ["", title, _section_underline(title)]
    if title == "Parameters":
        lines.extend(_format_parameter_body("\n".join(body_lines)))
        return lines

    for line in body_lines:
        stripped = line.rstrip()
        if not stripped:
            continue
        if stripped.startswith(" "):
            lines.append(stripped)
        else:
            lines.append(f"    {stripped}")
    return lines


def _normalize_google_param_lines(doc: str) -> str:
    """Convert Google-style ``name (type): desc`` lines to NumPy parameter lines.
"""

    def repl(match: re.Match[str]) -> str:
        indent, name, type_hint, desc = match.group(1), match.group(2), match.group(3), match.group(4)
        type_part = type_hint.strip() if type_hint else "Any"
        return f"{indent}{name} : {type_part}\n{indent}    {desc.strip()}"

    lines = doc.splitlines()
    out: list[str] = []
    in_google_args = False
    base_indent = ""

    for line in lines:
        if _GOOGLE_ARGS.match(line):
            in_google_args = True
            base_indent = _GOOGLE_ARGS.match(line).group(1)  # type: ignore[union-attr]
            out.append(f"{base_indent}Parameters")
            out.append(f"{base_indent}{_section_underline('Parameters')}")
            continue
        if in_google_args and re.match(rf"^{re.escape(base_indent)}(Returns|Raises|Yields|Notes|Examples):?\s*$", line):
            in_google_args = False
            out.append(line)
            continue
        if in_google_args:
            gmatch = _GOOGLE_PARAM_LINE.match(line)
            if gmatch:
                indent, name, type_hint, desc = gmatch.group(1), gmatch.group(2), gmatch.group(3), gmatch.group(4)
                type_part = type_hint.strip() if type_hint else "Any"
                out.append(f"{indent}{name} : {type_part}")
                out.append(f"{indent}    {desc.strip()}")
                continue
        out.append(line)

    return "\n".join(out)


def _dedent_section_headers(text: str) -> str:
    """Place NumPy section titles at column 0 for reliable parsing.
"""
    return re.sub(
        r"^[ \t]+(Parameters|Returns|Yields|Raises|Warns|Notes|References|Examples|Attributes)[ \t]*$",
        r"\1",
        text,
        flags=re.MULTILINE,
    )


def _fix_returns_indent(text: str) -> str:
    """Dedent over-indented return type lines in Returns/Yields sections.
"""
    lines = text.split("\n")
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if (
            line in ("Returns", "Yields")
            and i + 1 < len(lines)
            and re.match(r"^-+\s*$", lines[i + 1])
        ):
            out.append(line)
            out.append(lines[i + 1])
            i += 2
            section_lines: list[str] = []
            while i < len(lines):
                if (
                    lines[i] in _NUMPY_SECTIONS
                    and i + 1 < len(lines)
                    and re.match(r"^-+\s*$", lines[i + 1])
                ):
                    break
                section_lines.append(lines[i])
                i += 1
            first_content = next((ln for ln in section_lines if ln.strip()), None)
            if first_content and re.match(r"^    \S", first_content):
                for sl in section_lines:
                    if sl.strip() and sl.startswith("    "):
                        out.append(sl[4:])
                    else:
                        out.append(sl)
            else:
                out.extend(section_lines)
            continue
        out.append(line)
        i += 1
    return "\n".join(out)


def normalize_sections(doc: str) -> str:
    """Apply mechanical NumPy normalizations to a docstring.
"""
    if not doc:
        return doc

    text = doc.replace("\r\n", "\n")

    text = _GOOGLE_ARGS.sub(r"\1Parameters", text)
    text = _GOOGLE_RETURNS.sub(r"\1Returns", text)
    text = _GOOGLE_RAISES.sub(r"\1Raises", text)
    text = _GOOGLE_YIELDS.sub(r"\1Yields", text)

    text = _normalize_google_param_lines(text)

    def fix_colon_heading(match: re.Match[str]) -> str:
        indent, title = match.group(1), match.group(2)
        if title == "Args":
            title = "Parameters"
        return f"{indent}{title}"

    text = _SECTION_HEADING.sub(fix_colon_heading, text)

    def fix_underline(match: re.Match[str]) -> str:
        indent, title = match.group(1), match.group(2)
        return f"{indent}{title}\n{_section_underline(title)}"

    text = _SECTION_UNDERLINE.sub(fix_underline, text)

    text = _fix_section_header_gaps(text)
    text = _dedent_section_headers(text)
    text = _normalize_parameter_sections(text)
    text = _fix_returns_indent(text)
    text = re.sub(r", default=", ", default ", text)

    return text.strip() + ("\n" if text.strip() else "")


def _split_sections(doc: str) -> tuple[str, dict[str, str]]:
    """Split docstring into summary and section bodies.
"""
    if not doc:
        return "", {}

    text = normalize_sections(doc)
    section_pattern = re.compile(
        r"^([ \t]*)(Parameters|Returns|Yields|Raises|Warns|Notes|References|Examples|Attributes)\n[ \t]*-+\s*$",
        re.MULTILINE,
    )
    matches = list(section_pattern.finditer(text))
    if not matches:
        return text.strip(), {}

    summary = text[: matches[0].start()].strip()
    sections: dict[str, str] = {}

    for index, match in enumerate(matches):
        title = match.group(2)
        start = match.end()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        sections[title] = text[start:end].strip()

    return summary, sections


def merge_numpy_docstrings(wrapper_doc: str, impl_doc: str) -> str:
    """Merge wrapper and implementation docstrings without duplicating sections.
"""
    wrapper_doc = normalize_sections(wrapper_doc or "")
    impl_doc = normalize_sections(impl_doc or "")

    wrapper_summary, wrapper_sections = _split_sections(wrapper_doc)
    impl_summary, impl_sections = _split_sections(impl_doc)

    summary = wrapper_summary or impl_summary
    merged_sections = {**impl_sections, **wrapper_sections}

    lines: list[str] = []
    if summary:
        lines.append(summary)

    for title in ("Parameters", "Returns", "Yields", "Raises", "Warns", "Notes", "References", "Examples", "Attributes"):
        body = merged_sections.get(title, "").strip()
        if not body:
            continue
        lines.extend(format_section(title, body.splitlines()))

    return "\n".join(lines).strip() + "\n" if lines else ""


def _strip_prose_default(desc: str) -> str:
    text = _PROSE_DEFAULT.sub(". ", desc).strip()
    text = re.sub(r"\.\s*\.", ".", text)
    if text.endswith("."):
        return text
    return text + "." if text else ""


def _signature_param_map(func: Callable[..., Any]) -> dict[str, inspect.Parameter]:
    try:
        sig = inspect.signature(func)
    except (TypeError, ValueError):
        return {}
    return {
        name: param
        for name, param in sig.parameters.items()
        if name not in {"self", "cls"}
    }


def _resolve_type_annotation(
    name: str,
    doc_type: Optional[str],
    sig_params: dict[str, inspect.Parameter],
) -> str:
    doc_type = (doc_type or "Any").strip()
    param = sig_params.get(name)
    if param is None:
        return doc_type

    if param.default is inspect.Parameter.empty:
        if "optional" in doc_type.lower() or "default" in doc_type.lower():
            return doc_type
        return f"{doc_type.split(',')[0].strip()}, optional"

    default_val = param.default
    if isinstance(default_val, (str, int, float, bool, dict, list, tuple, type(None))):
        return infer_type_annotation(default_val, has_default=True)

    base = doc_type.split(",")[0].strip() or "Any"
    return f"{base}, default {repr(default_val)}"


def adapt_impl_doc_for_method(
    impl: Callable[..., Any],
    method: Callable[..., Any],
    *,
    drop_params: FrozenSet[str] = _DEFAULT_DROP_PARAMS,
    return_type: Optional[str] = "Sumstats",
    return_desc: Optional[str] = None,
    inplace: bool = True,
) -> str:
    """Adapt an implementation docstring for a Sumstats method wrapper.
"""
    if parse_docstring is None:
        return normalize_sections(impl.__doc__ or "")

    impl_doc = normalize_sections(impl.__doc__ or "")
    wrapper_doc = normalize_sections(method.__doc__ or "")

    wrapper_summary, _ = _split_sections(wrapper_doc)
    summary, sections = _split_sections(impl_doc)

    parsed = parse_docstring(impl_doc)
    sig_params = _signature_param_map(impl)

    meth_notes: list[str] = []
    param_lines: list[str] = []

    for param in parsed.params:
        if param.arg_name in drop_params:
            continue
        if not re.match(r"^[A-Za-z_]\w*$", param.arg_name):
            continue
        desc = _strip_prose_default(param.description or "")
        if _METH_NOTE.search(desc):
            meth_notes.append(desc.strip())
            desc = _strip_prose_default(re.sub(r"When called via :meth:`[^`]+`\(\),?\s*", "", desc, flags=re.IGNORECASE))
        type_ann = _resolve_type_annotation(param.arg_name, param.type_name, sig_params)
        default_val = ...
        sig_param = sig_params.get(param.arg_name)
        if sig_param is not None and sig_param.default is not inspect.Parameter.empty:
            default_val = sig_param.default
        param_lines.extend(format_parameter(param.arg_name, desc, default_val, type_annotation=type_ann))

    notes_parts: list[str] = []
    if parsed.meta and parsed.long_description:
        notes_parts.append(parsed.long_description.strip())

    for extra_section in ("Notes", "Raises", "Examples"):
        body = sections.get(extra_section, "").strip()
        if body:
            notes_parts.append(body)

    if parsed.returns and parsed.returns.description:
        ret_desc = parsed.returns.description.strip()
        for line in ret_desc.splitlines():
            line = line.strip()
            if _METH_NOTE.search(line):
                meth_notes.append(line)
        if parsed.returns.description and _METH_NOTE.search(parsed.returns.description):
            pass  # handled above

    notes_parts.extend(meth_notes)
    if sections.get("Notes"):
        notes_parts.append(sections["Notes"].strip())

    lines: list[str] = []
    final_summary = wrapper_summary or summary or (parsed.short_description or "").strip()
    if final_summary:
        lines.append(final_summary)
        if parsed.long_description and parsed.long_description.strip() not in final_summary:
            long_desc = parsed.long_description.strip()
            if long_desc and long_desc not in final_summary:
                lines.append("")
                lines.append(long_desc)

    if param_lines:
        lines.append("")
        lines.append("Parameters")
        lines.append("----------")
        lines.extend(param_lines)

    if return_type is not None:
        desc = return_desc
        if not desc:
            desc = (
                "The Sumstats object (``self``)."
                if inplace
                else "New Sumstats when ``inplace=False``; otherwise ``None``."
            )
        lines.append("")
        lines.append("Returns")
        lines.append("-------")
        lines.append(return_type)
        lines.append(f"    {desc.strip()}")

    unique_notes: list[str] = []
    seen: set[str] = set()
    for part in notes_parts:
        part = part.strip()
        if part and part not in seen:
            seen.add(part)
            unique_notes.append(part)
    if unique_notes:
        lines.append("")
        lines.append("Notes")
        lines.append("-----")
        for part in unique_notes:
            for note_line in part.splitlines():
                stripped = note_line.strip()
                if stripped:
                    lines.append(f"    {stripped}")

    return "\n".join(lines).strip() + "\n" if lines else ""

