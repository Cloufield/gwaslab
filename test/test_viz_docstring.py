"""Tests for registry-aligned plot docstrings (@add_viz_doc)."""

import re

import pytest

from gwaslab.g_Sumstats import Sumstats
from gwaslab.viz.viz_aux_doc import build_viz_docstring, get_doc_params_manager, get_doc_hidden_params
from gwaslab.viz.viz_plot_mqqplot import _mqqplot


_DOC_PM = get_doc_params_manager()
_PARAM_LINE = re.compile(r"^(\w+) : .+$", re.MULTILINE)
_DEFAULT_IS_PROSE = re.compile(r"\. Default is ", re.MULTILINE)


_PARAM_SECTION = re.compile(
    r"^\s*Parameters\n\s*-+\n(.*?)(?=^\s*(Returns|Yields|Raises|Notes|Examples|References|Attributes)\n\s*-+\n|\Z)",
    re.MULTILINE | re.DOTALL,
)


def _parameters_section(doc: str) -> str:
    match = _PARAM_SECTION.search(doc or "")
    return match.group(1) if match else ""


def _doc_param_names(doc: str) -> set[str]:
    return set(_PARAM_LINE.findall(_parameters_section(doc)))


@pytest.fixture(scope="module")
def plot_region_doc():
    return Sumstats.plot_region.__doc__ or ""


def test_build_viz_docstring_plot_region_has_core_params():
    doc = build_viz_docstring(_DOC_PM, "plot_region", "r", _mqqplot)
    names = _doc_param_names(doc)
    for required in ("region", "vcf_path", "ld_link_thr", "marker_size", "sig_level"):
        assert required in names, f"missing {required} in generated doc"


def test_plot_region_doc_excludes_qq_panel(plot_region_doc):
    names = _doc_param_names(plot_region_doc)
    for banned in ("qtitle", "qtitle_pad", "qq_scatter_kwargs", "mqqratio", "mtitle", "density_color", "pinpoint", "gc", "stratified"):
        assert banned not in names, f"{banned} should not appear in plot_region doc"


def test_plot_region_keeps_shared_params(plot_region_doc):
    names = _doc_param_names(plot_region_doc)
    for kept in ("chrpad", "skip", "windowsizekb", "suggestive_sig_line", "scatter_kwargs"):
        assert kept in names, f"{kept} should remain in plot_region doc"


def test_plot_manhattan_doc_excludes_regional_panel():
    doc = Sumstats.plot_manhattan.__doc__ or ""
    names = _doc_param_names(doc)
    for banned in ("region", "ld_block", "rr_path", "vcf_path", "track_n", "gtf_path"):
        assert banned not in names, f"{banned} should not appear in plot_manhattan doc"
    assert "qtitle" in names or "mqqratio" in names


def test_plot_manhattan_allowed_excludes_regional_panel():
    allowed = _DOC_PM.allowed("plot_manhattan", None) or set()
    for banned in ("region", "ld_block", "vcf_path"):
        assert banned not in allowed


def test_param_groups_no_overlap():
    groups = _DOC_PM._param_groups
    seen = {}
    for group_name, members in groups.items():
        member_set = set(members)
        for name in member_set:
            assert name not in seen, f"{name} in both {seen[name]} and {group_name}"
            seen[name] = group_name


def test_plot_region_doc_excludes_highlight_family(plot_region_doc):
    names = _doc_param_names(plot_region_doc)
    highlight_family = {n for n in names if n.startswith("highlight")}
    assert highlight_family == set()


def test_plot_region_allowed_excludes_highlight_family():
    allowed = _DOC_PM.allowed("plot_region", "r") or set()
    highlight_family = {n for n in allowed if n.startswith("highlight")}
    assert highlight_family == set()


def test_plot_region_doc_excludes_removed_sig_level_aliases(plot_region_doc):
    names = _doc_param_names(plot_region_doc)
    assert "sig_level_plot" not in names
    assert "sig_level_lead" not in names


def test_plot_region_doc_no_default_is_prose(plot_region_doc):
    assert not _DEFAULT_IS_PROSE.search(plot_region_doc)


def test_plot_region_marker_size_default(plot_region_doc):
    assert "marker_size : list, default [40, 65]" in plot_region_doc or "[40, 65]" in plot_region_doc


def test_plot_region_scatter_kwargs_banned_subkeys(plot_region_doc):
    assert "Do not set in ``scatter_kwargs``:" in plot_region_doc
    assert "s" in plot_region_doc


def test_plot_region_allowed_subset_of_doc(plot_region_doc):
    allowed = _DOC_PM.allowed("plot_region", "r") or set()
    doc_names = _doc_param_names(plot_region_doc)
    hidden = get_doc_hidden_params()
    for name in doc_names:
        assert name in allowed or name.startswith("_"), f"{name} not in allowed"
    for name in allowed:
        if name.startswith("_") or name in hidden:
            continue
        assert name in doc_names, f"{name} in allowed but missing from doc"


def test_plot_region_doc_hides_legacy_aliases(plot_region_doc):
    hidden = get_doc_hidden_params() - {"insumstats", "log", "object"}
    for name in hidden:
        assert name not in _doc_param_names(plot_region_doc), f"legacy alias {name} should be hidden"


_PARAM_DESC = re.compile(r"^(\w+) : .+\n    (.+)$", re.MULTILINE)


def _doc_param_descriptions(doc: str) -> dict[str, str]:
    return {name: desc.strip() for name, desc in _PARAM_DESC.findall(_parameters_section(doc))}


def test_plot_region_doc_param_count_reasonable(plot_region_doc):
    count = len(_doc_param_names(plot_region_doc))
    assert 145 <= count <= 165, f"unexpected param count: {count}"


def test_plot_region_doc_params_have_description(plot_region_doc):
    descriptions = _doc_param_descriptions(plot_region_doc)
    for name in _doc_param_names(plot_region_doc):
        assert descriptions.get(name), f"empty description line for {name} in doc"


def test_plot_region_doc_ld_params_precise_descriptions(plot_region_doc):
    descriptions = _doc_param_descriptions(plot_region_doc)
    ld_pairs = descriptions["region_ld_colors"].lower()
    assert "len(region_ld_threshold)" in ld_pairs or "len(thresholds)" in ld_pairs
    assert "region_ld_threshold" in ld_pairs
    shapes = descriptions["region_marker_shapes"].lower()
    assert "shape" in shapes
    assert "color" in shapes
    assert "region_ld_threshold" in shapes or "region_ld_colors" in shapes
    ld_link = descriptions["ld_link"].lower()
    assert "region_ld_threshold" in ld_link
    assert "region_ld_colors" in ld_link
    assert "region_ld_threshold" in descriptions["ld_link_thr"].lower()


def test_plot_region_doc_no_duplicate_impl_summary(plot_region_doc):
    assert plot_region_doc.count("Create an MQQ plot") == 0
    assert plot_region_doc.strip().startswith("Regional association plot")


def test_sumstats_plot_methods_have_parameters_section():
    for method in (
        Sumstats.plot_mqq,
        Sumstats.plot_manhattan,
        Sumstats.plot_qq,
        Sumstats.plot_region,
        Sumstats.plot_snp_density,
    ):
        doc = method.__doc__ or ""
        assert "Parameters" in doc, f"{method.__name__} missing Parameters section"
        assert _doc_param_names(doc), f"{method.__name__} has empty Parameters list"


@pytest.mark.parametrize(
    "method_name,key,mode",
    [
        ("plot_daf", "plot_daf", None),
        ("plot_gwheatmap", "plot_gwheatmap", None),
        ("plot_trumpet", "plot_trumpet", None),
        ("plot_phenogram", "plot_phenogram", None),
        ("plot_ld_block", "plot_ld_block", None),
        ("plot_effect", "plot_effect", None),
        ("plot_sankey", "plot_sankey", None),
        ("plot_associations", "plot_associations", None),
        ("plot_pipcs", "plot_pipcs", None),
    ],
)
def test_sumstats_viz_doc_has_parameters(method_name, key, mode):
    doc = getattr(Sumstats, method_name).__doc__ or ""
    assert "Parameters" in doc, f"{method_name} missing Parameters"
    assert len(_doc_param_names(doc)) >= 10, f"{method_name} doc too short"


def test_plot_trumpet_doc_not_empty():
    doc = Sumstats.plot_trumpet.__doc__ or ""
    names = _doc_param_names(doc)
    assert len(names) >= 50
    assert "prevalence" in names
    assert "mode" in names


def test_plot_pipcs_doc_includes_title():
    names = _doc_param_names(Sumstats.plot_pipcs.__doc__ or "")
    assert "title" in names
    assert "title_kwargs" in names
    assert "palette" not in names


def test_plot_qq_doc_excludes_manhattan_only():
    names = _doc_param_names(Sumstats.plot_qq.__doc__ or "")
    for banned in ("windowsizekb", "mqqratio", "region", "density_color"):
        assert banned not in names, f"{banned} in plot_qq doc"
    assert "qtitle" in names or "qq_scatter_kwargs" in names


def test_plot_snp_density_doc_includes_brisbane():
    names = _doc_param_names(Sumstats.plot_snp_density.__doc__ or "")
    for kept in ("bwindowsizekb", "density_color", "density_threshold"):
        assert kept in names, f"{kept} missing from plot_snp_density doc"
    assert "qtitle" not in names


def test_plot_manhattan_excludes_panel_layout():
    allowed = _DOC_PM.allowed("plot_manhattan", None) or set()
    for banned in ("height_ratios", "align_xaxis", "arc_height", "file_format"):
        assert banned not in allowed


def test_plot_panels_allowed_has_layout_params():
    allowed = _DOC_PM.allowed("plot_panels", None) or set()
    for kept in ("height_ratios", "hspace", "align_xaxis"):
        assert kept in allowed


def test_plot_forest_allowed_has_column_mapping():
    allowed = _DOC_PM.allowed("plot_forest", None) or set()
    for kept in ("beta_col", "se_col", "colors"):
        assert kept in allowed


def test_plot_chromatin_registry_exists():
    allowed = _DOC_PM.allowed("plot_chromatin", None) or set()
    assert "region_chromatin_labels" in allowed
    assert "xlim_i" in allowed


def test_plot_power_registry_modes():
    q_allowed = _DOC_PM.allowed("plot_power", "q") or set()
    b_allowed = _DOC_PM.allowed("plot_power", "b") or set()
    assert "prevalences" in b_allowed
    assert "prevalences" not in q_allowed


def test_gl_plot_forest_doc_from_registry():
    import gwaslab as gl
    doc = gl.plot_forest.__doc__ or ""
    assert "Parameters" in doc
    assert "beta_col" in _doc_param_names(doc)


@pytest.mark.parametrize(
    "method_name,key,mode,forbidden,required_min",
    [
        ("plot_daf", "plot_daf", None, (), 10),
        ("plot_gwheatmap", "plot_gwheatmap", None, (), 10),
        ("plot_trumpet", "plot_trumpet", None, ("marker_size", "style"), 50),
        ("plot_phenogram", "plot_phenogram", None, ("ancestry_col",), 50),
        ("plot_ld_block", "plot_ld_block", None, (), 10),
        ("plot_effect", "plot_effect", None, (), 10),
        ("plot_sankey", "plot_sankey", None, ("columns",), 10),
        ("plot_associations", "plot_associations", None, (), 10),
        ("plot_pipcs", "plot_pipcs", None, ("palette", "edgecolor"), 10),
        ("plot_qq", "plot_qq", "qq", ("windowsizekb", "region"), 10),
        ("plot_snp_density", "plot_snp_density", "b", ("qtitle",), 10),
    ],
)
def test_plot_doc_alignment(method_name, key, mode, forbidden, required_min):
    doc = getattr(Sumstats, method_name).__doc__ or ""
    names = _doc_param_names(doc)
    allowed = _DOC_PM.allowed(key, mode) or set()
    hidden = get_doc_hidden_params()
    assert len(names) >= required_min
    for ban in forbidden:
        assert ban not in names, f"{ban} should not be in {method_name} doc"
    for name in names:
        assert name in allowed or name.startswith("_"), f"{name} not allowed for {key}"
    for name in allowed:
        if name.startswith("_") or name in hidden:
            continue
        assert name in names, f"{name} missing from {method_name} doc"


def test_gl_wrappers_have_parameters_section():
    import gwaslab as gl
    for fn in (gl.plot_miami2, gl.compare_effect, gl.plot_lead_overlap, gl.plot_rg, gl.plot_stacked_mqq, gl.plot_power):
        doc = fn.__doc__ or ""
        assert "Parameters" in doc, f"{fn.__name__} missing Parameters"
        assert _doc_param_names(doc), f"{fn.__name__} empty param list"


def test_plot_region_any_optional_fraction(plot_region_doc):
    lines = [l for l in (plot_region_doc or "").splitlines() if re.match(r"^\w+ : ", l)]
    any_opt = sum(1 for l in lines if "Any, optional" in l)
    assert len(lines) > 0
    assert any_opt / len(lines) < 0.25, f"too many Any, optional: {any_opt}/{len(lines)}"


def test_gl_plot_panels_returns_excludes_examples():
    import gwaslab as gl
    from gwaslab.viz.viz_aux_doc import build_gl_viz_docstring
    from gwaslab.viz.viz_plot_stackedpanel import plot_panels

    doc = build_gl_viz_docstring(
        "plot_panels",
        plot_panels,
        None,
        summary="Stack multi-panel figure from Panel objects.",
        positional_params=(("panels", "list", "Panel layout objects from :class:`gwaslab.Panel`."),),
    )
    assert "Returns" in doc
    returns_body = doc.split("Returns\n-------\n", 1)[1]
    if "\n\n" in returns_body:
        returns_body = returns_body.split("\n\n", 1)[0]
    assert "Examples" not in returns_body
    assert ">>> import gwaslab" not in returns_body
    assert "Create panels" not in returns_body

    runtime_doc = gl.plot_panels.__doc__ or ""
    runtime_returns = runtime_doc.split("Returns\n-------\n", 1)[1].split("\n\n", 1)[0]
    assert ">>> import gwaslab" not in runtime_returns

