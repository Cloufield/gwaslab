"""Tests for Sumstats @add_doc adaptation (adapt_impl_doc_for_method)."""

from gwaslab.info.g_numpy_doc import adapt_impl_doc_for_method
from gwaslab.qc.qc_fix_sumstats import _fix_ID
from docstring_parser import parse


def _method_stub(name: str):
    def method(**_kwargs):
        return None

    method.__name__ = name
    return method


def test_adapt_fix_id_drops_sumstats_obj():
    doc = adapt_impl_doc_for_method(_fix_ID, _method_stub("fix_id"))
    params = {p.arg_name for p in parse(doc).params}
    assert "sumstats_obj" not in params
    assert "log" not in params
    assert "fixprefix" in params


def test_adapt_fix_id_returns_sumstats():
    doc = adapt_impl_doc_for_method(_fix_ID, _method_stub("fix_id"))
    assert "Returns" in doc
    assert "Sumstats" in doc
    ret = parse(doc).returns
    assert ret is not None
    assert ret.type_name == "Sumstats"


def test_adapt_strips_prose_default():
    def impl(x=None):
        pass

    impl.__doc__ = (
        "Do thing.\n\n"
        "Parameters\n"
        "----------\n"
        "x : str, optional\n"
        "    Path to file. Default is None.\n"
    )
    doc = adapt_impl_doc_for_method(impl, impl)
    assert "Default is" not in doc
    assert "Path to file." in doc


def test_sumstats_fix_id_runtime_doc():
    from gwaslab.g_Sumstats import Sumstats

    doc = Sumstats.fix_id.__doc__ or ""
    params = {p.arg_name for p in parse(doc).params}
    assert "sumstats_obj" not in params
    assert parse(doc).returns.type_name == "Sumstats"


def test_sumstats_get_lead_returns_dataframe():
    from gwaslab.g_Sumstats import Sumstats

    ret = parse(Sumstats.get_lead.__doc__ or "").returns
    assert ret is not None
    assert "DataFrame" in (ret.type_name or "")


def test_sumstats_filter_value_return_type():
    from gwaslab.g_Sumstats import Sumstats

    ret = parse(Sumstats.filter_value.__doc__ or "").returns
    assert ret is not None
    assert "Sumstats" in (ret.type_name or "")


def test_gl_compare_effect_documents_positional_paths():
    import gwaslab as gl

    doc = gl.compare_effect.__doc__ or ""
    assert "path1" in doc
    assert "path2" in doc


def test_plot_region_returns_excludes_log():
    from gwaslab.g_Sumstats import Sumstats

    doc = Sumstats.plot_region.__doc__ or ""
    if "Returns" in doc:
        returns_block = doc.split("Returns\n-------\n", 1)[1].split("\n\n")[0]
        assert "log :" not in returns_block.lower()


def test_gl_plot_sankey_documents_columns():
    import gwaslab as gl

    doc = gl.plot_sankey.__doc__ or ""
    assert "columns" in doc
    assert "data" in doc

    from gwaslab.g_Sumstats import Sumstats

    for name in ("fix_chr", "filter_in", "estimate_h2_by_ldsc", "get_lead"):
        params = {p.arg_name for p in parse(getattr(Sumstats, name).__doc__ or "").params}
        assert "sumstats_obj" not in params
        assert "insumstats_or_dataframe" not in params
