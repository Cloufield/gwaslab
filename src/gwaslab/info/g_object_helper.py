from typing import Callable, Any, Optional, FrozenSet, Dict
import functools
import matplotlib.pyplot as plt

from gwaslab.info.g_numpy_doc import adapt_impl_doc_for_method, merge_numpy_docstrings

_DEFAULT_DROP = frozenset({
    "sumstats_obj",
    "log",
    "insumstats",
    "insumstats_or_dataframe",
    "meta",
})

# Per-method docstring overrides keyed by Sumstats method name.
SUMSTATS_DOC_OVERRIDES: Dict[str, Dict[str, Any]] = {
    "summary": {
        "return_type": "dict",
        "return_desc": "QC summary dictionary from ``summarize()`` (also stored in ``self.meta``).",
        "inplace": False,
    },
    "lookup_status": {
        "return_type": "pandas.DataFrame",
        "return_desc": "Status-code explanations for each digit position in the status column.",
        "inplace": False,
    },
    "check_sumstats_qc_status": {
        "return_type": "dict",
        "return_desc": "QC status report from metadata and data checks.",
        "inplace": False,
    },
    "update_meta": {
        "return_type": "None",
        "return_desc": "Updates ``self.meta`` in place; returns None.",
        "inplace": True,
    },
    "set_build": {
        "return_type": "None",
        "return_desc": "Updates ``self.build`` and chromosome columns in place; returns None.",
        "inplace": True,
    },
    "infer_build": {
        "return_type": "None",
        "return_desc": "Infers build and updates ``self.build`` in place; returns None.",
        "inplace": True,
    },
    "liftover": {
        "return_type": "None",
        "return_desc": "Liftover coordinates in place; returns None.",
        "inplace": True,
    },
    "fill_data": {
        "return_type": "None",
        "return_desc": "Fills derived statistics columns in place; returns None.",
        "inplace": True,
    },
    "reload": {
        "return_type": "None",
        "return_desc": "Reloads offloaded data into ``self.data``; returns None.",
        "inplace": True,
    },
    "flip_snpid": {
        "return_type": "None",
        "return_desc": "Flips SNPID column in place; returns None.",
        "inplace": True,
    },
    "strip_snpid": {
        "return_type": "None",
        "return_desc": "Strips SNPID prefixes in place; returns None.",
        "inplace": True,
    },
    "get_gc": {
        "return_type": "float",
        "return_desc": "Genomic inflation factor (lambda GC); also stored in ``self.meta``.",
        "inplace": False,
    },
    "get_density": {
        "return_type": "None",
        "return_desc": "Adds density columns to ``self.data`` in place; returns None.",
        "inplace": True,
    },
    "get_per_snp_r2": {
        "return_type": "None",
        "return_desc": "Adds per-SNP R-squared column in place; returns None.",
        "inplace": True,
    },
    "get_ess": {
        "return_type": "None",
        "return_desc": "Adds effective sample size column in place; returns None.",
        "inplace": True,
    },
    "infer_ancestry": {
        "return_type": "Any",
        "return_desc": "Ancestry inference result from the reference panel comparison.",
        "inplace": False,
    },
    "get_lead": {
        "return_type": "pandas.DataFrame or Sumstats",
        "return_desc": "Lead variants as a DataFrame, or a filtered Sumstats when ``gls=True``.",
        "inplace": False,
    },
    "get_novel": {
        "return_type": "pandas.DataFrame or tuple",
        "return_desc": "Novel loci table, or ``(novel, known)`` when ``output_known=True``.",
        "inplace": False,
    },
    "get_top": {
        "return_type": "pandas.DataFrame or Sumstats",
        "return_desc": "Top variants as a DataFrame, or a filtered Sumstats when ``gls=True``.",
        "inplace": False,
    },
    "report": {
        "return_type": "str",
        "return_desc": "Path to the generated QC report HTML (or PDF when supported).",
        "inplace": False,
    },
    "view_sumstats": {
        "return_type": "pandas.DataFrame",
        "return_desc": "Filtered or full sumstats dataframe.",
        "inplace": False,
    },
    "estimate_h2_by_ldsc": {
        "return_type": "pandas.DataFrame",
        "return_desc": "LDSC h2 results; also stored on ``self.ldsc_h2`` / ``self.ldsc_h2_results``.",
        "inplace": False,
    },
    "estimate_h2_cts_by_ldsc": {
        "return_type": "pandas.DataFrame",
        "return_desc": "Cell-type-specific LDSC results; stored on ``self.ldsc_h2_cts``.",
        "inplace": False,
    },
    "estimate_partitioned_h2_by_ldsc": {
        "return_type": "pandas.DataFrame",
        "return_desc": "Partitioned LDSC results; stored on ``self.ldsc_partitioned_h2_results``.",
        "inplace": False,
    },
    "estimate_rg_by_ldsc": {
        "return_type": "pandas.DataFrame",
        "return_desc": "Genetic correlation results appended to ``self.ldsc_rg``.",
        "inplace": False,
    },
    "clump": {
        "return_type": "None",
        "return_desc": "Stores clumping results on ``self.clumps``; returns None.",
        "inplace": True,
    },
    "to_format": {
        "return_type": "str or None",
        "return_desc": "Output file path when writing to disk.",
        "inplace": False,
    },
    "to_pickle": {
        "return_type": "None",
        "return_desc": "Pickles the object to disk; returns None.",
        "inplace": False,
    },
    "to_gsf": {
        "return_type": "None",
        "return_desc": "Writes GWASLab standard format; returns None.",
        "inplace": False,
    },
    "filter_value": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None`` (updates ``self`` in place).",
        "inplace": False,
    },
    "filter_in": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
    "filter_out": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
    "filter_region_in": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
    "filter_region_out": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
    "filter_region": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
    "random_variants": {
        "return_type": "Sumstats or None",
        "return_desc": "New Sumstats when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
}

_PAIR_DEFAULT_DROP = _DEFAULT_DROP | frozenset({
    "gls",
    "glsp",
    "sumstatspair_object",
    "sumstats_pair",
    "sumstats_multi",
})

SUMSTATS_PAIR_DOC_OVERRIDES: Dict[str, Dict[str, Any]] = {
    "clump": {
        "return_type": "None",
        "return_desc": "Stores clumping results on ``self.clumps``; returns None.",
        "inplace": True,
    },
    "to_coloc": {
        "return_type": "None",
        "return_desc": "Stores coloc file list in ``self.coloc``; returns None.",
        "inplace": True,
    },
    "to_mesusie": {
        "return_type": "None",
        "return_desc": "Stores MESuSiE file list in ``self.mesusie``; returns None.",
        "inplace": True,
    },
    "run_mesusie": {
        "return_type": "None",
        "return_desc": "Stores finemapping results in ``self.mesusie_res``; returns None.",
        "inplace": True,
    },
    "run_multisusie_rss": {
        "return_type": "None",
        "return_desc": "Stores MultiSuSiE results in ``self.multisusie_res``; returns None.",
        "inplace": True,
    },
    "run_coloc_susie": {
        "return_type": "None",
        "return_desc": "Stores coloc.susie results in ``self.coloc_susie_res``; returns None.",
        "inplace": True,
    },
    "run_two_sample_mr": {
        "return_type": "None",
        "return_desc": "Stores MR results in ``self.mr``; returns None.",
        "inplace": True,
    },
    "run_ccgwas": {
        "return_type": "None",
        "return_desc": "Runs CC-GWAS from the merged pair; returns None.",
        "inplace": True,
    },
    "read_pipcs": {
        "return_type": "None",
        "return_desc": "Loads PIP and credible-set columns into ``self.mesusie_res``; returns None.",
        "inplace": True,
    },
    "filter_value": {
        "return_type": "SumstatsPair or None",
        "return_desc": "New SumstatsPair when ``inplace=False``; otherwise ``None``.",
        "inplace": False,
    },
    "extract_with_ld_proxy": {
        "return_type": "pandas.DataFrame",
        "return_desc": "Merged table with LD proxy replacements for missing variants.",
        "inplace": False,
    },
    "run_meta_analysis": {
        "return_type": "Sumstats",
        "return_desc": "Fixed- or random-effects meta-analysis across both studies.",
        "inplace": False,
    },
    "to_pickle": {
        "return_type": "None",
        "return_desc": "Pickles the SumstatsPair to disk; returns None.",
        "inplace": False,
    },
    "offload": {
        "return_type": "None",
        "return_desc": "Writes ``self.data`` to temp storage and frees memory; returns None.",
        "inplace": True,
    },
    "reload": {
        "return_type": "None",
        "return_desc": "Reloads offloaded data into ``self.data``; returns None.",
        "inplace": True,
    },
}


def add_sumstats_doc(
    impl: Callable[..., Any],
    *,
    drop_params: Optional[FrozenSet[str]] = None,
    return_type: Optional[str] = None,
    return_desc: Optional[str] = None,
    inplace: Optional[bool] = None,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Attach an impl docstring adapted for a Sumstats method wrapper.
"""

    def deco(func: Callable[..., Any]) -> Callable[..., Any]:
        overrides = dict(SUMSTATS_DOC_OVERRIDES.get(func.__name__, {}))
        if drop_params is not None:
            overrides["drop_params"] = drop_params
        if return_type is not None:
            overrides["return_type"] = return_type
        if return_desc is not None:
            overrides["return_desc"] = return_desc
        if inplace is not None:
            overrides["inplace"] = inplace

        opts: Dict[str, Any] = {
            "drop_params": drop_params or overrides.pop("drop_params", _DEFAULT_DROP),
            "return_type": overrides.pop("return_type", "Sumstats"),
            "return_desc": overrides.pop("return_desc", None),
            "inplace": overrides.pop("inplace", True),
        }

        wrapper_doc = func.__doc__ or ""
        adapted = adapt_impl_doc_for_method(impl, func, **opts)
        func.__doc__ = merge_numpy_docstrings(wrapper_doc, adapted) if wrapper_doc.strip() else adapted
        return func

    return deco


def add_doc(src: Callable[..., Any]) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Backward-compatible alias for :func:`add_sumstats_doc`.
"""
    return add_sumstats_doc(src)


def add_pair_doc(
    impl: Callable[..., Any],
    *,
    drop_params: Optional[FrozenSet[str]] = None,
    return_type: Optional[str] = None,
    return_desc: Optional[str] = None,
    inplace: Optional[bool] = None,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Attach an impl docstring adapted for a SumstatsPair method wrapper.
"""

    def deco(func: Callable[..., Any]) -> Callable[..., Any]:
        overrides = dict(SUMSTATS_PAIR_DOC_OVERRIDES.get(func.__name__, {}))
        if drop_params is not None:
            overrides["drop_params"] = drop_params
        if return_type is not None:
            overrides["return_type"] = return_type
        if return_desc is not None:
            overrides["return_desc"] = return_desc
        if inplace is not None:
            overrides["inplace"] = inplace

        opts: Dict[str, Any] = {
            "drop_params": drop_params or overrides.pop("drop_params", _PAIR_DEFAULT_DROP),
            "return_type": overrides.pop("return_type", "SumstatsPair"),
            "return_desc": overrides.pop("return_desc", None),
            "inplace": overrides.pop("inplace", True),
        }

        wrapper_doc = func.__doc__ or ""
        adapted = adapt_impl_doc_for_method(impl, func, **opts)
        func.__doc__ = merge_numpy_docstrings(wrapper_doc, adapted) if wrapper_doc.strip() else adapted
        return func

    return deco


def add_viz_doc(
    impl: Callable[..., Any],
    key: str,
    mode: Optional[str] = None,
    *,
    summary: str = "",
    hide_params: Optional[FrozenSet[str]] = None,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Attach a registry-aligned docstring for a Sumstats plot wrapper.
"""

    def deco(func: Callable[..., Any]) -> Callable[..., Any]:
        from gwaslab.viz.viz_aux_doc import _DOC_PM, build_viz_docstring

        doc = build_viz_docstring(
            _DOC_PM,
            key,
            mode,
            impl,
            wrapper_doc=summary or func.__doc__ or "",
            hide_params=hide_params,
        )
        func.__doc__ = doc
        wrapped = getattr(func, "__wrapped__", None)
        if wrapped is not None:
            wrapped.__doc__ = doc
        return func

    return deco


def suppress_display(func: Callable[..., Any]) -> Callable[..., Any]:
    """Close matplotlib figures returned by a plotting wrapper.

Parameters
----------
func : callable
    Plotting function that may return a ``Figure`` or ``(Figure, ...)`` tuple.
Returns
-------
callable
    Wrapped function with the same call signature as ``func``.
"""
    @functools.wraps(func)
    def _wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        try:
            from matplotlib.figure import Figure
            if isinstance(result, Figure):
                plt.close(result)
            elif isinstance(result, tuple) and len(result) > 0 and isinstance(result[0], Figure):
                plt.close(result[0])
            return result
        except Exception:
            return result
    return _wrapper
