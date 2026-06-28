from typing import TYPE_CHECKING, Optional, Dict, Any, Union, List
import pandas as pd

from gwaslab.g_Sumstats import Sumstats
from gwaslab.g_Sumstats_polars import Sumstatsp
from gwaslab.g_SumstatsT import SumstatsT
from gwaslab.g_SumstatsPair import SumstatsPair
from gwaslab.g_SumstatsMulti import SumstatsMulti
from gwaslab.g_SumstatsSet import SumstatsSet
from gwaslab.info.g_version import _show_version as show_version
from gwaslab.info.g_version import gwaslab_info

from gwaslab.util.util_in_convert_h2 import h2_obs_to_liab
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2
from gwaslab.util.util_in_convert_h2 import h2_se_to_p
from gwaslab.util.util_in_calculate_power import get_power as _get_power
from gwaslab.util.util_in_calculate_power import get_beta
from gwaslab.util.util_ex_gwascatalog import gwascatalog_trait
from gwaslab.util.util_ex_process_h5 import process_vcf_to_hfd5
from gwaslab.util.rwrapper.util_ex_run_susie import _run_susie_rss as _run_susie_rss_impl
from gwaslab.util.util_in_fill_data import rank_based_int
from gwaslab.util.util_in_simulate import simulate_sumstats_region
from gwaslab.util.util_in_simulate import simulate_sumstats_global

from gwaslab.io.io_read_ldsc import read_ldsc as _read_ldsc
from gwaslab.io.io_read_ldsc import read_popcorn
from gwaslab.io.io_read_ldsc import read_greml
from gwaslab.io.io_to_pickle import dump_pickle
from gwaslab.io.io_to_pickle import load_pickle as _load_pickle
from gwaslab.io.io_gwaslab_standard import load_gsf
from gwaslab.io.io_read_tabular import _read_tabular as read_tabular
from gwaslab.io.io_read_tabular import read_bim
from gwaslab.io.io_read_tabular import read_fam
from gwaslab.io.io_read_tabular import read_psam
from gwaslab.io.io_read_tabular import read_pvar
from gwaslab.io.io_read_tabular import read_bgen_sample
from gwaslab.io.io_gtf import read_gtf
from gwaslab.io.io_gtf import read_gtf_file
from gwaslab.io.io_bigwig_bigbed import read_bigwig
from gwaslab.io.io_bigwig_bigbed import read_bigwig_intervals
from gwaslab.io.io_bigwig_bigbed import read_bigwig_stats
from gwaslab.io.io_bigwig_bigbed import read_bigbed
from gwaslab.io.io_bedpe import read_bedpe
from gwaslab.io.io_ucsc_bed import read_bed

from gwaslab.viz.viz_plot_compare_effect import compare_effect as _compare_effect
from gwaslab.viz.viz_plot_forestplot import plot_forest as _plot_forest
from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2 as _plot_miami2
from gwaslab.viz.viz_plot_lead_overlap import plot_lead_overlap as _plot_lead_overlap
from gwaslab.viz.viz_plot_rg_heatmap import plot_rg as _plot_rg
from gwaslab.viz.viz_plot_stackedregional import plot_stacked_mqq as _plot_stacked_mqq
from gwaslab.viz.viz_plot_trumpetplot import plot_power as _plot_power
from gwaslab.viz.viz_plot_trumpetplot import plot_power_x as _plot_power_x
from gwaslab.viz.viz_plot_scatter_with_reg import scatter as _scatter
from gwaslab.viz.viz_plot_ld_block import plot_ld_block as _plot_ld_block
from gwaslab.viz.viz_plot_track import plot_track
from gwaslab.viz.viz_plot_arc import plot_arc
from gwaslab.viz.viz_aux_panel import Panel
from gwaslab.viz.viz_plot_stackedpanel import plot_panels as _plot_panels
from gwaslab.viz.viz_plot_sankey import plot_sankey as _plot_sankey

from gwaslab.bd.bd_common_data import get_NC_to_chr
from gwaslab.bd.bd_common_data import get_NC_to_number
from gwaslab.bd.bd_common_data import get_chr_to_NC
from gwaslab.bd.bd_common_data import get_number_to_NC
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_common_data import get_number_to_chr
from gwaslab.bd.bd_common_data import get_chr_to_number
from gwaslab.bd.bd_common_data import get_high_ld
from gwaslab.bd.bd_common_data import get_format_dict
from gwaslab.bd.bd_common_data import get_formats_list
from gwaslab.bd.bd_download import update_formatbook
from gwaslab.bd.bd_download import list_formats
from gwaslab.bd.bd_download import list_formats_with_descriptions
from gwaslab.bd.bd_download import check_format as _check_format
from gwaslab.bd.bd_download import check_available_ref
from gwaslab.bd.bd_download import update_available_ref
from gwaslab.bd.bd_download import check_downloaded_ref
from gwaslab.bd.bd_download import download_ref as _download_ref
from gwaslab.bd.bd_download import remove_file
from gwaslab.bd.bd_download import get_path as _get_path
from gwaslab.bd.bd_download import update_record
from gwaslab.bd.bd_download import scan_downloaded_files
from gwaslab.bd.bd_download import add_local_data
from gwaslab.bd.bd_download import remove_local_record
from gwaslab.bd.bd_download import filter_downloaded_registry
from gwaslab.bd.bd_download import infer_registry_kind
from gwaslab.bd.bd_download import set_default_directory
from gwaslab.bd.bd_download import get_default_directory
from gwaslab.bd.bd_download import check_and_download
from gwaslab.bd.bd_config import options
from gwaslab.bd.bd_config import ensure_user_layout
from gwaslab.extension.gwascatalog.sumstats_download import download_sumstats as _download_sumstats
from gwaslab.qc.qc_reserved_headers import researved_header
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config
from gwaslab.view.view_report import generate_qc_report

__version__ = gwaslab_info()["version"]


_viz_params = VizParamsManager()
load_viz_config(_viz_params)


def download_ref(
    name: str,
    directory: Optional[str] = None,
    local_filename: Optional[str] = None,
    overwrite: bool = False,
    log: Optional[Log] = None,
) -> Optional[str]:
    _log = log if log is not None else Log()
    return _download_ref(name, directory=directory, local_filename=local_filename, overwrite=overwrite, log=_log)


download_ref.__doc__ = _download_ref.__doc__


def get_path(
    name: str,
    log: Optional[Log] = None,
    verbose: bool = True,
    raise_on_missing: bool = False,
) -> Union[str, bool]:
    _log = log if log is not None else Log()
    return _get_path(name, log=_log, verbose=verbose, raise_on_missing=raise_on_missing)


get_path.__doc__ = _get_path.__doc__


def read_ldsc(filelist: Optional[List[str]] = None, mode: str = "h2") -> pd.DataFrame:
    return _read_ldsc(filelist=filelist, mode=mode)


read_ldsc.__doc__ = _read_ldsc.__doc__


def load_pickle(path: str) -> Any:
    return _load_pickle(path)


load_pickle.__doc__ = _load_pickle.__doc__


def download_sumstats(accession: str, **kwargs: Any) -> Optional[str]:
    return _download_sumstats(accession, **kwargs)


download_sumstats.__doc__ = _download_sumstats.__doc__


def check_format(path: str, **kwargs: Any) -> Any:
    return _check_format(path, **kwargs)


check_format.__doc__ = _check_format.__doc__


def get_power(*args: Any, **kwargs: Any) -> Any:
    return _get_power(*args, **kwargs)


get_power.__doc__ = _get_power.__doc__


def run_susie_rss(*args: Any, **kwargs: Any) -> Any:
    return _run_susie_rss_impl(*args, **kwargs)


run_susie_rss.__doc__ = _run_susie_rss_impl.__doc__


def compare_effect(path1: str, path2: str, **kwargs: Any) -> Any:
    params = _viz_params.merge("compare_effect", kwargs)
    params = _viz_params.filter(_compare_effect, params, key="compare_effect", log=Log(), verbose=params.get("verbose", True))
    return _compare_effect(path1, path2, **params)

def plot_forest(data: Union[pd.DataFrame, Sumstats], **kwargs: Any) -> Any:
    params = _viz_params.merge("plot_forest", kwargs)
    params = _viz_params.filter(_plot_forest, params, key="plot_forest", log=Log(), verbose=params.get("verbose", True))
    return _plot_forest(data, **params)

def plot_miami2(path1: Optional[str] = None, path2: Optional[str] = None, merged_sumstats: Optional[Union[pd.DataFrame, Sumstats]] = None, **kwargs: Any) -> Any:
    params = _viz_params.merge("plot_miami2", kwargs)
    params = _viz_params.filter(_plot_miami2, params, key="plot_miami2", log=Log(), verbose=params.get("verbose", True))
    return _plot_miami2(path1=path1, path2=path2, merged_sumstats=merged_sumstats, **params)

def plot_lead_overlap(objects: Any, **kwargs: Any) -> Any:
    params = _viz_params.merge("plot_lead_overlap", kwargs)
    params = _viz_params.filter(_plot_lead_overlap, params, key="plot_lead_overlap", log=Log(), verbose=params.get("verbose", True))
    return _plot_lead_overlap(objects=objects, **params)

def plot_rg(ldscrg: Union[str, pd.DataFrame], **kwargs: Any) -> Any:
    params = _viz_params.merge("plot_rg", kwargs)
    params = _viz_params.filter(_plot_rg, params, key="plot_rg", log=Log(), verbose=params.get("verbose", True))
    return _plot_rg(ldscrg, **params)

def plot_stacked_mqq(objects: Any, **kwargs: Any) -> Any:
    mode = kwargs.get("mode", "r")
    reg_mode = mode if mode in ("r", "m", "mqq") else None
    params = _viz_params.merge("plot_stacked_mqq", kwargs, mode=reg_mode)
    params = _viz_params.filter(
        _plot_stacked_mqq,
        params,
        key="plot_stacked_mqq",
        mode=reg_mode,
        log=Log(),
        verbose=params.get("verbose", True),
    )
    return _plot_stacked_mqq(objects, **params)

def plot_power(**kwargs: Any) -> Any:
    mode = "b" if kwargs.get("mode") == "b" else "q"
    params = _viz_params.merge("plot_power", kwargs, mode=mode)
    params = _viz_params.filter(_plot_power, params, key="plot_power", mode=mode, log=Log(), verbose=params.get("verbose", True))
    return _plot_power(**params)

def plot_power_x(**kwargs: Any) -> Any:
    params = _viz_params.merge("plot_power_x", kwargs)
    params = _viz_params.filter(_plot_power_x, params, key="plot_power_x", log=Log(), verbose=params.get("verbose", True))
    return _plot_power_x(**params)

def scatter(df: Union[pd.DataFrame, Sumstats], x: str, y: str, **kwargs: Any) -> Any:
    params = _viz_params.merge("plot_scatter", kwargs)
    params = _viz_params.filter(_scatter, params, key="plot_scatter", log=Log(), verbose=params.get("verbose", True))
    return _scatter(df=df, x=x, y=y, **params)

def plot_ld_block(**kwargs: Any) -> Any:
    params = _viz_params.merge("plot_ld_block", kwargs)
    params = _viz_params.filter(_plot_ld_block, params, key="plot_ld_block", log=Log(), verbose=params.get("verbose", True))
    return _plot_ld_block(**params)

def plot_panels(panels: Any, **kwargs: Any) -> Any:
    params = _viz_params.merge("plot_panels", kwargs)
    params = _viz_params.filter(_plot_panels, params, key="plot_panels", log=Log(), verbose=params.get("verbose", True))
    return _plot_panels(panels, **params)

def plot_sankey(data: Union[pd.DataFrame, Sumstats], **kwargs: Any) -> Any:
    if "columns" not in kwargs:
        raise TypeError("plot_sankey() missing required argument: columns")
    columns = kwargs.pop("columns")
    params = _viz_params.merge("plot_sankey", kwargs)
    params = _viz_params.filter(_plot_sankey, params, key="plot_sankey", log=Log(), verbose=params.get("verbose", True))
    return _plot_sankey(data, columns, **params)


from gwaslab.viz.viz_aux_doc import build_gl_viz_docstring

_GL_VIZ_DOCS = (
    ("compare_effect", _compare_effect, None, "Compare effect sizes between two GWAS summary-statistics files.", (
        ("path1", "str", "Path to the first summary-statistics file."),
        ("path2", "str", "Path to the second summary-statistics file."),
    )),
    ("plot_forest", _plot_forest, None, "Forest plot for meta-analysis study effects.", (
        ("data", "pandas.DataFrame or Sumstats", "Study-level effect summary table."),
    )),
    ("plot_miami2", _plot_miami2, None, "Mirrored Manhattan plot comparing two traits or studies.", (
        ("path1", "str, optional", "Path to the first summary-statistics file."),
        ("path2", "str, optional", "Path to the second summary-statistics file."),
        ("merged_sumstats", "pandas.DataFrame or Sumstats, optional", "Pre-merged sumstats instead of two files."),
    )),
    ("plot_lead_overlap", _plot_lead_overlap, None, "Venn/UpSet overlap of lead loci across studies.", (
        ("objects", "list", "Sumstats objects or lead-variant tables to compare."),
    )),
    ("plot_rg", _plot_rg, None, "Genetic correlation heatmap from LDSC results.", (
        ("ldscrg", "str or pandas.DataFrame", "LDSC genetic-correlation log or parsed results table."),
    )),
    ("plot_stacked_mqq", _plot_stacked_mqq, None, "Stacked Manhattan/QQ/regional panels from multiple Sumstats objects.", (
        ("objects", "list", "Sumstats objects to stack in one figure."),
    )),
    ("plot_power", _plot_power, None, "Theoretical GWAS power curves (mode='q' quantitative or 'b' binary).", ()),
    ("plot_power_x", _plot_power_x, None, "Extended power curves with custom MAF/beta grids.", ()),
    ("plot_panels", _plot_panels, None, "Stack multi-panel figure from Panel objects.", (
        ("panels", "list", "Panel layout objects from :class:`gwaslab.Panel`."),
    )),
    ("plot_ld_block", _plot_ld_block, None, "LD block as a 45-degree rotated inverted triangle from LD matrix or VCF.", ()),
    ("plot_sankey", _plot_sankey, None, "Sankey / alluvial diagram from categorical sumstats columns.", (
        ("data", "pandas.DataFrame or Sumstats", "Input sumstats with categorical columns."),
        ("columns", "list of str", "Column names defining Sankey stages (required)."),
    )),
)
for _key, _impl, _mode, _summary, _pos in _GL_VIZ_DOCS:
    _fn = globals()[_key]
    _fn.__doc__ = build_gl_viz_docstring(_key, _impl, _mode, summary=_summary, positional_params=_pos or None)

scatter.__doc__ = build_gl_viz_docstring(
    "plot_scatter",
    _scatter,
    None,
    summary="Scatter plot comparing two columns from sumstats.",
    positional_params=(
        ("df", "pandas.DataFrame or Sumstats", "Input data."),
        ("x", "str", "Column name for the x-axis."),
        ("y", "str", "Column name for the y-axis."),
    ),
)

