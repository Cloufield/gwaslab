from typing import TYPE_CHECKING, Optional, Dict, Any, Union
import pandas as pd

from gwaslab.g_Sumstats import Sumstats
from gwaslab.g_Sumstats_polars import Sumstatsp
from gwaslab.g_SumstatsT import SumstatsT
from gwaslab.g_SumstatsPair import SumstatsPair
from gwaslab.g_SumstatsMulti import SumstatsMulti
from gwaslab.g_SumstatsSet import SumstatsSet
from gwaslab.info.g_version import _show_version as show_version

from gwaslab.util.util_in_convert_h2 import h2_obs_to_liab
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2
from gwaslab.util.util_in_convert_h2 import h2_se_to_p
from gwaslab.util.util_in_calculate_power import get_power
from gwaslab.util.util_in_calculate_power import get_beta
from gwaslab.util.util_ex_gwascatalog import gwascatalog_trait
from gwaslab.util.util_ex_process_h5 import process_vcf_to_hfd5
from gwaslab.util.util_ex_run_susie import _run_susie_rss as run_susie_rss
from gwaslab.util.util_in_meta import meta_analyze
from gwaslab.util.util_in_fill_data import rank_based_int

from gwaslab.io.io_read_ldsc import read_ldsc
from gwaslab.io.io_read_ldsc import read_popcorn
from gwaslab.io.io_to_pickle import dump_pickle
from gwaslab.io.io_to_pickle import load_pickle
from gwaslab.io.io_gwaslab_standard import load_gsf
from gwaslab.io.io_read_tabular import _read_tabular as read_tabular

from gwaslab.viz.viz_plot_compare_effect import compare_effect as _compare_effect
from gwaslab.viz.viz_plot_forestplot import plot_forest as _plot_forest
from gwaslab.viz.viz_plot_miamiplot2 import plot_miami2 as _plot_miami2
from gwaslab.viz.viz_plot_rg_heatmap import plot_rg as _plot_rg
from gwaslab.viz.viz_plot_stackedregional import plot_stacked_mqq as _plot_stacked_mqq
from gwaslab.viz.viz_plot_trumpetplot import plot_power as _plot_power
from gwaslab.viz.viz_plot_trumpetplot import plot_power_x as _plot_power_x
from gwaslab.viz.viz_plot_scatter_with_reg import scatter as _scatter

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
from gwaslab.bd.bd_download import check_format
from gwaslab.bd.bd_download import check_available_ref
from gwaslab.bd.bd_download import update_available_ref
from gwaslab.bd.bd_download import check_downloaded_ref
from gwaslab.bd.bd_download import download_ref
from gwaslab.bd.bd_download import check_available_ref
from gwaslab.bd.bd_download import remove_file
from gwaslab.bd.bd_download import get_path
from gwaslab.bd.bd_download import update_record
from gwaslab.bd.bd_download import scan_downloaded_files
from gwaslab.bd.bd_download import add_local_data
from gwaslab.bd.bd_download import remove_local_record
from gwaslab.bd.bd_config import options
from gwaslab.qc.qc_reserved_headers import researved_header
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config
from gwaslab.view.view_report import generate_qc_report


_viz_params = VizParamsManager()
load_viz_config(_viz_params)

def compare_effect(path1: str, path2: str, **kwargs: Any) -> Any:
    """Compare effect sizes between two sumstats files."""
    params = _viz_params.merge("compare_effect", kwargs)
    params = _viz_params.filter(_compare_effect, params, key="compare_effect", log=Log(), verbose=params.get("verbose", True))
    return _compare_effect(path1, path2, **params)

def plot_forest(data: Union[pd.DataFrame, Sumstats], **kwargs: Any) -> Any:
    """Plot forest plot from sumstats data."""
    params = _viz_params.merge("plot_forest", kwargs)
    params = _viz_params.filter(_plot_forest, params, key="plot_forest", log=Log(), verbose=params.get("verbose", True))
    return _plot_forest(data, **params)

def plot_miami2(path1: Optional[str] = None, path2: Optional[str] = None, merged_sumstats: Optional[Union[pd.DataFrame, Sumstats]] = None, **kwargs: Any) -> Any:
    """Plot Miami plot (mirrored Manhattan plot) from sumstats."""
    params = _viz_params.merge("plot_miami2", kwargs)
    params = _viz_params.filter(_plot_miami2, params, key="plot_miami2", log=Log(), verbose=params.get("verbose", True))
    return _plot_miami2(path1=path1, path2=path2, merged_sumstats=merged_sumstats, **params)

def plot_rg(ldscrg: Union[str, pd.DataFrame], **kwargs: Any) -> Any:
    """Plot genetic correlation heatmap from LDSC results."""
    params = _viz_params.merge("plot_rg", kwargs)
    params = _viz_params.filter(_plot_rg, params, key="plot_rg", log=Log(), verbose=params.get("verbose", True))
    return _plot_rg(ldscrg, **params)

def plot_stacked_mqq(objects: Any, **kwargs: Any) -> Any:
    """Plot stacked Manhattan and QQ plots from multiple sumstats objects."""
    params = _viz_params.merge("plot_stacked_mqq", kwargs)
    params = _viz_params.filter(_plot_stacked_mqq, params, key="plot_stacked_mqq", log=Log(), verbose=params.get("verbose", True))
    return _plot_stacked_mqq(objects, **params)

def plot_power(**kwargs: Any) -> Any:
    """Plot power analysis (trumpet plot)."""
    params = _viz_params.merge("plot_power", kwargs)
    params = _viz_params.filter(_plot_power, params, key="plot_power", log=Log(), verbose=params.get("verbose", True))
    return _plot_power(**params)

def plot_power_x(**kwargs: Any) -> Any:
    """Plot power analysis with extended features."""
    params = _viz_params.merge("plot_power_x", kwargs)
    params = _viz_params.filter(_plot_power_x, params, key="plot_power_x", log=Log(), verbose=params.get("verbose", True))
    return _plot_power_x(**params)

def scatter(df: Union[pd.DataFrame, Sumstats], x: str, y: str, **kwargs: Any) -> Any:
    """Create scatter plot comparing two columns from sumstats."""
    params = _viz_params.merge("plot_scatter", kwargs)
    params = _viz_params.filter(_scatter, params, key="plot_scatter", log=Log(), verbose=params.get("verbose", True))
    return _scatter(df=df, x=x, y=y, **params)
