import pandas as pd
from typing import TYPE_CHECKING, Optional, Dict, Any, Callable

if TYPE_CHECKING:
    from gwaslab.info.g_Log import Log

from gwaslab.info.g_version import _show_version
from gwaslab.info.g_version import gwaslab_info
from gwaslab.info.g_meta import _init_meta
from gwaslab.qc.qc_fix_sumstats import _process_build
from gwaslab.util.util_in_merge import _extract_variant
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config
from gwaslab.viz.viz_plot_effect import _plot_effect
from gwaslab.info.g_Log import Log
from gwaslab.g_Sumstats import Sumstats


#20250215
class SumstatsSet(Sumstats):
    def __init__(
        self,
        sumstats_dic: Dict[str, Any],
        variant_set: Optional[Any] = None,
        build: str = "99",
        species: str = "homo sapiens",
        build_infer: bool = False,
        set: str = "set1",
        verbose: bool = True,
        **readargs: Any
    ) -> None:

        # basic attributes
        self.data = pd.DataFrame()
        self.log = Log()
        # meta information

        self.meta = _init_meta() 
        self.meta["gwaslab"]["genome_build"] = _process_build(build, log=self.log, verbose=False, species=species)
        self._build = self.meta["gwaslab"]["genome_build"]
        self.meta["gwaslab"]["set_name"] =  set
        self.meta["gwaslab"]["species"] = species

        # print gwaslab version information
        _show_version(self.log,  verbose=verbose)

        self.data = _extract_variant(variant_set, sumstats_dic,log=self.log, verbose=verbose)
        
        self.viz_params = VizParamsManager()
        load_viz_config(self.viz_params)
 
    def _apply_viz_params(self, func: Callable[..., Any], kwargs: Dict[str, Any], key: Optional[str] = None, mode: Optional[str] = None) -> Dict[str, Any]:
        params = self.viz_params.merge(key or func.__name__, kwargs, mode=mode)
        return self.viz_params.filter(func, params, key=key or func.__name__, mode=mode, log=self.log, verbose=kwargs.get("verbose", True))

    def plot_effect(self, **kwargs: Any) -> None:
        _plot_effect(self.data,**self._apply_viz_params(_plot_effect, kwargs, key="plot_effect"))


 
