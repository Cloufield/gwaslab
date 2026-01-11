"""
Panel class for storing panel information for stacked plots.

This module provides a Panel class that stores configuration information
for individual panels in a stacked figure, which can be used with plot_panels
to create multi-panel visualizations.
"""

from typing import Dict, Any, Optional
from gwaslab.viz.viz_aux_params import VizParamsManager, load_viz_config
from gwaslab.viz.viz_plot_track import plot_track
from gwaslab.viz.viz_plot_arc import plot_arc
from gwaslab.viz.viz_plot_ld_block import plot_ld_block
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.viz.viz_aux_chromatin import _plot_chromatin_state
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.info.g_Log import Log


class Panel:
    """
    A class to store panel configuration information for stacked plots.
    
    This class stores the panel type and all parameters needed to create
    a panel in a stacked figure. Panels can be of different types such as
    "track" (for genomic tracks) or "arc" (for BEDPE arc plots).
    
    Parameters
    ----------
    panel_type : str
        Type of panel to create. Supported types:
        - "track": For genomic tracks (GTF, BED, bigWig, bigBed files)
        - "arc": For BEDPE arc plots
        - "ld_block": For LD block plots (45° rotated inverted triangle)
        - "region": For regional plots with LD information and gene track (requires 2 axes)
        - "chromatin": For chromatin state tracks (Roadmap 15-state model)
        - "pipcs": For PIP (Posterior Inclusion Probability) and Credible Sets plots
    **kwargs
        All keyword arguments specific to the panel type. These will be
        processed through the parameter manager (similar to Sumstats.plot_* methods)
        to merge defaults, filter allowed parameters, and validate arguments before
        being stored. The processed parameters will be passed to the appropriate
        plotting function when plot_panels is called.
    
    Attributes
    ----------
    panel_type : str
        The type of panel (e.g., "track", "arc")
    kwargs : dict
        Dictionary storing all panel-specific parameters
    
    Examples
    --------
    >>> import gwaslab as gl
    >>> 
    >>> # Create a track panel
    >>> panel1 = gl.Panel(
    ...     "track",
    ...     track_path="genes.gtf",
    ...     region=(1, 1000000, 2000000),
    ...     color="#020080"
    ... )
    >>> 
    >>> # Create an arc panel
    >>> panel2 = gl.Panel(
    ...     "arc",
    ...     bedpe_path="contacts.bedpe.gz",
    ...     region=(1, 1000000, 2000000),
    ...     color="#FF0000",
    ...     alpha=0.3
    ... )
    >>> 
    >>> # Create an LD block panel
    >>> panel3 = gl.Panel(
    ...     "ld_block",
    ...     vcf_path="ld_ref.vcf.gz",
    ...     region=(1, 1000000, 2000000),
    ...     sumstats=sumstats
    ... )
    >>> 
    >>> # Create a region panel
    >>> panel4 = gl.Panel(
    ...     "region",
    ...     sumstats=sumstats,
    ...     region=(1, 1000000, 2000000),
    ...     vcf_path="ld_ref.vcf.gz",
    ...     build="38"
    ... )
    >>> 
    >>> # Create a chromatin panel
    >>> panel5 = gl.Panel(
    ...     "chromatin",
    ...     region_chromatin_files=["E098_15_coreMarks_mnemonics.bed.gz"],
    ...     region_chromatin_labels=["E098"],
    ...     region=(1, 1000000, 2000000)
    ... )
    >>> 
    >>> # Create a pipcs panel
    >>> panel6 = gl.Panel(
    ...     "pipcs",
    ...     pipcs_raw=pipcs_dataframe,
    ...     region=(1, 1000000, 2000000)
    ... )
    >>> 
    >>> # Plot panels together
    >>> gl.plot_panels([panel1, panel2, panel3, panel4, panel5, panel6])
    """
    
    def __init__(self, panel_type: str, verbose: bool = True, log: Optional[Log] = None, **kwargs):
        """
        Initialize a Panel object.
        
        Parameters
        ----------
        panel_type : str
            Type of panel ("track", "arc", etc.)
        verbose : bool, default=True
            Whether to show verbose messages during parameter processing
        log : Log, optional
            Logger instance. If None, creates a new Log instance.
        **kwargs
            Panel-specific parameters. These will be processed through the
            parameter manager before being stored.
        """
        if not isinstance(panel_type, str):
            raise TypeError(f"panel_type must be a string, got {type(panel_type)}")
        
        self.panel_type = panel_type.lower()
        
        # Apply parameter manager to kwargs
        if log is None:
            log = Log()
        
        # Map panel types to parameter manager keys, plotting functions, and modes
        param_config = {
            "track": ("plot_track", plot_track, None),
            "arc": ("plot_arc", plot_arc, None),
            "ld_block": ("plot_ld_block", plot_ld_block, None),
            "region": ("plot_region", _mqqplot, "r"),  # Use _mqqplot with mode="r" for regional plots
            "pipcs": ("plot_pipcs", _plot_cs, None)
        }
        
        if self.panel_type in param_config:
            param_key, plot_func, mode = param_config[self.panel_type]
            
            # Preserve data parameters that might be filtered out (e.g., sumstats, track_path, bedpe_path)
            # These are typically not in function signatures but are required for plotting
            data_params = {}
            data_param_keys = ["insumstats", "sumstats", "track_path", "bedpe_path", "region_chromatin_files", "region_chromatin_labels", "pipcs_raw"]
            for key in data_param_keys:
                if key in kwargs:
                    data_params[key] = kwargs[key]
            
            # Create parameter manager instance and load config
            viz_params = VizParamsManager()
            load_viz_config(viz_params)
            
            # Merge and filter parameters
            params = viz_params.merge(param_key, kwargs, mode=mode)
            self.kwargs = viz_params.filter(plot_func, params, key=param_key, mode=mode, log=log, verbose=verbose)
            
            # Add back data parameters that were filtered out
            self.kwargs.update(data_params)
        else:
            # If panel type doesn't have parameter management, store kwargs as-is
            self.kwargs = kwargs.copy()
    
    def __repr__(self) -> str:
        """String representation of the Panel."""
        kwargs_str = ", ".join(f"{k}={v!r}" for k, v in list(self.kwargs.items())[:3])
        if len(self.kwargs) > 3:
            kwargs_str += f", ... ({len(self.kwargs) - 3} more)"
        return f"Panel(panel_type='{self.panel_type}', {kwargs_str})"
    
    def __str__(self) -> str:
        """Human-readable string representation."""
        return f"Panel(type='{self.panel_type}', params={len(self.kwargs)} kwargs)"
    
    def get_type(self) -> str:
        """
        Get the panel type.
        
        Returns
        -------
        str
            The panel type
        """
        return self.panel_type
    
    def get_kwargs(self) -> Dict[str, Any]:
        """
        Get all stored keyword arguments.
        
        Returns
        -------
        dict
            Dictionary of all panel parameters
        """
        return self.kwargs.copy()
    
    def update_kwargs(self, **kwargs):
        """
        Update stored keyword arguments.
        
        Parameters
        ----------
        **kwargs
            Keyword arguments to update or add
        """
        self.kwargs.update(kwargs)
    
    def set_kwarg(self, key: str, value: Any):
        """
        Set a single keyword argument.
        
        Parameters
        ----------
        key : str
            Parameter name
        value : Any
            Parameter value
        """
        self.kwargs[key] = value
    
    def get_kwarg(self, key: str, default: Any = None) -> Any:
        """
        Get a single keyword argument.
        
        Parameters
        ----------
        key : str
            Parameter name
        default : Any, optional
            Default value if key is not found
        
        Returns
        -------
        Any
            Parameter value or default
        """
        return self.kwargs.get(key, default)
