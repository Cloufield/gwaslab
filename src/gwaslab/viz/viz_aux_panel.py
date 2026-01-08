"""
Panel class for storing panel information for stacked plots.

This module provides a Panel class that stores configuration information
for individual panels in a stacked figure, which can be used with plot_panels
to create multi-panel visualizations.
"""

from typing import Dict, Any


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
    **kwargs
        All keyword arguments specific to the panel type. These will be
        stored and passed to the appropriate plotting function when
        plot_panels is called.
    
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
    >>> # Plot panels together
    >>> gl.plot_panels([panel1, panel2, panel3])
    """
    
    def __init__(self, panel_type: str, **kwargs):
        """
        Initialize a Panel object.
        
        Parameters
        ----------
        panel_type : str
            Type of panel ("track", "arc", etc.)
        **kwargs
            Panel-specific parameters
        """
        if not isinstance(panel_type, str):
            raise TypeError(f"panel_type must be a string, got {type(panel_type)}")
        
        self.panel_type = panel_type.lower()
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
