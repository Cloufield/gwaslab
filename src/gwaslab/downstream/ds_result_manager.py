"""Downstream analysis result manager for GWASLab Sumstats objects.

This module provides a centralized manager for storing and accessing downstream
analysis results, including LDSC, finemapping, clumping, and other analysis outputs.

Key Features:
- Centralized storage of all downstream analysis results
- Property-based accessors for backward compatibility
- Methods to get, set, clear, and check result availability
- Support for LDSC, finemapping, clumping, and PIPCS results
"""

import pandas as pd
from typing import Optional, Dict, Any, Union


class DownstreamResultManager:
    """Manager for downstream analysis results in Sumstats objects.
    
    This class centralizes the storage and management of all downstream analysis
    results, providing a clean interface for accessing and manipulating results
    from various analyses like LDSC, finemapping, clumping, etc.
    
    Attributes
    ----------
    ldsc_h2 : Optional[Any]
        LDSC heritability estimate (scalar value).
    ldsc_h2_results : Optional[pd.DataFrame]
        Detailed LDSC heritability results.
    ldsc_rg : pd.DataFrame
        LDSC genetic correlation results.
    ldsc_h2_cts : Optional[Any]
        LDSC cell-type-specific heritability results.
    ldsc_partitioned_h2_summary : Optional[Any]
        LDSC partitioned heritability summary.
    ldsc_partitioned_h2_results : Optional[pd.DataFrame]
        Detailed LDSC partitioned heritability results.
    finemapping : Dict[str, Any]
        Finemapping results dictionary (keys: "path", "file", "plink_log", "pipcs").
    clumps : Dict[str, Any]
        Clumping results dictionary (keys: "clumps", "clumps_raw", "plink_log").
    pipcs : pd.DataFrame
        PIPCS (Posterior Inclusion Probabilities for Credible Sets) results.
    
    Examples
    --------
    >>> manager = DownstreamResultManager()
    >>> manager.set_ldsc_h2(0.5)
    >>> manager.get_ldsc_h2()
    0.5
    >>> manager.has_ldsc_h2()
    True
    >>> manager.clear_ldsc()
    """
    
    def __init__(self):
        """Initialize the downstream result manager with default empty values."""
        # LDSC results
        self._ldsc_h2 = None
        self._ldsc_h2_results = None
        self._ldsc_rg = pd.DataFrame()
        self._ldsc_h2_cts = None
        self._ldsc_partitioned_h2_summary = None
        self._ldsc_partitioned_h2_results = None
        
        # Finemapping results
        self._finemapping = {}
        
        # Clumping results
        self._clumps = {}
        
        # PIPCS results
        self._pipcs = pd.DataFrame()
    
    # ============================================================================
    # LDSC Properties and Methods
    # ============================================================================
    
    @property
    def ldsc_h2(self) -> Optional[Any]:
        """LDSC heritability estimate."""
        return self._ldsc_h2
    
    @ldsc_h2.setter
    def ldsc_h2(self, value: Optional[Any]):
        """Set LDSC heritability estimate."""
        self._ldsc_h2 = value
    
    @property
    def ldsc_h2_results(self) -> Optional[pd.DataFrame]:
        """Detailed LDSC heritability results."""
        return self._ldsc_h2_results
    
    @ldsc_h2_results.setter
    def ldsc_h2_results(self, value: Optional[pd.DataFrame]):
        """Set detailed LDSC heritability results."""
        self._ldsc_h2_results = value
    
    @property
    def ldsc_rg(self) -> pd.DataFrame:
        """LDSC genetic correlation results."""
        return self._ldsc_rg
    
    @ldsc_rg.setter
    def ldsc_rg(self, value: Union[pd.DataFrame, None]):
        """Set LDSC genetic correlation results."""
        if value is None:
            self._ldsc_rg = pd.DataFrame()
        else:
            self._ldsc_rg = value
    
    @property
    def ldsc_h2_cts(self) -> Optional[Any]:
        """LDSC cell-type-specific heritability results."""
        return self._ldsc_h2_cts
    
    @ldsc_h2_cts.setter
    def ldsc_h2_cts(self, value: Optional[Any]):
        """Set LDSC cell-type-specific heritability results."""
        self._ldsc_h2_cts = value
    
    @property
    def ldsc_partitioned_h2_summary(self) -> Optional[Any]:
        """LDSC partitioned heritability summary."""
        return self._ldsc_partitioned_h2_summary
    
    @ldsc_partitioned_h2_summary.setter
    def ldsc_partitioned_h2_summary(self, value: Optional[Any]):
        """Set LDSC partitioned heritability summary."""
        self._ldsc_partitioned_h2_summary = value
    
    @property
    def ldsc_partitioned_h2_results(self) -> Optional[pd.DataFrame]:
        """Detailed LDSC partitioned heritability results."""
        return self._ldsc_partitioned_h2_results
    
    @ldsc_partitioned_h2_results.setter
    def ldsc_partitioned_h2_results(self, value: Optional[pd.DataFrame]):
        """Set detailed LDSC partitioned heritability results."""
        self._ldsc_partitioned_h2_results = value
    
    def has_ldsc_h2(self) -> bool:
        """Check if LDSC heritability estimate exists."""
        return self._ldsc_h2 is not None
    
    def has_ldsc_h2_results(self) -> bool:
        """Check if detailed LDSC heritability results exist."""
        return self._ldsc_h2_results is not None and not self._ldsc_h2_results.empty
    
    def has_ldsc_rg(self) -> bool:
        """Check if LDSC genetic correlation results exist."""
        return not self._ldsc_rg.empty
    
    def has_ldsc_h2_cts(self) -> bool:
        """Check if LDSC cell-type-specific heritability results exist."""
        return self._ldsc_h2_cts is not None
    
    def has_ldsc_partitioned_h2(self) -> bool:
        """Check if LDSC partitioned heritability results exist."""
        return (self._ldsc_partitioned_h2_summary is not None or 
                (self._ldsc_partitioned_h2_results is not None and 
                 not self._ldsc_partitioned_h2_results.empty))
    
    def clear_ldsc(self):
        """Clear all LDSC results."""
        self._ldsc_h2 = None
        self._ldsc_h2_results = None
        self._ldsc_rg = pd.DataFrame()
        self._ldsc_h2_cts = None
        self._ldsc_partitioned_h2_summary = None
        self._ldsc_partitioned_h2_results = None
    
    # ============================================================================
    # Finemapping Properties and Methods
    # ============================================================================
    
    @property
    def finemapping(self) -> Dict[str, Any]:
        """Finemapping results dictionary."""
        return self._finemapping
    
    @finemapping.setter
    def finemapping(self, value: Dict[str, Any]):
        """Set finemapping results dictionary."""
        if value is None:
            self._finemapping = {}
        else:
            self._finemapping = value
    
    def get_finemapping(self, key: str, default: Any = None) -> Any:
        """Get a specific finemapping result by key.
        
        Parameters
        ----------
        key : str
            Key to retrieve (e.g., "path", "file", "plink_log", "pipcs").
        default : Any, optional
            Default value if key doesn't exist.
        
        Returns
        -------
        Any
            The value associated with the key, or default if not found.
        """
        return self._finemapping.get(key, default)
    
    def set_finemapping(self, key: str, value: Any):
        """Set a specific finemapping result.
        
        Parameters
        ----------
        key : str
            Key to set (e.g., "path", "file", "plink_log", "pipcs").
        value : Any
            Value to store.
        """
        self._finemapping[key] = value
    
    def has_finemapping(self, key: Optional[str] = None) -> bool:
        """Check if finemapping results exist.
        
        Parameters
        ----------
        key : str, optional
            Specific key to check. If None, checks if any finemapping data exists.
        
        Returns
        -------
        bool
            True if finemapping results exist (for the specified key or any key).
        """
        if key is None:
            return len(self._finemapping) > 0
        return key in self._finemapping
    
    def clear_finemapping(self):
        """Clear all finemapping results."""
        self._finemapping = {}
    
    # ============================================================================
    # Clumping Properties and Methods
    # ============================================================================
    
    @property
    def clumps(self) -> Dict[str, Any]:
        """Clumping results dictionary."""
        return self._clumps
    
    @clumps.setter
    def clumps(self, value: Dict[str, Any]):
        """Set clumping results dictionary."""
        if value is None:
            self._clumps = {}
        else:
            self._clumps = value
    
    def get_clumps(self, key: str, default: Any = None) -> Any:
        """Get a specific clumping result by key.
        
        Parameters
        ----------
        key : str
            Key to retrieve (e.g., "clumps", "clumps_raw", "plink_log").
        default : Any, optional
            Default value if key doesn't exist.
        
        Returns
        -------
        Any
            The value associated with the key, or default if not found.
        """
        return self._clumps.get(key, default)
    
    def set_clumps(self, key: str, value: Any):
        """Set a specific clumping result.
        
        Parameters
        ----------
        key : str
            Key to set (e.g., "clumps", "clumps_raw", "plink_log").
        value : Any
            Value to store.
        """
        self._clumps[key] = value
    
    def has_clumps(self, key: Optional[str] = None) -> bool:
        """Check if clumping results exist.
        
        Parameters
        ----------
        key : str, optional
            Specific key to check. If None, checks if any clumping data exists.
        
        Returns
        -------
        bool
            True if clumping results exist (for the specified key or any key).
        """
        if key is None:
            return len(self._clumps) > 0
        return key in self._clumps
    
    def clear_clumps(self):
        """Clear all clumping results."""
        self._clumps = {}
    
    # ============================================================================
    # PIPCS Properties and Methods
    # ============================================================================
    
    @property
    def pipcs(self) -> pd.DataFrame:
        """PIPCS (Posterior Inclusion Probabilities for Credible Sets) results."""
        return self._pipcs
    
    @pipcs.setter
    def pipcs(self, value: Union[pd.DataFrame, None]):
        """Set PIPCS results."""
        if value is None:
            self._pipcs = pd.DataFrame()
        else:
            self._pipcs = value
    
    def has_pipcs(self) -> bool:
        """Check if PIPCS results exist."""
        return not self._pipcs.empty
    
    def clear_pipcs(self):
        """Clear PIPCS results."""
        self._pipcs = pd.DataFrame()
    
    # ============================================================================
    # General Management Methods
    # ============================================================================
    
    def clear_all(self):
        """Clear all downstream analysis results."""
        self.clear_ldsc()
        self.clear_finemapping()
        self.clear_clumps()
        self.clear_pipcs()
    
    def get_summary(self) -> Dict[str, Any]:
        """Get a summary of available downstream analysis results.
        
        Returns
        -------
        Dict[str, Any]
            Dictionary with keys indicating which results are available.
        """
        return {
            "ldsc_h2": self.has_ldsc_h2(),
            "ldsc_h2_results": self.has_ldsc_h2_results(),
            "ldsc_rg": self.has_ldsc_rg(),
            "ldsc_h2_cts": self.has_ldsc_h2_cts(),
            "ldsc_partitioned_h2": self.has_ldsc_partitioned_h2(),
            "finemapping": self.has_finemapping(),
            "clumps": self.has_clumps(),
            "pipcs": self.has_pipcs(),
        }
    
    def __repr__(self) -> str:
        """String representation of the manager."""
        summary = self.get_summary()
        available = [k for k, v in summary.items() if v]
        if available:
            return f"DownstreamResultManager(available: {', '.join(available)})"
        else:
            return "DownstreamResultManager(no results)"


