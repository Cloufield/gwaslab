import pandas as pd
import glob
import os
from typing import TYPE_CHECKING, Optional, Dict, Any, Callable, Union

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
    """
    A class for working with multiple GWAS summary statistics datasets.

    SumstatsSet allows you to load and combine variants from multiple studies,
    either by providing a dictionary of Sumstats objects or by using a glob pattern
    to automatically discover and load multiple files.

    Parameters
    ----------
    sumstats_dic : Union[Dict[str, Any], str]
        Either:
        - A dictionary mapping study names to Sumstats objects or file paths
        - A glob pattern string (e.g., "./study_*.txt") to auto-discover files
    variant_set : Optional[Any], default=None
        Collection of variants to extract. See `_extract_variant` for details.
    build : str, default="99"
        Genome build version.
    species : str, default="homo sapiens"
        Species name.
    build_infer : bool, default=False
        Whether to infer genome build.
    set : str, default="set1"
        Name for this set of summary statistics.
    verbose : bool, default=True
        Whether to print log messages.
    **readargs : Any
        Additional arguments passed to Sumstats initialization when loading from files.

    Examples
    --------
    >>> # Load from a dictionary of Sumstats objects
    >>> ss = SumstatsSet({"study1": sumstats1, "study2": sumstats2}, variant_set=["rs12345"])

    >>> # Load from a glob pattern
    >>> ss = SumstatsSet("./data/study_*.txt", variant_set=["rs12345"], fmt="auto")
    """

    def __init__(
        self,
        sumstats_dic: Union[Dict[str, Any], str],
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

        # Handle glob pattern input
        if isinstance(sumstats_dic, str):
            sumstats_dic = self._load_from_glob_pattern(
                sumstats_dic, 
                build=build, 
                species=species, 
                verbose=verbose, 
                **readargs
            )

        # If variant_set is None, concatenate all data; otherwise extract specific variants
        if variant_set is None:
            self.data = self._concat_all_variants(sumstats_dic, log=self.log, verbose=verbose)
        else:
            self.data = _extract_variant(variant_set, sumstats_dic, log=self.log, verbose=verbose)
        
        self.viz_params = VizParamsManager()
        load_viz_config(self.viz_params)

    def _concat_all_variants(
        self,
        sumstats_dic: Dict[str, 'Sumstats'],
        log: 'Log' = None,
        verbose: bool = True
    ) -> pd.DataFrame:
        """
        Concatenate all variants from multiple Sumstats objects.

        Used when variant_set is None to load all data from all studies.

        Parameters
        ----------
        sumstats_dic : Dict[str, Sumstats]
            Dictionary mapping study names to Sumstats objects.
        log : Log, optional
            Log object for messages.
        verbose : bool, default=True
            Whether to print log messages.

        Returns
        -------
        pd.DataFrame
            Combined DataFrame with all variants from all studies,
            with a "STUDY" column indicating the source.
        """
        if log is None:
            log = Log()

        combined = pd.DataFrame()
        
        for key, sumstats_gls in sumstats_dic.items():
            log.write(f" -{key} : {sumstats_gls}", verbose=verbose)

        for key, sumstats_gls in sumstats_dic.items():
            sumstats_single = sumstats_gls.data.copy()
            sumstats_single["STUDY"] = key
            
            log.write(f" -Loaded {len(sumstats_single)} variants from {key}", verbose=verbose)
            
            # Reorder columns to put STUDY first
            cols = ["STUDY"] + [c for c in sumstats_single.columns if c != "STUDY"]
            sumstats_single = sumstats_single[cols]
            
            combined = pd.concat([combined, sumstats_single], ignore_index=True)
        
        log.write(f" -Total: {len(combined)} variants from {len(sumstats_dic)} studies", verbose=verbose)
        return combined

    def _load_from_glob_pattern(
        self,
        pattern: str,
        build: str = "99",
        species: str = "homo sapiens",
        verbose: bool = True,
        **readargs: Any
    ) -> Dict[str, 'Sumstats']:
        """
        Load multiple Sumstats objects from files matching a glob pattern.

        Parameters
        ----------
        pattern : str
            A glob pattern to match files (e.g., "./study_*.txt", "./data/trait_?.tsv").
            Supports standard glob wildcards: * (any characters), ? (single character).
        build : str, default="99"
            Genome build version passed to each Sumstats object.
        species : str, default="homo sapiens"
            Species name passed to each Sumstats object.
        verbose : bool, default=True
            Whether to print log messages.
        **readargs : Any
            Additional arguments passed to Sumstats initialization.

        Returns
        -------
        Dict[str, Sumstats]
            A dictionary mapping derived study names to Sumstats objects.
            Study names are derived from filenames by removing directory path
            and common extensions (.txt, .tsv, .gz, etc.).

        Raises
        ------
        FileNotFoundError
            If no files match the given pattern.

        Notes
        -----
        The study name is derived from the filename using the following logic:
        1. Extract the basename (remove directory path)
        2. Remove common extensions: .gz, .txt, .tsv, .csv, .sumstats
        
        For example:
        - "./data/study_EUR.txt" -> "study_EUR"
        - "./trait_1.sumstats.gz" -> "trait_1"
        """
        # Expand glob pattern
        matched_files = sorted(glob.glob(pattern))
        
        if not matched_files:
            raise FileNotFoundError(f"No files match pattern: {pattern}")
        
        self.log.write(f" -Detected glob pattern: {pattern}", verbose=verbose)
        self.log.write(f" -Found {len(matched_files)} matching file(s)", verbose=verbose)
        
        sumstats_dict: Dict[str, 'Sumstats'] = {}
        
        for filepath in matched_files:
            # Derive study name from filename
            study_name = self._derive_study_name(filepath)
            
            self.log.write(f" -Loading: {filepath} as '{study_name}'", verbose=verbose)
            
            # Create Sumstats object for this file
            sumstats_obj = Sumstats(
                filepath,
                build=build,
                species=species,
                verbose=verbose,
                **readargs
            )
            
            sumstats_dict[study_name] = sumstats_obj
        
        self.log.write(f" -Successfully loaded {len(sumstats_dict)} Sumstats object(s)", verbose=verbose)
        
        return sumstats_dict

    @staticmethod
    def _derive_study_name(filepath: str) -> str:
        """
        Derive a study name from a file path.

        Parameters
        ----------
        filepath : str
            Full path to the file.

        Returns
        -------
        str
            Derived study name based on the filename.

        Examples
        --------
        >>> SumstatsSet._derive_study_name("./data/study_EUR.txt")
        'study_EUR'
        >>> SumstatsSet._derive_study_name("./trait_1.sumstats.gz")
        'trait_1'
        """
        # Get basename
        name = os.path.basename(filepath)
        
        # Remove common extensions (order matters - .gz first, then others)
        extensions_to_remove = ['.gz', '.bgz', '.zst', '.sumstats', '.txt', '.tsv', '.csv']
        for ext in extensions_to_remove:
            if name.lower().endswith(ext):
                name = name[:-len(ext)]
        
        return name

    def _apply_viz_params(self, func: Callable[..., Any], kwargs: Dict[str, Any], key: Optional[str] = None, mode: Optional[str] = None) -> Dict[str, Any]:
        params = self.viz_params.merge(key or func.__name__, kwargs, mode=mode)
        return self.viz_params.filter(func, params, key=key or func.__name__, mode=mode, log=self.log, verbose=kwargs.get("verbose", True))

    def plot_effect(self, **kwargs: Any) -> None:
        _plot_effect(self.data,**self._apply_viz_params(_plot_effect, kwargs, key="plot_effect"))


 
