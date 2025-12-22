"""
Forest plot visualization for meta-analysis results.

This module provides functionality to create forest plots showing
individual study results and combined meta-analysis estimates.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.axes import Axes
from typing import Union, Optional, Dict, List, Tuple
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_aux_save_figure import save_figure

try:
    from statsmodels.stats.meta_analysis import effectsize_smd, combine_effects
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    from scipy import stats


def plot_forest(
    data: Union[pd.DataFrame, str],
    study_col: str,
    group_col: Optional[Union[str, bool]] = None,
    beta_col: str = "beta",
    se_col: str = "se",
    compact_factor: float = 1.0,
    width_ratios: List[float] = [2, 6, 2],
    sharex: str = "col",
    meta: bool = True,
    combine_effects_kwargs: Optional[Dict] = None,
    fig_kwargs: Optional[Dict] = None,
    save: Union[bool, str, None] = None,
    save_kwargs: Optional[Dict] = None,
    fontsize: int = 12,
    font_family: str = "Arial",
    colors: Optional[List[str]] = None,
    verbose: bool = True,
    log: Log = Log(),
) -> Tuple[plt.Figure, List]:
    """
    Create a forest plot for meta-analysis results.
    
    A forest plot displays individual study effect estimates with confidence
    intervals, along with a combined meta-analysis estimate (fixed or random effects).
    
    Parameters
    ----------
    data : pd.DataFrame or str
        DataFrame containing study data, or path to a whitespace-separated file.
        Required columns: study_col, beta_col, se_col
    study_col : str
        Column name containing study identifiers
    group_col : str, bool, or None, optional
        Column name for grouping studies. If False or None, all studies are
        placed in a single group "Group1". Default: None
    beta_col : str, optional
        Column name for effect estimates (beta). Default: "beta"
    se_col : str, optional
        Column name for standard errors. Default: "se"
    compact_factor : float, optional
        Factor to adjust figure height (higher = more compact). Default: 1.0
    width_ratios : list of float, optional
        Width ratios for the three subplot columns [group, plot, text].
        Default: [2, 6, 2]
    sharex : str, optional
        Whether to share x-axis across rows. Default: "col"
    meta : bool, optional
        Whether to perform and display meta-analysis statistics. Default: True
    combine_effects_kwargs : dict, optional
        Additional arguments passed to combine_effects function.
        Default: None
    fig_kwargs : dict, optional
        Additional arguments passed to plt.subplots(). Default: None
    save : str, bool, or None, optional
        If str: file path to save figure
        If True: save to default path
        If None/False: don't save. Default: None
    save_kwargs : dict, optional
        Additional arguments for saving figure. Default: None
    fontsize : int, optional
        Font size for labels and text. Default: 12
    font_family : str, optional
        Font family for text. Default: "Arial"
    colors : list of str, optional
        List of colors for alternating studies. If None, uses default
        colors from plot_mqq: ["#597FBD", "#74BAD3"]. Default: None
    verbose : bool, optional
        Whether to print progress messages. Default: True
    log : Log, optional
        Logging object. Default: Log()
    
    Returns
    -------
    tuple
        (fig, axes) where fig is matplotlib Figure and axes is list of axes
    
    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.DataFrame({
    ...     'study': ['Study1', 'Study2', 'Study3'],
    ...     'beta': [0.5, 0.7, 0.6],
    ...     'se': [0.2, 0.15, 0.18]
    ... })
    >>> fig, axes = plot_forest(data, study_col='study')
    """
    log.write("Start to create forest plot...", verbose=verbose)
    
    # Load data
    if isinstance(data, str):
        log.write(f" -Loading data from file: {data}", verbose=verbose)
        meta_df_all = pd.read_csv(data, sep="\s+")
    elif isinstance(data, pd.DataFrame):
        if data.empty:
            log.warning("Input DataFrame is empty!")
            return None, None
        meta_df_all = data.copy()
    else:
        raise TypeError(f"data must be pd.DataFrame or str, got {type(data)}")
    
    # Validate required columns
    required_cols = [study_col, beta_col, se_col]
    missing_cols = [col for col in required_cols if col not in meta_df_all.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Handle grouping
    if group_col is False or group_col is None:
        group_col = "_GROUP"
        meta_df_all[group_col] = "Group1"
        log.write(" -No grouping column specified, using single group", verbose=verbose)
    elif group_col not in meta_df_all.columns:
        log.warning(f"Group column '{group_col}' not found, using single group")
        group_col = "_GROUP"
        meta_df_all[group_col] = "Group1"
    
    # Prepare combine_effects kwargs
    if combine_effects_kwargs is None:
        combine_effects_kwargs = {}
    
    # Process each group
    df_array = []
    df_array_name = []
    df_array_het = []
    
    unique_groups = meta_df_all[group_col].unique()
    log.write(f" -Processing {len(unique_groups)} group(s)...", verbose=verbose)
    
    for group_name in unique_groups:
        meta_df = meta_df_all.loc[meta_df_all[group_col] == group_name, :].copy()
        log.write(f"  -Group '{group_name}': {len(meta_df)} studies", verbose=verbose)
        
        # Check for missing values
        if meta_df[beta_col].isna().any() or meta_df[se_col].isna().any():
            log.warning(f"  -Removing {meta_df[[beta_col, se_col]].isna().any(axis=1).sum()} studies with missing beta/SE")
            meta_df = meta_df.dropna(subset=[beta_col, se_col])
        
        if len(meta_df) == 0:
            log.warning(f"  -No valid studies in group '{group_name}', skipping")
            continue
        
        # Perform meta-analysis
        try:
            if HAS_STATSMODELS:
                combined_results = combine_effects(
                    meta_df[beta_col].values,
                    meta_df[se_col].values ** 2,  # combine_effects expects variance
                    row_names=meta_df[study_col].values,
                    **combine_effects_kwargs
                )
            else:
                # Fallback: simple fixed-effects meta-analysis
                log.write("  -statsmodels not available, using simple fixed-effects", verbose=verbose)
                combined_results = _simple_meta_analysis(
                    meta_df[beta_col].values,
                    meta_df[se_col].values,
                    meta_df[study_col].values
                )
        except Exception as e:
            log.warning(f"  -Error in meta-analysis for group '{group_name}': {str(e)}")
            continue
        
        # Extract heterogeneity statistics if requested
        het_stats = None
        if meta and HAS_STATSMODELS:
            try:
                het_test = combined_results.test_homogeneity()
                Q = het_test.statistic
                Qdf = het_test.df
                Qp = het_test.pvalue
                I2 = max(0, (Q - Qdf) / Q * 100) if Q > 0 else 0
                het_stats = (Q, Qdf, Qp, I2)
            except Exception as e:
                log.warning(f"  -Could not compute heterogeneity statistics: {str(e)}")
                het_stats = None
        
        # Prepare data for plotting
        if HAS_STATSMODELS:
            df_to_plot = combined_results.summary_frame().iloc[:-3, :].copy()
        else:
            df_to_plot = combined_results.copy()
        
        df_to_plot["YORDER"] = range(len(df_to_plot), 0, -1)
        df_to_plot["OR"] = np.exp(df_to_plot["eff"])
        df_to_plot["OR_95L"] = np.exp(df_to_plot["ci_low"])
        df_to_plot["OR_95U"] = np.exp(df_to_plot["ci_upp"])
        df_to_plot["OR_95L-OR"] = df_to_plot["OR"] - df_to_plot["OR_95L"]
        df_to_plot["OR_95U-OR"] = df_to_plot["OR_95U"] - df_to_plot["OR"]
        df_to_plot["MARKERSIZE"] = 1 / (df_to_plot["sd_eff"] ** 2)
        
        df_array.append(df_to_plot)
        df_array_name.append(group_name)
        df_array_het.append(het_stats if het_stats else [])
    
    if len(df_array) == 0:
        log.warning("No valid groups to plot!")
        return None, None
    
    # Set up figure
    nrows = len(df_array)
    style = set_plot_style(
        plot="plot_forest",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        save=save,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    
    fig_kwargs_style = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    fontsize = style.get("fontsize", fontsize)
    font_family = style.get("font_family", font_family)
    
    # Calculate figure size
    base_width = fig_kwargs_style.get("figsize", (15, 5))[0]
    base_height = fig_kwargs_style.get("figsize", (15, 5))[1]
    total_height = sum(len(df) for df in df_array) / compact_factor
    fig_kwargs_style["figsize"] = (base_width, max(base_height, total_height))
    
    # Create subplots
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=3,
        sharex=sharex,
        gridspec_kw={"width_ratios": width_ratios},
        **fig_kwargs_style,
    )
    
    # Handle single row case
    if nrows == 1:
        axes = axes.reshape(1, -1)
    
    # Set up colors (default from plot_mqq)
    if colors is None:
        colors = ["#597FBD", "#74BAD3"]  # Default colors from plot_mqq
    
    # Plot each group
    for i in range(nrows):
        plot_row(
            axes=axes[i],
            df_to_plot=df_array[i],
            group_name=df_array_name[i],
            het_stats=df_array_het[i],
            meta=meta,
            is_first_row=(i == 0),
            fontsize=fontsize,
            font_family=font_family,
            colors=colors,
        )
    
    fig.tight_layout()
    
    # Save figure
    save_figure(fig, save, keyword="forest", save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    log.write("Finished creating forest plot.", verbose=verbose)
    return fig, axes


def plot_row(
    axes: Union[np.ndarray, List[Axes]],
    df_to_plot: pd.DataFrame,
    group_name: str,
    het_stats: Optional[Tuple[float, float, float, float]],
    meta: bool,
    is_first_row: bool,
    fontsize: int = 12,
    font_family: str = "Arial",
    colors: Optional[List[str]] = None,
) -> None:
    """
    Plot a single row of the forest plot.
    
    Parameters
    ----------
    axes : np.ndarray
        Array of 3 axes [group_label, plot, text]
    df_to_plot : pd.DataFrame
        DataFrame with study results and meta-analysis
    group_name : str
        Name of the group
    het_stats : tuple or None
        Heterogeneity statistics (Q, Qdf, Qp, I2) or None
    meta : bool
        Whether to show meta-analysis statistics
    is_first_row : bool
        Whether this is the first row (for titles)
    fontsize : int, optional
        Font size. Default: 12
    font_family : str, optional
        Font family. Default: "Arial"
    colors : list of str, optional
        List of colors for alternating studies. Default: None
    """
    # Set up colors
    if colors is None:
        colors = ["#597FBD", "#74BAD3"]  # Default colors from plot_mqq
    
    # Separate individual studies from fixed effect
    studies = df_to_plot.loc[df_to_plot.index != "fixed effect", :].copy()
    fixed_effect = df_to_plot.loc[df_to_plot.index == "fixed effect", :]
    
    if len(studies) == 0:
        return
    
    # Calculate marker sizes (proportional to precision)
    studies["MARKERSIZE_N"] = studies["MARKERSIZE"] / np.nanmax(studies["MARKERSIZE"]) * 0.8
    
    # Convert to pixels for proper scaling
    pixels = axes[1].transData.transform((0, len(df_to_plot) + 0.5))[1] - axes[1].transData.transform((0, 0))[1]
    studies["MARKERSIZE_P"] = studies["MARKERSIZE_N"] * (pixels / (len(df_to_plot) + 0.5))
    
    # Assign color to studies (use lighter blue #74BAD3)
    study_color = colors[1] if len(colors) > 1 else colors[0]  # Use second color (lighter blue)
    
    # Plot individual study points with color
    axes[1].scatter(
        x=studies["OR"],
        y=studies["YORDER"],
        s=studies["MARKERSIZE_P"] * 2,
        marker="s",
        color=study_color,
        edgecolor="black",
        linewidth=0.5,
        zorder=3
    )
    
    # Plot error bars with matching color
    axes[1].errorbar(
        x=studies["OR"],
        y=studies["YORDER"],
        xerr=studies[["OR_95L-OR", "OR_95U-OR"]].T.values,
        elinewidth=2,
        lw=0,
        marker="",
        color=study_color,
        zorder=2
    )
    
    # Plot fixed effect diamond
    if len(fixed_effect) > 0:
        xm = fixed_effect["OR"].values[0]
        xl = fixed_effect["OR_95L"].values[0]
        xr = fixed_effect["OR_95U"].values[0]
        ym = fixed_effect["YORDER"].values[0]
        yu = ym + 0.25
        yd = ym - 0.25
        
        xy = [
            (xl, ym),  # left
            (xm, yu),  # top
            (xr, ym),  # right
            (xm, yd),  # bottom
        ]
        
        patches = [Polygon(xy, closed=True)]
        # Use darker blue (#597FBD) for fixed effect
        fixed_effect_color = colors[0] if len(colors) > 0 else "#597FBD"
        p = PatchCollection(patches, color=fixed_effect_color, edgecolor="black", alpha=0.8, zorder=4)
        axes[1].add_collection(p)
        
        # Add vertical line for fixed effect (matching color)
        axes[1].axvline(x=xm, color=fixed_effect_color, ls="--", zorder=1)
    
    # Add reference line at OR=1
    axes[1].axvline(x=1, color="black", zorder=1)
    
    # Set y-axis limits
    y_min, y_max = 0.25, len(df_to_plot) + 0.75
    axes[1].set_ylim(y_min, y_max)
    axes[2].set_ylim(y_min, y_max)
    axes[2].set_xlim([0, 1])
    
    # Hide spines
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    # Hide ticks for side panels
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[2].set_xticks([])
    axes[2].set_yticks([])
    
    # Add titles (only on first row)
    if is_first_row:
        axes[0].set_title("Groups", size=fontsize, loc="left", fontweight="bold", family=font_family)
        axes[1].set_title("Odds Ratio", size=fontsize, fontweight="bold", family=font_family)
        axes[1].set_ylabel("Studies", fontsize=fontsize, family=font_family)
        axes[2].set_title("Odds Ratio, [95% CI]", size=fontsize, loc="right", fontweight="bold", family=font_family)
    
    # Add heterogeneity statistics
    if meta and het_stats is not None and len(het_stats) == 4:
        Q, Qdf, Qp, I2 = het_stats
        
        # Format p-value nicely (similar to plot_compare_effect style)
        try:
            if Qp > 1e-300:
                p_str = f"{Qp:.2e}"
                if 'e' in p_str:
                    mantissa, exponent = p_str.split('e')
                    # Remove trailing zeros from mantissa
                    mantissa = f"{float(mantissa):.2f}".rstrip('0').rstrip('.')
                    exponent = int(exponent)
                    if exponent == 0:
                        p_formatted = f"${mantissa}$"
                    else:
                        p_formatted = rf"${mantissa} \times 10^{{{exponent}}}$"
                else:
                    p_formatted = f"${p_str}$"
            else:
                p_formatted = r"$< 1 \times 10^{-300}$"
        except (ValueError, OverflowError):
            # Fallback to simple format
            if Qp < 0.001:
                p_formatted = f"${Qp:.2e}$"
            elif Qp < 0.01:
                p_formatted = f"${Qp:.3f}$"
            else:
                p_formatted = f"${Qp:.2f}$"
        
        het_string1 = rf'$Q={Q:.2f}, Qdf={Qdf:.0f}$'
        het_string2 = rf'$p={p_formatted}, I^2={I2:.1f}\%$'
        axes[0].text(0, 0, s=het_string1 + "\n" + het_string2, size=fontsize, family=font_family)
    
    # Set y-axis ticks and labels
    axes[1].set_yticks(df_to_plot["YORDER"])
    axes[1].set_yticklabels(df_to_plot.index, size=fontsize, family=font_family)
    
    # Bold the fixed effect label
    yticklabels = axes[1].get_yticklabels()
    for label in yticklabels:
        if label.get_text() == "fixed effect":
            label.set_fontweight("bold")
    
    # Add group name
    axes[0].text(x=0, y=0.5, s=group_name, size=fontsize, family=font_family,
                 transform=axes[0].transAxes, va="center", ha="left")
    
    # Add text annotations with OR and CI
    for index, row in df_to_plot.iterrows():
        fontweight = "bold" if index == "fixed effect" else "normal"
        axes[2].text(
            x=1,
            y=row["YORDER"],
            s=f"{row['OR']:.2f} [{row['OR_95L']:.2f}, {row['OR_95U']:.2f}]",
            size=fontsize,
            fontweight=fontweight,
            family=font_family,
            ha="right",
            va="center",
        )


def _simple_meta_analysis(
    beta: Union[np.ndarray, List[float]],
    se: Union[np.ndarray, List[float]],
    row_names: Union[np.ndarray, List[str]]
) -> pd.DataFrame:
    """
    Simple fixed-effects meta-analysis (fallback when statsmodels unavailable).
    
    Parameters
    ----------
    beta : np.ndarray
        Effect estimates
    se : np.ndarray
        Standard errors
    row_names : np.ndarray
        Study names
    
    Returns
    -------
    pd.DataFrame
        Summary frame with individual studies and fixed effect
    """
    # Calculate weights (inverse variance)
    weights = 1 / (se ** 2)
    total_weight = weights.sum()
    
    # Fixed effect estimate
    beta_fixed = (beta * weights).sum() / total_weight
    se_fixed = np.sqrt(1 / total_weight)
    
    # Create summary DataFrame
    results = []
    
    # Individual studies
    for i, (b, s, name) in enumerate(zip(beta, se, row_names)):
        results.append({
            "eff": b,
            "sd_eff": s,
            "ci_low": b - 1.96 * s,
            "ci_upp": b + 1.96 * s,
        })
    
    # Fixed effect
    results.append({
        "eff": beta_fixed,
        "sd_eff": se_fixed,
        "ci_low": beta_fixed - 1.96 * se_fixed,
        "ci_upp": beta_fixed + 1.96 * se_fixed,
    })
    
    df = pd.DataFrame(results, index=list(row_names) + ["fixed effect"])
    return df
