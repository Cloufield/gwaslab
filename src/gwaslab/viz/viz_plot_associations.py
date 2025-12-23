import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from gwaslab.io.io_process_kwargs import _update_kwargs
from gwaslab.io.io_process_kwargs import _update_arg
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_style_options import set_plot_style

def _plot_associations(associations, 
                       values="BETA_GCV2",
                        fig_kwargs=None,
                        log=Log(),
                        verbose=True,
                        sort="P_GCV2",
                        fontsize=12,
                        font_family="Arial",
                        cmap=None,
                        ylabel="Traits",
                        xlabel="rsID - Gene Name",
                        save=None,
                        save_kwargs=None,
                        title=None,
                        heatmap_kwargs=None,
                        annot=True,
                        fmt=".2f",
                        vmin=None,
                        vmax=None,
                        center=None,
                        cbar=True,
                        cbar_kws=None,
                        xticklabel_kwargs=None,
                        yticklabel_kwargs=None,
                        **args               
                       ):
    """
    Plot associations as a heatmap using GCV2 format.
    
    Parameters
    ----------
    associations : pd.DataFrame
        Associations DataFrame in GCV2 format (with _GCV2 suffix columns)
    values : str, optional
        Column name for values to plot. Default: "BETA_GCV2"
    sort : str, optional
        Column name to sort by. Default: "P_GCV2"
    save : str, optional
        File path to save the figure. Default: None
    save_kwargs : dict, optional
        Additional arguments for saving the figure. Default: None
    title : str, optional
        Title for the plot. Default: None
    heatmap_kwargs : dict, optional
        Additional arguments passed to seaborn.heatmap. Default: None
    annot : bool, optional
        Whether to annotate cells with values. Default: True
    fmt : str, optional
        Format string for annotations. Default: ".2f"
    vmin, vmax : float, optional
        Values to anchor the colormap. Default: None (auto)
    center : float, optional
        Value to center the colormap. Default: None
    cbar : bool, optional
        Whether to draw a colorbar. Default: True
    cbar_kws : dict, optional
        Additional arguments for colorbar. Default: None
    xticklabel_kwargs : dict, optional
        Additional arguments for x-axis tick labels. Default: None
    yticklabel_kwargs : dict, optional
        Additional arguments for y-axis tick labels. Default: None
    """
    
    # Extract dataframe if Sumstats object is passed
    if hasattr(associations, 'data') and not isinstance(associations, pd.DataFrame):
        associations = associations.data
    
    log.write("Start to create heatmap for associations...", verbose=verbose)
    
    cmap = _update_arg(cmap, "RdBu")
    
    # GCV2 format column names
    rsid_col = 'rsID_GCV2' if 'rsID_GCV2' in associations.columns else 'rsID'
    trait_col = 'efo_traits_GCV2' if 'efo_traits_GCV2' in associations.columns else 'reported_trait_GCV2'
    gene_col = 'GENENAME_GCV2' if 'GENENAME_GCV2' in associations.columns else None
    
    # Use defaults if provided values don't exist
    if values not in associations.columns:
        values = 'BETA_GCV2' if 'BETA_GCV2' in associations.columns else values
    if sort not in associations.columns:
        sort = 'P_GCV2' if 'P_GCV2' in associations.columns else sort
    
    log.write(f" -Value to plot: {values}", verbose=verbose)
    log.write(f" -Total number of associations :{len(associations)} ", verbose=verbose)
    log.write(f" -Sorting associations by :{sort} ", verbose=verbose)
    
    # Check required columns
    if rsid_col not in associations.columns:
        log.warning(f"Required column '{rsid_col}' not found. Available columns: {list(associations.columns)}", verbose=verbose)
        return None, None
    
    if trait_col not in associations.columns:
        log.warning(f"Required column '{trait_col}' not found. Available columns: {list(associations.columns)}", verbose=verbose)
        return None, None
    
    # Sort and deduplicate
    if sort in associations.columns:
        associations = associations.sort_values(by=[sort])
    associations = associations.drop_duplicates(subset=[rsid_col, trait_col])
    
    log.write(f" -Keeping unique rsID-traits pairs :{len(associations)} ", verbose=verbose)
    
    # Drop NA values
    if values in associations.columns:
        associations = associations.dropna(subset=[values])
        log.write(f" -Keeping associations without NA in {values} :{len(associations)} ", verbose=verbose)
    else:
        log.warning(f"Values column '{values}' not found. Available columns: {list(associations.columns)}", verbose=verbose)
        return None, None
    
    log.write(f" -Total number of unique variants for plotting :{associations[rsid_col].nunique()} ", verbose=verbose)
    log.write(f" -Total number of unique traits for plotting:{associations[trait_col].nunique()} ", verbose=verbose)
    
    # Create pivot table
    if gene_col and gene_col in associations.columns:
        # Include gene name in index
        index_cols = [rsid_col, gene_col]
    else:
        # Just use rsID
        index_cols = [rsid_col]
        log.write(f" -Gene column '{gene_col}' not found, using rsID only", verbose=verbose)
    
    matrix_beta = associations.pivot_table(
        index=index_cols,
        columns=trait_col,
        values=values,
        aggfunc='first'  # Take first value if duplicates
    )
    
    matrix_beta = matrix_beta.astype("float64")
    
    # Drop column level if it exists (from pivot_table with values parameter)
    if matrix_beta.columns.nlevels > 1:
        matrix_beta.columns = matrix_beta.columns.droplevel(0)
    
    height = max(2, len(matrix_beta)//2)
    width = max(2, len(matrix_beta.columns)//2)
    
    # Update fig_kwargs with calculated size
    if fig_kwargs is None:
        fig_kwargs = {}
    if "figsize" not in fig_kwargs:
        fig_kwargs["figsize"] = (width, height)
    
    style = set_plot_style(
        plot="plot_associations",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig, ax = plt.subplots(**(style.get("fig_kwargs", {})))
    save_kwargs = style.get("save_kwargs", {})
    
    # Prepare heatmap arguments
    heatmap_defaults = {
        "annot": annot,
        "fmt": fmt,
        "cmap": cmap,
        "ax": ax,
        "cbar": cbar,
    }
    if vmin is not None:
        heatmap_defaults["vmin"] = vmin
    if vmax is not None:
        heatmap_defaults["vmax"] = vmax
    if center is not None:
        heatmap_defaults["center"] = center
    if cbar_kws is not None:
        heatmap_defaults["cbar_kws"] = cbar_kws
    
    # Merge with user-provided heatmap_kwargs
    if heatmap_kwargs is not None:
        heatmap_defaults.update(heatmap_kwargs)
    
    sns.heatmap(matrix_beta.T, **heatmap_defaults)
    
    # Set labels
    ax.set_xlabel(xlabel, fontsize=fontsize, fontfamily=font_family)
    ax.set_ylabel(ylabel, fontsize=fontsize, fontfamily=font_family)
    
    # Set title if provided
    if title is not None:
        ax.set_title(title, fontsize=fontsize, fontfamily=font_family)
    
    # Apply tick label kwargs if provided
    if xticklabel_kwargs is not None:
        ax.set_xticklabels(ax.get_xticklabels(), **xticklabel_kwargs)
    if yticklabel_kwargs is not None:
        ax.set_yticklabels(ax.get_yticklabels(), **yticklabel_kwargs)
    
    # Save figure if requested
    if save is not None:
        from gwaslab.viz.viz_aux_save_figure import save_figure
        save_figure(fig, save, keyword="associations", save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    log.write("Finished creating heatmap for associations.", verbose=verbose)
    return fig, ax
