from typing import TYPE_CHECKING, Optional, Union, Tuple, Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
    from gwaslab.g_Sumstats import Sumstats
    from gwaslab.bd.bd_chromosome_mapper import ChromosomeMapper


def _edges_from_centers(x):
    """Convert center coordinates (length n) into edge coordinates (length n+1)."""
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        # arbitrary edge padding
        return np.array([x[0] - 0.5, x[0] + 0.5], dtype=float)

    mid = (x[:-1] + x[1:]) / 2.0
    first = x[0] - (mid[0] - x[0])
    last = x[-1] + (x[-1] - mid[-1])
    return np.concatenate([[first], mid, [last]])




def _prepare_ld_data_from_vcf(
    vcf_path: str,
    region: Tuple[Union[int, str], int, int],
    sumstats: Union[pd.DataFrame, 'Sumstats'],
    pos_col: str,
    nea_col: str,
    ea_col: str,
    tabix: Optional[bool],
    mapper: Optional[Any],
    log: Log,
    verbose: bool
) -> Tuple[np.ndarray, pd.DataFrame, pd.DataFrame, bool, Optional[np.ndarray]]:
    """
    Extract LD matrix from VCF and prepare data for plotting.
    
    Returns
    -------
    ld_matrix : np.ndarray
        LD matrix from VCF (only matched variants).
    matched_sumstats : pd.DataFrame
        Matched sumstats from VCF.
    region_sumstats : pd.DataFrame
        All variants in region from sumstats.
    use_i_coordinate : bool
        Whether "i" coordinate system is available.
    all_i_positions : np.ndarray or None
        "i" coordinates for all variants in region, if available.
    """
    if sumstats is None:
        raise ValueError("sumstats must be provided when using vcf_path and region")
    
    log.write("Extracting LD matrix from VCF file...", verbose=verbose)
    from gwaslab.io.io_vcf import _get_ld_matrix_from_vcf
    
    matched_sumstats, ld_matrix = _get_ld_matrix_from_vcf(
        sumstats_or_dataframe=sumstats,
        vcf_path=vcf_path,
        region=region,
        log=log,
        verbose=verbose,
        pos=pos_col,
        nea=nea_col,
        ea=ea_col,
        mapper=mapper,
        tabix=tabix
    )
    
    if ld_matrix.size == 0:
        raise ValueError("No valid variants found in VCF for the specified region")
    
    # Get original sumstats filtered to region (all variants in region, not just matched)
    if hasattr(sumstats, 'data'):
        original_sumstats = sumstats.data
    else:
        original_sumstats = sumstats
    
    # Filter original sumstats to region to get all variant positions
    from gwaslab.util.util_in_filter_value import _filter_region
    region_sumstats = _filter_region(
        sumstats_or_dataframe=original_sumstats,
        region=region,
        chrom="CHR" if "CHR" in original_sumstats.columns else None,
        pos=pos_col,
        log=log,
        verbose=verbose
    )
    
    if len(region_sumstats) == 0:
        raise ValueError("No variants found in sumstats for the specified region")
    
    # Check if we should use "i" coordinate for alignment
    use_i_coordinate = False
    all_i_positions = None
    if "i" in region_sumstats.columns:
        all_i_positions = region_sumstats.sort_values(pos_col)["i"].values.astype(float)
        use_i_coordinate = True
        log.write("Using 'i' coordinate system for x-axis alignment with regional plot", verbose=verbose)
    
    return ld_matrix, matched_sumstats, region_sumstats, use_i_coordinate, all_i_positions


def _prepare_ld_matrix(
    ld_matrix: np.ndarray,
    matched_sumstats: pd.DataFrame,
    region_sumstats: pd.DataFrame,
    pos_col: str,
    log: Log,
    verbose: bool
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Prepare full LD matrix with all variant positions, filling in LD values for matched variants.
    
    Parameters
    ----------
    ld_matrix : np.ndarray
        LD matrix from VCF (only matched variants).
    matched_sumstats : pd.DataFrame
        Matched sumstats from VCF.
    region_sumstats : pd.DataFrame
        All variants in region from sumstats.
    pos_col : str
        Position column name.
    log : Log
        Logger instance.
    verbose : bool
        Whether to print log messages.
    
    Returns
    -------
    full_ld_matrix : np.ndarray
        Full LD matrix (n_all x n_all) with NaN for empty positions.
    all_positions : np.ndarray
        All variant positions in region, sorted.
    """
    # Get all positions from region_sumstats (all variants in region, sorted by position)
    region_sumstats = region_sumstats.sort_values(pos_col)
    all_positions = region_sumstats[pos_col].values.astype(float)
    n_all = len(all_positions)
    
    # Create mapping from position to index in all_positions
    pos_to_idx = {int(p): idx for idx, p in enumerate(all_positions)}
    
    # Get matched positions from VCF
    if pos_col not in matched_sumstats.columns:
        raise ValueError(f"{pos_col} column not found in matched sumstats")
    matched_positions = matched_sumstats[pos_col].values.astype(int)
    
    # Find indices of matched variants in the full position grid
    matched_indices = []
    valid_matched = []
    for i, mp in enumerate(matched_positions):
        if mp in pos_to_idx:
            matched_indices.append(pos_to_idx[mp])
            valid_matched.append(i)
    
    matched_indices = np.array(matched_indices, dtype=int)
    valid_matched = np.array(valid_matched, dtype=int)
    n_matched = len(matched_indices)
    
    if n_matched != ld_matrix.shape[0]:
        log.write(f"Warning: Number of matched variants ({n_matched}) doesn't match LD matrix size ({ld_matrix.shape[0]})", verbose=verbose)
        # Use only valid matched variants
        n_matched = min(n_matched, ld_matrix.shape[0], len(valid_matched))
        matched_indices = matched_indices[:n_matched]
        valid_matched = valid_matched[:n_matched]
    
    # Create full LD matrix with white (NaN) for empty positions
    full_ld_matrix = np.full((n_all, n_all), np.nan, dtype=float)
    
    # Fill in LD values for matched variants
    for i, valid_i in enumerate(valid_matched[:n_matched]):
        idx_i = matched_indices[i]
        for j, valid_j in enumerate(valid_matched[:n_matched]):
            idx_j = matched_indices[j]
            if idx_j >= idx_i:  # Upper triangle (including diagonal)
                if valid_i < ld_matrix.shape[0] and valid_j < ld_matrix.shape[1]:
                    full_ld_matrix[idx_i, idx_j] = ld_matrix[valid_i, valid_j]
                    if idx_i != idx_j:  # Symmetric
                        full_ld_matrix[idx_j, idx_i] = ld_matrix[valid_i, valid_j]
    
    log.write(f"Created LD matrix: {n_all} variant positions in region, {n_matched} with LD data, {n_all - n_matched} empty (white)", verbose=verbose)
    
    return full_ld_matrix, all_positions


def _calculate_regional_xlim(
    region: Tuple[Union[int, str], int, int],
    sumstats: Union[pd.DataFrame, 'Sumstats'],
    pos_col: str
) -> Tuple[float, float]:
    """
    Calculate xlim for regional plot mode to align with regional plot's x-axis.
    
    Parameters
    ----------
    region : tuple
        Region specification (chromosome, start, end).
    sumstats : DataFrame or Sumstats
        Sumstats filtered to region.
    pos_col : str
        Position column name.
    
    Returns
    -------
    xlim : tuple
        (xmin, xmax) for x-axis limits.
    """
    if hasattr(sumstats, 'data'):
        region_sumstats = sumstats.data
    else:
        region_sumstats = sumstats
    
    if "i" not in region_sumstats.columns or pos_col not in region_sumstats.columns:
        raise ValueError("Sumstats must contain 'i' and position columns for regional mode")
    
    # Calculate gene_track_start_i the same way as regional plot
    most_left_snp = region_sumstats["i"].idxmin()
    gene_track_offset = region_sumstats.loc[most_left_snp, pos_col] - region[1]
    gene_track_start_i = region_sumstats.loc[most_left_snp, "i"] - gene_track_offset - region[1]
    
    # The regional plot x-axis shows: [gene_track_start_i + region[1], gene_track_start_i + region[2]]
    regional_xlim = [gene_track_start_i + region[1], gene_track_start_i + region[2]]
    
    return regional_xlim[0], regional_xlim[1]


def _set_axis_limits(
    ax: 'Axes',
    ranks: np.ndarray,
    ld_tri: np.ma.MaskedArray,
    mode: str
) -> None:
    """
    Set xlim and ylim for the plot based on mode.
    
    Uses edge coordinates (not center coordinates) to prevent clipping of cell corners.
    Now uses ranks for uniform cell sizes.
    
    Parameters
    ----------
    ax : Axes
        Axes object to set limits on.
    ranks : np.ndarray
        Variant rank array (length n), should be [0, 1, 2, ..., n-1].
    ld_tri : np.ma.MaskedArray
        Masked LD matrix (upper triangle).
    mode : str
        'regional' or 'standalone'.
    """
    # Calculate edge coordinates (same as used in pcolormesh)
    # These are the actual coordinates used by pcolormesh, so we need to include them
    e = _edges_from_centers(ranks)  # (n+1,)
    X_edges, Y_edges = np.meshgrid(e, e, indexing="ij")  # (n+1, n+1)
    U_edges = (X_edges + Y_edges) / 2.0
    V_edges = (Y_edges - X_edges) / 2.0
    V_edges_negated = -V_edges
    
    # Get the full range of edge coordinates (these define the cell boundaries)
    # For upper triangle, we only need the lower half (V_negated <= 0)
    u_min_edge = np.min(U_edges)
    u_max_edge = np.max(U_edges)
    v_neg_min_edge = np.min(V_edges_negated)
    # For upper triangle, maximum V_negated is 0 (at the diagonal)
    v_neg_max_edge = 0.0
    
    # Set xlim: use U edge coordinate range with padding
    u_range = u_max_edge - u_min_edge
    padding = max(0.005 * u_range, abs(u_min_edge) * 1e-6) if u_range > 0 else 0.01
    # Only set xlim if not already set by shared axis
    current_xlim = ax.get_xlim()
    if current_xlim == (0.0, 1.0):  # Default matplotlib xlim
        ax.set_xlim(u_min_edge - padding, u_max_edge + padding)
    
    # Set ylim using edge coordinates with padding
    # Only show lower half: from v_neg_min_edge to 0
    v_neg_range = v_neg_max_edge - v_neg_min_edge  # This is 0 - v_neg_min_edge
    if v_neg_range > 0:
        padding = max(0.005 * v_neg_range, abs(v_neg_min_edge) * 1e-6)
        ax.set_ylim(v_neg_min_edge - padding, v_neg_max_edge + padding)
    else:
        padding = 0.1
        ax.set_ylim(v_neg_min_edge - padding, v_neg_max_edge + padding)


def _plot_ld_triangle(
    ax: 'Axes',
    ld: np.ndarray,
    ranks: np.ndarray,
    cmap: Union[str, matplotlib.colors.Colormap],
    vmin: float,
    vmax: float,
    i_coords: Optional[np.ndarray] = None,
    **kwargs
) -> Tuple[matplotlib.collections.QuadMesh, np.ma.MaskedArray]:
    """
    Plot the LD triangle using pcolormesh with 45° rotation.
    Uses variant ranks (0, 1, 2, ..., n-1) for uniform cell sizes.
    In regional mode, can transform ranks to i coordinates for x-axis alignment.
    
    Parameters
    ----------
    ax : Axes
        Axes object to plot on.
    ld : np.ndarray
        LD matrix (n x n).
    ranks : np.ndarray
        Variant rank array (length n), should be [0, 1, 2, ..., n-1] for uniform cells.
    cmap : str or Colormap
        Colormap for LD values.
    vmin : float
        Minimum value for colormap.
    vmax : float
        Maximum value for colormap.
    i_coords : np.ndarray, optional
        Deprecated: Always uses ranks for uniform cells. This parameter is ignored.
    **kwargs
        Additional arguments passed to pcolormesh.
    
    Returns
    -------
    QuadMesh
        The pcolormesh object.
    ld_tri : np.ma.MaskedArray
        Masked LD matrix (upper triangle).
    """
    # Mask lower triangle (keep diagonal + upper triangle)
    mask_lower = np.tril(np.ones((ld.shape[0], ld.shape[0]), dtype=bool), k=-1)
    # Also mask NaN values (empty positions) - these will appear white
    mask_nan = np.isnan(ld)
    # Combine masks: mask lower triangle OR NaN values
    combined_mask = mask_lower | mask_nan
    ld_tri = np.ma.array(ld, mask=combined_mask)
    
    # Build corner grid using ranks (uniform spacing)
    # Always use ranks for the grid to maintain uniform cell sizes
    # In regional mode, ax_ld_block does NOT share x-axis, so it uses rank-based coordinates
    e = _edges_from_centers(ranks)  # (n+1,)
    
    X, Y = np.meshgrid(e, e, indexing="ij")  # corners for pcolormesh
    
    # 45°-style transform (scaled; visually nice for LD triangles)
    U = (X + Y) / 2.0
    V = (Y - X) / 2.0
    
    # Negate V to create inverted triangle (without using invert_yaxis)
    V_negated = -V
    
    # Plot using pcolormesh with negated V to create inverted triangle
    m = ax.pcolormesh(
        U, V_negated, ld_tri,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        **kwargs
    )
    
    return m, ld_tri


def _plot_position_bar(
    ax_pos: 'Axes',
    actual_positions: np.ndarray,
    ranks: np.ndarray,
    fontsize: float = 10,
    font_family: str = "Arial",
    show_xticks: bool = True,
    region_step: int = 21,
    xlabel: Optional[str] = None,
    region: Optional[Tuple[Union[int, str], int, int]] = None
) -> None:
    """
    Plot a position bar showing actual chromosome positions.
    Points are positioned according to their actual basepair positions.
    
    Parameters
    ----------
    ax_pos : Axes
        Axes object for the position bar (typically above the LD block).
    actual_positions : np.ndarray
        Actual chromosome positions (length n).
    ranks : np.ndarray
        Variant ranks (length n), should be [0, 1, 2, ..., n-1].
    fontsize : float
        Font size for labels. Default: 10.
    font_family : str
        Font family. Default: "Arial".
    show_xticks : bool
        Whether to show x-axis ticks. Default: True.
    region_step : int
        Number of x-axis tick positions. Default: 21.
    xlabel : str, optional
        Label for x-axis. If None, will be set based on region if provided. Default: None.
    region : tuple, optional
        Region specification (chromosome, start, end). If provided, xlabel will include chromosome number. Default: None.
    """
    n = len(ranks)
    
    # Use actual positions for x-axis (not ranks)
    pos_min = actual_positions.min()
    pos_max = actual_positions.max()
    pos_range = pos_max - pos_min
    
    # Plot the position bar as a horizontal line/bar spanning the actual position range
    bar_height = 0.1
    ax_pos.barh(0, width=pos_range, left=pos_min, height=bar_height, color='gray', alpha=0.3)
    
    # Add vertical lines at ALL variant positions (using actual positions) that fill ylim(0, 0.05)
    ax_pos.vlines(actual_positions, 0, 0.05, colors='black', linewidth=1, zorder=3)
    
    # Set x-axis limits based on actual positions (only if not shared)
    # In regional mode, xlim is set by shared x-axis
    if not (hasattr(ax_pos, '_shared_axes') and 'x' in ax_pos._shared_axes):
        padding = pos_range * 0.02 if pos_range > 0 else (pos_max * 0.01 if pos_max > 0 else 1)
        ax_pos.set_xlim(pos_min - padding, pos_max + padding)
    ax_pos.set_ylim(0, 0.05)

    if show_xticks:
        # Set xticks using np.linspace like in viz_plot_regional2.py (line 533)
        tick_positions = np.linspace(pos_min, pos_max, num=region_step)
        
        # Format position labels like region_ticks in viz_plot_regional2.py (line 523)
        # Format as MB with 3 decimal places
        tick_labels = list(map('{:.3f}'.format, tick_positions.astype("int")/1000000))
        
        # Set xticks and xticklabels like in viz_plot_regional2.py (lines 533-534)
        ax_pos.set_xticks(tick_positions)
        ax_pos.set_xticklabels(tick_labels, rotation=45)
        
        # Configure x-axis to show ticks on top
        ax_pos.xaxis.set_ticks_position('top')
        ax_pos.xaxis.set_label_position('top')
        ax_pos.tick_params(axis='x', which='major', top=True, bottom=False, labeltop=True, labelbottom=False)

        # Set xlabel like in viz_plot_manhattan_like.py (lines 388-394)
        if region is not None:
            xlabel = "Chromosome " + str(region[0]) + " (MB)"
        else:
            xlabel = "Genomic position"
        ax_pos.set_xlabel(xlabel, fontsize=fontsize, family=font_family)
    
    # Remove y-axis
    ax_pos.set_yticks([])
    ax_pos.spines['left'].set_visible(True)
    ax_pos.spines['right'].set_visible(True)
    ax_pos.spines['top'].set_visible(True)  # Show top spine since bar is on top
    ax_pos.spines['bottom'].set_visible(True)  # Hide bottom spine


def _add_annotation_lines(
    ax_ld: 'Axes',
    ax_pos: 'Axes',
    ranks: np.ndarray,
    actual_positions: np.ndarray,
    n_lines: Optional[int] = None,
    line_color: str = 'gray',
    line_alpha: float = 0.3,
    line_style: str = '--',
    line_width: float = 0.5,
    lead_snp_is: Optional[list] = None,
    lead_snp_is_color: Optional[list] = None
) -> None:
    """
    Add annotation lines connecting position bar to LD block cell right top corners.
    Annotates all variants if n_lines is None, otherwise annotates n_lines variants.
    
    Parameters
    ----------
    ax_ld : Axes
        Axes object for the LD block plot.
    ax_pos : Axes
        Axes object for the position bar.
    ranks : np.ndarray
        Variant ranks (length n), should be [0, 1, 2, ..., n-1].
    actual_positions : np.ndarray
        Actual chromosome positions (length n).
    n_lines : int, optional
        Number of annotation lines to draw. If None, draws lines for all variants.
        Default: None (all variants).
    line_color : str
        Color for annotation lines. Default: 'gray'.
    line_alpha : float
        Alpha (transparency) for annotation lines. Default: 0.3.
    line_style : str
        Line style for annotation lines. Default: '--'.
    line_width : float
        Line width for annotation lines. Default: 0.5.
    lead_snp_is : list, optional
        List of lead SNP i-coordinates (x-axis positions). Default: None.
    lead_snp_is_color : list, optional
        List of colors corresponding to lead SNPs. Default: None.
    """
    n = len(ranks)
    
    # Create a mapping from actual positions to lead SNP colors
    lead_snp_color_map = {}
    if lead_snp_is is not None and lead_snp_is_color is not None:
        for lead_pos, lead_color in zip(lead_snp_is, lead_snp_is_color):
            lead_snp_color_map[lead_pos] = lead_color
    
    # Select positions to annotate - all if n_lines is None, otherwise subset
    if n_lines is None:
        indices = np.arange(n)  # All variants
    elif n_lines >= n:
        indices = np.arange(n)
    else:
        indices = np.linspace(0, n-1, n_lines, dtype=int)
    
    # Get position bar ylim
    ylim_pos = [0]*len(actual_positions) #ax_pos.get_ylim()
    
    # For each selected variant, draw a line from position bar to LD block
    for idx in indices:
        rank = ranks[idx]
        actual_pos = actual_positions[idx]
        
        # Determine color and line width for this line: use lead SNP color and thicker width if available
        current_color = line_color
        current_linewidth = line_width
        is_lead_snp = actual_pos in lead_snp_color_map
        if is_lead_snp:
            current_color = lead_snp_color_map[actual_pos]
            # Make lead SNP lines more visible with thicker width
            current_linewidth = max(line_width * 2.5, 2.0)  # At least 2.0, or 2.5x the default width
        
        # Start point: position bar (actual position, bottom of bar)
        x_start = actual_pos  # Use actual position, not rank
        y_start = ylim_pos[0]  # Bottom of position bar (since it's on top)
        
        # End point: LD block cell right top corner
        # For a cell at rank i, the right top corner in the rotated coordinate system:
        # - The cell spans from rank i to rank i+1 in both X and Y
        # - The right edge is at rank i+1
        # - The top edge is at the diagonal (where X = Y)
        # - In transformed coordinates: U = (X + Y) / 2, V_negated = -(Y - X) / 2
        # - For the right top corner of cell (i, i): X = i+1, Y = i+1
        # - So U = (i+1 + i+1) / 2 = i+1, V_negated = -(i+1 - i+1) / 2 = 0
        
        # Calculate edge coordinates
        e = _edges_from_centers(ranks)
        if int(rank) + 1 < len(e):
            cell_right = e[int(rank) + 1]
        else:
            cell_right = e[-1]
        
        # Right top corner: U = cell_right, V_negated = 0
        x_end = cell_right - 0.5
        y_end = 0.0 + 0.5
        
        # Draw line using ConnectionPatch
        con = ConnectionPatch(
            xyA=(x_start, y_start), xyB=(x_end, y_end),
            coordsA=ax_pos.transData, coordsB=ax_ld.transData,
            axesA=ax_pos, axesB=ax_ld,
            color=current_color, alpha=line_alpha, linestyle=line_style, linewidth=current_linewidth,
            zorder=1
        )
        ax_ld.add_patch(con)


def _draw_ld_block_grid(
    ax: 'Axes',
    ranks: np.ndarray,
    grid_kwargs: Optional[dict] = None
) -> None:
    """
    Draw grid lines on the LD block plot.
    
    Parameters
    ----------
    ax : Axes
        Axes object for the LD block plot.
    ranks : np.ndarray
        Variant rank array (length n).
    grid_kwargs : dict, optional
        Additional keyword arguments for grid lines (e.g., {'color': 'gray', 'linewidth': 0.5, 'alpha': 0.5}).
        Default: None.
    """
    n = len(ranks)
    
    # Get edge coordinates (same as used in _plot_ld_triangle)
    e = _edges_from_centers(ranks)  # (n+1,)
    
    # Default grid styling
    default_grid_kwargs = {
        'color': 'gray',
        'linewidth': 0.5,
        'alpha': 0.3,
        'linestyle': '-',
        'zorder': 5
    }
    
    # Update with user-provided kwargs if any
    if grid_kwargs is not None:
        default_grid_kwargs.update(grid_kwargs)
    
    # Draw grid lines along cell boundaries in the transformed coordinate system
    # For the upper triangle, we draw lines along the edges of cells
    
    # Draw lines along constant X (from diagonal to top boundary)
    for i in range(len(e)):
        x_edge = e[i]
        y_max = e[-1]
        
        # Only draw if this X edge is within the triangle (x_edge <= y_max)
        if x_edge <= y_max:
            # Start at diagonal where Y = X: U = x_edge, V_negated = 0
            u_start = x_edge
            v_neg_start = 0.0
            
            # End at top boundary where Y = y_max
            u_end = (x_edge + y_max) / 2.0
            v_neg_end = -(y_max - x_edge) / 2.0
            
            # Only draw if end point is in upper triangle (V_negated <= 0)
            if v_neg_end <= 0:
                ax.plot([u_start, u_end], [v_neg_start, v_neg_end], **default_grid_kwargs)
    
    # Draw lines along constant Y (from diagonal to right boundary)
    for j in range(len(e)):
        y_edge = e[j]
        x_max = e[-1]
        
        # Only draw if this Y edge is within the triangle (y_edge <= x_max for upper triangle)
        # Actually, for upper triangle, we want y_edge >= 0 and the line should go from diagonal to right
        if y_edge >= 0:
            # Start at diagonal where X = Y: U = y_edge, V_negated = 0
            u_start = y_edge
            v_neg_start = 0.0
            
            # End at right boundary where X = x_max (only if y_edge <= x_max)
            if y_edge <= x_max:
                u_end = (x_max + y_edge) / 2.0
                v_neg_end = -(y_edge - x_max) / 2.0
                
                # Only draw if end point is in upper triangle (V_negated <= 0)
                if v_neg_end <= 0:
                    ax.plot([u_start, u_end], [v_neg_start, v_neg_end], **default_grid_kwargs)


def _annotate_ld_block_left_side(
    ax: 'Axes',
    ranks: np.ndarray,
    actual_positions: np.ndarray,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    region_sumstats: Optional[pd.DataFrame] = None,
    pos_col: str = "POS",
    snpid_col: str = "SNPID",
    anno_col: Optional[str] = None,
    anno_set: Optional[Union[list, set, np.ndarray]] = None,
    fontsize: float = 8,
    font_family: str = "Arial",
    color: str = "black",
    tick_length: float = 0.02,
    offset: float = 0.05,
    ld_block_anno_kwargs: Optional[dict] = None
) -> None:
    """
    Add annotations on the left side of the LD block triangle.
    Annotations are placed along the left vertical edge of the triangle and styled like axis ticks and labels.
    Font size is automatically calculated to prevent overlapping between annotations.
    
    Parameters
    ----------
    ax : Axes
        Axes object for the LD block plot.
    ranks : np.ndarray
        Variant rank array (length n).
    actual_positions : np.ndarray
        Actual chromosome positions (length n).
    region : tuple, optional
        Region specification (chromosome, start, end). Used to get chromosome if not in sumstats.
    region_sumstats : pd.DataFrame, optional
        Sumstats DataFrame that may contain CHR, SNPID, and annotation columns.
    pos_col : str
        Position column name. Default: "POS".
    snpid_col : str
        SNPID column name for filtering. Default: "SNPID".
    anno_col : str, optional
        Column name to use for annotation text. If None, uses "chr:pos" format. Default: None.
    anno_set : list, set, or np.ndarray, optional
        List of SNPIDs to annotate. If None, annotates all variants. Default: None.
    fontsize : float
        Initial font size for annotations. Will be auto-reduced if needed to prevent overlap.
        Default: 8.
    font_family : str
        Font family for annotations. Default: "Arial".
    color : str
        Text and tick color. Default: "black".
    tick_length : float
        Length of tick marks as fraction of y-axis range. Default: 0.02.
    offset : float
        Horizontal offset from the left edge as fraction of x-axis range. Default: 0.05.
    ld_block_anno_kwargs : dict, optional
        Additional keyword arguments for text annotations. If 'fontsize' is explicitly set here,
        it will override the auto-calculation. Default: None.
    """
    n = len(ranks)
    
    # Get annotation text for each variant
    annotation_texts = []
    snpids = []
    
    if region_sumstats is not None:
        # Sort region_sumstats by position to match ranks order
        if pos_col in region_sumstats.columns:
            region_sumstats_sorted = region_sumstats.sort_values(pos_col).reset_index(drop=True)
        else:
            region_sumstats_sorted = region_sumstats.copy()
        
        # Get SNPIDs for filtering
        if snpid_col in region_sumstats_sorted.columns:
            snpids = region_sumstats_sorted[snpid_col].values.tolist()
        else:
            snpids = [None] * len(region_sumstats_sorted)
        
        # Get annotation text from anno_col if provided, otherwise use chr:pos
        if anno_col is not None and anno_col in region_sumstats_sorted.columns:
            # Use specified column for annotation
            annotation_texts = region_sumstats_sorted[anno_col].astype(str).values.tolist()
        else:
            # Use chr:pos format
            chromosomes = []
            if "CHR" in region_sumstats_sorted.columns and pos_col in region_sumstats_sorted.columns:
                chromosomes = region_sumstats_sorted["CHR"].astype(str).values.tolist()
            elif region is not None:
                chromosomes = [str(region[0])] * len(region_sumstats_sorted)
            else:
                chromosomes = ["?"] * len(region_sumstats_sorted)
            
            positions = region_sumstats_sorted[pos_col].values if pos_col in region_sumstats_sorted.columns else actual_positions
            annotation_texts = [f"{chr}:{int(pos)}" for chr, pos in zip(chromosomes, positions)]
        
        # Ensure we have the same length as ranks
        if len(annotation_texts) != n:
            # If lengths don't match, try to match by position
            annotation_texts_matched = []
            snpids_matched = []
            pos_to_text = {}
            pos_to_snpid = {}
            
            for idx, row in region_sumstats_sorted.iterrows():
                if pos_col in row and pd.notna(row[pos_col]):
                    pos_key = int(row[pos_col])
                    if anno_col is not None and anno_col in row:
                        pos_to_text[pos_key] = str(row[anno_col])
                    else:
                        chr_val = str(int(row["CHR"])) if "CHR" in row and pd.notna(row["CHR"]) else (str(region[0]) if region is not None else "?")
                        pos_to_text[pos_key] = f"{chr_val}:{int(pos_key)}"
                    if snpid_col in row:
                        pos_to_snpid[pos_key] = row[snpid_col]
            
            for pos in actual_positions:
                pos_int = int(pos)
                if pos_int in pos_to_text:
                    annotation_texts_matched.append(pos_to_text[pos_int])
                    snpids_matched.append(pos_to_snpid.get(pos_int, None))
                else:
                    # Fallback to chr:pos
                    chr_str = str(region[0]) if region is not None else "?"
                    annotation_texts_matched.append(f"{chr_str}:{pos_int}")
                    snpids_matched.append(None)
            
            annotation_texts = annotation_texts_matched
            snpids = snpids_matched
    else:
        # No region_sumstats, use chr:pos format
        if region is not None:
            chr_str = str(region[0])
            chromosomes = [chr_str] * n
        else:
            chromosomes = ["?"] * n
        
        annotation_texts = [f"{chr}:{int(pos)}" for chr, pos in zip(chromosomes, actual_positions)]
        snpids = [None] * n
    
    # Filter variants based on anno_set if provided
    if anno_set is not None:
        # Convert anno_set to set for fast lookup
        anno_set_set = set(anno_set) if not isinstance(anno_set, set) else anno_set
        # Create boolean mask for variants to annotate
        annotate_mask = []
        for snpid in snpids:
            if snpid is not None and snpid in anno_set_set:
                annotate_mask.append(True)
            elif snpid is None:
                # If SNPID is not available, don't annotate (can't match)
                annotate_mask.append(False)
            else:
                annotate_mask.append(False)
    else:
        # Annotate all variants
        annotate_mask = [True] * n
    
    # Get current axis limits to calculate offsets
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]
    
    # Calculate offset from left edge (in U coordinate)
    # The left edge of the triangle is at X = 0 (first variant)
    # In transformed coordinates: U = (X + Y) / 2, V_negated = -(Y - X) / 2
    # For the left edge (X = 0): U = Y / 2, V_negated = -Y / 2
    # So along left edge: U = -V_negated
    
    # For each variant at rank i, find its position on the left edge
    # The left edge point for variant i would be at X = 0, Y = i (if it exists)
    # But since we're in upper triangle, for variant i, the leftmost point is at X = i, Y = i
    # However, for the left side annotation, we want to place it along the left boundary
    # which is the line X = 0. For variant i, we project to this line.
    
    # Actually, for the left side, we want to place annotations along the left vertical edge
    # of the triangle. The left edge is where X = 0 (the first column).
    # For variant at rank i, we find where the left edge (X=0) intersects with its row.
    # But in upper triangle, row i starts at X = i, so we need to extrapolate.
    
    # Better approach: Place annotations along the left boundary line X = 0
    # For variant at rank i, the left boundary point would be at X = 0, Y = i
    # U = (0 + i) / 2 = i / 2
    # V_negated = -(i - 0) / 2 = -i / 2
    
    # Calculate offset perpendicular to the left edge
    # The left edge has slope: U = -V_negated, so a perpendicular offset
    # For simplicity, offset in U direction (to the left)
    u_offset = -offset * x_range
    
    # Calculate tick length in data coordinates
    tick_len_data = tick_length * y_range
    
    # Auto-calculate font size to prevent overlapping
    # Calculate minimum spacing between consecutive annotations along the left edge
    # Only consider variants that will actually be annotated
    annotated_indices = [i for i in range(n) if annotate_mask[i]]
    
    if len(annotated_indices) > 1:
        # Calculate positions for annotations along the left edge
        # For rank i: V_negated = -i/2
        v_negated_positions = [-rank / 2.0 for rank in ranks]
        
        # Find minimum spacing between consecutive annotated variants in data coordinates
        min_spacing_data = float('inf')
        for j in range(len(annotated_indices) - 1):
            i1 = annotated_indices[j]
            i2 = annotated_indices[j + 1]
            spacing = abs(v_negated_positions[i2] - v_negated_positions[i1])
            if spacing < min_spacing_data:
                min_spacing_data = spacing
        
        # Convert minimum spacing to display coordinates (pixels)
        # Get transform from data to display coordinates
        transform = ax.transData
        # Use actual annotation positions to calculate spacing in pixels
        # Find the pair with minimum spacing and convert to pixels
        min_spacing_pixels = float('inf')
        for j in range(len(annotated_indices) - 1):
            i1 = annotated_indices[j]
            i2 = annotated_indices[j + 1]
            # Get positions for consecutive annotations
            u1 = ranks[i1] / 2.0
            v1 = v_negated_positions[i1]
            u2 = ranks[i2] / 2.0
            v2 = v_negated_positions[i2]
            
            # Transform to display coordinates
            point1 = transform.transform((u1, v1))
            point2 = transform.transform((u2, v2))
            
            # Calculate distance in pixels (vertical distance along the left edge)
            spacing_pix = abs(point2[1] - point1[1])
            if spacing_pix < min_spacing_pixels:
                min_spacing_pixels = spacing_pix
        
        spacing_pixels = min_spacing_pixels
        
        # Estimate text height in pixels for a given font size
        # Font size in points: 1 point ≈ 1.33 pixels at 96 DPI
        # Text height is approximately fontsize * 1.2 (line height factor) * 1.33
        dpi = ax.figure.get_dpi() if hasattr(ax.figure, 'get_dpi') else 100
        points_to_pixels = dpi / 72.0
        text_height_factor = 1.3  # Factor to account for line height
        
        # Calculate maximum font size that fits in the spacing
        # We want: fontsize * text_height_factor * points_to_pixels < spacing_pixels * overlap_factor
        overlap_factor = 0.85  # Use 85% of spacing to prevent overlap
        max_fontsize_points = (spacing_pixels * overlap_factor) / (text_height_factor * points_to_pixels)
        
        # Clamp fontsize between reasonable bounds
        min_fontsize = 4.0
        max_fontsize = 20.0
        auto_fontsize = max(min_fontsize, min(max_fontsize, max_fontsize_points))
        
        # Use the smaller of user-specified fontsize or auto-calculated fontsize
        # But allow user override if explicitly set in kwargs
        if ld_block_anno_kwargs is not None and 'fontsize' in ld_block_anno_kwargs:
            # User explicitly set fontsize, use it
            final_fontsize = ld_block_anno_kwargs['fontsize']
        else:
            # Use auto-calculated fontsize, but don't exceed user's initial fontsize
            final_fontsize = min(fontsize, auto_fontsize)
    else:
        final_fontsize = fontsize
    
    # Prepare default text kwargs
    default_text_kwargs = {
        'ha': 'right',  # Right alignment so right end is aligned
        'va': 'center',
        'fontsize': final_fontsize,
        'fontfamily': font_family,
        'color': color,
        'zorder': 10  # High zorder to appear on top
    }
    
    # Update with user-provided kwargs if any
    if ld_block_anno_kwargs is not None:
        default_text_kwargs.update(ld_block_anno_kwargs)
    
    # Annotate each variant along the left edge of the triangle (only those in annotate_mask)
    for i in range(n):
        if not annotate_mask[i]:
            continue  # Skip variants not in anno_set
        
        rank = ranks[i]
        text = annotation_texts[i] if i < len(annotation_texts) else f"?:{int(actual_positions[i])}"
        
        # Position on the left edge: X = 0, Y = rank
        # In transformed coordinates:
        # U = (0 + rank) / 2 = rank / 2
        # V_negated = -(rank - 0) / 2 = -rank / 2
        u_pos = rank / 2.0
        v_negated_pos = -rank / 2.0
        
        # Add offset to place text to the left of the edge
        # Offset perpendicular to the left edge (which has direction along U = -V_negated)
        # For a perpendicular offset, we offset in both U and V_negated
        # A vector perpendicular to (1, 1) is (1, -1) normalized
        # So offset: delta_U = -offset * cos(45°), delta_V = offset * sin(45°)
        # But simpler: just offset in U direction (to the left)
        u_text = u_pos + u_offset
        
        # Draw tick mark extending from left edge to the left (perpendicular to edge)
        # Tick extends along the perpendicular direction
        tick_end_u = u_pos + u_offset * 0.3  # Tick extends 30% of the way to text
        # Keep V_negated the same for horizontal tick (perpendicular to left edge)
        ax.plot([u_pos, tick_end_u], [v_negated_pos, v_negated_pos], 
                color=color, linewidth=0.8, zorder=9, clip_on=False)
        
        # Add text annotation (right-aligned so right end is aligned)
        ax.text(u_text, v_negated_pos, text, **default_text_kwargs)


def _annotate_ld_cells(
    ax: 'Axes',
    ld: np.ndarray,
    ranks: np.ndarray,
    ld_tri: np.ma.MaskedArray,
    anno_cell_fmt: str = "{:.2f}",
    fontsize: float = 8,
    font_family: str = "Arial",
    color: str = "black",
    anno_cell_kwargs: Optional[dict] = None
) -> None:
    """
    Add text annotations to LD cells showing the r² values.
    
    Parameters
    ----------
    ax : Axes
        Axes object for the LD block plot.
    ld : np.ndarray
        LD matrix (n x n).
    ranks : np.ndarray
        Variant rank array (length n).
    ld_tri : np.ma.MaskedArray
        Masked LD matrix (upper triangle).
    anno_cell_fmt : str
        Format string for cell annotations. Default: "{:.2f}".
    fontsize : float
        Font size for annotations. Default: 8.
    font_family : str
        Font family for annotations. Default: "Arial".
    color : str
        Text color. Default: "black".
    anno_cell_kwargs : dict, optional
        Additional keyword arguments for text annotations. Default: None.
        These will override default text styling (fontsize, fontfamily, color, etc.).
    """
    n = len(ranks)
    
    # Get edge coordinates (same as used in _plot_ld_triangle)
    e = _edges_from_centers(ranks)  # (n+1,)
    
    # Iterate over upper triangle (j >= i)
    for i in range(n):
        for j in range(i, n):
            # Skip if masked (lower triangle or NaN)
            if ld_tri.mask[i, j]:
                continue
            
            # Get LD value from the masked array
            ld_value = ld_tri[i, j]
            # Check if value is masked or NaN
            if isinstance(ld_value, np.ma.core.MaskedConstant) or np.isnan(ld_value):
                continue
            
            # Calculate cell center in rank coordinates
            # Cell spans from e[i] to e[i+1] in X, and e[j] to e[j+1] in Y
            x_center = (e[i] + e[i+1]) / 2.0
            y_center = (e[j] + e[j+1]) / 2.0
            
            # Transform to U, V_negated coordinates (same as in _plot_ld_triangle)
            u_center = (x_center + y_center) / 2.0
            v_negated_center = -(y_center - x_center) / 2.0
            
            # Format the text using the format string
            try:
                text = anno_cell_fmt.format(float(ld_value))
            except (ValueError, TypeError):
                # Fallback to default format if custom format fails
                text = f"{float(ld_value):.2f}"
            
            # Prepare default text kwargs
            default_text_kwargs = {
                'ha': 'center',
                'va': 'center',
                'fontsize': fontsize,
                'fontfamily': font_family,
                'color': color,
                'zorder': 10  # High zorder to appear on top
            }
            
            # Update with user-provided kwargs if any
            if anno_cell_kwargs is not None:
                default_text_kwargs.update(anno_cell_kwargs)
            
            # Add text annotation
            ax.text(u_center, v_negated_center, text, **default_text_kwargs)


def plot_ld_block(
    ld: Optional[np.ndarray] = None,
    pos: Optional[np.ndarray] = None,
    vcf_path: Optional[str] = None,
    region: Optional[Tuple[Union[int, str], int, int]] = None,
    sumstats: Optional[Union[pd.DataFrame, 'Sumstats']] = None,
    pos_col: str = "POS",
    nea_col: str = "NEA",
    ea_col: str = "EA",
    tabix: Optional[bool] = None,
    mapper: Optional[Any] = None,
    ax: Optional['Axes'] = None,
    ax_pos: Optional['Axes'] = None,
    mode: Optional[str] = None,
    cmap: Union[str, matplotlib.colors.Colormap, None] = None,
    vmin: float = 0.0,
    vmax: float = 1.0,
    xlabel: str = "Genomic position",
    title: Optional[str] = None,
    cbar: bool = True,
    cbar_label: str = "LD (r²)",
    cbar_kwargs: Optional[dict] = None,
    fig_kwargs: Optional[dict] = None,
    save: Optional[Union[bool, str]] = None,
    save_kwargs: Optional[dict] = None,
    fontsize: float = 10,
    font_family: str = "Arial",
    region_step: int = 21,
    lead_snp_is: Optional[list] = None,
    lead_snp_is_color: Optional[list] = None,
    anno_cell: bool = False,
    anno_cell_fmt: str = "{:.2f}",
    anno_cell_kwargs: Optional[dict] = None,
    ld_block_grid: bool = False,
    ld_block_grid_kwargs: Optional[dict] = None,
    ld_block_anno: Union[bool, str] = False,
    ld_block_anno_kwargs: Optional[dict] = None,
    ld_block_anno_set: Optional[Union[list, set, np.ndarray]] = None,
    ld_block_anno_max_rows: int = 100,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs
) -> Tuple['Figure', 'Axes']:
    """
    Plot the upper triangle of an LD matrix as a 45°-rotated inverted triangle.
    
    This function supports two modes:
    1. **Standalone mode**: Creates its own figure and plots LD block independently.
    2. **Regional mode**: Plots on provided axes (typically from plot_mqq) and aligns
       x-axis with regional plot using "i" coordinate system.
    
    Parameters
    ----------
    mode : str, optional
        Plot mode: 'regional' or 'standalone'. If None, automatically determines
        based on whether ax is provided and i-coordinates are available.
        Default: None (auto-detect).
    ld : np.ndarray, optional
        (n, n) array (e.g., r^2), symmetric LD matrix. If None, will be extracted
        from VCF file when vcf_path and region are provided.
    pos : np.ndarray, optional
        Length n positions for x-axis (bp or index). If None, uses np.arange(n) or
        extracts from sumstats when using VCF.
    vcf_path : str, optional
        Path to VCF file for LD calculation. If provided, region must also be provided.
    region : tuple, optional
        Region specification (chromosome, start, end) for extracting LD from VCF.
        Example: (1, 1000000, 2000000) for chr1:1000000-2000000.
    sumstats : pd.DataFrame or Sumstats, optional
        Sumstats object or DataFrame for matching variants with VCF. Required when
        using vcf_path. Should contain POS, NEA, and EA columns.
    pos_col : str, optional
        Column name for position in sumstats. Default: "POS".
    nea_col : str, optional
        Column name for non-effect allele in sumstats. Default: "NEA".
    ea_col : str, optional
        Column name for effect allele in sumstats. Default: "EA".
    tabix : bool, optional
        Whether to use tabix indexing for VCF. If None, auto-detects.
    mapper : ChromosomeMapper, optional
        ChromosomeMapper instance for chromosome name conversion. If None, creates
        a default mapper with automatic format detection.
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. If None, creates a new figure (standalone mode).
        If provided with region/sumstats/i-coordinates, uses regional mode.
    cmap : str or matplotlib.colors.Colormap, optional
        Colormap for the LD values. Default: continuous gradient matching region plot
        LD categories (dark blue -> light blue -> green -> orange -> red).
    vmin : float, optional
        Minimum value for colormap scaling. Default: 0.0.
    vmax : float, optional
        Maximum value for colormap scaling. Default: 1.0.
    xlabel : str, optional
        Label for x-axis. Default: "Genomic position".
    title : str, optional
        Title for the plot. Default: None (no title).
    cbar : bool, optional
        Whether to draw a colorbar. Default: True.
    cbar_label : str, optional
        Label for the colorbar. Default: "LD (r²)".
    cbar_kwargs : dict, optional
        Additional arguments for colorbar. Default: None.
    fig_kwargs : dict, optional
        Additional arguments for figure creation. Default: None.
    save : bool or str, optional
        Whether to save the figure. If str, path to save file. Default: None.
    save_kwargs : dict, optional
        Additional arguments for saving the figure. Default: None.
    fontsize : float, optional
        Font size for labels. Default: 12.
    font_family : str, optional
        Font family. Default: "Arial".
    region_step : int, optional
        Number of x-axis tick positions for regional mode. Default: 21.
    lead_snp_is : list, optional
        List of lead SNP i-coordinates for colored annotation lines. Default: None.
    lead_snp_is_color : list, optional
        List of colors corresponding to lead SNPs. Default: None.
    anno_cell : bool, optional
        Whether to annotate cells with LD r² values. Default: False.
    anno_cell_fmt : str, optional
        Format string for cell annotations (e.g., "{:.2f}" for 2 decimal places).
        Default: "{:.2f}".
    anno_cell_kwargs : dict, optional
        Additional keyword arguments for text annotations (e.g., {'fontsize': 10, 'weight': 'bold'}).
        Default: None.
    ld_block_anno : bool or str, optional
        Whether to add annotations on the left side of the LD block triangle.
        If True, uses 'chr:pos' format. If a string (e.g., 'SNPID', 'rsID'), uses that column for annotation text.
        Annotations are styled like axis ticks and labels, with right-aligned text.
        Default: False.
    ld_block_anno_kwargs : dict, optional
        Additional keyword arguments for left-side annotations (e.g., {'fontsize': 8, 'color': 'black'}).
        Default: None.
    ld_block_anno_set : list, set, or np.ndarray, optional
        List of SNPIDs to annotate. If None, annotates all variants. Variants are matched by SNPID column.
        Default: None.
    ld_block_anno_max_rows : int, optional
        Maximum number of variants to annotate. If the number of variants to annotate exceeds this limit,
        annotations will be skipped. Default: 100.
    log : Log, optional
        Logger instance. Default: Log().
    verbose : bool, optional
        Whether to print log messages. Default: True.
    **kwargs
        Additional arguments passed to pcolormesh.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes.Axes
        Axes object.
    """
    log.write("Start to create LD block plot (45° rotated inverted triangle)...", verbose=verbose)
 
    # ============================================================================
    # STEP 1: Prepare LD data
    # ============================================================================
    use_i_coordinate = False
    all_i_positions = None
    region_sumstats = None
    
    if vcf_path is not None and region is not None:
        # Extract LD from VCF
        ld_matrix, matched_sumstats, region_sumstats, use_i_coordinate, all_i_positions = \
            _prepare_ld_data_from_vcf(
                vcf_path=vcf_path,
                region=region,
                sumstats=sumstats,
                pos_col=pos_col,
                nea_col=nea_col,
                ea_col=ea_col,
                tabix=tabix,
                mapper=mapper,
                log=log,
                verbose=verbose
            )
        
        # Prepare full LD matrix with all variant positions
        ld, all_positions = _prepare_ld_matrix(
            ld_matrix=ld_matrix,
            matched_sumstats=matched_sumstats,
            region_sumstats=region_sumstats,
            pos_col=pos_col,
            log=log,
            verbose=verbose
        )
        
        # Store actual positions for position bar
        actual_positions = all_positions.copy()
    else:
        # If not using VCF, we need to get actual positions from pos parameter
        actual_positions = None
    
    # Validate that ld is provided
    if ld is None:
        raise ValueError("ld must be provided, or vcf_path and region must be provided")
    
    # Convert to numpy array
    ld = np.asarray(ld, dtype=float)
    n = ld.shape[0]
    
    if ld.shape != (n, n):
        raise ValueError(f"ld must be square, got shape {ld.shape}")
    
    # Use variant ranks (0, 1, 2, ..., n-1) for uniform cell sizes
    ranks = np.arange(n, dtype=float)
    log.write(f"Using variant ranks (0 to {n-1}) for uniform cell sizes", verbose=verbose)
    
    # Store actual positions if provided, otherwise use ranks as positions
    if actual_positions is None:
        if pos is not None:
            actual_positions = np.asarray(pos, dtype=float)
            if actual_positions.shape != (n,):
                raise ValueError(f"pos must have length n={n}, got length {len(actual_positions)}")
            log.write(f"Using provided positions for position bar (range: {actual_positions.min():.0f} to {actual_positions.max():.0f})", verbose=verbose)
        else:
            actual_positions = ranks.copy()
            log.write(f"Using ranks as positions for position bar", verbose=verbose)
    
    # ============================================================================
    # STEP 2: Determine plot mode
    # ============================================================================
    if mode is None:
        # Auto-detect mode: regional if ax is provided AND we have region/sumstats with i-coordinates
        if ax is not None and region is not None and sumstats is not None and use_i_coordinate:
            mode = 'regional'
        else:
            mode = 'standalone'
    elif mode not in ['regional', 'standalone']:
        raise ValueError(f"mode must be 'regional' or 'standalone', got '{mode}'")
    
    # ============================================================================
    # STEP 3: Set up figure and axes
    # ============================================================================
    style = set_plot_style(
        plot="plot_ld_block",
        fig_kwargs=fig_kwargs or {"figsize": (8, 4), "dpi": 200},  # Increased height for position bar
        save_kwargs=save_kwargs,
        fontsize=fontsize,
        fontfamily=font_family,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    fontsize = style["fontsize"]
    
    # Create figure/axes if not provided
    if ax is None:
        # Create figure with subplots: position bar on top, LD block on bottom
        fig = plt.figure(**fig_kwargs)
        # Create grid for subplots: 2 rows, 1 column
        # Top: position bar (height ratio ~0.2), Bottom: LD block (height ratio ~0.8)
        gs = fig.add_gridspec(2, 1, height_ratios=[1,14], hspace=.10)
        ax_pos = fig.add_subplot(gs[0])  # Position bar axes (top)
        ax = fig.add_subplot(gs[1])  # LD block axes (bottom)
        # Note: We don't share x-axis because position bar uses actual positions, LD block uses ranks
    else:
        fig = ax.figure
        # If ax is provided (regional mode), ax_pos should also be provided
        # Both axes are created by _process_layout in plot_mqq
    
    # Set default colormap to match region plot LD categories
    if cmap is None:
        # Create continuous gradient from region plot LD colors
        # Colors: dark blue -> light blue -> green -> orange -> red
        ld_colors = ["#020080", "#86CEF9", "#24FF02", "#FDA400", "#FF0000"]
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('ld_block', ld_colors, N=256)
    
    # ============================================================================
    # STEP 4: Plot LD triangle using ranks
    # ============================================================================
    # Always use ranks for the grid (uniform cells)
    # In regional mode, ax_ld_block does NOT share x-axis, so it uses rank-based coordinates
    m, ld_tri = _plot_ld_triangle(
        ax=ax,
        ld=ld,
        ranks=ranks,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        i_coords=None,  # Always use ranks, not i coordinates
        **kwargs
    )
    
    # ============================================================================
    # STEP 5: Set axis limits
    # ============================================================================
    if mode == 'regional':
        # In regional mode, ax_ld_block does NOT share x-axis
        # Grid uses ranks (uniform cells), x-axis is in rank space
        # Calculate edge coordinates from ranks
        e = _edges_from_centers(ranks)  # (n+1,) in rank space
        X_edges, Y_edges = np.meshgrid(e, e, indexing="ij")  # (n+1, n+1)
        U_edges = (X_edges + Y_edges) / 2.0
        V_edges = (Y_edges - X_edges) / 2.0
        V_edges_negated = -V_edges
        
        # For upper triangle, we only need the lower half (V_negated <= 0)
        u_min_edge = np.min(U_edges)
        u_max_edge = np.max(U_edges)
        v_neg_min_edge = np.min(V_edges_negated)
        v_neg_max_edge = 0.0
        
        # Set xlim based on rank-based U edges (no padding, as requested)
        ax.set_xlim(u_min_edge, u_max_edge)
        
        # Set ylim considering the 4/n part (upper half of a cell when rotated 45°)
        # The aspect ratio is 2:(1 + 4/n), so y-range = x-range / 2 * (1 + 4/n)
        x_range = u_max_edge - u_min_edge
        y_range = (x_range / 2.0) * (1.0 + 4.0 / n)
        
        # Center ylim around the V_negated range, but ensure it covers from v_neg_min_edge to 0
        y_center = v_neg_max_edge  # Center at 0 (the diagonal)
        ylim_top = y_center + y_range / 2.0
        ylim_bottom = y_center - y_range / 2.0
        
        # Ensure we cover the full range from v_neg_min_edge to 0
        if ylim_bottom > v_neg_min_edge:
            ylim_bottom = v_neg_min_edge
            # Adjust top to maintain y_range
            ylim_top = ylim_bottom + y_range
        
        ax.set_ylim(ylim_bottom, ylim_top)
    else:
        # Standalone mode: use _set_axis_limits
        _set_axis_limits(
            ax=ax,
            ranks=ranks,
            ld_tri=ld_tri,
            mode=mode
        )
    
    # ============================================================================
    # STEP 6: Draw grid if requested
    # ============================================================================
    if ld_block_grid:
        _draw_ld_block_grid(
            ax=ax,
            ranks=ranks,
            grid_kwargs=ld_block_grid_kwargs
        )
    
    # ============================================================================
    # STEP 7: Add cell annotations if requested
    # ============================================================================
    if anno_cell:
        _annotate_ld_cells(
            ax=ax,
            ld=ld,
            ranks=ranks,
            ld_tri=ld_tri,
            anno_cell_fmt=anno_cell_fmt,
            fontsize=fontsize * 0.7,  # Slightly smaller than default fontsize
            font_family=font_family,
            color="black",
            anno_cell_kwargs=anno_cell_kwargs
        )
    
    # ============================================================================
    # STEP 7.5: Add left-side annotations if requested
    # ============================================================================
    if ld_block_anno:
        # Count how many variants would be annotated
        n_variants_to_annotate = n  # Default: all variants
        if ld_block_anno_set is not None:
            # Count variants that match the anno_set
            if region_sumstats is not None and "SNPID" in region_sumstats.columns:
                anno_set_set = set(ld_block_anno_set) if not isinstance(ld_block_anno_set, set) else ld_block_anno_set
                # Match by position to get SNPIDs
                if pos_col in region_sumstats.columns:
                    pos_to_snpid = {}
                    for idx, row in region_sumstats.iterrows():
                        if pos_col in row and "SNPID" in row and pd.notna(row[pos_col]) and pd.notna(row["SNPID"]):
                            pos_to_snpid[int(row[pos_col])] = row["SNPID"]
                    
                    n_variants_to_annotate = sum(1 for pos in actual_positions 
                                                  if int(pos) in pos_to_snpid 
                                                  and pos_to_snpid[int(pos)] in anno_set_set)
                else:
                    n_variants_to_annotate = 0
            else:
                n_variants_to_annotate = 0
        
        # Check if number of variants exceeds the limit
        if n_variants_to_annotate > ld_block_anno_max_rows:
            log.write(f"Skipping LD block annotations: {n_variants_to_annotate} variants exceed maximum limit of {ld_block_anno_max_rows}", verbose=verbose)
        else:
            # Extract column name from ld_block_anno if it's a string, otherwise use None (chr:pos format)
            anno_col = ld_block_anno if isinstance(ld_block_anno, str) else None
            
            _annotate_ld_block_left_side(
                ax=ax,
                ranks=ranks,
                actual_positions=actual_positions,
                region=region,
                region_sumstats=region_sumstats,
                pos_col=pos_col,
                snpid_col="SNPID",  # Default SNPID column name
                anno_col=anno_col,
                anno_set=ld_block_anno_set,
                fontsize=fontsize * 0.8,  # Slightly smaller than default fontsize
                font_family=font_family,
                color="black",
                ld_block_anno_kwargs=ld_block_anno_kwargs
            )
    
    # ============================================================================
    # STEP 8: Apply cosmetics
    # ============================================================================
    # Remove spine lines
    for spine in ax.spines.values():
        spine.set_visible(False)
    

    # For standalone mode, adjust figure to fill properly
    if mode == 'standalone':
        # Calculate edge coordinates from ranks (same as used in pcolormesh)
        # This ensures we use the actual converted data coordinates
        e = _edges_from_centers(ranks)  # (n+1,)
        X_edges, Y_edges = np.meshgrid(e, e, indexing="ij")  # (n+1, n+1)
        U_edges = (X_edges + Y_edges) / 2.0
        V_edges = (Y_edges - X_edges) / 2.0
        V_edges_negated = -V_edges
        
        # Get the full range of edge coordinates (these define the cell boundaries)
        # Using edge coordinates ensures complete edge cells
        # For upper triangle, we only need the lower half (V_negated <= 0)
        u_min_edge = np.min(U_edges)
        u_max_edge = np.max(U_edges)
        v_neg_min_edge = np.min(V_edges_negated)
        # For upper triangle, maximum V_negated is 0 (at the diagonal)
        v_neg_max_edge = 0.0
        
        # Set figure size: aspect ratio for ax = 2:(1 + 4/n)
        # Account for position bar (height_ratios=[1, 14])
        ax_aspect_ratio = 2.0 / (1.0 + 4.0 / n)
        current_size = fig.get_size_inches()
        fig_width = current_size[0]
        # LD block gets 14/15 of height, so adjust for that
        fig_height = (fig_width / ax_aspect_ratio) * (15.0 / 14.0) * 1.1  # 1.1 for colorbar/labels
        fig.set_size_inches(fig_width, fig_height)
        
        # Use tight_layout
        fig.tight_layout()
        
        # Set xlim from edge coordinates (no padding)
        ax.set_xlim(u_min_edge, u_max_edge)
        
        # Set ylim considering the 4/n part (upper half of a cell when rotated 45°)
        # The aspect ratio is 2:(1 + 4/n), so y-range = x-range / 2 * (1 + 4/n)
        xlim = ax.get_xlim()
        x_range = xlim[1] - xlim[0]
        y_range = (x_range / 2.0) * (1.0 + 4.0 / n)
        
        # Center ylim around the V_negated range, but ensure it covers from v_neg_min_edge to 0
        # The 4/n part represents the upper half, so we extend upward from v_neg_max_edge (0)
        y_center = v_neg_max_edge  # Center at 0 (the diagonal)
        ylim_top = y_center + y_range / 2.0
        ylim_bottom = y_center - y_range / 2.0
        
        # Ensure we cover the full range from v_neg_min_edge to 0
        if ylim_bottom > v_neg_min_edge:
            ylim_bottom = v_neg_min_edge
            # Adjust top to maintain y_range
            ylim_top = ylim_bottom + y_range
        
        ax.set_ylim(ylim_bottom, ylim_top)
    
    if title is not None:
        ax.set_title(title, fontsize=fontsize, fontfamily=font_family)
    
    # ============================================================================
    # STEP 6.5: Add position bar and annotation lines
    # ============================================================================
    if ax_pos is not None:
        if mode == 'standalone':
            log.write("Adding position bar and annotation lines (standalone mode)...", verbose=verbose)
            # Plot position bar (on top, using actual positions)
            _plot_position_bar(
                ax_pos=ax_pos,
                actual_positions=actual_positions,
                ranks=ranks,
                fontsize=fontsize,
                font_family=font_family,
                show_xticks=True,
                region_step=region_step,
                xlabel=xlabel,
                region=region
            )
            
            # Add annotation lines for all variants
            _add_annotation_lines(
                ax_ld=ax,
                ax_pos=ax_pos,
                ranks=ranks,
                actual_positions=actual_positions,
                n_lines=None,  # None means annotate all variants
                line_color='gray',
                line_alpha=0.3,
                line_style='--',
                line_width=0.5,
                lead_snp_is=lead_snp_is,
                lead_snp_is_color=lead_snp_is_color
            )
        elif mode == 'regional':
            log.write("Adding position bar (regional mode, using i coordinates)...", verbose=verbose)
            # In regional mode, position bar should use i coordinates to align with regional plot
            # Use i coordinates if available, otherwise use actual positions
            if use_i_coordinate and all_i_positions is not None:
                pos_for_bar = all_i_positions
            else:
                pos_for_bar = actual_positions
            
            # Plot position bar using i coordinates
            _plot_position_bar(
                ax_pos=ax_pos,
                actual_positions=pos_for_bar,
                ranks=ranks,
                fontsize=fontsize,
                font_family=font_family,
                show_xticks=False,
                region_step=region_step,
                xlabel=xlabel,
                region=region
            )
            
            # Add annotation lines
            _add_annotation_lines(
                ax_ld=ax,
                ax_pos=ax_pos,
                ranks=ranks,
                actual_positions=pos_for_bar,
                n_lines=None,
                line_color='gray',
                line_alpha=0.3,
                line_style='--',
                line_width=0.5,
                lead_snp_is=lead_snp_is,
                lead_snp_is_color=lead_snp_is_color
            )
    
    # ============================================================================
    # STEP 7: Add colorbar
    # ============================================================================
    if cbar:
        if cbar_kwargs is None:
            cbar_kwargs = {}
        # Use inset_axes to place colorbar at lower right corner
        # Position: lower right, inside axes coordinates (0-1 range)
        # For inverted triangle, the lower right area has space
        # bbox_to_anchor must be 4-tuple when using relative units for width/height
        cax = inset_axes(ax, width="25%", height="5%", loc='lower right', 
                        bbox_to_anchor=(0.0, 0.05, 1, 1), bbox_transform=ax.transAxes,
                        borderpad=0)
        cbar_obj = plt.colorbar(m, cax=cax, orientation='horizontal', **cbar_kwargs)
        # Set label and ticks on top of horizontal colorbar
        cbar_obj.ax.xaxis.set_label_position('top')
        cbar_obj.ax.xaxis.set_ticks_position('top')
        cbar_obj.ax.tick_params(labeltop=True, labelbottom=False)
        cbar_obj.set_label(cbar_label, fontsize=fontsize, fontfamily=font_family)
    
    log.write("Finished creating LD block plot!", verbose=verbose)
    
    # Remove y-axis ticks
    ax.set_yticks([])
    
    # Set x ticks and labels based on mode
    if mode == 'regional':
        # In regional mode, ax_ld_block does NOT share x-axis, no ticks or labels
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labeltop=False)
    else:
        # Standalone mode: remove x-axis labels (position bar will show them)
        ax.set_xticks([])
        ax.set_xlabel("")
    
    # ============================================================================
    # STEP 8: Save figure (only in standalone mode)
    # ============================================================================
    if mode == 'standalone':
        save_figure(fig, save, keyword="ld_block", save_kwargs=save_kwargs, log=log, verbose=verbose)
    
    return fig, ax
