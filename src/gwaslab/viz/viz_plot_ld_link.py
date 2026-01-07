from typing import TYPE_CHECKING, Optional, Union, Tuple, Any
import numpy as np
import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_plot_ld_block import (
    _prepare_ld_data_from_vcf,
    _prepare_ld_matrix
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from gwaslab.g_Sumstats import Sumstats


def _plot_ld_link(
    ax,
    vcf_path,
    region,
    sumstats,
    pos_col="POS",
    region_ld_threshold=[0.2, 0.4, 0.6, 0.8],
    region_ld_colors=None,
    palette=None,
    link_alpha_scale=0.2,
    link_linewidth=1.0,
    sig_level=None,
    log=Log(),
    verbose=True
):
    """
    Internal function to plot LD lines on regional plot axes.
    Draws straight lines connecting variant pairs with LD above minimum threshold.
    Lines connect the actual variant positions (i, scaled_P).
    Line colors are assigned based on LD categories using region_ld_threshold and region_ld_colors.
    Only draws lines, does not modify axes labels or other properties.
    
    Parameters
    ----------
    region_ld_threshold : list
        LD r² thresholds for color categories, e.g., [0.2, 0.4, 0.6, 0.8]
    region_ld_colors : list
        Colors for each LD category. Should have length len(region_ld_threshold) + 3:
        [no_data_color, 0<thr[0], thr[0]<thr[1], ..., thr[-1]<1.0, lead_color]
    """
    log.write("Adding LD link plot...", verbose=verbose)
    
    # Use region_ld_colors if provided, otherwise extract from palette, otherwise use default
    if region_ld_colors is not None and isinstance(region_ld_colors, list):
        region_ld_colors_for_link = region_ld_colors
    elif palette is not None and isinstance(palette, dict):
        # Extract colors from palette dictionary
        # For single reference: palette keys are 100+i, values are region_ld_colors[i]
        # We need to extract colors in order: palette[100], palette[101], palette[102], ...
        sorted_keys = sorted([k for k in palette.keys() if isinstance(k, (int, float)) and k >= 100 and k < 200])
        if sorted_keys and len(sorted_keys) >= len(region_ld_threshold) + 3:
            # Extract colors in order, mapping palette keys 100+i to list index i
            region_ld_colors_for_link = [palette[k] for k in sorted_keys]
        else:
            region_ld_colors_for_link = ["#E4E4E4", "#020080", "#86CEF9", "#24FF02", "#FDA400", "#FF0000", "#FF0000"]
    else:
        region_ld_colors_for_link = ["#E4E4E4", "#020080", "#86CEF9", "#24FF02", "#FDA400", "#FF0000", "#FF0000"]
    
    # Ensure we have enough colors (at least len(region_ld_threshold) + 3)
    min_colors_needed = len(region_ld_threshold) + 3
    if len(region_ld_colors_for_link) < min_colors_needed:
        # Pad with last color if needed
        region_ld_colors_for_link = region_ld_colors_for_link + [region_ld_colors_for_link[-1]] * (min_colors_needed - len(region_ld_colors_for_link))
    
    # Use the processed colors
    region_ld_colors = region_ld_colors_for_link
    
    # Extract LD matrix from VCF
    ld_matrix, matched_sumstats, region_sumstats, use_i_coord, all_i_pos = _prepare_ld_data_from_vcf(
        vcf_path=vcf_path,
        region=region,
        sumstats=sumstats,
        pos_col=pos_col,
        nea_col="NEA",
        ea_col="EA",
        tabix=None,
        mapper=None,
        log=log,
        verbose=verbose
    )
    
    # Prepare full LD matrix
    full_ld_matrix, all_positions = _prepare_ld_matrix(
        ld_matrix=ld_matrix,
        matched_sumstats=matched_sumstats,
        region_sumstats=region_sumstats,
        pos_col=pos_col,
        log=log,
        verbose=verbose
    )
    
    # Use i coordinates and scaled_P for line coordinates
    region_sumstats_sorted = region_sumstats.sort_values(pos_col)
    
    # Get x coordinates (i if available, otherwise positions)
    if use_i_coord and all_i_pos is not None:
        if "i" in region_sumstats_sorted.columns:
            x_pos = region_sumstats_sorted["i"].values.astype(float)
        else:
            x_pos = all_positions
    else:
        x_pos = all_positions
    
    # Get y coordinates (scaled_P)
    if "scaled_P" in region_sumstats_sorted.columns:
        y_pos = region_sumstats_sorted["scaled_P"].values.astype(float)
    else:
        log.warning("scaled_P column not found. Cannot draw lines connecting variant positions.")
        return
    
    ld = full_ld_matrix
    ld = np.asarray(ld)
    x_pos = np.asarray(x_pos)
    y_pos = np.asarray(y_pos)
    
    if ld.shape[0] != ld.shape[1]:
        raise ValueError(f"LD matrix must be square, got shape {ld.shape}")
    if len(x_pos) != ld.shape[0]:
        raise ValueError(f"Position array length ({len(x_pos)}) must match LD matrix size ({ld.shape[0]})")
    if len(y_pos) != ld.shape[0]:
        raise ValueError(f"scaled_P array length ({len(y_pos)}) must match LD matrix size ({ld.shape[0]})")
    
    # Get significance information if sig_level is provided
    sig_mask = None
    if sig_level is not None:
        # Check if we have scaled_P or P column
        if "scaled_P" in region_sumstats_sorted.columns:
            # scaled_P is -log10(P), so check if scaled_P >= -log10(sig_level)
            sig_threshold = -np.log10(sig_level)
            sig_mask = region_sumstats_sorted["scaled_P"].values >= sig_threshold
            log.write(f"Using significance threshold: P <= {sig_level} (scaled_P >= {sig_threshold:.2f})", verbose=verbose)
        elif "P" in region_sumstats_sorted.columns:
            # Use P column directly
            sig_mask = region_sumstats_sorted["P"].values <= sig_level
            log.write(f"Using significance threshold: P <= {sig_level}", verbose=verbose)
        else:
            log.warning("Significance threshold specified but no P or scaled_P column found. Ignoring sig_level.")
            sig_mask = None
    
    n = len(x_pos)
    pos_range = x_pos[-1] - x_pos[0] if n > 1 else 1
    
    # Minimum threshold: use the first threshold as minimum (typically 0.2)
    min_thr = region_ld_threshold[0] if len(region_ld_threshold) > 0 else 0.0
    
    # Function to get color index based on LD value
    # Matches the logic in process_vcf:
    # LD = 0 or NaN: index 0
    # 0 < LD <= 0.2: index 1
    # 0.2 < LD <= 0.4: index 2
    # 0.4 < LD <= 0.6: index 3
    # 0.6 < LD <= 0.8: index 4
    # 0.8 < LD <= 1.0: index 5
    # Lead variant: index 6 (not used for links)
    def get_ld_color_index(ld_value):
        """Get color index for LD value based on thresholds.
        Returns index into region_ld_colors array.
        Index 0: no data/LD=0
        Index 1: 0 < LD <= threshold[0]
        Index 2: threshold[0] < LD <= threshold[1]
        ...
        Index len(thresholds)+1: threshold[-1] < LD <= 1.0
        """
        if np.isnan(ld_value) or ld_value <= 0:
            return 0
        # Check thresholds in reverse order to find the highest threshold that LD exceeds
        # If LD > threshold[i], it belongs to category i+2
        for idx in range(len(region_ld_threshold) - 1, -1, -1):
            threshold = region_ld_threshold[idx]
            if ld_value > threshold:
                return idx + 2
        # LD <= first threshold, so it belongs to category 1
        return 1
    
    # Plot lines for pairs with LD >= minimum threshold
    line_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not np.isnan(ld[i, j]) and ld[i, j] >= min_thr:
                # Check significance threshold if provided
                if sig_mask is not None:
                    # Only draw line if at least one variant in the pair is significant
                    if not (sig_mask[i] or sig_mask[j]):
                        continue
                x1, x2 = x_pos[i], x_pos[j]
                y1, y2 = y_pos[i], y_pos[j]
                
                # Use fixed line width
                lw = link_linewidth
                
                # Alpha based on LD value
                alpha = min(ld[i, j] * link_alpha_scale, 1.0)
                
                # Get color based on LD category
                color_idx = get_ld_color_index(ld[i, j])
                if color_idx < len(region_ld_colors):
                    line_color = region_ld_colors[color_idx]
                else:
                    # Fallback to last color if index out of range
                    line_color = region_ld_colors[-1]
                
                # Draw straight line connecting the two variant positions
                ax.plot([x1, x2], [y1, y2], color=line_color, alpha=alpha, linewidth=lw, zorder=1)
                line_count += 1
    
    if sig_level is not None:
        log.write(f"Plotted {line_count} lines with LD >= {min_thr} and at least one variant with P <= {sig_level}", verbose=verbose)
    else:
        log.write(f"Plotted {line_count} lines with LD >= {min_thr}", verbose=verbose)

