"""
X-Axis Alignment Manager for aligning x-axes across multiple panels.

This module provides the XAxisManager class for aligning x-axes in terms of
data range (xlim) and tick positions/labels without using matplotlib's sharex.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional, Tuple, Union, Any
from gwaslab.info.g_Log import Log

# Try to import for inset axes detection
try:
    from mpl_toolkits.axes_grid1.parasite_axes import AxesHostAxes
except ImportError:
    AxesHostAxes = None


class XAxisManager:
    """
    Manager for aligning x-axes across multiple matplotlib axes.
    
    This class takes all axes in the figure and applies:
    - X-tick and label alignment to specified axes only
    - Vertical spacing adjustment to all axes
    
    Parameters
    ----------
    all_axes : list of matplotlib.axes.Axes
        List of ALL axes in the figure (preserves order from top to bottom).
        all_axes[0] is the top panel, all_axes[-1] is the bottom panel.
        Spacing adjustment is applied to all axes in this list.
    xlim : tuple of (float, float), optional
        Explicit x-axis limits (xmin, xmax). If None, will be calculated
        from region parameter.
    xticks : array-like, optional
        Explicit tick positions. If None, will be calculated from region
        and region_step.
    xticklabels : list of str, optional
        Explicit tick labels. If None, will be auto-generated or calculated
        from region.
    region : tuple of (int, int, int), optional
        Genomic region as (chromosome, start, end) in 1-based coordinates.
        Used to calculate xlim and ticks if not explicitly provided.
    region_step : int, optional
        Number of ticks to generate. Used with region to calculate tick
        positions using np.linspace.
    track_start_i : float, default=0.0
        Offset for track-based plots (track, arc). The xlim will be calculated
        as (track_start_i + region[1], track_start_i + region[2]).
    gene_track_start_i : float, optional
        Offset for gene track in regional plots. If None, uses track_start_i.
    xlabel : str, optional
        X-axis label to set on aligned axes. If None, preserves existing labels.
    fontsize : float, optional
        Font size for tick labels and xlabel.
    font_family : str, optional
        Font family for tick labels and xlabel.
    fig : matplotlib.figure.Figure, optional
        Figure object containing the axes. If None, will be extracted from
        the first axis in all_axes. Required for adjusting vertical spacing.
    adjust_spacing : bool, default=True
        Whether to automatically adjust vertical spacing to prevent
        overlapping tick labels and ensure proper spacing.
    min_hspace : float, default=0.05
        Minimum vertical spacing between subplots (as fraction of figure height).
    verbose : bool, default=True
        Whether to show progress messages.
    log : Log, optional
        Logger instance for messages.
    
    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from gwaslab.viz.viz_aux_xaxis_manager import XAxisManager
    >>> 
    >>> fig, axes = plt.subplots(3, 1)
    >>> # ... plot on axes ...
    >>> 
    >>> # Align specified axes, apply spacing to all axes
    >>> xm = XAxisManager(
    ...     all_axes=list(fig.get_axes()),  # All axes for spacing
    ...     region=(1, 1000000, 2000000),
    ...     region_step=21,
    ...     track_start_i=0.0
    ... )
    >>> xm.register_align(axes[:2])  # Only align first 2 axes
    >>> xm.align()
    """
    
    def __init__(
        self,
        all_axes: List[plt.Axes],
        xlim: Optional[Tuple[float, float]] = None,
        xticks: Optional[np.ndarray] = None,
        xticklabels: Optional[List[str]] = None,
        region: Optional[Tuple[int, int, int]] = None,
        region_step: Optional[int] = None,
        track_start_i: float = 0.0,
        gene_track_start_i: Optional[float] = None,
        xlabel: Optional[str] = None,
        fontsize: Optional[float] = None,
        font_family: Optional[str] = None,
        fig: Optional[plt.Figure] = None,
        adjust_spacing: bool = True,
        min_hspace: float = 0.05,
        verbose: bool = True,
        log: Optional[Log] = None
    ):
        # All axes in the figure (for spacing adjustment)
        # all_axes[0] is top panel, all_axes[-1] is bottom panel
        self.all_axes = list(all_axes) if all_axes is not None else []
        # X-axis alignment parameters
        self.xlim = xlim
        self.xticks = xticks
        self.xticklabels = xticklabels
        self.region = region
        self.region_step = region_step
        self.track_start_i = track_start_i
        self.gene_track_start_i = gene_track_start_i if gene_track_start_i is not None else track_start_i
        self.xlabel = xlabel
        self.fontsize = fontsize
        self.font_family = font_family
        
        # Figure and spacing parameters
        self.fig = fig
        if self.fig is None and len(self.all_axes) > 0:
            self.fig = self.all_axes[0].figure
        self.adjust_spacing = adjust_spacing
        self.min_hspace = min_hspace
        self.verbose = verbose
        self.log = log if log is not None else Log()
        
        # Axes to align (x-tick and label alignment only)
        # Spacing is applied to all_axes, but x-tick/label alignment is only for axes_to_align
        self.axes_to_align: List[plt.Axes] = []
    
    def register_align(self, ax: Union[List[plt.Axes], plt.Axes]) -> None:
        """
        Register axis(es) for x-tick and label alignment.
        
        These axes will have their xlim, xticks, and xticklabels aligned.
        Spacing is still applied to all axes in all_axes.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes or list of matplotlib.axes.Axes
            Axis(es) to register for x-tick/label alignment.
        """
        if isinstance(ax, plt.Axes):
            ax = [ax]
        
        for a in ax:
            if a not in self.axes_to_align:
                self.axes_to_align.append(a)
    
    def register_align_many(self, axes: Union[List[plt.Axes], plt.Axes]) -> None:
        """
        Register multiple axes for x-tick and label alignment.
        
        Alias for register_align for backward compatibility.
        
        Parameters
        ----------
        axes : list of matplotlib.axes.Axes or matplotlib.axes.Axes
            Axes to register for x-tick/label alignment.
        """
        self.register_align(axes)
    
    def _calculate_xlim(self) -> Tuple[float, float]:
        """Calculate xlim from region or use provided xlim."""
        if self.xlim is not None:
            return self.xlim
        
        if self.region is None:
            raise ValueError(
                "Either xlim or region must be provided to calculate x-axis limits."
            )
        
        # Calculate xlim based on region
        # For most panels, use track_start_i offset
        xmin = self.track_start_i + self.region[1]
        xmax = self.track_start_i + self.region[2]
        
        return (xmin, xmax)
    
    def _calculate_xticks(self) -> np.ndarray:
        """Calculate tick positions from region or use provided xticks."""
        if self.xticks is not None:
            return np.asarray(self.xticks)
        
        if self.region is None or self.region_step is None:
            # If no region/region_step, return empty array (let matplotlib auto-generate)
            return np.array([])
        
        # Calculate ticks based on region
        xmin, xmax = self._calculate_xlim()
        ticks = np.linspace(xmin, xmax, num=self.region_step)
        
        return ticks
    
    def _calculate_xticklabels(self) -> Optional[List[str]]:
        """Calculate tick labels from region or use provided labels."""
        if self.xticklabels is not None:
            return self.xticklabels
        
        if self.region is None or self.region_step is None:
            return None
        
        # Generate tick labels in MB format (similar to track and regional plots)
        # Calculate positions in base pairs
        tick_positions_bp = np.linspace(self.region[1], self.region[2], num=self.region_step).astype(int)
        # Format as MB with 3 decimal places
        tick_labels = [f'{pos/1000000:.3f}' for pos in tick_positions_bp]
        
        return tick_labels
    
    def _calculate_xlabel(self) -> Optional[str]:
        """Calculate xlabel from region or use provided label."""
        if self.xlabel is not None:
            return self.xlabel
        
        if self.region is None:
            return None
        
        # Generate xlabel in format "Chromosome X (MB)" (similar to track plots)
        xlabel = f"Chromosome {self.region[0]} (MB)"
        return xlabel
    
    def _get_tick_position(self, ax: plt.Axes) -> str:
        """
        Determine if x-axis ticks are on top or bottom.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis to check.
        
        Returns
        -------
        str
            'top' if ticks are on top, 'bottom' if ticks are on bottom.
        """
        # Check tick_params to see where labels are shown
        try:
            # Check if labeltop is True (labels on top)
            # Access tick_params through the axis
            tick_params = ax.xaxis._major_tick_kw
            if tick_params.get('labeltop', False):
                return 'top'
        except:
            pass
        
        # Check using get_ticks_position (returns 'top', 'bottom', 'both', or 'default')
        try:
            tick_position = ax.xaxis.get_ticks_position()
            if tick_position == 'top':
                return 'top'
            elif tick_position == 'both':
                # If both, check which labels are actually visible
                # Check if labeltop is True
                try:
                    if ax.xaxis._major_tick_kw.get('labeltop', False):
                        return 'top'
                except:
                    pass
        except:
            pass
        
        # Check if any visible tick labels are positioned above the axis
        try:
            tick_labels = ax.get_xticklabels()
            ylim = ax.get_ylim()
            axis_center = (ylim[0] + ylim[1]) / 2
            for label in tick_labels:
                if label.get_visible() and label.get_text().strip():
                    try:
                        label_pos = label.get_position()
                        # If label y position is above axis center, it's on top
                        if label_pos[1] > axis_center:
                            return 'top'
                    except:
                        pass
        except:
            pass
        
        # Default to bottom
        return 'bottom'
    
    def _calculate_required_space_for_labels(self, ax: plt.Axes, direction: str = 'below') -> float:
        """
        Calculate the space required for tick labels and xlabel.
        
        This method accounts for:
        - Rotation angle of tick labels
        - Actual bounding boxes of visible labels
        - X-axis label position
        - Tick position (top or bottom)
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis to check for label space requirements.
        direction : str, default='below'
            'below' to calculate space below axis, 'above' to calculate space above axis.
        
        Returns
        -------
        float
            Required space in figure coordinates (0-1 range).
        """
        if self.fig is None:
            self.fig = ax.figure
        
        fig_height = self.fig.get_figheight()
        pos = ax.get_position()
        required_space = 0.0
        
        # Determine actual tick position
        tick_position = self._get_tick_position(ax)
        
        # Get tick labels and check if they're actually visible and have text
        tick_labels = ax.get_xticklabels()
        has_visible_labels = False
        
        # Check if there are any visible labels with text
        for label in tick_labels:
            if label.get_visible() and label.get_text().strip():
                has_visible_labels = True
                break
        
        if has_visible_labels:
            # Get rotation angle (default is 45 degrees for our use case)
            rotation = 45  # Default rotation
            if len(tick_labels) > 0:
                # Try to get actual rotation from first label
                try:
                    rotation = tick_labels[0].get_rotation()
                except:
                    rotation = 45
            
            # Try to get actual bounding boxes if renderer is available
            renderer = self.fig.canvas.get_renderer()
            if renderer is not None:
                if direction == 'below' or (direction == 'auto' and tick_position == 'bottom'):
                    # Check space below axis
                    max_bottom = pos.y0  # Start with axis bottom
                    
                    for label in tick_labels:
                        if label.get_visible() and label.get_text().strip():
                            # Get bounding box in figure coordinates
                            bbox = label.get_window_extent(renderer=renderer)
                            bbox_fig = bbox.transformed(self.fig.transFigure.inverted())
                            # Check if this label extends below the axis
                            if bbox_fig.y0 < pos.y0:
                                max_bottom = min(max_bottom, bbox_fig.y0)
                    
                    # Calculate required space: distance from axis bottom to lowest label
                    if max_bottom < pos.y0:
                        required_space = (pos.y0 - max_bottom) * 1.3 + 0.01  # Add padding
                else:
                    # Check space above axis
                    min_top = pos.y1  # Start with axis top
                    
                    for label in tick_labels:
                        if label.get_visible() and label.get_text().strip():
                            # Get bounding box in figure coordinates
                            bbox = label.get_window_extent(renderer=renderer)
                            bbox_fig = bbox.transformed(self.fig.transFigure.inverted())
                            # Check if this label extends above the axis
                            if bbox_fig.y1 > pos.y1:
                                min_top = max(min_top, bbox_fig.y1)
                    
                    # Calculate required space: distance from axis top to highest label
                    if min_top > pos.y1:
                        required_space = (min_top - pos.y1) * 1.3 + 0.01  # Add padding
            else:
                # Fallback: estimate based on font size and rotation
                if self.fontsize is not None:
                    label_fontsize = self.fontsize
                else:
                    # Try to get fontsize from label
                    try:
                        label_fontsize = tick_labels[0].get_fontsize()
                    except:
                        label_fontsize = 10  # Default
                
                # Calculate vertical extent based on rotation
                # For rotated labels, the vertical extent depends on the rotation angle
                label_height_points = label_fontsize * 1.2
                label_height_fig = (label_height_points / 72.0) / fig_height
                
                # Convert rotation to radians
                rotation_rad = np.deg2rad(rotation)
                # Vertical extent = label_height * sin(rotation) for labels rotated
                # For 45 degrees: sin(45) ≈ 0.707
                vertical_extent = label_height_fig * abs(np.sin(rotation_rad))
                
                # Add padding
                padding = 0.02
                required_space = vertical_extent * 1.3 + padding
        
        # Add space for xlabel if present
        xlabel = ax.get_xlabel()
        if xlabel and xlabel.strip():
            # Check xlabel position
            try:
                xlabel_position = ax.xaxis.label.get_position()[1]
                # If xlabel is above axis center, it's on top
                xlabel_on_top = xlabel_position > (pos.y0 + pos.y1) / 2
            except:
                # Default based on tick position
                xlabel_on_top = (tick_position == 'top')
            
            # Get xlabel fontsize
            if self.fontsize is not None:
                xlabel_fontsize = self.fontsize
            else:
                # Try to get from axis
                try:
                    xlabel_obj = ax.xaxis.label
                    xlabel_fontsize = xlabel_obj.get_fontsize()
                except:
                    xlabel_fontsize = 10
            
            # Calculate xlabel height
            xlabel_height_points = xlabel_fontsize * 1.2
            xlabel_height_fig = (xlabel_height_points / 72.0) / fig_height
            
            # Add padding between tick labels and xlabel
            xlabel_padding = 0.01
            
            # Add xlabel space in the correct direction
            if (direction == 'below' or (direction == 'auto' and tick_position == 'bottom')) and not xlabel_on_top:
                required_space += xlabel_height_fig + xlabel_padding
            elif (direction == 'above' or (direction == 'auto' and tick_position == 'top')) and xlabel_on_top:
                required_space += xlabel_height_fig + xlabel_padding
        
        return required_space
    
    def _has_labels(self, ax: plt.Axes) -> bool:
        """
        Check if an axis has visible tick labels or xlabel.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axis to check.
        
        Returns
        -------
        bool
            True if axis has visible labels, False otherwise.
        """
        # Check for tick labels
        tick_labels = ax.get_xticklabels()
        for label in tick_labels:
            if label.get_visible() and label.get_text().strip():
                return True
        
        # Check for xlabel
        xlabel = ax.get_xlabel()
        if xlabel and xlabel.strip():
            return True
        
        return False
    
    def _separate_overlapping_axes(self, axes_for_spacing: List[plt.Axes], compact_spacing: float) -> None:
        """
        Separate overlapping adjacent axes by moving them apart.
        
        Parameters
        ----------
        axes_for_spacing : List[plt.Axes]
            List of axes to process (preserves order from top to bottom).
        compact_spacing : float
            Minimum spacing between axes.
        """
        max_iterations = 10
        for _ in range(max_iterations):
            has_overlaps = False
            for i in range(len(axes_for_spacing) - 1):
                top_ax = axes_for_spacing[i]
                bottom_ax = axes_for_spacing[i + 1]
                pos_top = top_ax.get_position()
                pos_bottom = bottom_ax.get_position()
                
                gap = pos_top.y1 - pos_bottom.y0
                
                if gap < compact_spacing:
                    has_overlaps = True
                    target_bottom_y0 = max(0.0, pos_top.y1 - compact_spacing - pos_bottom.height)
                    movement = target_bottom_y0 - pos_bottom.y0
                    
                    # Move bottom axis and all axes below it down
                    for j in range(i + 1, len(axes_for_spacing)):
                        ax_to_move = axes_for_spacing[j]
                        pos = ax_to_move.get_position()
                        new_y0 = max(0.0, min(1.0 - pos.height, pos.y0 + movement))
                        ax_to_move.set_position([pos.x0, new_y0, pos.width, pos.height])
                    break
            
            if not has_overlaps:
                break
    
    def _calculate_required_gap(self, top_ax: plt.Axes, bottom_ax: plt.Axes, 
                                 axes_with_labels_set: set, compact_spacing: float) -> float:
        """
        Calculate the required gap between two axes based on their tick positions.
        
        Parameters
        ----------
        top_ax : plt.Axes
            Top axis in the gap.
        bottom_ax : plt.Axes
            Bottom axis in the gap.
        axes_with_labels_set : set
            Set of axes that have labels (for spacing calculations).
        compact_spacing : float
            Minimum spacing when no labels are present.
        
        Returns
        -------
        float
            Required gap between the two axes.
        """
        required_gap_below_top = compact_spacing
        required_gap_above_bottom = compact_spacing
        
        # Check top axis: if it has ticks on bottom, need space below it
        if top_ax in axes_with_labels_set and self._has_labels(top_ax):
            if self._get_tick_position(top_ax) == 'bottom':
                required_gap_below_top = self._calculate_required_space_for_labels(top_ax, direction='below')
        
        # Check bottom axis: if it has ticks on top, need space above it
        if bottom_ax in axes_with_labels_set and self._has_labels(bottom_ax):
            if self._get_tick_position(bottom_ax) == 'top':
                required_gap_above_bottom = self._calculate_required_space_for_labels(bottom_ax, direction='above')
        
        return max(required_gap_below_top, required_gap_above_bottom, compact_spacing)
    
    def _adjust_figure_spacing(self) -> None:
        """
        Adjust vertical spacing by moving axes positions (moving approach).
        
        This method:
        1. Separates overlapping axes
        2. Calculates target positions based on label requirements
        3. Applies movements to all axes simultaneously
        
        Spacing is applied to ALL axes in self.all_axes.
        Label calculations use self.axes_to_align (to determine spacing requirements).
        """
        if not self.adjust_spacing:
            return
        
        if self.fig is None:
            if len(self.all_axes) > 0:
                self.fig = self.all_axes[0].figure
            else:
                return
        
        self.log.write("Adjusting vertical spacing using moving approach...", verbose=self.verbose)
        
        # Use all_axes for spacing (preserve original order from top to bottom)
        axes_for_spacing = list(self.all_axes)
        
        if len(axes_for_spacing) <= 1:
            return
        
        # Create a set of axes that have labels (for spacing calculations)
        axes_with_labels_set = set(self.axes_to_align) if self.axes_to_align else set()
        compact_spacing = 0.01  # Minimal spacing when no labels
        
        # Step 1: Separate overlapping adjacent axes
        self._separate_overlapping_axes(axes_for_spacing, compact_spacing)
        
        # Step 2: Store original positions (after separation)
        original_positions = []
        for ax in axes_for_spacing:
            pos = ax.get_position()
            original_positions.append({
                'y0': max(0.0, min(1.0, pos.y0)),
                'y1': max(0.0, min(1.0, pos.y1)),
                'height': max(0.0, min(1.0, pos.height)),
                'x0': pos.x0,
                'x1': pos.x1,
                'width': pos.width
            })
        
        num_gaps = len(axes_for_spacing) - 1
        if num_gaps <= 0:
            return
        
        # Step 3: Calculate target positions from top to bottom
        # In matplotlib: y=1.0 is TOP, y=0.0 is BOTTOM
        # We stack from TOP (y=1.0) downward (toward y=0.0)
        fig_height = 1.0
        target_y0 = [fig_height - original_positions[0]['height']]  # Top axis at top of figure
        
        for gap_num in range(num_gaps):
            top_ax_idx = gap_num
            bottom_ax_idx = gap_num + 1
            
            top_y0_target = target_y0[top_ax_idx]
            top_y1_target = top_y0_target + original_positions[top_ax_idx]['height']
            
            # Calculate required gap based on tick positions
            top_ax = axes_for_spacing[top_ax_idx]
            bottom_ax = axes_for_spacing[bottom_ax_idx]
            required_gap = self._calculate_required_gap(top_ax, bottom_ax, axes_with_labels_set, compact_spacing)
            
            # Calculate target position for bottom axis (below top axis)
            target_bottom_y0 = top_y0_target - required_gap - original_positions[bottom_ax_idx]['height']
            
            # Clamp to valid range
            bottom_ax_height = original_positions[bottom_ax_idx]['height']
            max_allowed_y0 = top_y0_target - compact_spacing - bottom_ax_height
            
            if target_bottom_y0 < 0.0:
                target_bottom_y0 = max(0.0, max_allowed_y0)
            elif target_bottom_y0 > max_allowed_y0:
                target_bottom_y0 = max_allowed_y0
            
            target_y0.append(target_bottom_y0)
        
        # Verify target_y0 has correct length
        if len(target_y0) != len(axes_for_spacing):
            while len(target_y0) < len(axes_for_spacing):
                target_y0.append(original_positions[len(target_y0)]['y0'])
            target_y0 = target_y0[:len(axes_for_spacing)]
        
        # Step 4: Calculate and apply movements simultaneously
        fig_height = 1.0
        for i, ax in enumerate(axes_for_spacing):
            movement = target_y0[i] - original_positions[i]['y0']
            original_pos = original_positions[i]
            pos = ax.get_position()
            
            new_y0 = original_pos['y0'] + movement
            new_y1 = new_y0 + original_pos['height']
            
            # Clamp to valid range (0-1)
            if new_y1 > fig_height:
                new_y0 = max(0.0, fig_height - original_pos['height'])
            elif new_y0 < 0.0:
                new_y0 = 0.0
            
            # Update position if changed
            if abs(new_y0 - pos.y0) > 0.001:
                ax.set_position([pos.x0, new_y0, pos.width, pos.height])
            
        
        self.log.write("Finished adjusting vertical spacing.", verbose=self.verbose)
    
    def align(self) -> None:
        """
        Apply x-tick and label alignment to registered axes, and spacing to all axes.
        
        This method:
        1. Calculates xlim, xticks, and xticklabels (if not explicitly provided)
        2. Applies x-tick and label alignment to axes_to_align only
        3. Applies vertical spacing adjustment to all axes in all_axes
        """
        # Step 1: Calculate xlim, ticks, and labels
        try:
            xlim = self._calculate_xlim()
            xticks = self._calculate_xticks()
            xticklabels = self._calculate_xticklabels()
            xlabel = self._calculate_xlabel()
        except ValueError as e:
            self.log.warning(f"Could not calculate axis parameters: {e}", verbose=self.verbose)
            return
        
        # Step 2: Apply x-tick and label alignment to registered axes only
        if len(self.axes_to_align) > 0:
            self.log.write(
                f"Aligning {len(self.axes_to_align)} axes: xlim={xlim}, num_ticks={len(xticks) if len(xticks) > 0 else 'auto'}",
                verbose=self.verbose
            )
            
            # Apply alignment to each registered axis
            # Since all axes will have the same ticks/labels after alignment,
            # show tick labels only on the bottom-most axis to avoid redundancy
            # Upper axes show ticks but no labels
            # IMPORTANT: Determine bottom axis based on position in all_axes, not order in axes_to_align
            # Find which axis in axes_to_align is at the bottom (lowest y position in all_axes)
            bottom_axis_in_align = None
            bottom_y0 = 1.0  # Start with top (y=1.0 is top in matplotlib)
            for ax in self.axes_to_align:
                if ax in self.all_axes:
                    pos = ax.get_position()
                    if pos.y0 < bottom_y0:  # Lower y0 means it's closer to bottom
                        bottom_y0 = pos.y0
                        bottom_axis_in_align = ax
            
            for i, ax in enumerate(self.axes_to_align):
                is_bottom_axis = (ax == bottom_axis_in_align)
                
                # Set xlim
                ax.set_xlim(xlim)
                
                # Set ticks
                if len(xticks) > 0:
                    ax.set_xticks(xticks)
                    
                    # Set tick labels only on bottom axis
                    if is_bottom_axis:
                        # Bottom axis: show tick labels
                        if xticklabels is not None:
                            ax.set_xticklabels(xticklabels, rotation=45)
                        elif self.fontsize is not None or self.font_family is not None:
                            # Update tick label properties if font settings provided
                            tick_kwargs = {}
                            if self.fontsize is not None:
                                tick_kwargs['fontsize'] = self.fontsize
                            if self.font_family is not None:
                                tick_kwargs['family'] = self.font_family
                            if tick_kwargs:
                                ax.tick_params(axis='x', **tick_kwargs)
                    else:
                        # Upper axes: show ticks but hide labels (since next axis has same ticks/labels)
                        ax.set_xticklabels([])
                
                # Set xlabel only on bottom axis
                if is_bottom_axis and xlabel is not None:
                    label_kwargs = {}
                    if self.fontsize is not None:
                        label_kwargs['fontsize'] = self.fontsize
                    if self.font_family is not None:
                        label_kwargs['family'] = self.font_family
                    ax.set_xlabel(xlabel, **label_kwargs)
                elif not is_bottom_axis:
                    # Upper axes: remove xlabel (bottom axis will show it)
                    ax.set_xlabel("")
        
        # Step 3: Adjust vertical spacing for ALL axes in all_axes
        # Uses self.axes_to_align internally to determine which axes have labels
        self._adjust_figure_spacing()
        
        if len(self.axes_to_align) > 0:
            self.log.write(f"Successfully aligned {len(self.axes_to_align)} axes.", verbose=self.verbose)
