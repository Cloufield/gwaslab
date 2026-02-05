import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log


def _get_text_rotated_height(text, renderer, fig=None):
    """
    Calculate the maximum y-coordinate of a text object accounting for rotation.
    
    For rotated text (especially 90 degrees), the standard get_window_extent()
    may not always report accurate bounds. This function explicitly calculates
    the vertical extent based on the text's position, actual text length, and rotation.
    
    Parameters
    ----------
    text : matplotlib.text.Text
        The text object to measure.
    renderer : matplotlib.backend_bases.RendererBase
        The renderer used to calculate text extents.
    fig : matplotlib.figure.Figure, optional
        The figure object, used for DPI calculation.
    
    Returns
    -------
    float
        The highest y-coordinate in pixels that the text reaches.
    """
    try:
        # Get the text string
        text_string = text.get_text()
        if not text_string or not text_string.strip():
            return 0
        
        # Get rotation angle in degrees
        # Handle string rotation values ('vertical', 'horizontal')
        rotation = text.get_rotation()
        if isinstance(rotation, str):
            rotation = 90.0 if rotation == "vertical" else 0.0
        rotation = float(rotation)
        
        # Get the text position in display coordinates
        position = text.get_transform().transform(text.get_position())
        x_pos, y_pos = position
        
        # Get the UNROTATED text dimensions using font metrics
        # This gives us the actual width based on the text string length
        fontsize = text.get_fontsize()
        
        # Determine ismath parameter for get_text_width_height_descent
        # ismath should be "TeX" if usetex is True, otherwise check for mathtext
        if text.get_usetex():
            ismath = "TeX"
        elif "$" in text_string:
            # Simple heuristic for mathtext detection
            ismath = True
        else:
            ismath = False
        
        # Use renderer to get actual text dimensions
        # get_text_width_height_descent returns (width, height, descent) in display coords
        descent = 0
        try:
            # Use public API get_fontproperties() instead of private _fontproperties
            width, height, descent = renderer.get_text_width_height_descent(
                text_string, 
                text.get_fontproperties(), 
                ismath=ismath
            )
        except Exception:
            # Fallback: estimate based on fontsize and string length
            # Rough estimate: each character is about 0.6 * fontsize wide
            dpi = fig.dpi if fig is not None else 100
            width = len(text_string) * fontsize * 0.6 * dpi / 72
            height = fontsize * dpi / 72
        
        # If no significant rotation, calculate simple top position
        if abs(rotation) < 1:
            bbox = text.get_window_extent(renderer)
            if bbox is not None:
                return bbox.y1
            return y_pos + height
        
        # Convert rotation to radians
        rotation_rad = np.radians(rotation)
        
        # Get alignment to determine anchor point
        ha = text.get_ha()
        va = text.get_va()
        
        # Check rotation_mode - affects how rotation is applied relative to anchor
        # 'default': rotate around anchor point (our assumption)
        # 'anchor': text is positioned then rotated around its anchor
        rotation_mode = text.get_rotation_mode()
        
        # Determine anchor offsets (where is the anchor relative to bbox corner)
        if ha == 'left':
            x_offset = 0
        elif ha == 'right':
            x_offset = -width
        else:  # center
            x_offset = -width / 2
            
        if va == 'bottom':
            y_offset = 0
        elif va == 'top':
            y_offset = -height
        elif va == 'baseline':
            # For baseline, use descent to approximate baseline position
            # This is an approximation as exact behavior depends on font metrics
            y_offset = -descent if descent > 0 else -height * 0.2
        elif va == 'center_baseline':
            y_offset = -height / 2
        else:  # center
            y_offset = -height / 2
        
        # Four corners of unrotated text bbox (relative to anchor)
        corners = [
            (x_offset, y_offset),
            (x_offset + width, y_offset),
            (x_offset + width, y_offset + height),
            (x_offset, y_offset + height),
        ]
        
        # Rotate corners and find max y
        cos_r = np.cos(rotation_rad)
        sin_r = np.sin(rotation_rad)
        
        max_y = y_pos
        for cx, cy in corners:
            # Apply rotation around anchor point
            # Note: rotation_mode='anchor' vs 'default' can affect this,
            # but for title positioning purposes, we use a consistent approach
            rotated_y = cx * sin_r + cy * cos_r
            # Translate to absolute position
            abs_y = y_pos + rotated_y
            max_y = max(max_y, abs_y)
        
        return max_y
    
    except Exception:
        return 0


def _get_highest_y_pixels(target_ax, renderer, fig=None):
    """
    Get the highest y pixel of all content including rotated text.
    
    This function finds the maximum y-coordinate (in pixels) of all visible
    content in the given axes, including rotated text annotations which may
    extend beyond the axes tight bounding box.
    
    Parameters
    ----------
    target_ax : matplotlib.axes.Axes
        The axes to check for content bounds.
    renderer : matplotlib.backend_bases.RendererBase
        The renderer used to calculate text extents.
    fig : matplotlib.figure.Figure, optional
        The figure object, used for DPI calculation. If not provided,
        will attempt to get from target_ax.
    
    Returns
    -------
    float
        The highest y-coordinate in pixels of any content in the axes.
    """
    highest_y = 0
    
    # Get figure from axes if not provided
    if fig is None:
        try:
            fig = target_ax.get_figure()
        except Exception:
            pass
    
    # Check axes tight bbox
    try:
        bbox = target_ax.get_tightbbox(renderer)
        if bbox is not None:
            highest_y = max(highest_y, bbox.y1)
    except Exception:
        pass
    
    # Explicitly check all text artists with rotation handling
    for text in target_ax.texts:
        try:
            # Use rotation-aware height calculation
            text_max_y = _get_text_rotated_height(text, renderer, fig)
            highest_y = max(highest_y, text_max_y)
            
            # Also check standard window extent as fallback
            text_bbox = text.get_window_extent(renderer)
            if text_bbox is not None:
                highest_y = max(highest_y, text_bbox.y1)
        except Exception:
            pass
    
    # Check annotations (they have text + arrow)
    # Only check actual text/annotation objects, NOT all children (spines, patches, etc. can report inflated bounds)
    from matplotlib.text import Text, Annotation
    annotations = [c for c in target_ax.get_children() 
                   if isinstance(c, (Text, Annotation))]
    for child in annotations:
        try:
            # Skip if already processed in ax.texts
            if child in target_ax.texts:
                continue
            
            # Handle text/annotation objects
            if hasattr(child, 'get_text') and hasattr(child, 'get_rotation'):
                child_max_y = _get_text_rotated_height(child, renderer, fig)
                highest_y = max(highest_y, child_max_y)
            
            # Also check window extent
            child_bbox = child.get_window_extent(renderer)
            if child_bbox is not None:
                highest_y = max(highest_y, child_bbox.y1)
        except Exception:
            pass
    
    return highest_y


def adjust_text_position(positions, yspan, repel_force=0.01, max_iter=100,amode="int",log=Log(),verbose=True, min_factor=None):
    # check the number of variants to annotate 
    #if repel_force>0:
    #    if 1/(repel_force*2 +0.01) < len(positions):
    #        log.write(" -Too many variants to annotate; maybe it is better to reduce the number of variants")
    #else:
    if len(positions)>30:
        log.write(" -Too many variants to annotate; maybe it is better to reduce the number of variants",verbose=verbose)

    # calculate the steps
    if amode=="int":
        step = int(yspan*repel_force) 
    elif amode=="log":
        if min_factor is None:
            min_factor = np.min(positions)
        #(1, max) -> (0, log(max)))
        positions = np.log2(positions/min_factor)
        step = max(positions)*repel_force 

    else:
        step = yspan*repel_force 

    # start iteration
    for i in range(max_iter):
        # check overlap
        index = check_overlap(positions, step)

        # fix overlap if needed
        if index == len(positions)+1:
            # if no overlap
            if amode=="int":
                return  np.floor(pd.to_numeric(positions, errors='coerce')).astype('Int64').copy()
            elif amode=="log":
                
                return  np.power(2, pd.to_numeric(positions, errors='coerce'))* min_factor
            else:
                return  pd.to_numeric(positions, errors='coerce')
        else:
            # if overlap exists
            if amode=="int":
                move_position_from_center_int(positions, index, step)
            else:
                move_position_from_center_float(positions, index, step)
    
    # when reaching maximum iteration, return anyway
    log.write(" -Reaching maximum iteration: {}; Skipping...".format(max_iter),verbose=verbose)
    if amode=="int":
        return np.floor(pd.to_numeric(positions, errors='coerce')).astype('Int64').copy()
    elif amode=="log":
        return np.exp(pd.to_numeric(positions, errors='coerce')) * min_factor
    else:
        return  pd.to_numeric(positions, errors='coerce')

def check_overlap(positions,step):
    # check overlap by walkthrough
    for i in range(1,len(positions)):
        if positions[i] - positions[i-1]< step:
            # if overlap (distance between two positions < step), return the index
            return i
    #if no overlap , return maxmum length + 1
    return len(positions)+1

def move_position_from_center_int(positions, center_i, step):
    # left side
    #print("center_i",center_i)
    #print("before",positions)
    
    # if not pass left bound
    if positions[center_i-1] - step//2 > 0:
        positions[center_i-1] = positions[center_i-1] - step//2
    
    # if not the second element
    if positions[center_i-2] - step//2 > 0:
        for i in range(center_i-2, 0, -1):
            if abs(positions[i] - positions[i+1]) < step:
                 positions[i] = positions[i] - step//2
            else:
                break
    
    # right side            
    positions[center_i] = positions[center_i] + step//2
    for i in range(center_i+1, len(positions)):
        if abs(positions[i] - positions[i-1]) < step:
             positions[i] = positions[i] + step//2
        else:
            break
    #print("after",positions)
    return positions
    
def move_position_from_center_float(positions, center_i, step):
    # left side
    #print("center_i",center_i)
    #print("before",positions)
    
    # if not pass left bound
    if positions[center_i-1] - step/2 > 0:
        positions[center_i-1] = positions[center_i-1] - step/2
    
    # if not the second element
    if positions[center_i-2] - step/2 > 0:
        for i in range(center_i-2, 0, -1):
            if abs(positions[i] - positions[i+1]) < step:
                 positions[i] = positions[i] - step/2
            else:
                break
    
    # right side            
    positions[center_i] = positions[center_i] + step/2
    for i in range(center_i+1, len(positions)):
        if abs(positions[i] - positions[i-1]) < step:
             positions[i] = positions[i] + step/2
        else:
            break
    #print("after",positions)
    return positions

