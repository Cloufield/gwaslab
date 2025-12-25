import matplotlib.pyplot as plt
from gwaslab.info.g_Log import Log
import time
import os
import functools

def save_figure(fig, save, keyword, save_kwargs=None, log = Log(), verbose=True):
    """
    Save matplotlib figure to file.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save
    save : str, bool, or None
        If str: file path to save to
        If True: save to default path
        If None/False: skip saving
    keyword : str
        Keyword for default filename generation
    save_kwargs : dict, optional
        Additional arguments passed to fig.savefig()
    log : Log, optional
        Logging object
    verbose : bool, optional
        Whether to print verbose messages
    """
    # Set verbose to False when save is None
    if save is None:
        verbose = False
    
    log.write("Start to save figure..." ,verbose=verbose)
    if save_kwargs is None:
        save_kwargs = {}
    
    if not save:
        log.write(" -Skip saving figure!" ,verbose=verbose)
        log.write("Finished saving figure..." ,verbose=verbose)
        return
    
    # Determine file path
    if save is True:
        file_path = get_default_path(keyword)
    else:
        file_path = str(save)
    
    # Get file extension for format-specific handling
    _, ext = os.path.splitext(file_path)
    ext = ext.lower().lstrip('.')
    
    # Vector formats (PDF, SVG, EPS) don't use bbox_inches="tight" by default
    # unless explicitly specified in save_kwargs
    vector_formats = {'pdf', 'svg', 'eps'}
    is_vector_format = ext in vector_formats
    
    # Prepare save arguments
    save_params = save_kwargs.copy()
    
    # For non-vector formats, use bbox_inches="tight" unless overridden
    if not is_vector_format and 'bbox_inches' not in save_params:
        save_params['bbox_inches'] = 'tight'
    
    # Check if file exists
    file_exists = os.path.exists(file_path)
    
    try:
        # Save the figure
        fig.savefig(file_path, **save_params)
        
        # Log success message
        if file_exists:
            overwrite_msg = f" (overwrite)" if ext else ""
            format_msg = f" ({ext})" if ext else ""
            log.write(f" -Saved to {file_path} successfully!{format_msg}{overwrite_msg}" ,verbose=verbose)
        else:
            format_msg = f" ({ext})" if ext else ""
            log.write(f" -Saved to {file_path} successfully!{format_msg}" ,verbose=verbose)
            
    except Exception as e:
        log.warning(f"Failed to save figure to {file_path}: {str(e)}")
        raise
    
    log.write("Finished saving figure..." ,verbose=verbose)

def get_default_path(keyword, fmt="png"):
    """
    Generate default file path for saving figures.
    
    Parameters
    ----------
    keyword : str
        Keyword to look up in path_dictionary
    fmt : str, optional
        File format extension (default: "png")
    
    Returns
    -------
    str
        Generated file path with timestamp and counter
    """
    path_dictionary = { 
        "m": "manhattan",
        "qq": "qq",
        "mqq": "mqq",
        "qqm": "qqm",
        "b": "brisbane",
        "r": "regional",
        "stacked_r": "stacked_regional",
        "trumpet_b": "trumpet_binary",
        "trumpet_q": "trumpet_quant",
        "power_b": "power_binary",
        "power_q": "power_quant",
        "power_xb": "power_x_binary",
        "power_xq": "power_x_quant",
        "ldscrg": "ldscrg_heatmap",
        "miami": "miami",
        "esc": "effect_size_comparision",
        "afc": "allele_frequency_comparision",
        "gwheatmap": "genome_wide_heatmap",
        "scatter": "scatter",
        "forest": "forest",
        "associations": "associations",
        "plot_associations": "associations",
        "gwaslab": "gwaslab"
    }
    
    prefix = path_dictionary.get(keyword, "gwaslab")
    timestamp = time.strftime('%Y%m%d')
    count = 1
    
    # Find first available filename (up to 10000 attempts)
    for _ in range(10000):
        file_path = f"./gwaslab_{prefix}_{timestamp}_{count}.{fmt}"
        if not os.path.exists(file_path):
            return file_path
        count += 1
    
    # Fallback if all attempts fail
    return f"./gwaslab_{prefix}_{timestamp}_{count}.{fmt}"

def safefig(func):
    """
    Decorator to safely handle matplotlib figures in case of exceptions.
    
    Closes all matplotlib figures if an exception occurs during function execution,
    preventing empty figure windows from being displayed.
    
    Parameters
    ----------
    func : callable
        Function to wrap
    
    Returns
    -------
    callable
        Wrapped function with exception handling
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            # Close all figures to prevent empty figure output on error
            plt.close('all')
            raise e
    return wrapper