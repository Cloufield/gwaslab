import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from typing import Optional, Dict, Any, Union, List
from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_save_figure import save_figure

#STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
#1	TssA	Active TSS	Red	255,0,0
#2	TssAFlnk	Flanking Active TSS	Orange Red	255,69,0
#3	TxFlnk	Transcr. at gene 5' and 3'	LimeGreen	50,205,50
#4	Tx	Strong transcription	Green	0,128,0
#5	TxWk	Weak transcription	DarkGreen	0,100,0
#6	EnhG	Genic enhancers	GreenYellow	194,225,5
#7	Enh	Enhancers	Yellow	255,255,0
#8	ZNF/Rpts	ZNF genes & repeats	Medium Aquamarine	102,205,170
#9	Het	Heterochromatin	PaleTurquoise	138,145,208
#10	TssBiv	Bivalent/Poised TSS	IndianRed	205,92,92
#11	BivFlnk	Flanking Bivalent TSS/Enh	DarkSalmon	233,150,122
#12	EnhBiv	Bivalent Enhancer	DarkKhaki	189,183,107
#13	ReprPC	Repressed PolyComb	Silver	128,128,128
#14	ReprPCWk	Weak Repressed PolyComb	Gainsboro	192,192,192
#15	Quies	Quiescent/Low	White	255,255,255

color_dict={
    "E1": np.array([255,0,0]),
    "E2": np.array([255,69,0]),
    "E3": np.array([50,205,50]),
    "E4": np.array([0,128,0]),
    "E5": np.array([0,100,0]),
    "E6": np.array([194,225,5]),
    "E7": np.array([255,255,0]),
    "E8": np.array([102,205,170]),
    "E9": np.array([138,145,208]),
    "E10":np.array([205,92,92]),
    "E11":np.array([233,150,122]),
    "E12":np.array([189,183,107]),
    "E13":np.array([128,128,128]),
    "E14":np.array([192,192,192]),
    "E15":np.array([255,255,255])
}

color_dict_i={
    1: np.array([255,0,0]),
    2: np.array([255,69,0]),
    3: np.array([50,205,50]),
    4: np.array([0,128,0]),
    5: np.array([0,100,0]),
    6: np.array([194,225,5]),
    7: np.array([255,255,0]),
    8: np.array([102,205,170]),
    9: np.array([138,145,208]),
    10:np.array([205,92,92]),
    11:np.array([233,150,122]),
    12:np.array([189,183,107]),
    13:np.array([128,128,128]),
    14:np.array([192,192,192]),
    15:np.array([255,255,255])
}


def _load_chromatin_file(file_path: str, target_chr: int, target_start: int, target_end: int, 
                         log: Log, verbose: bool) -> pd.DataFrame:
    """
    Load and filter chromatin state data from a BED file.
    
    Parameters
    ----------
    file_path : str
        Path to chromatin state BED file
    target_chr : int
        Target chromosome number
    target_start : int
        Region start position
    target_end : int
        Region end position
    log : Log
        Logger instance
    verbose : bool
        Whether to show progress messages
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: STATE_i, START, END
        Contains only records overlapping the target region
    """
    log.write(f" -Loading: {file_path}", verbose=verbose)
    
    # Read BED file
    df = pd.read_csv(file_path, sep="\t", header=None)
    df.columns = ["ID", "START", "END", "STATE"]
    
    # Extract chromosome number from ID column (format: "chr1_...")
    df["CHR"] = df["ID"].str.extract(r"chr([0-9]+)").astype("float").astype("Int64")
    
    # Extract state number from STATE column (format: "1_TssA" or "1_")
    df["STATE_i"] = df["STATE"].str.extract(r"([0-9]+)_*").astype("float").astype("Int64")
    
    # Filter for records in the target region
    # Record overlaps if: same chromosome AND (end > region_start AND start < region_end)
    in_region_mask = (
        (df["CHR"] == target_chr) & 
        (df["END"] > target_start) & 
        (df["START"] < target_end)
    )
    df_filtered = df.loc[in_region_mask, ["STATE_i", "START", "END"]].copy()
    
    # Sort by state number (descending) for consistent plotting order
    df_filtered = df_filtered.sort_values("STATE_i", ascending=False)
    
    log.write(f"  -Number of records in specified region: {len(df_filtered)}", verbose=verbose)
    
    return df_filtered


def _calculate_line_width(fig, ax, row_height: float = 0.1) -> float:
    """
    Calculate appropriate line width for chromatin state bars.
    
    The line width is calculated to match the visual height of each row,
    ensuring bars appear as solid rectangles rather than thin lines.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object for DPI calculation
    ax : matplotlib.axes.Axes
        Axes object for coordinate transformation
    row_height : float, default=0.1
        Height of each row in data coordinates
    
    Returns
    -------
    float
        Line width in points
    """
    # Get pixel height for row_height in data coordinates
    pixel_height = abs(
        ax.transData.transform([0, 0])[1] - 
        ax.transData.transform([0, row_height])[1]
    )
    
    # Convert pixels to points (1 point = 1/72 inch)
    # points = pixels * (72 / dpi)
    line_width_points = pixel_height * 72 / fig.dpi
    
    return line_width_points


def _plot_chromatin_state(
    region_chromatin_files, 
    region, 
    region_chromatin_labels: Optional[List[str]] = None,
    fig: Optional[plt.Figure] = None, 
    ax: Optional[plt.Axes] = None,
    xlim_i: Optional[list] = None,
    fontsize: int = 12,
    font_family: str = "Arial",
    fig_kwargs: Optional[Dict[str, Any]] = None,
    save: Optional[Union[bool, str]] = None,
    save_kwargs: Optional[Dict[str, Any]] = None,
    log: Log = Log(),
    verbose: bool = True
):
    """
    Plot chromatin state tracks from Roadmap 15-state model files.
    
    This function visualizes chromatin states from Roadmap Epigenomics Project
    15-state core marks files. Each state is displayed as a colored horizontal
    bar, with different tissues/cell types shown as separate rows.
    
    Parameters
    ----------
    region_chromatin_files : list of str
        List of paths to Roadmap 15_coreMarks_mnemonics.bed.gz files.
        Each file should contain chromatin state annotations in BED format.
    region : tuple of (int, int, int)
        Genomic region as (chromosome, start, end) in 1-based coordinates.
        Example: (1, 1000000, 2000000) for chr1:1000000-2000000
    region_chromatin_labels : list of str, optional
        List of labels for each chromatin file (one label per file).
        Should match the length of region_chromatin_files.
        If None, labels are automatically extracted from filenames by splitting
        on "_" and using the first part (e.g., "E098_15_coreMarks_mnemonics.bed.gz" -> "E098").
    fig : matplotlib.figure.Figure, optional
        Figure object to plot on. If None, a new figure will be created.
        If ax is provided but fig is None, fig will be extracted from ax.
    ax : matplotlib.axes.Axes, optional
        Axes object to plot on. If None, a new figure and axes will be created.
    xlim_i : list of float, optional
        X-axis limit offset. If None, defaults to [0].
        Used for aligning with other plots in stacked panels.
    fontsize : int, default=12
        Font size for labels.
    font_family : str, default="Arial"
        Font family for text labels.
    fig_kwargs : dict, optional
        Additional keyword arguments for figure creation (e.g., {'figsize': (10, 4), 'dpi': 200}).
        Only used when creating a new figure. Default figsize is (10, 2 + 0.3 * n_tracks).
    save : str, bool, or None, optional
        If str: file path to save figure.
        If True: save to default path.
        If None/False: skip saving.
    save_kwargs : dict, optional
        Additional arguments for saving (e.g., {'dpi': 300, 'bbox_inches': 'tight'}).
    log : Log, default=Log()
        Logger instance for messages.
    verbose : bool, default=True
        Whether to show progress messages.
    
    Returns
    -------
    matplotlib.figure.Figure
        The figure object.
    """
    # Calculate number of tracks for default figsize
    n_tracks = len(region_chromatin_files)
    
    # Auto-generate labels from filenames if not provided
    if region_chromatin_labels is None:
        region_chromatin_labels = []
        for file_path in region_chromatin_files:
            # Extract filename from path
            filename = os.path.basename(file_path)
            # Split by "_" and use the first part as label
            # Example: "E098_15_coreMarks_mnemonics.bed.gz" -> "E098"
            label = filename.split("_")[0]
            region_chromatin_labels.append(label)
        log.write(f" -Auto-generated labels from filenames: {region_chromatin_labels}", verbose=verbose)
    
    # Handle figure and axes creation
    if fig is None:
        if ax is None:
            # Create new figure and axes
            if fig_kwargs is None:
                fig_kwargs = {}
            # Set default figsize if not provided
            if "figsize" not in fig_kwargs:
                # Default: width=10, height based on number of tracks
                fig_kwargs["figsize"] = (10, 2 + 0.3 * n_tracks)
            fig, ax = plt.subplots(**fig_kwargs)
        else:
            # Extract figure from axes
            fig = ax.figure
    elif ax is None:
        # Get axes from figure, or create one if needed
        ax = fig.gca()
        if ax is None:
            ax = fig.add_subplot(111)
    
    # Parse region coordinates
    target_chr, target_start, target_end = region[0], region[1], region[2]
    
    # Calculate x-axis offset for alignment with other panels
    if xlim_i is None:
        xlim_i = [0]
    x_offset = xlim_i[0] - target_start
    
    # Set up axes limits
    row_height = 0.1
    y_min = -0.05
    y_max = row_height * n_tracks - 0.05
    ax.set_ylim([y_min, y_max])
    ax.set_xlim([x_offset + target_start, x_offset + target_end])
    
    # Calculate line width for chromatin state bars
    line_width = _calculate_line_width(fig, ax, row_height)
    
    # Plot chromatin states for each tissue/cell type
    for track_index, file_path in enumerate(region_chromatin_files):
        # Load and filter chromatin data for this file
        df_chromatin = _load_chromatin_file(
            file_path, target_chr, target_start, target_end, log, verbose
        )
        
        # Calculate y-position for this track (bottom of row)
        y_position = track_index * row_height
        
        # Plot each chromatin state segment as a horizontal bar
        for _, row in df_chromatin.iterrows():
            state_number = int(row["STATE_i"])
            segment_start = x_offset + row["START"]
            segment_end = x_offset + row["END"]
            
            # Get color for this state (RGB values 0-255, convert to 0-1 range)
            state_color_rgb = color_dict_i[state_number]
            state_color_normalized = state_color_rgb / 255.0
            
            # Plot horizontal line representing the chromatin state segment
            ax.plot(
                [segment_start, segment_end],
                [y_position, y_position],
                c=state_color_normalized,
                linewidth=line_width,
                solid_capstyle="butt",
                rasterized=True
            )
    
    # Set y-axis labels for each track
    if len(region_chromatin_labels) == n_tracks:
        y_tick_positions = [i * row_height for i in range(n_tracks)]
        ax.set_yticks(y_tick_positions, region_chromatin_labels, 
                     fontsize=fontsize, family=font_family)
    else:
        ax.set_yticks(ticks=[])
    
    # Invert y-axis so first track appears at top
    ax.invert_yaxis()
    
    # Save figure if requested (save_figure handles None/False internally)
    save_figure(fig=fig, save=save, keyword="chromatin", save_kwargs=save_kwargs, 
               log=log, verbose=verbose)
    
    return fig
