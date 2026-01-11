import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.collections import PatchCollection
from pathlib import Path
from typing import Optional, Union, Dict, Any, List
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.bd.bd_common_data import get_chr_to_number, get_number_to_chr
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style
from gwaslab.viz.viz_aux_reposition_text import adjust_text_position
from adjustText import adjust_text


def _plot_phenogram(
    insumstats,
    snpid: str = "SNPID",
    chrom: str = "CHR",
    pos: str = "POS",
    p: str = "P",
    mlog10p: str = "MLOG10P",
    windowsizekb: int = 500,
    sig_level: float = 5e-8,
    cytoband_path: Optional[str] = None,
    build: str = "19",
    ncols: int = 11,
    figsize: tuple = (20, 40),
    dpi: int = 100,
    annotate_snps: bool = True,
    annotation_kwargs: Optional[Dict[str, Any]] = None,
    anno_style: str = "expand",
    repel_force: float = 0.5,
    anno_max_iter: int = 100,
    save: Union[bool, str] = False,
    save_kwargs: Optional[Dict[str, Any]] = None,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    verbose: bool = True,
    log: Log = Log(),
) -> plt.Figure:
    """
    Create a phenogram plot showing cytobands for each chromosome with lead SNPs annotated.
    
    This function generates a karyotype-like visualization (phenogram) displaying:
    - Cytogenetic bands for each chromosome
    - Lead SNPs from the sumstats data
    - Optional annotations for SNP IDs
    
    Parameters
    ----------
    mysumstats : Sumstats
        Sumstats object containing the GWAS data
    snpid : str, default="SNPID"
        Column name for SNP identifier
    chrom : str, default="CHR"
        Column name for chromosome
    pos : str, default="POS"
        Column name for position
    p : str, default="P"
        Column name for p-value
    mlog10p : str, default="MLOG10P"
        Column name for -log10(p-value)
    windowsizekb : int, default=500
        Window size in kilobases for lead variant identification
    sig_level : float, default=5e-8
        Significance threshold for variant selection
    cytoband_path : str, optional
        Path to cytoband file. If None, uses default hg19 cytoband file
    build : str, default="19"
        Genome build version ("19" for hg19, "38" for hg38)
    ncols : int, default=11
        Number of columns for arranging chromosomes
    figsize : tuple, default=(20, 40)
        Figure size in inches
    dpi : int, default=100
        Resolution of the figure
    annotate_snps : bool, default=True
        If True, annotate SNP IDs on the plot
    annotation_kwargs : dict, optional
        Additional keyword arguments for text annotations
    save : bool or str, default=False
        If True or str, save the figure
    save_kwargs : dict, optional
        Additional keyword arguments for saving
    fig_kwargs : dict, optional
        Additional keyword arguments for figure creation
    verbose : bool, default=True
        If True, print progress messages
    log : Log, default=Log()
        Logging object
        
    Returns
    -------
    matplotlib.figure.Figure
        The created matplotlib figure object
    """
    
    log.write("Start to create phenogram plot...", verbose=verbose)
    
    # Extract dataframe if Sumstats object is passed
    if hasattr(insumstats, 'data') and not isinstance(insumstats, pd.DataFrame):
        insumstats = insumstats.data
    
    # Create working copy to preserve original (plotting functions should not modify input)
    sumstats = insumstats.copy()
    
    # Get leads from sumstats
    log.write(" -Extracting lead variants...", verbose=verbose)
    leads = _get_sig(
        insumstats_or_dataframe=sumstats,
        variant_id=snpid,
        chrom=chrom,
        pos=pos,
        p=p,
        mlog10p=mlog10p,
        windowsizekb=windowsizekb,
        sig_level=sig_level,
        log=log,
        verbose=verbose
    )
    
    if leads is None or len(leads) == 0:
        log.write(" -No lead variants found. Plotting chromosomes without annotations.", verbose=verbose)
        leads = pd.DataFrame()
    else:
        log.write(" -Found {} lead variants to annotate.".format(len(leads)), verbose=verbose)
    
    # Load cytoband data
    if cytoband_path is None:
        # Use default cytoband file path
        data_dir = Path(__file__).parent.parent / "data" / "cytoband"
        if build == "19":
            cytoband_path = data_dir / "cytoBand_hg19.txt.gz"
        elif build == "38":
            cytoband_path = data_dir / "cytoBand_hg38.txt.gz"
        else:
            # Default to hg19
            cytoband_path = data_dir / "cytoBand_hg19.txt.gz"
            log.write(" -Unknown build '{}', using hg19 cytoband data.".format(build), verbose=verbose)
    
    log.write(" -Loading cytoband data from: {}".format(cytoband_path), verbose=verbose)
    
    try:
        cytobands = pd.read_csv(cytoband_path, sep="\s+", header=None, compression='gzip')
        cytobands.columns = ["CHR", "START", "END", "ARM", "STAIN"]
    except Exception as e:
        log.write(" -Error loading cytoband data: {}".format(e), verbose=verbose)
        raise ValueError("Could not load cytoband data from {}".format(cytoband_path))
    
    # Color dictionary for cytobands
    color_dict = {
        "gpos100": (100/255, 100/255, 100/255),
        "gpos": (100/255, 100/255, 100/255),
        "gpos75": (110/255, 110/255, 110/255),
        "gpos66": (130/255, 130/255, 130/255),
        "gpos50": (160/255, 160/255, 160/255),
        "gpos33": (180/255, 180/255, 180/255),
        "gpos25": (200/255, 200/255, 200/255),
        "gvar": (220/255, 220/255, 220/255),
        "gneg": (255/255, 255/255, 255/255),
        "acen": (217/255, 217/255, 217/255),
        "stalk": (100/255, 127/255, 164/255)
    }
    cytobands["COLOR"] = cytobands["STAIN"].map(color_dict)
    
    # Get chromosome sizes from cytoband data
    chr_sizes = {}
    for chr_name in cytobands["CHR"].unique():
        chr_bands = cytobands[cytobands["CHR"] == chr_name]
        chr_sizes[chr_name] = chr_bands["END"].max()
    
    # Convert chromosome names to numbers for sorting
    chr_to_num = get_chr_to_number(out_chr=True, xymt=["X", "Y", "MT"])
    num_to_chr = get_number_to_chr()
    
    # Get numeric chromosome list (1-22, X, Y, MT)
    chr_list = []
    for chr_name in sorted(cytobands["CHR"].unique()):
        # Extract number from chr name (e.g., "chr1" -> 1)
        if chr_name.startswith("chr"):
            chr_num_str = chr_name[3:]
            try:
                if chr_num_str in ["X", "Y", "MT"]:
                    chr_list.append(chr_name)
                else:
                    chr_num = int(chr_num_str)
                    if 1 <= chr_num <= 22:
                        chr_list.append(chr_name)
            except ValueError:
                continue
    
    # Sort chromosomes numerically (1-22, then X, Y, MT)
    def chr_sort_key(chr_name):
        if not chr_name.startswith("chr"):
            return (999, chr_name)
        chr_num_str = chr_name[3:]
        if chr_num_str.isdigit():
            return (int(chr_num_str), chr_name)
        elif chr_num_str == "X":
            return (23, chr_name)
        elif chr_num_str == "Y":
            return (24, chr_name)
        elif chr_num_str == "MT":
            return (25, chr_name)
        else:
            return (999, chr_name)
    
    chr_list = sorted(chr_list, key=chr_sort_key)
    
    # Limit to autosomes (1-22) for now - can be extended later
    chr_list = [c for c in chr_list if c.startswith("chr") and c[3:].isdigit() and 1 <= int(c[3:]) <= 22]
    n_chr = len(chr_list)
    
    log.write(" -Plotting {} chromosomes...".format(n_chr), verbose=verbose)
    
    # Set up style
    style = set_plot_style(
        plot="plot_phenogram",
        fig_kwargs=fig_kwargs if fig_kwargs is not None else {},
        save_kwargs=save_kwargs if save_kwargs is not None else {},
        save=save,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    
    # Create figure
    if "figsize" not in fig_kwargs:
        fig_kwargs["figsize"] = figsize
    if "dpi" not in fig_kwargs:
        fig_kwargs["dpi"] = dpi
    
    fig, axes = plt.subplots(nrows=1, ncols=ncols, figsize=fig_kwargs["figsize"], 
                             dpi=fig_kwargs["dpi"])
    
    # Handle single axis case
    if ncols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Calculate offsets for each chromosome
    max_chr_size = max(chr_sizes.values())
    offset = [0 for _ in range(ncols)]
    max_row_offset = 0
    row_n = 0
    
    # Prepare leads data for annotation
    leads_dict = {}
    if len(leads) > 0:
        # Convert chromosome format in leads to match cytoband format
        leads_copy = leads.copy()
        if chrom in leads_copy.columns and pos in leads_copy.columns:
            # Convert numeric chromosomes to chr format
            for idx, row in leads_copy.iterrows():
                chr_val = row[chrom]
                pos_val = row[pos]
                if pd.notna(chr_val) and pd.notna(pos_val):
                    # Convert to chr format
                    if isinstance(chr_val, (int, float)):
                        chr_str = "chr{}".format(int(chr_val))
                    else:
                        chr_str = str(chr_val)
                        if not chr_str.startswith("chr"):
                            # Try to extract number
                            try:
                                chr_num = int(chr_str)
                                chr_str = "chr{}".format(chr_num)
                            except (ValueError, TypeError):
                                # If it's already a string like "X", "Y", etc., add "chr" prefix
                                chr_str = "chr{}".format(chr_str)
                    
                    # Only include if it's a valid chromosome (1-22, X, Y, MT)
                    if chr_str in chr_list or (chr_str.startswith("chr") and 
                                               (chr_str[3:].isdigit() or chr_str[3:] in ["X", "Y", "MT"])):
                        leads_dict.setdefault(chr_str, []).append({
                            'pos': int(pos_val) if pd.notna(pos_val) else None,
                            'snpid': row[snpid] if snpid in leads_copy.columns and pd.notna(row[snpid]) else None,
                            'index': idx
                        })
    
    # Plot each chromosome
    chr_bottom_positions = {}  # Store bottom position for each chromosome for label placement
    for i, chr_name in enumerate(chr_list):
        chr_cytobands = cytobands.loc[cytobands["CHR"] == chr_name, :].copy()
        
        if i // ncols > row_n:
            row_n += 1
            max_row_offset += 1.2
        
        chr_size = chr_sizes[chr_name]
        
        # Get the offset for this specific chromosome
        chr_offset = max_row_offset
        
        # Calculate bottom position of chromosome for label placement
        # The chromosome extends from offset to offset + chr_size/max_chr_size
        # Plus telomere at bottom (0.02)
        # Since y-axis is inverted, bottom is at higher y value
        telemere_full_length = 0.02
        # Calculate the actual bottom of the chromosome including telomere
        # The telomere is centered at offset + chr_size/max_chr_size
        # and extends telemere_full_length/2 below that center
        # So the bottom edge is at: offset + chr_size/max_chr_size + telemere_full_length/2
        chr_bottom_with_telomere = chr_size / max_chr_size + chr_offset + telemere_full_length / 2
        chr_bottom_positions[i] = chr_bottom_with_telomere
        
        # Get centromere boundaries
        acen_bands = chr_cytobands.loc[chr_cytobands["STAIN"] == "acen", :]
        if len(acen_bands) > 0:
            chr_centromere_u = acen_bands["START"].min()
            chr_centromere_l = acen_bands["END"].max()
        else:
            # Default centromere position if not found
            chr_centromere_u = chr_size * 0.44
            chr_centromere_l = chr_size * 0.46
        
        # Plot chromosome
        _plot_chr(
            axes[i % ncols],
            chr_size=chr_size,
            chr_centromere_u=chr_centromere_u,
            chr_centromere_l=chr_centromere_l,
            max_chr_size=max_chr_size,
            offset=max_row_offset,
            chr_cytobands=chr_cytobands,
            chr_name=chr_name,
            leads=leads_dict.get(chr_name, []),
            annotate_snps=annotate_snps,
            annotation_kwargs=annotation_kwargs if annotation_kwargs else {},
            anno_style=anno_style,
            repel_force=repel_force,
            anno_max_iter=anno_max_iter,
            log=log,
            verbose=verbose
        )
    
    # Set axis properties
    for i in range(n_chr):
        axes[i % ncols].set_ylim(-0.1, max_row_offset + 1)
        axes[i % ncols].invert_yaxis()
        axes[i % ncols].spines['top'].set_visible(False)
        axes[i % ncols].spines['right'].set_visible(False)
        axes[i % ncols].spines['bottom'].set_visible(False)
        axes[i % ncols].spines['left'].set_visible(False)
        # Add chromosome label - position right below the chromosome (avoiding telomere overlap)
        chr_num = chr_list[i][3:] if chr_list[i].startswith("chr") else chr_list[i]
        # Use stored bottom position for this chromosome
        if i in chr_bottom_positions:
            chr_bottom_with_telomere = chr_bottom_positions[i]
            # Position label below the telomere to avoid overlap
            # Since y-axis is inverted, "below" means higher y value
            # chr_bottom_with_telomere is already the bottom edge of telomere
            # Add padding (0.03) to ensure no overlap
            label_y = chr_bottom_with_telomere + 0.03
            axes[i % ncols].text(0.1, label_y, chr_num, ha='center', va='bottom', 
                                fontsize=10, zorder=300)
    
    # Hide unused axes
    for i in range(n_chr, len(axes)):
        axes[i].set_visible(False)
    
    fig.tight_layout()
    
    # Save figure
    save_figure(fig=fig, save=save, keyword="phenogram", save_kwargs=save_kwargs, 
                log=log, verbose=verbose)
    
    log.write("Finished creating phenogram plot.", verbose=verbose)
    
    return fig


def _plot_chr(
    ax,
    chr_size: float,
    chr_centromere_u: float,
    chr_centromere_l: float,
    max_chr_size: float,
    offset: float,
    chr_cytobands: pd.DataFrame,
    chr_name: str,
    leads: List[Dict[str, Any]],
    annotate_snps: bool = True,
    annotation_kwargs: Dict[str, Any] = None,
    anno_style: str = "expand",
    repel_force: float = 0.02,
    anno_max_iter: int = 100,
    log: Log = Log(),
    verbose: bool = True
):
    """
    Plot a single chromosome with cytobands and optional SNP annotations.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to plot on
    chr_size : float
        Size of the chromosome
    chr_centromere_u : float
        Upper boundary of centromere
    chr_centromere_l : float
        Lower boundary of centromere
    max_chr_size : float
        Maximum chromosome size (for normalization)
    offset : float
        Y-axis offset for this chromosome
    chr_cytobands : pd.DataFrame
        Cytoband data for this chromosome
    chr_name : str
        Chromosome name
    leads : list
        List of lead SNP dictionaries with 'pos', 'snpid', 'index' keys
    annotate_snps : bool
        If True, annotate SNP positions
    annotation_kwargs : dict
        Additional keyword arguments for annotations
    log : Log
        Logging object
    verbose : bool
        If True, print progress messages
    """
    
    if annotation_kwargs is None:
        annotation_kwargs = {}
    
    positions = [
        0 + offset,
        chr_centromere_u / max_chr_size + offset,
        chr_centromere_l / max_chr_size + offset,
        chr_size / max_chr_size + offset
    ]
    
    centromere_full_length = (chr_centromere_l - chr_centromere_u) / max_chr_size
    telemere_full_length = 0.02
    
    height_for_arm1 = positions[1] - positions[0]
    height_for_arm2 = positions[3] - positions[2]
    
    full_length = height_for_arm1 + height_for_arm2 + centromere_full_length
    
    # Draw chromosome arms
    arms = [
        Rectangle((0, positions[0]), width=0.21, height=height_for_arm1),
        Rectangle((0, positions[1] + centromere_full_length), width=0.21, height=height_for_arm2)
    ]
    
    # Draw centromeres
    centromeres = [
        Ellipse(xy=(0.1, positions[1] + centromere_full_length), width=0.2, height=centromere_full_length),
        Ellipse(xy=(0.1, positions[1]), width=0.2, height=centromere_full_length)
    ]
    
    # Draw telomeres
    telemeres = [
        Ellipse(xy=(0.105, positions[0]), width=0.21, height=telemere_full_length),
        Ellipse(xy=(0.105, positions[0] + full_length), width=0.21, height=telemere_full_length)
    ]
    
    # Create patch collections
    arms_pc = PatchCollection(arms, facecolor="white", edgecolor="grey", zorder=100)
    centromeres_pc = PatchCollection(centromeres, facecolor="grey", edgecolor="black", zorder=99)
    telemeres_pc = PatchCollection(telemeres, facecolor="white", edgecolor="grey", zorder=1)
    
    # Add collections to axes
    ax.add_collection(arms_pc)
    ax.add_collection(centromeres_pc)
    ax.add_collection(telemeres_pc)
    
    # Draw cytobands
    for index, row in chr_cytobands.iterrows():
        if row["END"] <= chr_centromere_u:
            # Upper arm
            band_start = row["START"] / max_chr_size + offset
            band_end = row["END"] / max_chr_size + offset
        elif row["START"] >= chr_centromere_l:
            # Lower arm
            band_start = (row["START"] - chr_centromere_l) / max_chr_size + height_for_arm1 + offset + centromere_full_length
            band_end = (row["END"] - chr_centromere_l) / max_chr_size + height_for_arm1 + offset + centromere_full_length
        else:
            # Skip centromere bands (already drawn)
            continue
        
        band_height = band_end - band_start
        
        if band_height > 0:
            band = Rectangle((0.01, band_start), width=0.2, height=band_height)
            facecolor = row["COLOR"]
            
            if row["STAIN"] == "stalk":
                bands_pc = PatchCollection([band], facecolor=facecolor, edgecolor=None, 
                                           linewidths=0, zorder=102)
            else:
                bands_pc = PatchCollection([band], facecolor=facecolor, edgecolor=None, 
                                         linewidths=0, zorder=101)
            ax.add_collection(bands_pc)
    
    # Annotate lead SNPs with arrows using expand style (horizontal layout)
    if annotate_snps and len(leads) > 0:
        # Filter and sort leads by position
        leads_sorted = sorted([l for l in leads if l['pos'] is not None and pd.notna(l['pos'])], 
                             key=lambda x: x['pos'])
        
        if len(leads_sorted) == 0:
            return
        
        # Calculate y positions on chromosome for all leads
        leads_with_y = []
        for lead in leads_sorted:
            pos = lead['pos']
            
            # Calculate y position on chromosome based on genomic position
            if pos <= chr_centromere_u:
                # Upper arm
                y_pos_chr = pos / max_chr_size + offset
            elif pos >= chr_centromere_l:
                # Lower arm
                y_pos_chr = (pos - chr_centromere_l) / max_chr_size + height_for_arm1 + offset + centromere_full_length
            else:
                # In centromere region, skip annotation
                continue
            
            leads_with_y.append({
                'lead': lead,
                'y_pos_chr': y_pos_chr,
                'original_y': y_pos_chr
            })
        
        if len(leads_with_y) == 0:
            return
        
        # Extract y positions for adjustment
        y_positions = np.array([l['y_pos_chr'] for l in leads_with_y])
        
        # Always adjust positions to avoid overlap (using expand style logic)
        # Calculate y span for this chromosome
        if len(y_positions) > 1:
            y_span = max(y_positions) - min(y_positions)
            # If variants are very close together, use a minimum span
            if y_span < 0.01:
                y_span = 0.1
        else:
            y_span = 0.1
        
        # Adjust positions using adjust_text_position to prevent overlap
        adjusted_y = adjust_text_position(
            y_positions.copy(),
            y_span,
            repel_force=repel_force,
            max_iter=anno_max_iter,
            amode="float",
            log=log,
            verbose=verbose
        )
        
        # Update y positions with adjusted values
        for i, l in enumerate(leads_with_y):
            l['y_pos_chr'] = adjusted_y[i]
        
        # Base x position for annotations (on the right side)
        annotation_x_base = 0.5
        
        # Calculate fixed armB for this chromosome (vertical distance from annotation area to chromosome)
        # Use a representative y position (middle of annotation area) to calculate the distance
        if len(leads_with_y) > 0:
            # Use the first lead's adjusted y position as reference
            ref_y = leads_with_y[0]['y_pos_chr']
            ref_original_y = leads_with_y[0]['original_y']
        else:
            ref_y = offset + chr_size / max_chr_size / 2
            ref_original_y = ref_y
        
        # Calculate fixed armB in pixels (horizontal distance)
        transform = ax.transData
        annotation_pixel = transform.transform((annotation_x_base, ref_y))
        # Arrow points to right edge of chromosome (x=0.21)
        chr_right_edge_pixel = transform.transform((0.21, ref_original_y))
        # Fixed horizontal distance from annotation area to right edge of chromosome
        armB_length_pixel = abs(annotation_pixel[0] - chr_right_edge_pixel[0])
        
        # Create annotation objects
        anno_objects = []
        
        for l in leads_with_y:
            lead = l['lead']
            y_pos_chr = l['y_pos_chr']
            original_y = l['original_y']
            
            # Plot marker on chromosome at original position (as a horizontal line)
            # Draw a horizontal line across the full width of the chromosome
            # Chromosome rectangle: Rectangle((0, y), width=0.21) extends from x=0 to x=0.21
            # The line should exactly match the rectangle edges
            # Use exact coordinates to align with rectangle boundaries
            ax.plot([0.0, 0.21], [original_y, original_y], 
                   color="red", linewidth=2, clip_on=False, zorder=200)
            
            if lead['snpid'] is not None:
                snp_text = str(lead['snpid'])
                # Truncate long SNP IDs
                if len(snp_text) > 20:
                    snp_text = snp_text[:17] + "..."
                
                # Point on chromosome (right edge where arrow points)
                chr_point_x = 0.21  # Right edge of chromosome
                chr_point_y = original_y
                
                # Text position (on the right, using adjusted y)
                text_x = annotation_x_base
                text_y = y_pos_chr
                
                # Set default annotation style
                anno_default = {
                    "fontsize": 8,
                    "ha": "left",
                    "va": "center",
                    "fontweight": "normal"
                }
                if annotation_kwargs:
                    anno_default.update(annotation_kwargs)
                
                # Create annotation with arrow
                # armA = 0 (no horizontal arm, go directly horizontal)
                # armB = fixed for this chromosome (horizontal distance to chromosome)
                anno_obj = ax.annotate(
                    snp_text,
                    xy=(chr_point_x, chr_point_y),  # Point on chromosome (arrow target) - FIXED
                    xytext=(text_x, text_y),  # Text position (arrow start)
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color="gray",
                        lw=0.8,
                        # No initial arm (armA=0), fixed horizontal distance (armB)
                        # angleA=180: point left from text, angleB=0: point horizontally to chromosome
                        connectionstyle=f"arc,angleA=180,armA=0,angleB=0,armB={armB_length_pixel},rad=0"
                    ),
                    clip_on=False,
                    zorder=201,
                    **anno_default
                )
                anno_objects.append(anno_obj)
        
        # Use adjust_text for final fine-tuning (vertical movement only for horizontal layout)
        # Note: We skip adjust_text for now as it interferes with arrow positioning
        # The initial adjust_text_position should be sufficient for spacing
        # If needed, we can enable this but would need to recalculate arrows after adjustment
        # if anno_style == "expand" and len(anno_objects) > 1:
        #     log.write(" -Auto-adjusting annotation positions...", verbose=verbose)
        #     adjust_text(
        #         texts=anno_objects,
        #         autoalign=False,
        #         only_move={'points': 'y', 'text': 'y', 'objects': 'y'},  # Only move vertically
        #         ax=ax,
        #         precision=0.02,
        #         force_text=(repel_force, repel_force),
        #         expand_text=(1, 1),
        #         expand_objects=(0, 0),
        #         expand_points=(0, 0),
        #         va="center",
        #         ha='left',
        #         avoid_points=False,
        #         lim=100
        #     )
    
    ax.set_xticks(ticks=[])
    ax.set_yticks(ticks=[])
    # Extend xlim to accommodate annotations on the right
    # Need space for chromosome (0-0.21), gap, and annotations (up to 0.5+)
    ax.set_xlim([0, 0.7])
