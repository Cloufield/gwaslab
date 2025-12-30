"""
Report generation module for GWASLab.

This module provides functions to generate comprehensive QC reports
including basic QC, lead variant extraction, and visualization.
"""

from typing import TYPE_CHECKING, Optional, Dict, Any, List
import os
import pandas as pd
from pathlib import Path
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_filter_value import _get_region_start_and_end

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

# Try to import PDF generation library
try:
    import weasyprint
    HAS_WEASYPRINT = True
except ImportError:
    HAS_WEASYPRINT = False
    weasyprint = None


def generate_qc_report(
    sumstats: 'Sumstats',
    output_path: str = "gwas_qc_report.html",
    basic_check_kwargs: Optional[Dict[str, Any]] = None,
    harmonize_kwargs: Optional[Dict[str, Any]] = None,
    get_lead_kwargs: Optional[Dict[str, Any]] = None,
    mqq_plot_kwargs: Optional[Dict[str, Any]] = None,
    regional_plot_kwargs: Optional[Dict[str, Any]] = None,
    output_kwargs: Optional[Dict[str, Any]] = None,
    report_title: str = "GWAS Quality Control Report",
    verbose: bool = True
) -> str:
    """
    Generate a comprehensive QC report including basic QC, harmonization (optional),
    lead variants, MQQ plots, regional plots, and output (optional).
    
    This function performs:
    1. Basic QC using basic_check()
    2. Harmonization using harmonize() (optional, if harmonize_kwargs provided)
    3. Lead variant extraction using get_lead()
    4. MQQ plot generation
    5. Regional plots for each lead variant locus
    6. Output to specified format using to_format() (optional, if output_kwargs provided)
    7. HTML/PDF report generation with all results
    
    Parameters
    ----------
    sumstats : gwaslab.Sumstats
        Sumstats object to analyze
    output_path : str, optional
        Path to save the report. Supports both HTML (.html) and PDF (.pdf) formats.
        Default: "gwas_qc_report.html"
    basic_check_kwargs : dict, optional
        Keyword arguments passed to basic_check(). Default: None
    harmonize_kwargs : dict, optional
        Keyword arguments passed to harmonize(). If provided, harmonization will be performed
        after basic_check. Set to {} to use default harmonization settings.
        Default: None (no harmonization)
    get_lead_kwargs : dict, optional
        Keyword arguments passed to get_lead(). Default: None
    mqq_plot_kwargs : dict, optional
        Keyword arguments passed to plot_mqq() for MQQ plot. Default: None
    regional_plot_kwargs : dict, optional
        Keyword arguments passed to plot_mqq() for regional plots. Default: None
    output_kwargs : dict, optional
        Keyword arguments passed to to_format() for outputting sumstats.
        Must include 'path' key. If provided, sumstats will be saved to the specified format.
        Example: {"path": "output_sumstats", "fmt": "ldsc", "gzip": True}
        Default: None (no output)
    report_title : str, optional
        Title for the HTML report. Default: "GWAS Quality Control Report"
    verbose : bool, optional
        Whether to print progress messages. Default: True
    
    Returns
    -------
    str
        Path to the generated report (HTML or PDF)
    
    Notes
    -----
    For PDF output, `weasyprint` must be installed.
    If PDF format is requested but weasyprint is not available, the function will
    generate HTML instead and issue a warning.
    
    Examples
    --------
    >>> import gwaslab as gl
    >>> mysumstats = gl.Sumstats("sumstats.txt.gz")
    >>> # Basic report
    >>> gl.generate_qc_report(
    ...     mysumstats,
    ...     output_path="my_report.html",
    ...     get_lead_kwargs={"sig_level": 5e-8, "windowsizekb": 500}
    ... )
    >>> # With harmonization
    >>> gl.generate_qc_report(
    ...     mysumstats,
    ...     output_path="my_report.html",
    ...     harmonize_kwargs={"ref_seq": "ref.fa", "ref_infer": "ref.vcf.gz"}
    ... )
    >>> # With output
    >>> gl.generate_qc_report(
    ...     mysumstats,
    ...     output_path="my_report.html",
    ...     output_kwargs={"path": "clean_sumstats", "fmt": "ldsc", "gzip": True}
    ... )
    """
    log = Log()
    
    if basic_check_kwargs is None:
        basic_check_kwargs = {}
    if harmonize_kwargs is None:
        harmonize_kwargs = None  # Keep as None to indicate no harmonization
    elif harmonize_kwargs == {}:
        harmonize_kwargs = {}  # Empty dict means use defaults
    if get_lead_kwargs is None:
        get_lead_kwargs = {}
    if mqq_plot_kwargs is None:
        mqq_plot_kwargs = {}
    if regional_plot_kwargs is None:
        regional_plot_kwargs = {}
    if output_kwargs is None:
        output_kwargs = None  # Keep as None to indicate no output
    
    # Create output directory if needed
    output_path_obj = Path(output_path)
    output_dir = output_path_obj.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a directory for plot images
    plots_dir = output_dir / f"{output_path_obj.stem}_plots"
    plots_dir.mkdir(exist_ok=True)
    
    log.write("=" * 80, verbose=verbose)
    log.write(f"Generating QC Report: {report_title}", verbose=verbose)
    log.write("=" * 80, verbose=verbose)
    
    # Track processing steps for report
    processing_steps = []
    output_files = []
    
    # Step 1: Basic QC
    step_num = 1
    total_steps = 4 + (1 if harmonize_kwargs is not None else 0) + (1 if output_kwargs is not None else 0)
    log.write(f"\n[Step {step_num}/{total_steps}] Running basic QC...", verbose=verbose)
    sumstats.basic_check(**basic_check_kwargs)
    log.write("Basic QC completed.", verbose=verbose)
    processing_steps.append("Basic QC")
    step_num += 1
    
    # Step 2: Harmonization (optional)
    harmonization_performed = False
    if harmonize_kwargs is not None:
        log.write(f"\n[Step {step_num}/{total_steps}] Running harmonization...", verbose=verbose)
        # If basic_check was already done, skip it in harmonize
        if "basic_check" not in harmonize_kwargs:
            harmonize_kwargs["basic_check"] = False
        sumstats.harmonize(**harmonize_kwargs)
        log.write("Harmonization completed.", verbose=verbose)
        processing_steps.append("Harmonization")
        harmonization_performed = True
        step_num += 1
    
    # Step 3: Get lead variants
    log.write(f"\n[Step {step_num}/{total_steps}] Extracting lead variants...", verbose=verbose)
    lead_variants = sumstats.get_lead(**get_lead_kwargs)
    
    if lead_variants is None or len(lead_variants) == 0:
        log.write("WARNING: No lead variants found. Report will contain only MQQ plot.", verbose=verbose)
        lead_variants = pd.DataFrame()
    else:
        log.write(f"Found {len(lead_variants)} lead variant(s).", verbose=verbose)
    
    # Step 4: Create MQQ plot
    step_num += 1
    log.write(f"\n[Step {step_num}/{total_steps}] Creating MQQ plot...", verbose=verbose)
    mqq_plot_path = plots_dir / "mqq_plot.png"
    
    # Set default save path for MQQ plot
    mqq_kwargs = mqq_plot_kwargs.copy()
    if "save" not in mqq_kwargs:
        mqq_kwargs["save"] = str(mqq_plot_path)
    if "save_kwargs" not in mqq_kwargs:
        mqq_kwargs["save_kwargs"] = {"dpi": 300, "bbox_inches": "tight", "facecolor": "white"}
    
    # Create MQQ plot
    try:
        sumstats.plot_mqq(**mqq_kwargs)
        # Verify file was created
        if mqq_plot_path.exists():
            log.write(f"MQQ plot saved to {mqq_plot_path}", verbose=verbose)
        else:
            log.write(f"WARNING: MQQ plot file not found at {mqq_plot_path}", verbose=verbose)
            mqq_plot_path = None
    except Exception as e:
        log.write(f"WARNING: Failed to create MQQ plot: {e}", verbose=verbose)
        mqq_plot_path = None
    
    # Step 5: Create regional plots for each lead variant
    step_num += 1
    log.write(f"\n[Step {step_num}/{total_steps}] Creating regional plots for lead variants...", verbose=verbose)
    regional_plots = []
    
    if len(lead_variants) > 0:
        # Get default window size for regional plots
        windowsizekb = get_lead_kwargs.get("windowsizekb", 500)
        
        for idx, row in lead_variants.iterrows():
            try:
                # Get chromosome and position
                chrom = row.get("CHR", None)
                pos = row.get("POS", None)
                
                if chrom is None or pd.isna(chrom) or pos is None or pd.isna(pos):
                    log.write(f"WARNING: Skipping lead variant {idx} - missing CHR or POS", verbose=verbose)
                    continue
                
                # Get region coordinates
                region = _get_region_start_and_end(
                    chrom=chrom,
                    pos=pos,
                    windowsizekb=windowsizekb,
                    verbose=False,
                    log=log
                )
                
                # Create regional plot
                region_plot_path = plots_dir / f"regional_plot_chr{chrom}_{pos}.png"
                
                # Prepare regional plot kwargs
                reg_kwargs = regional_plot_kwargs.copy()
                reg_kwargs["mode"] = "r"
                reg_kwargs["region"] = region
                if "save" not in reg_kwargs:
                    reg_kwargs["save"] = str(region_plot_path)
                if "save_kwargs" not in reg_kwargs:
                    reg_kwargs["save_kwargs"] = {"dpi": 300, "bbox_inches": "tight", "facecolor": "white"}
                
                # Create regional plot
                sumstats.plot_mqq(**reg_kwargs)
                
                # Verify file was created
                if not region_plot_path.exists():
                    log.write(f"WARNING: Regional plot file not found at {region_plot_path}", verbose=verbose)
                    continue
                
                # Store plot info
                variant_id = row.get("SNPID", f"chr{chrom}:{pos}")
                regional_plots.append({
                    "variant_id": variant_id,
                    "chrom": chrom,
                    "pos": pos,
                    "p": row.get("P", None),
                    "mlog10p": row.get("MLOG10P", None),
                    "beta": row.get("BETA", None),
                    "se": row.get("SE", None),
                    "plot_path": region_plot_path,
                    "region": region
                })
                
                log.write(f"Created regional plot for {variant_id} (chr{chrom}:{pos})", verbose=verbose)
                
            except Exception as e:
                log.write(f"WARNING: Failed to create regional plot for variant {idx}: {e}", verbose=verbose)
                continue
    
    log.write(f"\nCreated {len(regional_plots)} regional plot(s).", verbose=verbose)
    
    # Step 6: Generate summary
    log.write("\nGenerating summary statistics...", verbose=verbose)
    summary_dict = sumstats.summary()
    log.write("Summary generated.", verbose=verbose)
    
    # Step 7: Output sumstats (optional)
    if output_kwargs is not None:
        step_num += 1
        log.write(f"\n[Step {step_num}/{total_steps}] Outputting sumstats...", verbose=verbose)
        if "path" not in output_kwargs:
            log.write("WARNING: output_kwargs must include 'path' key. Skipping output.", verbose=verbose)
        else:
            # Create a copy of output_kwargs to avoid modifying the original
            output_kwargs_copy = output_kwargs.copy()
            output_path_sumstats = output_kwargs_copy.pop("path")
            sumstats.to_format(path=output_path_sumstats, **output_kwargs_copy)
            log.write(f"Sumstats saved to: {output_path_sumstats}", verbose=verbose)
            processing_steps.append("Output")
            output_files.append(output_path_sumstats)
    
    # Determine output format from file extension
    output_format = output_path_obj.suffix.lower()
    
    # Generate HTML report (always generate HTML first, then convert if needed)
    log.write("\nGenerating HTML report...", verbose=verbose)
    
    # Get log text
    log_text = sumstats.log.log_text if hasattr(sumstats, 'log') and hasattr(sumstats.log, 'log_text') else ""
    
    html_content = _generate_html_report(
        sumstats=sumstats,
        lead_variants=lead_variants,
        mqq_plot_path=mqq_plot_path,
        regional_plots=regional_plots,
        report_title=report_title,
        output_dir=output_dir,
        plots_dir=plots_dir,
        use_absolute_paths=False,
        processing_steps=processing_steps,
        harmonization_performed=harmonization_performed,
        output_files=output_files,
        summary_dict=summary_dict,
        log_text=log_text
    )
    
    # Save HTML report (always save HTML as intermediate)
    html_path = output_path_obj.with_suffix('.html')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    # Convert to PDF if requested
    if output_format == '.pdf':
        pdf_path = output_path_obj
        log.write("\nConverting HTML to PDF...", verbose=verbose)
        
        if HAS_WEASYPRINT:
            try:
                # Convert HTML to PDF using weasyprint
                # base_url allows weasyprint to resolve relative paths for images
                weasyprint.HTML(string=html_content, base_url=str(output_dir)).write_pdf(pdf_path)
                log.write(f"PDF report saved to: {pdf_path}", verbose=verbose)
                log.write("=" * 80, verbose=verbose)
                return str(pdf_path)
            except Exception as e:
                log.write(f"WARNING: Failed to generate PDF with weasyprint: {e}", verbose=verbose)
                log.write("  Note: weasyprint may require additional system dependencies.", verbose=verbose)
                log.write("  On Ubuntu/Debian: sudo apt-get install python3-weasyprint", verbose=verbose)
                log.write("  Or install via pip: pip install weasyprint", verbose=verbose)
                log.write("HTML report saved instead.", verbose=verbose)
                log.write(f"HTML report saved to: {html_path}", verbose=verbose)
                log.write("=" * 80, verbose=verbose)
                return str(html_path)
        else:
            log.write("WARNING: weasyprint is not installed. HTML report saved instead.", verbose=verbose)
            log.write("To generate PDF reports, install weasyprint:", verbose=verbose)
            log.write("  pip install weasyprint", verbose=verbose)
            log.write(f"HTML report saved to: {html_path}", verbose=verbose)
            log.write("=" * 80, verbose=verbose)
            return str(html_path)
    else:
        # HTML format
        log.write(f"\nReport saved to: {output_path}", verbose=verbose)
        log.write("=" * 80, verbose=verbose)
        return str(output_path)


def _generate_html_report(
    sumstats: 'Sumstats',
    lead_variants: pd.DataFrame,
    mqq_plot_path: Optional[Path],
    regional_plots: Dict[str, Path],
    report_title: str,
    output_dir: Path,
    plots_dir: Path,
    use_absolute_paths: bool = False,
    processing_steps: Optional[List[str]] = None,
    harmonization_performed: bool = False,
    output_files: Optional[List[str]] = None,
    summary_dict: Optional[Dict[str, Any]] = None,
    log_text: str = ""
) -> str:
    """
    Generate HTML content for the QC report.
    
    Parameters
    ----------
    sumstats : gwaslab.Sumstats
        Sumstats object
    lead_variants : pd.DataFrame
        DataFrame containing lead variants
    mqq_plot_path : Path or None
        Path to MQQ plot image
    regional_plots : list
        List of dictionaries containing regional plot information
    report_title : str
        Title for the report
    output_dir : Path
        Output directory
    plots_dir : Path
        Directory containing plot images
    use_absolute_paths : bool, optional
        Whether to use absolute paths for images (deprecated, kept for compatibility)
    processing_steps : list, optional
        List of processing steps performed
    harmonization_performed : bool, optional
        Whether harmonization was performed
    output_files : list, optional
        List of output file paths
    
    Returns
    -------
    str
        HTML content as string
    """
    if processing_steps is None:
        processing_steps = []
    if output_files is None:
        output_files = []
    html_parts = []
    
    # HTML header
    html_parts.append("""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #7f8c8d;
            margin-top: 20px;
        }}
        .section {{
            margin: 20px 0;
        }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 4px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        .summary {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .summary-item {{
            margin: 5px 0;
        }}
        .variant-info {{
            background-color: #f8f9fa;
            padding: 15px;
            border-left: 4px solid #3498db;
            margin: 15px 0;
        }}
        .timestamp {{
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 30px;
        }}
        .log-container {{
            background-color: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
            font-family: 'Courier New', monospace;
            font-size: 0.85em;
            max-height: 500px;
            overflow-y: auto;
            white-space: pre-wrap;
            word-wrap: break-word;
        }}
        .summary-table {{
            margin: 20px 0;
        }}
        .summary-table table {{
            width: 100%;
        }}
        .summary-category {{
            margin: 20px 0;
        }}
        .summary-category h3 {{
            color: #3498db;
            border-bottom: 1px solid #3498db;
            padding-bottom: 5px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{}</h1>
        <div class="timestamp">Generated: {}</div>
    """.format(
        report_title,
        report_title,
        pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    ))
    
    # Summary section
    html_parts.append("""
        <div class="section">
            <h2>Summary</h2>
            <div class="summary">
    """)
    
    # Add basic statistics
    if hasattr(sumstats, 'data') and sumstats.data is not None:
        n_variants = len(sumstats.data)
        html_parts.append(f'<div class="summary-item"><strong>Total Variants:</strong> {n_variants:,}</div>')
    
    if lead_variants is not None and len(lead_variants) > 0:
        html_parts.append(f'<div class="summary-item"><strong>Lead Variants:</strong> {len(lead_variants)}</div>')
        if "P" in lead_variants.columns:
            min_p = lead_variants["P"].min()
            html_parts.append(f'<div class="summary-item"><strong>Minimum P-value:</strong> {min_p:.2e}</div>')
    else:
        html_parts.append('<div class="summary-item"><strong>Lead Variants:</strong> 0</div>')
    
    # Add processing steps
    if processing_steps:
        html_parts.append(f'<div class="summary-item"><strong>Processing Steps:</strong> {", ".join(processing_steps)}</div>')
    
    # Add harmonization info
    if harmonization_performed:
        html_parts.append('<div class="summary-item"><strong>Harmonization:</strong> Performed</div>')
    
    # Add output files
    if output_files:
        html_parts.append(f'<div class="summary-item"><strong>Output Files:</strong> {", ".join(output_files)}</div>')
    
    html_parts.append("""
            </div>
        </div>
    """)
    
    # MQQ Plot section
    if mqq_plot_path and mqq_plot_path.exists():
        html_parts.append("""
        <div class="section">
            <h2>Manhattan-QQ Plot</h2>
            <div class="plot-container">
        """)
        
        # Use relative path for images (weasyprint can resolve relative paths via base_url)
        img_path = os.path.relpath(mqq_plot_path, output_dir)
        html_parts.append(f'<img src="{img_path}" alt="MQQ Plot">')
        
        html_parts.append("""
            </div>
        </div>
        """)
    
    # Lead Variants Table
    if lead_variants is not None and len(lead_variants) > 0:
        html_parts.append("""
        <div class="section">
            <h2>Lead Variants</h2>
            <table>
        """)
        
        # Table header
        columns_to_show = ["SNPID", "CHR", "POS", "EA", "NEA", "P", "MLOG10P", "BETA", "SE"]
        available_columns = [col for col in columns_to_show if col in lead_variants.columns]
        
        if not available_columns:
            # Fallback to all columns
            available_columns = lead_variants.columns.tolist()[:10]  # Limit to first 10 columns
        
        html_parts.append("<thead><tr>")
        for col in available_columns:
            html_parts.append(f"<th>{col}</th>")
        html_parts.append("</tr></thead><tbody>")
        
        # Table rows
        for idx, row in lead_variants.iterrows():
            html_parts.append("<tr>")
            for col in available_columns:
                value = row.get(col, "")
                if pd.isna(value):
                    value = ""
                elif isinstance(value, float):
                    if col in ["P"]:
                        value = f"{value:.2e}"
                    elif col in ["MLOG10P", "BETA", "SE"]:
                        value = f"{value:.4f}"
                    else:
                        value = f"{value:.4f}"
                else:
                    value = str(value)
                html_parts.append(f"<td>{value}</td>")
            html_parts.append("</tr>")
        
        html_parts.append("</tbody></table></div>")
    
    # Regional Plots section
    if regional_plots:
        html_parts.append("""
        <div class="section">
            <h2>Regional Plots</h2>
        """)
        
        for i, plot_info in enumerate(regional_plots, 1):
            html_parts.append(f"""
            <div class="variant-info">
                <h3>Locus {i}: {plot_info['variant_id']}</h3>
                <p><strong>Location:</strong> chr{plot_info['chrom']}:{plot_info['pos']:,}</p>
            """)
            
            if plot_info.get('p') is not None:
                html_parts.append(f"<p><strong>P-value:</strong> {plot_info['p']:.2e}</p>")
            if plot_info.get('mlog10p') is not None:
                html_parts.append(f"<p><strong>-log10(P):</strong> {plot_info['mlog10p']:.4f}</p>")
            if plot_info.get('beta') is not None:
                html_parts.append(f"<p><strong>BETA:</strong> {plot_info['beta']:.4f}</p>")
            if plot_info.get('se') is not None:
                html_parts.append(f"<p><strong>SE:</strong> {plot_info['se']:.4f}</p>")
            
            if plot_info['plot_path'].exists():
                # Use relative path for images (weasyprint can resolve relative paths via base_url)
                img_path = os.path.relpath(plot_info['plot_path'], output_dir)
                html_parts.append(f"""
                <div class="plot-container">
                    <img src="{img_path}" alt="Regional Plot for {plot_info['variant_id']}">
                </div>
                """)
            
            html_parts.append("</div>")
        
        html_parts.append("</div>")
    
    # Summary Statistics section
    if summary_dict:
        html_parts.append("""
        <div class="section">
            <h2>Summary Statistics</h2>
        """)
        
        # Overview
        if "overview" in summary_dict or "META" in summary_dict:
            overview = summary_dict.get("overview", summary_dict.get("META", {}))
            html_parts.append("""
            <div class="summary-category">
                <h3>Overview</h3>
                <table>
                    <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                    <tbody>
            """)
            for key, value in overview.items():
                html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        # Chromosomes
        if "chromosomes" in summary_dict or "CHR" in summary_dict:
            chr_info = summary_dict.get("chromosomes", summary_dict.get("CHR", {}))
            html_parts.append("""
            <div class="summary-category">
                <h3>Chromosome Distribution</h3>
                <table>
                    <thead><tr><th>Chromosome</th><th>Count</th></tr></thead>
                    <tbody>
            """)
            for key, value in chr_info.items():
                if key.startswith("chr") or key == "Chromosomes_numbers":
                    html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        # Missing values
        if "missing_values" in summary_dict or "MISSING" in summary_dict:
            missing_info = summary_dict.get("missing_values", summary_dict.get("MISSING", {}))
            if missing_info and missing_info.get("Missing_total", 0) > 0:
                html_parts.append("""
                <div class="summary-category">
                    <h3>Missing Values</h3>
                    <table>
                        <thead><tr><th>Column</th><th>Missing Count</th></tr></thead>
                        <tbody>
                """)
                for key, value in missing_info.items():
                    if key.startswith("Missing_"):
                        col_name = key.replace("Missing_", "")
                        html_parts.append(f"<tr><td><strong>{col_name}</strong></td><td>{value:,}</td></tr>")
                html_parts.append("</tbody></table></div>")
        
        # MAF distribution
        if "MAF" in summary_dict:
            maf_info = summary_dict["MAF"]
            html_parts.append("""
            <div class="summary-category">
                <h3>Minor Allele Frequency (MAF) Distribution</h3>
                <table>
                    <thead><tr><th>Category</th><th>Count</th></tr></thead>
                    <tbody>
            """)
            for key, value in maf_info.items():
                html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:,}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        # P-value statistics
        if "p_values" in summary_dict or "P" in summary_dict:
            p_info = summary_dict.get("p_values", summary_dict.get("P", {}))
            html_parts.append("""
            <div class="summary-category">
                <h3>P-value Statistics</h3>
                <table>
                    <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                    <tbody>
            """)
            for key, value in p_info.items():
                if isinstance(value, float):
                    if "Minimum" in key:
                        html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:.2e}</td></tr>")
                    else:
                        html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:,}</td></tr>")
                else:
                    html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        # Variant status
        if "variant_status" in summary_dict:
            status_info = summary_dict["variant_status"]
            html_parts.append("""
            <div class="summary-category">
                <h3>Variant Status Summary</h3>
                <table>
                    <thead><tr><th>Status Code</th><th>Description</th><th>Count</th></tr></thead>
                    <tbody>
            """)
            for key, value in status_info.items():
                if isinstance(value, dict):
                    desc = value.get("explanation", value.get("description", ""))
                    count = value.get("count", 0)
                    html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{desc}</td><td>{count:,}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        # Variants metadata
        if "variants" in summary_dict:
            variants_info = summary_dict["variants"]
            html_parts.append("""
            <div class="summary-category">
                <h3>Variant Metadata</h3>
                <table>
                    <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                    <tbody>
            """)
            for key, value in variants_info.items():
                if isinstance(value, float):
                    if "P" in key:
                        html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:.2e}</td></tr>")
                    else:
                        html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:.6f}</td></tr>")
                else:
                    html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:,}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        # Samples metadata
        if "samples" in summary_dict:
            samples_info = summary_dict["samples"]
            html_parts.append("""
            <div class="summary-category">
                <h3>Sample Metadata</h3>
                <table>
                    <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                    <tbody>
            """)
            for key, value in samples_info.items():
                if isinstance(value, float):
                    html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:,.0f}</td></tr>")
                else:
                    html_parts.append(f"<tr><td><strong>{key}</strong></td><td>{value:,}</td></tr>")
            html_parts.append("</tbody></table></div>")
        
        html_parts.append("</div>")
    
    # Log section
    if log_text:
        html_parts.append("""
        <div class="section">
            <h2>Processing Log</h2>
            <div class="log-container">
        """)
        # Escape HTML special characters in log text
        import html as html_module
        escaped_log = html_module.escape(log_text)
        html_parts.append(escaped_log)
        html_parts.append("""
            </div>
        </div>
        """)
    
    # Footer
    html_parts.append("""
    </div>
</body>
</html>
    """)
    
    return "\n".join(html_parts)

