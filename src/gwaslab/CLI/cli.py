from typing import Optional, Any
import argparse
import sys
import gwaslab as gl

def _load_sumstats(path: str, fmt: str, nrows: Optional[int]) -> gl.Sumstats:
    return gl.Sumstats(path, fmt=fmt, nrows=nrows) if nrows is not None else gl.Sumstats(path, fmt=fmt)

def _cmd_version(_args: argparse.Namespace) -> None:
    gl.show_version()

def _process_sumstats(args: argparse.Namespace) -> gl.Sumstats:
    """Main processing function that handles QC, harmonization, and formatting"""
    # Load sumstats
    s = _load_sumstats(args.input, args.fmt, args.nrows)
    
    # Determine what processing to do based on options
    do_qc = args.qc or args.remove or args.remove_dup or args.normalize
    do_harmonize = (args.ref_seq is not None or 
                   args.ref_rsid_tsv is not None or 
                   args.ref_rsid_vcf is not None or 
                   args.ref_infer is not None or
                   args.harmonize)
    
    # Perform QC if requested
    if do_qc:
        s.basic_check(
            remove=args.remove,
            remove_dup=args.remove_dup,
            threads=args.threads,
            normalize=args.normalize,
            verbose=not args.quiet
        )
    
    # Perform harmonization if requested
    if do_harmonize:
        harmonize_kwargs = {
            "basic_check": args.basic_check if not do_qc else False,  # Skip if QC already done
            "ref_seq": args.ref_seq,
            "ref_rsid_tsv": args.ref_rsid_tsv,
            "ref_rsid_vcf": args.ref_rsid_vcf,
            "ref_infer": args.ref_infer,
            "ref_alt_freq": args.ref_alt_freq,
            "ref_maf_threshold": args.ref_maf_threshold,
            "maf_threshold": args.maf_threshold,
            "threads": args.threads,
            "remove": args.remove,
            "verbose": not args.quiet,
            "sweep_mode": args.sweep_mode
        }
        # Remove None values
        harmonize_kwargs = {k: v for k, v in harmonize_kwargs.items() if v is not None}
        s.harmonize(**harmonize_kwargs)
    
    # Handle assign-rsid if specified
    if args.assign_rsid:
        s.assign_rsid(
            ref_rsid_tsv=args.ref_rsid_tsv,
            ref_rsid_vcf=args.ref_rsid_vcf,
            threads=args.threads,
            overwrite=args.overwrite
        )
    
    # Handle rsid-to-chrpos if specified
    if args.rsid_to_chrpos:
        # Use HDF5-based parallel processing (faster than old TSV approach)
        # Accepts VCF file (will auto-generate HDF5 path) or direct HDF5 path
        s.rsid_to_chrpos(
            ref_rsid_to_chrpos_vcf=args.ref_rsid_vcf if args.ref_rsid_vcf else None,
            ref_rsid_to_chrpos_hdf5=args.ref_rsid_tsv if args.ref_rsid_tsv else None,  # Reuse TSV arg for HDF5 path
            build=args.build,
            threads=args.threads if args.threads else 4,
            verbose=not args.quiet
        )
    
    # Output formatting
    if args.out is not None:
        to_format_kwargs = {
            "fmt": args.to_fmt,
            "tab_fmt": args.tab_fmt,
            "gzip": not args.no_gzip,
            "bgzip": args.bgzip,
            "tabix": args.tabix,
            "hapmap3": args.hapmap3,
            "exclude_hla": args.exclude_hla,
            "hla_range": (args.hla_lower, args.hla_upper) if args.exclude_hla else None,
            "n": args.n,
            "chr_prefix": args.chr_prefix,
            "xymt_number": args.xymt_number,
            "verbose": not args.quiet
        }
        # Remove None values
        to_format_kwargs = {k: v for k, v in to_format_kwargs.items() if v is not None}
        s.to_format(args.out, **to_format_kwargs)
    
    return s

def main(argv: Optional[list] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    
    # Check if version command
    if "version" in argv:
        parser = argparse.ArgumentParser(prog="gwaslab", add_help=True)
        parser.add_argument("version", nargs="?", help="Show version")
        args = parser.parse_args(argv)
        _cmd_version(args)
        return
    
    parser = argparse.ArgumentParser(
        prog="gwaslab",
        description="GWASLab: A Python package for processing GWAS summary statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic QC and output
  gwaslab --input sumstats.tsv --fmt auto --qc --out cleaned.tsv --to-fmt gwaslab
  
  # Harmonization with reference
  gwaslab --input sumstats.tsv --fmt auto --ref-seq ref.fasta --harmonize --out harmonized.tsv --to-fmt gwaslab
  
  # Format conversion only
  gwaslab --input sumstats.tsv --fmt gwaslab --out sumstats.ldsc --to-fmt ldsc
        """
    )
    
    # Input/Output arguments
    parser.add_argument("--input", required=True, help="Input sumstats file path")
    parser.add_argument("--fmt", default="auto", help="Input format (default: auto)")
    parser.add_argument("--nrows", type=int, help="Number of rows to read (for testing)")
    parser.add_argument("--out", help="Output file path")
    parser.add_argument("--to-fmt", default="gwaslab", help="Output format (default: gwaslab)")
    parser.add_argument("--tab-fmt", default="tsv", choices=["tsv", "csv", "parquet"], help="Tabular format (default: tsv)")
    
    # Processing mode flags
    processing_group = parser.add_argument_group("Processing Options")
    processing_group.add_argument("--qc", action="store_true", help="Perform quality control (basic_check)")
    processing_group.add_argument("--harmonize", action="store_true", help="Perform harmonization")
    processing_group.add_argument("--assign-rsid", action="store_true", help="Assign rsID to variants")
    processing_group.add_argument("--rsid-to-chrpos", action="store_true", help="Convert rsID to CHR:POS")
    
    # QC options
    qc_group = parser.add_argument_group("QC Options")
    qc_group.add_argument("--remove", action="store_true", help="Remove bad quality variants")
    qc_group.add_argument("--remove-dup", action="store_true", help="Remove duplicated variants")
    qc_group.add_argument("--normalize", action="store_true", help="Normalize indels")
    
    # Harmonization options
    harm_group = parser.add_argument_group("Harmonization Options")
    harm_group.add_argument("--basic-check", action="store_true", default=True, help="Run basic QC in harmonization (default: True)")
    harm_group.add_argument("--no-basic-check", dest="basic_check", action="store_false", help="Skip basic QC in harmonization")
    harm_group.add_argument("--ref-seq", help="Reference sequence file (FASTA) for allele flipping")
    harm_group.add_argument("--ref-rsid-tsv", help="Reference rsID HDF5 file (legacy name, accepts HDF5 path)")
    harm_group.add_argument("--ref-rsid-vcf", help="Reference rsID VCF/BCF file")
    harm_group.add_argument("--ref-infer", help="Reference VCF/BCF file for strand inference")
    harm_group.add_argument("--ref-alt-freq", help="Allele frequency field name in VCF INFO (default: AF)")
    harm_group.add_argument("--ref-maf-threshold", type=float, default=0.4, help="MAF threshold for reference (default: 0.4)")
    harm_group.add_argument("--maf-threshold", type=float, default=0.40, help="MAF threshold for sumstats (default: 0.40)")
    harm_group.add_argument("--sweep-mode", action="store_true", help="Use sweep mode for large datasets")
    
    # Assign rsID options
    assign_group = parser.add_argument_group("Assign rsID Options")
    assign_group.add_argument("--overwrite", choices=["all", "invalid", "empty"], default="empty", help="Overwrite mode for rsID assignment (default: empty)")
    
    # rsID to CHR:POS options
    rtc_group = parser.add_argument_group("rsID to CHR:POS Options")
    rtc_group.add_argument("--build", default="19", help="Genome build (default: 19)")
    rtc_group.add_argument("--overwrite-rtc", action="store_true", help="Overwrite existing CHR:POS")
    rtc_group.add_argument("--chunksize", type=int, default=5000000, help="Chunk size for processing (default: 5000000)")
    
    # Output formatting options
    format_group = parser.add_argument_group("Output Formatting Options")
    format_group.add_argument("--no-gzip", action="store_true", help="Disable gzip compression")
    format_group.add_argument("--bgzip", action="store_true", help="Use bgzip compression")
    format_group.add_argument("--tabix", action="store_true", help="Create tabix index")
    format_group.add_argument("--hapmap3", action="store_true", help="Extract HapMap3 variants only")
    format_group.add_argument("--exclude-hla", action="store_true", help="Exclude HLA region")
    format_group.add_argument("--hla-lower", type=int, default=25, help="HLA region lower bound (default: 25)")
    format_group.add_argument("--hla-upper", type=int, default=34, help="HLA region upper bound (default: 34)")
    format_group.add_argument("--n", type=float, help="Add N column with specified value")
    format_group.add_argument("--chr-prefix", default="", help="Prefix for chromosome column (e.g., 'chr')")
    format_group.add_argument("--xymt-number", action="store_true", help="Use numeric notation for X, Y, MT (23, 24, 25)")
    
    # General options
    general_group = parser.add_argument_group("General Options")
    general_group.add_argument("--threads", type=int, default=1, help="Number of threads for parallel processing (default: 1)")
    general_group.add_argument("--quiet", action="store_true", help="Suppress output messages")
    
    args = parser.parse_args(argv)
    
    # If no processing is specified but output is requested, do basic check
    # This ensures output is properly formatted
    if not (args.qc or args.harmonize or args.assign_rsid or args.rsid_to_chrpos):
        if args.out is not None:
            args.qc = True
        else:
            # If no output and no processing, just load and return (no-op)
            return
    
    # Process sumstats
    _process_sumstats(args)

if __name__ == "__main__":
    main()
