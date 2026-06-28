"""
GWASLab CLI - Main entry point.

Simple flat CLI structure with flags for operations.
"""

import argparse
import os
import sys
from typing import List, Optional

import numpy as np
import pandas as pd

from gwaslab.util.util_in_filter_value import _exclude, _extract

from gwaslab_cli.cli_banner import emit_cli_mode_banner
from gwaslab_cli.pair import run_pair
from gwaslab_cli.handlers import (
    load_sumstats,
    run_config,
    run_config_show,
    run_config_set,
    run_path,
    run_fb_list,
    run_fb_show,
    run_fb_update,
    run_version,
    run_download,
    run_download_ref,
    run_list_ref,
    run_init,
    run_ref_add,
    run_ref_remove,
    run_report,
)


def _apply_cli_fix_steps(s, args) -> bool:
    """Run optional Sumstats fix_* steps from CLI flags.

    Returns True if both fix_chr and fix_pos were applied (so coordinates match
    ``ensure_coords_ready()``); False otherwise.
    """
    verbose = not args.quiet
    need_chr = args.fix_chr or args.fix_chr_pos or args.fix_chr_pos_allele
    need_pos = args.fix_pos or args.fix_chr_pos or args.fix_chr_pos_allele
    need_allele = args.fix_allele or args.fix_chr_pos_allele
    if not (need_chr or need_pos or need_allele or args.fix_id):
        return False
    if need_chr:
        s.fix_chr(verbose=verbose)
    if need_pos:
        s.fix_pos(verbose=verbose)
    if need_allele:
        s.fix_allele(verbose=verbose)
    if args.fix_id:
        s.fix_id(verbose=verbose)
    return bool(need_chr and need_pos)


def _read_variant_id_list(path: str) -> List[str]:
    """Read one variant ID per line (first column if whitespace-separated); skip blanks and # comments."""
    out: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.append(line.split()[0])
    return out


def _variant_id_column(df: pd.DataFrame) -> Optional[str]:
    if "SNPID" in df.columns:
        return "SNPID"
    if "rsID" in df.columns:
        return "rsID"
    return None


def _strip_chr_token(u: str) -> str:
    u = str(u).strip().upper()
    return u[3:] if u.startswith("CHR") else u


def _apply_chr_filter(s, chr_filter: List[str], verbose: bool) -> None:
    allowed = {_strip_chr_token(t) for t in chr_filter}
    ch = s.data["CHR"]
    norm = ch.map(lambda x: _strip_chr_token(x) if pd.notna(x) else None)
    mask = norm.isin(allowed)
    before = len(s.data)
    s.data = s.data.loc[mask].copy()
    s.log.write(
        f" -Chromosome filter: kept {len(s.data)} of {before} variants", verbose=verbose
    )


def _maf_series(df: pd.DataFrame) -> Optional[pd.Series]:
    if "MAF" in df.columns:
        return pd.to_numeric(df["MAF"], errors="coerce")
    if "EAF" in df.columns:
        e = pd.to_numeric(df["EAF"], errors="coerce")
        return pd.Series(np.minimum(e, 1.0 - e), index=df.index)
    if "FRQ" in df.columns:
        e = pd.to_numeric(df["FRQ"], errors="coerce")
        return pd.Series(np.minimum(e, 1.0 - e), index=df.index)
    return None


def _apply_cli_variant_filters(s, args, parser: argparse.ArgumentParser) -> None:
    """PLINK-style variant filters (list, BED, chr, MAF/MAC, SNPs-only, INFO)."""
    verbose = not args.quiet
    vpath = getattr(args, "exclude", None)
    if vpath:
        ids = _read_variant_id_list(vpath)
        idcol = _variant_id_column(s.data)
        if idcol is None:
            parser.error("--exclude requires SNPID or rsID column in sumstats")
        s.data = _exclude(s.data, ids, id_use=idcol, log=s.log, verbose=verbose)
    if args.extract_list_file:
        ids = _read_variant_id_list(args.extract_list_file)
        idcol = _variant_id_column(s.data)
        if idcol is None:
            parser.error("--extract list requires SNPID or rsID column in sumstats")
        s.data = _extract(s.data, ids, id_use=idcol, log=s.log, verbose=verbose)
    if getattr(args, "extract_bed", None):
        p = os.path.expanduser(args.extract_bed)
        if not os.path.isfile(p):
            parser.error(f"--extract-bed: file not found: {args.extract_bed}")
        s.filter_region_in(inplace=True, path=p, build=args.build, verbose=verbose)
    if getattr(args, "exclude_bed", None):
        p = os.path.expanduser(args.exclude_bed)
        if not os.path.isfile(p):
            parser.error(f"--exclude-bed: file not found: {args.exclude_bed}")
        s.filter_region_out(inplace=True, path=p, build=args.build, verbose=verbose)
    cf = getattr(args, "chr_filter", None)
    if cf:
        _apply_chr_filter(s, cf, verbose)
    maf_s = None
    if args.maf is not None or args.max_maf is not None or args.mac is not None:
        maf_s = _maf_series(s.data)
        if maf_s is None:
            parser.error("MAF filters require MAF, EAF, or FRQ column")
    if args.maf is not None:
        before = len(s.data)
        s.data = s.data.loc[maf_s >= float(args.maf)].copy()
        s.log.write(
            f" -MAF >= {args.maf}: {len(s.data)} of {before} variants", verbose=verbose
        )
        maf_s = _maf_series(s.data)
    if args.max_maf is not None:
        before = len(s.data)
        s.data = s.data.loc[maf_s <= float(args.max_maf)].copy()
        s.log.write(
            f" -MAF <= {args.max_maf}: {len(s.data)} of {before} variants",
            verbose=verbose,
        )
        maf_s = _maf_series(s.data)
    if args.mac is not None:
        if "MAC" in s.data.columns:
            mac_s = pd.to_numeric(s.data["MAC"], errors="coerce")
        elif "N" in s.data.columns and maf_s is not None:
            n = pd.to_numeric(s.data["N"], errors="coerce")
            mac_s = 2.0 * n * maf_s
        else:
            parser.error("--mac requires a MAC column or N plus MAF/EAF/FRQ")
        before = len(s.data)
        s.data = s.data.loc[mac_s >= float(args.mac)].copy()
        s.log.write(
            f" -MAC >= {args.mac}: {len(s.data)} of {before} variants", verbose=verbose
        )
    if getattr(args, "snps_only", False):
        s.filter_snp(inplace=True, mode="in", verbose=verbose)
    if args.min_info is not None:
        info_col = None
        for c in s.data.columns:
            if str(c).upper() == "INFO":
                info_col = c
                break
        if info_col is None:
            parser.error("--min-info requires an INFO column")
        before = len(s.data)
        vals = pd.to_numeric(s.data[info_col], errors="coerce")
        s.data = s.data.loc[vals >= float(args.min_info)].copy()
        s.log.write(
            f" -INFO >= {args.min_info}: {len(s.data)} of {before} variants",
            verbose=verbose,
        )


def _add_download_ref_args(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "keyword",
        help=(
            "Reference keyword (see: gwaslab list ref --available; gwaslab config show reference)"
        ),
    )
    p.add_argument(
        "--directory",
        "-d",
        metavar="DIR",
        help="Directory for downloaded file (default: GWASLab data directory)",
    )
    p.add_argument(
        "--local-filename",
        metavar="NAME",
        help="Optional local filename (default: from download URL)",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing local file if present",
    )


def _add_download_sumstats_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("gcst_id", help="GWAS Catalog accession (e.g. GCST90270926)")
    p.add_argument(
        "-o",
        "--output-dir",
        "--directory",
        "-d",
        dest="output_dir",
        metavar="DIR",
        help="Output directory for downloaded files",
    )


def main(argv: Optional[list] = None) -> None:
    """Main entry point for GWASLab CLI."""
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="gwaslab",
        description=(
            "GWASLab: A Python package for processing GWAS summary statistics\n"
            "\n"
            "Processing order (when multiple flags are given):\n"
            "  1. Optional fixes --fix-chr / --fix-pos / --fix-chr-pos / --fix-chr-pos-allele / --fix-allele / --fix-id\n"
            "  2. QC            --qc / --remove / --remove-dup / --normalize\n"
            "  3. Filter region --filter-region CHR START END (auto fix_chr+fix_pos if QC not run)\n"
            "  4. Variant filters --extract/--exclude lists, --extract-bed/--exclude-bed, --chr, --maf/--max-maf/--mac, ...\n"
            "  5. Harmonize     --harmonize [--ref-seq ...]\n"
            "  6. Assign rsID   --assign-rsid  (auto fix_chr+fix_pos if QC not run)\n"
            "  7. rsID→CHR:POS  --rsid-to-chrpos\n"
            "  8. Infer build   --infer-build\n"
            "  9. Liftover      --liftover FROM TO  (auto fix_chr+fix_pos if QC not run)\n"
            " 10. Plot          --plot TYPE         (auto fix_chr+fix_pos if QC not run; writes to --output, then exits)\n"
            " 11. Get            --get lead|novel|proxy (writes to --output, then exits)\n"
            " 12. Save          --output FILE [--to-fmt FORMAT]"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process sumstats
  gwaslab --input sumstats.tsv --qc --output cleaned.tsv
  gwaslab --input sumstats.tsv --harmonize --ref-seq ref.fa --output harmonized.tsv
  gwaslab --input sumstats.tsv --liftover 19 38 --output lifted_hg38.tsv
  gwaslab --input sumstats.tsv --filter-region 7 126253550 128253550 --output chr7_region.tsv.gz
  gwaslab --input sumstats.tsv --infer-build --liftover 19 38 --output lifted_hg38.tsv
  gwaslab --input sumstats.tsv --output converted.ldsc --to-fmt ldsc

  # Plot
  gwaslab --input sumstats.tsv --plot manhattan --output manhattan.png
  gwaslab --input sumstats.tsv --plot qq --output qq.png

  # Extract
  gwaslab --input sumstats.tsv --get lead --output leads.tsv
  gwaslab --input sumstats.tsv --get novel --efo EFO_0004330 --output novel.tsv

  # Other commands
  gwaslab version
  gwaslab config
  gwaslab config show config
  gwaslab config show reference
  gwaslab path config
  gwaslab formatbook list
  gwaslab list ref --available
  gwaslab download ref 1kg_eas_hg19
  gwaslab download sumstats GCST90270926 --directory ./downloads
  gwaslab download-sumstats GCST90270926
  gwaslab download-ref 1kg_eas_hg19
  gwaslab init
  gwaslab init --directory /path/to/reference/cache

  # QC report (HTML or PDF)
  gwaslab report --input sumstats.tsv --output qc_report.html
  gwaslab report -i sumstats.tsv -o report.pdf --vcf ref.vcf.gz --quiet
        """,
    )

    # Input/Output
    parser.add_argument("--input", "-i", help="Input sumstats file path")
    parser.add_argument("--out", "--output", "-o", dest="output", help="Output file path")
    parser.add_argument("--fmt", "-f", default="auto", help="Input format (default: auto)")
    parser.add_argument("--to-fmt", default="gwaslab", help="Output format (default: gwaslab)")
    parser.add_argument("--nrows", type=int, help="Number of rows to read (for testing)")

    # Processing flags
    parser.add_argument("--qc", action="store_true", help="Perform quality control")
    parser.add_argument("--harmonize", action="store_true", help="Perform harmonization")
    parser.add_argument(
        "--fix-chr",
        action="store_true",
        help="Run fix_chr() only (chromosome notation)",
    )
    parser.add_argument(
        "--fix-pos",
        action="store_true",
        help="Run fix_pos() only (position dtype/range)",
    )
    parser.add_argument(
        "--fix-chr-pos",
        action="store_true",
        help="Run fix_chr() then fix_pos() (same as --fix-chr --fix-pos)",
    )
    parser.add_argument(
        "--fix-chr-pos-allele",
        action="store_true",
        help="Run fix_chr(), fix_pos(), and fix_allele()",
    )
    parser.add_argument(
        "--fix-allele",
        action="store_true",
        help="Run fix_allele() only (allele notation)",
    )
    parser.add_argument(
        "--fix-id",
        action="store_true",
        help="Run fix_id() only (SNPID / rsID column)",
    )
    parser.add_argument(
        "--assign-rsid",
        action="store_true",
        help="Assign rsID to variants (runs fix_chr+fix_pos if basic_check was not run)",
    )
    parser.add_argument("--rsid-to-chrpos", action="store_true", help="Convert rsID to CHR:POS")
    parser.add_argument(
        "--infer-build",
        action="store_true",
        help="Infer genome build from HapMap3 SNP coordinates (hg19/hg38)",
    )
    parser.add_argument(
        "--liftover",
        nargs=2,
        metavar=("FROM_BUILD", "TO_BUILD"),
        help="Liftover from one build to another (runs fix_chr+fix_pos if basic_check was not run)",
    )
    parser.add_argument(
        "--filter-region",
        nargs=3,
        metavar=("CHR", "START", "END"),
        help="Filter variants to a genomic region (runs fix_chr+fix_pos if basic_check was not run)",
    )
    parser.add_argument(
        "--extract",
        metavar="FILE",
        default=None,
        help="Variant ID list file (one ID per line); keep only these variants in the pipeline",
    )
    parser.add_argument(
        "--get",
        choices=["lead", "novel", "proxy"],
        dest="get_mode",
        default=None,
        help="Write lead variants, novel associations, or proxies to --output, then exit",
    )
    parser.add_argument(
        "--exclude",
        metavar="FILE",
        default=None,
        help="Exclude variant IDs from this file (one per line; matches SNPID or rsID column)",
    )
    parser.add_argument(
        "--extract-bed",
        metavar="FILE",
        default=None,
        help="Keep variants overlapping intervals in a BED file (0-based, half-open)",
    )
    parser.add_argument(
        "--exclude-bed",
        metavar="FILE",
        default=None,
        help="Remove variants overlapping intervals in a BED file",
    )
    parser.add_argument(
        "--chr",
        nargs="+",
        metavar="CHR",
        dest="chr_filter",
        default=None,
        help="Keep variants on these chromosomes only (PLINK 2 --chr analog)",
    )
    parser.add_argument("--maf", type=float, default=None, help="Minimum MAF (uses MAF, EAF, or FRQ column)")
    parser.add_argument("--max-maf", type=float, default=None, dest="max_maf", help="Maximum MAF")
    parser.add_argument(
        "--mac",
        type=float,
        default=None,
        help="Minimum MAC (uses MAC column, or 2*N*MAF if N and MAF/EAF/FRQ present)",
    )
    parser.add_argument(
        "--snps-only",
        action="store_true",
        help="Keep only SNPs (single-nucleotide EA and NEA)",
    )
    parser.add_argument(
        "--min-info",
        type=float,
        default=None,
        dest="min_info",
        help="Minimum INFO score (requires INFO column)",
    )

    # Plot flag
    parser.add_argument(
        "--plot",
        choices=["manhattan", "qq", "mqq", "regional", "miami"],
        help="Generate a plot (runs fix_chr+fix_pos if basic_check was not run)",
    )

    # Common options
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress output")

    # QC options
    parser.add_argument("--remove", action="store_true", help="Remove bad quality variants")
    parser.add_argument("--remove-dup", action="store_true", help="Remove duplicated variants")
    parser.add_argument("--normalize", action="store_true", help="Normalize indels")

    # Harmonization options
    parser.add_argument("--ref-seq", help="Reference FASTA file")
    parser.add_argument("--ref-rsid-tsv", help="Reference rsID HDF5 file")
    parser.add_argument("--ref-rsid-vcf", help="Reference rsID VCF file")
    parser.add_argument("--ref-infer", help="Reference VCF for strand inference")
    parser.add_argument(
        "--ref-alt-freq",
        default=None,
        help="VCF INFO field for ALT allele frequency when using --ref-infer (default: AF)",
    )
    parser.add_argument("--ref-maf-threshold", type=float, default=0.4)
    parser.add_argument("--maf-threshold", type=float, default=0.40)
    parser.add_argument("--sweep-mode", action="store_true")
    parser.add_argument("--basic-check", action="store_true", help="Enable basic check during harmonization")
    parser.add_argument(
        "--no-basic-check", action="store_true", help="Disable basic check during harmonization"
    )

    # Output options
    parser.add_argument("--no-gzip", action="store_true")
    parser.add_argument("--bgzip", action="store_true")
    parser.add_argument("--tabix", action="store_true")
    parser.add_argument("--hapmap3", action="store_true")
    parser.add_argument("--exclude-hla", action="store_true")
    parser.add_argument("--hla-lower", type=int, help="Lower bound for HLA exclusion (in Mb)")
    parser.add_argument("--hla-upper", type=int, help="Upper bound for HLA exclusion (in Mb)")
    parser.add_argument("--build", default="19")
    parser.add_argument("--chr-prefix", default="")
    parser.add_argument("--xymt-number", action="store_true", help="Use numeric encoding for XYMT chromosomes")
    parser.add_argument("--n", type=int, help="Sample size to add to output")
    parser.add_argument("--tab-fmt", default="tsv", help="Output table format (tsv, csv, parquet)")
    parser.add_argument("--overwrite", default="empty", help="Overwrite mode for assign-rsid (empty, all, etc.)")

    # Plot options
    parser.add_argument("--sig-level", type=float, default=5e-8)
    parser.add_argument("--ylim", type=float, nargs=2)
    parser.add_argument("--highlight", nargs="+")
    parser.add_argument(
        "--plot-chr",
        metavar="CHR",
        default=None,
        help="Chromosome for --plot regional (or pass a single --chr with --start/--end)",
    )
    parser.add_argument("--start", type=int, help="Start position for regional plot")
    parser.add_argument("--end", type=int, help="End position for regional plot")

    # Extract options
    parser.add_argument("--sig-level-extract", type=float, default=5e-8, dest="sig_level_extract")
    parser.add_argument("--windowsizekb", type=int, default=500)
    parser.add_argument("--efo", nargs="+", help="EFO IDs for novel extraction")
    parser.add_argument("--only-novel", action="store_true")

    # Subcommands (version, config, etc.)
    subparsers = parser.add_subparsers(dest="command")

    # Version
    version_parser = subparsers.add_parser("version", help="Show version")
    version_parser.set_defaults(func=run_version)

    # Config
    config_parser = subparsers.add_parser("config", help="Show configuration")
    config_parser.add_argument("--json", action="store_true")
    config_parser.set_defaults(func=run_config)
    config_sub = config_parser.add_subparsers(dest="config_action")
    config_show = config_sub.add_parser("show", help="Show one configured path by keyword")
    config_show.add_argument("keyword", help="Path keyword (e.g., config, reference)")
    config_show.set_defaults(func=run_config_show)
    config_set = config_sub.add_parser("set", help="Set a configuration path (persisted to settings.json)")
    config_set.add_argument("key", choices=["data_directory", "config"], help="Configuration key")
    config_set.add_argument("path", help="New path value")
    config_set.set_defaults(func=run_config_set)

    # Path
    path_parser = subparsers.add_parser("path", help="Resolve path")
    path_parser.add_argument("keyword")
    path_parser.set_defaults(func=run_path)

    # Formatbook
    fb_parser = subparsers.add_parser("formatbook", help="Manage formats")
    fb_sub = fb_parser.add_subparsers(dest="fb_action")
    fb_list = fb_sub.add_parser("list", help="List formats")
    fb_list.add_argument("--json", action="store_true")
    fb_list.set_defaults(func=run_fb_list)
    fb_show = fb_sub.add_parser("show", help="Show format")
    fb_show.add_argument("format")
    fb_show.set_defaults(func=run_fb_show)
    fb_update = fb_sub.add_parser("update", help="Update formatbook")
    fb_update.set_defaults(func=run_fb_update)

    # download ref | sumstats (unified entry; same flags as legacy flat commands)
    download_group = subparsers.add_parser(
        "download",
        help="Download reference data or GWAS Catalog sumstats",
    )
    download_sub = download_group.add_subparsers(dest="download_action", required=True)
    dl_ref_cmd = download_sub.add_parser(
        "ref",
        help="Download a reference file by keyword (same as download-ref)",
    )
    _add_download_ref_args(dl_ref_cmd)
    dl_ref_cmd.set_defaults(func=run_download_ref)
    dl_sum_cmd = download_sub.add_parser(
        "sumstats",
        help="Download sumstats from GWAS Catalog (same as download-sumstats)",
    )
    _add_download_sumstats_args(dl_sum_cmd)
    dl_sum_cmd.set_defaults(func=run_download)

    # Download sumstats (legacy flat command)
    dl_parser = subparsers.add_parser(
        "download-sumstats",
        help="Download from GWAS Catalog (alias of: gwaslab download sumstats)",
    )
    _add_download_sumstats_args(dl_parser)
    dl_parser.set_defaults(func=run_download)

    # Download reference (legacy flat command)
    dl_ref_parser = subparsers.add_parser(
        "download-ref",
        help="Download a reference by keyword (alias of: gwaslab download ref)",
    )
    _add_download_ref_args(dl_ref_parser)
    dl_ref_parser.set_defaults(func=run_download_ref)

    # List reference catalogs
    list_parser = subparsers.add_parser(
        "list",
        help="List items (e.g. reference keywords available or already registered)",
    )
    list_sub = list_parser.add_subparsers(dest="list_action", required=True)
    list_ref_parser = list_sub.add_parser(
        "ref",
        help="List reference keywords from the catalog and/or local download registry",
    )
    list_ref_parser.add_argument(
        "--available",
        action="store_true",
        help="Show keywords from the remote reference catalog (reference.json)",
    )
    list_ref_parser.add_argument(
        "--downloaded",
        action="store_true",
        help="Show keywords registered in local config (downloaded references)",
    )
    list_ref_parser.add_argument(
        "--source",
        choices=["catalog", "gwas_catalog", "local"],
        help="Filter --downloaded entries by source metadata",
    )
    list_ref_parser.add_argument(
        "--kind",
        choices=["ref", "sumstats"],
        help="Filter --downloaded entries by registry kind",
    )
    list_ref_parser.add_argument("--json", action="store_true", help="Print JSON")
    list_ref_parser.add_argument("--quiet", "-q", action="store_true", help="Less logging")
    list_ref_parser.set_defaults(func=run_list_ref)

    # Init — scan data directory for downloaded reference files
    init_parser = subparsers.add_parser(
        "init",
        help="Prepare local reference registry (create dirs, migrate config, scan files)",
    )
    init_parser.add_argument(
        "--dir",
        "-d",
        "--directory",
        dest="scan_dir",
        metavar="DIR",
        help="Directory to scan (default: GWASLab data directory, see gwaslab path data_directory)",
    )
    init_parser.add_argument(
        "--recursive",
        "-r",
        action="store_true",
        help="Scan subdirectories recursively when registering files",
    )
    init_parser.add_argument("--quiet", "-q", action="store_true", help="Less logging")
    init_parser.set_defaults(func=run_init)

    ref_parser = subparsers.add_parser(
        "ref",
        help="Manage local reference registry entries",
    )
    ref_sub = ref_parser.add_subparsers(dest="ref_action", required=True)
    ref_add = ref_sub.add_parser("add", help="Register a local file by keyword")
    ref_add.add_argument("keyword", help="Registry keyword (e.g. 1kg_eas_hg19)")
    ref_add.add_argument("path", help="Absolute or relative path to the local file")
    ref_add.add_argument("--tbi", help="Path to .tbi index (VCF)")
    ref_add.add_argument("--description", help="Human-readable description")
    ref_add.set_defaults(func=run_ref_add)
    ref_remove = ref_sub.add_parser("remove", help="Remove a registry entry")
    ref_remove.add_argument("keyword", help="Registry keyword to remove")
    ref_remove.add_argument(
        "--delete-file",
        action="store_true",
        help="Also delete the registered file (and .tbi if recorded)",
    )
    ref_remove.set_defaults(func=run_ref_remove)

    pair_parser = subparsers.add_parser(
        "pair",
        help="Two-sumstats workflows (operation token can appear among flags)",
    )
    pair_parser.add_argument(
        "pair_argv",
        nargs=argparse.REMAINDER,
        help=(
            "Pair command arguments. Example: --input1 a.tsv --input2 b.tsv "
            "--fmt1 auto --fmt2 auto --build 19 --sync-alleles miami --output out.png"
        ),
    )
    pair_parser.set_defaults(func=lambda a: run_pair(a.pair_argv))

    report_parser = subparsers.add_parser(
        "report",
        help="Generate QC HTML/PDF report (QC, lead variants, MQQ, optional regional plots)",
    )
    report_parser.add_argument("--input", "-i", required=True, help="Input sumstats file")
    report_parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Report file path (.html or .pdf; PDF needs weasyprint)",
    )
    report_parser.add_argument("--fmt", "-f", default="auto", help="Input format (default: auto)")
    report_parser.add_argument("--build", default="19", help="Genome build of input (default: 19)")
    report_parser.add_argument(
        "--filter-region",
        nargs=3,
        metavar=("CHR", "START", "END"),
        help="Filter input to a genomic region before report generation",
    )
    report_parser.add_argument("--nrows", type=int, help="Max rows to read (for testing)")
    report_parser.add_argument("--threads", "-t", type=int, default=1, help="Thread count")
    report_parser.add_argument("--quiet", "-q", action="store_true", help="Less logging")
    report_parser.add_argument(
        "--title",
        default="GWAS Quality Control Report",
        dest="report_title",
        help="Report title",
    )
    report_parser.add_argument(
        "--remove",
        action="store_true",
        help="Remove low-quality variants in basic_check (same as main --remove)",
    )
    report_parser.add_argument("--remove-dup", action="store_true", help="Remove duplicates in basic_check")
    report_parser.add_argument(
        "--no-normalize",
        action="store_true",
        help="Disable indel normalization in basic_check",
    )
    report_parser.add_argument(
        "--sig-level",
        type=float,
        default=5e-8,
        help="Genome-wide significance for get_lead and MQQ (default: 5e-8)",
    )
    report_parser.add_argument(
        "--windowsizekb",
        type=int,
        default=500,
        help="Lead clumping window (kb) and regional window basis (default: 500)",
    )
    report_parser.add_argument(
        "--lead-anno",
        action="store_true",
        help="Annotate lead variants (slower)",
    )
    report_parser.add_argument("--mqq-dpi", type=int, default=300, help="MQQ figure DPI (default: 300)")
    report_parser.add_argument(
        "--vcf",
        dest="regional_vcf",
        metavar="PATH",
        help="Reference VCF (tabix-indexed) for LD in regional plots",
    )
    report_parser.add_argument(
        "--harmonize",
        action="store_true",
        help="Run harmonization after basic_check (supply reference options)",
    )
    report_parser.add_argument("--ref-seq", help="Reference FASTA for harmonization")
    report_parser.add_argument("--ref-rsid-tsv", help="Reference rsID HDF5 for harmonization")
    report_parser.add_argument("--ref-rsid-vcf", help="Reference rsID VCF for harmonization")
    report_parser.add_argument("--ref-infer", help="Reference VCF for strand inference (harmonization)")
    report_parser.add_argument(
        "--ref-alt-freq",
        default=None,
        help="VCF INFO field for ALT allele frequency with --ref-infer (default: AF)",
    )
    report_parser.add_argument("--ref-maf-threshold", type=float, default=0.4)
    report_parser.add_argument("--maf-threshold", type=float, default=0.40)
    report_parser.add_argument("--sweep-mode", action="store_true", help="Harmonization sweep mode")
    report_parser.add_argument(
        "--save-sumstats",
        metavar="PATH",
        dest="save_sumstats",
        help="Also write processed sumstats with to_format() (base path)",
    )
    report_parser.add_argument("--save-fmt", default="gwaslab", help="Format for --save-sumstats (default: gwaslab)")
    report_parser.add_argument(
        "--save-tab-fmt",
        default="tsv",
        help="Table format for --save-sumstats (default: tsv)",
    )
    report_parser.add_argument(
        "--save-no-gzip",
        action="store_true",
        help="Disable gzip for --save-sumstats",
    )
    report_parser.set_defaults(func=run_report)

    args = parser.parse_args(argv)
    emit_cli_mode_banner(parser.prog, argv, quiet=getattr(args, "quiet", False))

    # Handle subcommands
    if args.command:
        if hasattr(args, "func"):
            args.func(args)
        return

    # Require --input for main operations
    if not args.input:
        parser.error("--input is required (or use a subcommand like 'version', 'config')")

    extract_list_file: Optional[str] = None
    if args.extract is not None:
        pe = os.path.expanduser(args.extract)
        if not os.path.isfile(pe):
            parser.error(f"--extract: file not found: {args.extract!r}")
        extract_list_file = pe
    args.extract_list_file = extract_list_file
    args.extract_mode = args.get_mode

    if args.exclude:
        pe = os.path.expanduser(args.exclude)
        if not os.path.isfile(pe):
            parser.error(f"--exclude: file not found: {args.exclude!r}")
        args.exclude = pe

    if (
        getattr(args, "plot", None) == "regional"
        and args.plot_chr is None
        and args.chr_filter
        and len(args.chr_filter) == 1
        and args.start is not None
        and args.end is not None
    ):
        args.plot_chr = args.chr_filter[0]
        args.chr_filter = None

    # Derive input build from --build or --liftover FROM build.
    input_build = args.build
    if args.liftover and args.build == "19" and not args.infer_build:
        input_build = args.liftover[0]

    # Load sumstats
    s = load_sumstats(args.input, args.fmt, args.nrows, build=input_build)

    if args.infer_build:
        s.infer_build(verbose=not args.quiet)

    ran_basic_check = False
    coords_ready = False

    def ensure_coords_ready() -> None:
        """Ensure CHR/POS are standardized for operations that require them."""
        nonlocal coords_ready
        if not coords_ready:
            s.fix_chr(verbose=not args.quiet)
            s.fix_pos(verbose=not args.quiet)
            coords_ready = True

    if _apply_cli_fix_steps(s, args):
        coords_ready = True

    # Processing
    if args.qc or args.remove or args.remove_dup or args.normalize:
        normalize = args.normalize if args.normalize else (True if args.qc else False)
        s.basic_check(
            remove=args.remove,
            remove_dup=args.remove_dup,
            threads=args.threads,
            normalize=normalize,
            verbose=not args.quiet,
        )
        ran_basic_check = True
        coords_ready = True

    if args.filter_region:
        ensure_coords_ready()
        region_chr, region_start, region_end = args.filter_region
        try:
            region_start = int(region_start)
            region_end = int(region_end)
        except ValueError:
            parser.error("--filter-region START and END must be integers")
        if region_start > region_end:
            parser.error("--filter-region requires START <= END")
        s.filter_region(inplace=True, region=(region_chr, region_start, region_end), verbose=not args.quiet)

    need_coords_for_variant_filters = bool(
        args.extract_bed or args.exclude_bed or args.chr_filter
    )
    if need_coords_for_variant_filters:
        ensure_coords_ready()
    if (
        args.exclude
        or args.extract_list_file
        or args.extract_bed
        or args.exclude_bed
        or args.chr_filter
        or args.maf is not None
        or args.max_maf is not None
        or args.mac is not None
        or args.snps_only
        or args.min_info is not None
    ):
        _apply_cli_variant_filters(s, args, parser)

    if args.harmonize or args.ref_seq:
        if args.basic_check:
            basic_check = True
        elif args.no_basic_check:
            basic_check = False
        else:
            basic_check = not args.qc
        ref_alt_freq_kw = None
        if args.ref_infer:
            ref_alt_freq_kw = args.ref_alt_freq or "AF"
        s.harmonize(
            basic_check=basic_check,
            ref_seq=args.ref_seq,
            ref_rsid_tsv=args.ref_rsid_tsv,
            ref_rsid_vcf=args.ref_rsid_vcf,
            ref_infer=args.ref_infer,
            ref_alt_freq=ref_alt_freq_kw,
            ref_maf_threshold=args.ref_maf_threshold,
            maf_threshold=args.maf_threshold,
            threads=args.threads,
            sweep_mode=args.sweep_mode,
            verbose=not args.quiet,
        )
        ran_basic_check = ran_basic_check or basic_check
        coords_ready = coords_ready or basic_check

    if args.assign_rsid:
        ensure_coords_ready()
        s.assign_rsid(
            ref_rsid_tsv=args.ref_rsid_tsv,
            ref_rsid_vcf=args.ref_rsid_vcf,
            overwrite=args.overwrite,
            threads=args.threads,
            verbose=not args.quiet,
        )

    if args.rsid_to_chrpos:
        s.rsid_to_chrpos(
            ref_rsid_to_chrpos_vcf=args.ref_rsid_vcf,
            ref_rsid_to_chrpos_hdf5=args.ref_rsid_tsv,
            build=args.build,
            threads=args.threads if args.threads > 1 else 4,
            verbose=not args.quiet,
        )

    if args.liftover:
        ensure_coords_ready()
        from_build, to_build = args.liftover
        s.liftover(
            from_build=from_build,
            to_build=to_build,
            verbose=not args.quiet,
        )

    # Plotting
    if args.plot:
        if not ran_basic_check:
            # Plotting needs sanitized CHR/POS. If basic_check was not run,
            # apply the lightweight fallback of fix_chr + fix_pos.
            ensure_coords_ready()
        plot_common_kwargs = {
            "save": args.output,
            "verbose": not args.quiet,
            "sig_level": args.sig_level,
        }
        if args.ylim:
            plot_common_kwargs["ylim"] = tuple(args.ylim)
        if args.highlight:
            plot_common_kwargs["highlight"] = args.highlight
        if args.plot == "manhattan":
            s.plot_mqq(mode="m", **plot_common_kwargs)
        elif args.plot == "qq":
            s.plot_mqq(mode="qq", **plot_common_kwargs)
        elif args.plot == "mqq":
            s.plot_mqq(mode="mqq", **plot_common_kwargs)
        elif args.plot == "regional":
            if args.plot_chr is None or args.start is None or args.end is None:
                parser.error(
                    "--plot regional requires --plot-chr, --start, and --end "
                    "(or a single --chr CHROM together with --start and --end)"
                )
            reg_chr = args.plot_chr
            try:
                reg_chr = int(str(reg_chr))
            except (TypeError, ValueError):
                pass
            s.plot_mqq(
                mode="r",
                region=(reg_chr, args.start, args.end),
                **plot_common_kwargs,
            )
        elif args.plot == "miami":
            parser.error("Miami plot requires two inputs. Use Python API: gl.plot_miami2()")
        return

    # --get lead|novel|proxy (variant lists use --extract FILE in the filter step)
    if args.extract_mode:
        if args.extract_mode == "lead":
            result = s.get_lead(
                sig_level=args.sig_level_extract,
                windowsizekb=args.windowsizekb,
                verbose=not args.quiet,
            )
            if args.output:
                result.to_csv(args.output, sep="\t", index=False)
        elif args.extract_mode == "novel":
            result = s.get_novel(
                efo=args.efo or False,
                only_novel=args.only_novel,
                sig_level=args.sig_level_extract,
                verbose=not args.quiet,
            )
            if args.output:
                result.to_csv(args.output, sep="\t", index=False)
        elif args.extract_mode == "proxy":
            parser.error("Proxy extraction not yet implemented in CLI")
        return

    # Output (if not plot/extract which handle their own output)
    if args.output:
        to_format_kwargs = {
            "fmt": args.to_fmt,
            "tab_fmt": args.tab_fmt,
            "gzip": not args.no_gzip,
            "bgzip": args.bgzip,
            "tabix": args.tabix,
            "hapmap3": args.hapmap3,
            "exclude_hla": args.exclude_hla,
            "chr_prefix": args.chr_prefix,
            "xymt_number": args.xymt_number,
            "n": args.n,
            "verbose": not args.quiet,
        }

        if args.hla_lower is not None or args.hla_upper is not None:
            to_format_kwargs["hla_range"] = (args.hla_lower or 25, args.hla_upper or 34)

        s.to_format(args.output, **to_format_kwargs)


if __name__ == "__main__":
    main()

