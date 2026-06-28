"""
GWASLab CLI - Handler functions.

All helper functions used by the CLI live here so that package import remains
lightweight until actual commands require loading ``gwaslab``.
"""

from __future__ import annotations

import sys
import os
import json
from pathlib import Path
from typing import Optional


def _resolve_path_keyword(gl, keyword: str) -> Optional[str]:
    """Resolve a path keyword from built-in paths or downloaded records."""
    if keyword in gl.options.paths:
        return gl.options.paths[keyword]
    return gl.get_path(keyword, verbose=False)


def load_sumstats(path: str, fmt: str, nrows: Optional[int] = None, build: str = "19"):
    """Load summary statistics from file.

    Parameters
    ----------
    path : str
        Path to input file.
    fmt : str
        Input format.
    nrows : int, optional
        Number of rows to read (for testing).
    build : str, optional
        Genome build of the input file (default: "19"). Encoded into the
        STATUS column at load time; must match the actual build of the input
        data.

    Returns
    -------
    gl.Sumstats
        Loaded Sumstats object.
    """
    import gwaslab as gl

    kwargs: dict = dict(fmt=fmt, build=build)
    if nrows is not None:
        kwargs["nrows"] = nrows
    return gl.Sumstats(path, **kwargs)


def run_config(args) -> None:
    """Handle the ``config`` subcommand."""
    import json
    import gwaslab as gl

    if args.json:
        print(json.dumps({"paths": gl.options.paths}, indent=2))
    else:
        for key, value in gl.options.paths.items():
            print(f"{key}: {value}")


def run_config_show(args) -> None:
    """Handle the ``config show`` subcommand."""
    import gwaslab as gl

    path = _resolve_path_keyword(gl, args.keyword)
    if not path:
        print(f"Path not found: {args.keyword}", file=sys.stderr)
        sys.exit(1)

    # For built-in JSON-backed config keys, display file content directly.
    if args.keyword in {"config", "reference", "formatbook"}:
        if not os.path.exists(path):
            print(f"Config file not found: {path}", file=sys.stderr)
            sys.exit(1)
        try:
            with open(path, "r", encoding="utf-8") as f:
                payload = json.load(f)
            print(json.dumps(payload, indent=2))
            return
        except Exception as exc:
            print(f"Failed to load JSON from {path}: {exc}", file=sys.stderr)
            sys.exit(1)

    # For non-JSON keys (e.g. data_directory or downloaded keyword), print resolved path.
    print(path)


def run_path(args) -> None:
    """Handle the ``path`` subcommand."""
    import gwaslab as gl

    path = _resolve_path_keyword(gl, args.keyword)
    if path:
        print(path)
    else:
        print(f"Path not found: {args.keyword}", file=sys.stderr)
        sys.exit(1)


def run_fb_list(args) -> None:
    """Handle the ``formatbook list`` subcommand."""
    import json
    import gwaslab as gl

    rows = gl.list_formats_with_descriptions(silent=True)
    if args.json:
        payload = {
            "formats": [
                {"index": i, "name": name, "description": desc}
                for i, (name, desc) in enumerate(rows, start=1)
            ]
        }
        print(json.dumps(payload, indent=2))
    else:
        w_idx = len(str(len(rows))) if rows else 1
        w_name = max((len(name) for name, _ in rows), default=0)
        for i, (name, desc) in enumerate(rows, start=1):
            print(f"{i:>{w_idx}}  {name:<{w_name}}  {desc}")


def run_fb_show(args) -> None:
    """Handle the ``formatbook show`` subcommand."""
    import json
    import gwaslab as gl

    mapping = gl.check_format(args.format)
    print(json.dumps({args.format: mapping}, indent=2))


def run_fb_update(args) -> None:
    """Handle the ``formatbook update`` subcommand."""
    import gwaslab as gl

    gl.update_formatbook()


def run_version(args) -> None:
    """Handle the ``version`` subcommand."""
    import gwaslab as gl

    gl.show_version()


def run_download(args) -> None:
    """Handle ``download-sumstats`` and ``download sumstats`` subcommands."""
    import gwaslab as gl

    gl.download_sumstats(
        gcst_id=args.gcst_id,
        output_dir=args.output_dir,
        verbose=True,
    )


def run_download_ref(args) -> None:
    """Handle the ``download-ref`` and ``download ref`` subcommands."""
    import gwaslab as gl

    gl.download_ref(
        name=args.keyword,
        directory=args.directory,
        local_filename=args.local_filename,
        overwrite=args.overwrite,
    )


def run_list_ref(args) -> None:
    """List reference keywords: remote catalog (--available) and/or local registry (--downloaded)."""
    import json
    import gwaslab as gl

    show_a = args.available
    show_d = args.downloaded
    if not show_a and not show_d:
        show_a = True
        show_d = True

    verbose = not getattr(args, "quiet", False)
    payload: dict = {}
    if show_a:
        payload["available"] = gl.check_available_ref(verbose=verbose, show_all=False)
    if show_d:
        payload["downloaded"] = gl.check_downloaded_ref(verbose=verbose)
        from gwaslab.bd.bd_download import filter_downloaded_registry

        payload["downloaded"] = filter_downloaded_registry(
            payload["downloaded"],
            source=getattr(args, "source", None),
            kind=getattr(args, "kind", None),
        )

    if args.json:
        print(json.dumps(payload, indent=2))
        return

    if show_a:
        print("Available reference keywords (install with: gwaslab download ref KEY):")
        for key in sorted(payload.get("available", {}).keys()):
            meta = payload["available"].get(key) or {}
            line = f"  {key}"
            desc = meta.get("description")
            if desc:
                line += f" — {desc}"
            print(line)

    if show_d:
        if show_a:
            print()
        print("Downloaded / registered references (local paths):")
        dl = payload.get("downloaded") or {}
        for key in sorted(dl.keys()):
            meta = dl[key] if isinstance(dl[key], dict) else {}
            path = meta.get("local_path", "")
            extra = f"\t{path}" if path else ""
            print(f"  {key}{extra}")


def run_config_set(args) -> None:
    """Handle the ``config set`` subcommand."""
    import gwaslab as gl

    key = args.key
    value = os.path.abspath(os.path.expanduser(args.path))
    allowed = {"data_directory", "config"}
    if key not in allowed:
        print(f"Unsupported config key: {key}. Allowed: {', '.join(sorted(allowed))}", file=sys.stderr)
        sys.exit(1)

    if key == "data_directory":
        gl.set_default_directory(value, persist=True)
    else:
        gl.options.set_option(key, value, persist=True)
        if key == "config":
            parent = os.path.dirname(value)
            if parent:
                os.makedirs(parent, exist_ok=True)
            if not os.path.exists(value):
                with open(value, "w", encoding="utf-8") as handle:
                    json.dump({"downloaded": {}}, handle, indent=4)

    print(f"{key}: {gl.options.paths[key]}")


def run_init(args) -> None:
    """Prepare local reference registry: dirs, config migration, optional data_directory, scan."""
    import gwaslab as gl
    from gwaslab.bd.bd_config import ensure_user_layout

    ensure_user_layout()

    scan_dir = None
    if getattr(args, "scan_dir", None):
        scan_dir = os.path.abspath(os.path.expanduser(args.scan_dir))
        if not os.path.isdir(scan_dir):
            os.makedirs(scan_dir, exist_ok=True)
        gl.set_default_directory(scan_dir, persist=True)
    else:
        scan_dir = gl.get_default_directory()

    ok = gl.scan_downloaded_files(
        verbose=not args.quiet,
        directory=scan_dir,
        recursive=getattr(args, "recursive", False),
    )
    if not ok:
        sys.exit(1)


def run_ref_add(args) -> None:
    """Register a local reference file in the registry."""
    import gwaslab as gl

    local_path = os.path.abspath(os.path.expanduser(args.path))
    tbi = os.path.abspath(os.path.expanduser(args.tbi)) if getattr(args, "tbi", None) else None
    ok = gl.add_local_data(
        args.keyword,
        local_path,
        description=getattr(args, "description", None),
        tbi=tbi,
    )
    if not ok:
        sys.exit(1)
    resolved = gl.get_path(args.keyword, verbose=False)
    print(resolved if resolved else local_path)


def run_ref_remove(args) -> None:
    """Remove a registry entry; optionally delete files from disk."""
    import gwaslab as gl

    ok = gl.remove_local_record(
        args.keyword,
        delete_file=getattr(args, "delete_file", False),
    )
    if not ok:
        sys.exit(1)


def run_report(args) -> None:
    """Generate QC HTML/PDF report via ``Sumstats.report()``."""
    suffix = Path(args.output).suffix.lower()
    if suffix not in {".html", ".pdf"}:
        print(
            "report: --output must end with .html or .pdf",
            file=sys.stderr,
        )
        sys.exit(1)

    s = load_sumstats(
        args.input,
        args.fmt,
        nrows=getattr(args, "nrows", None),
        build=args.build,
    )
    if getattr(args, "filter_region", None):
        region_chr, region_start, region_end = args.filter_region
        try:
            region_start = int(region_start)
            region_end = int(region_end)
        except ValueError:
            print("report: --filter-region START and END must be integers", file=sys.stderr)
            sys.exit(1)
        if region_start > region_end:
            print("report: --filter-region requires START <= END", file=sys.stderr)
            sys.exit(1)
        # Ensure coordinate columns are normalized before region filtering.
        s.fix_chr(verbose=not args.quiet)
        s.fix_pos(verbose=not args.quiet)
        s.filter_region(
            inplace=True,
            region=(region_chr, region_start, region_end),
            verbose=not args.quiet,
        )

    basic_check_kwargs = {
        "remove": args.remove,
        "remove_dup": args.remove_dup,
        "normalize": not args.no_normalize,
        "threads": args.threads,
        "verbose": not args.quiet,
    }

    get_lead_kwargs = {
        "sig_level": args.sig_level,
        "windowsizekb": args.windowsizekb,
        "anno": args.lead_anno,
        "verbose": not args.quiet,
    }

    mqq_plot_kwargs = {
        "mode": "mqq",
        "sig_level": args.sig_level,
        "dpi": args.mqq_dpi,
        "verbose": not args.quiet,
    }

    regional_plot_kwargs: dict = {"verbose": not args.quiet}
    if getattr(args, "regional_vcf", None):
        regional_plot_kwargs["vcf_path"] = args.regional_vcf

    harmonize_kwargs = None
    if args.harmonize:
        ref_alt_freq_kw = None
        if args.ref_infer:
            ref_alt_freq_kw = args.ref_alt_freq or "AF"
        harmonize_kwargs = {
            "basic_check": False,
            "ref_seq": args.ref_seq,
            "ref_rsid_tsv": args.ref_rsid_tsv,
            "ref_rsid_vcf": args.ref_rsid_vcf,
            "ref_infer": args.ref_infer,
            "ref_alt_freq": ref_alt_freq_kw,
            "ref_maf_threshold": args.ref_maf_threshold,
            "maf_threshold": args.maf_threshold,
            "threads": args.threads,
            "sweep_mode": args.sweep_mode,
            "verbose": not args.quiet,
        }

    output_kwargs = None
    save_path = getattr(args, "save_sumstats", None)
    if save_path:
        output_kwargs = {
            "path": save_path,
            "fmt": args.save_fmt,
            "tab_fmt": args.save_tab_fmt,
            "gzip": not args.save_no_gzip,
            "verbose": not args.quiet,
        }

    out = s.report(
        output_path=args.output,
        basic_check_kwargs=basic_check_kwargs,
        harmonize_kwargs=harmonize_kwargs,
        get_lead_kwargs=get_lead_kwargs,
        mqq_plot_kwargs=mqq_plot_kwargs,
        regional_plot_kwargs=regional_plot_kwargs,
        output_kwargs=output_kwargs,
        report_title=args.report_title,
        verbose=not args.quiet,
    )
    if not args.quiet:
        print(out)

