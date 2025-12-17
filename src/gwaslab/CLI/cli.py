import argparse
import sys
import gwaslab as gl

def _load_sumstats(path, fmt, nrows):
    return gl.Sumstats(path, fmt=fmt, nrows=nrows) if nrows is not None else gl.Sumstats(path, fmt=fmt)

def _cmd_version(_args):
    gl.show_version()

def _cmd_check(args):
    s = _load_sumstats(args.input, args.fmt, args.nrows)
    s.basic_check(remove=args.remove, remove_dup=args.remove_dup, n_cores=args.n_cores, normalize=args.normalize, verbose=not args.quiet)
    if args.out is not None:
        s.to_format(args.out, fmt="gwaslab", verbose=not args.quiet)

def _cmd_assign_rsid(args):
    s = _load_sumstats(args.input, args.fmt, args.nrows)
    s.assign_rsid(ref_rsid_tsv=args.ref_tsv, ref_rsid_vcf=args.ref_vcf, n_cores=args.n_cores, overwrite=args.overwrite)
    if args.out is not None:
        s.to_format(args.out, fmt="gwaslab", verbose=not args.quiet)

def _cmd_rsid_to_chrpos(args):
    s = _load_sumstats(args.input, args.fmt, args.nrows)
    s.rsid_to_chrpos(ref_rsid_to_chrpos_tsv=args.ref_tsv, build=args.build, overwrite=args.overwrite, remove=args.remove, chunksize=args.chunksize, verbose=not args.quiet)
    if args.out is not None:
        s.to_format(args.out, fmt="gwaslab", verbose=not args.quiet)

def _cmd_to_format(args):
    s = _load_sumstats(args.input, args.fmt, args.nrows)
    s.to_format(args.out, fmt=args.out_fmt, tab_fmt=args.tab_fmt, gzip=not args.no_gzip, bgzip=args.bgzip, tabix=args.tabix, hapmap3=args.hapmap3, exclude_hla=args.exclude_hla, hla_range=(args.hla_lower, args.hla_upper), n=args.n, chr_prefix=args.chr_prefix, xymt_number=args.xymt_number, verbose=not args.quiet)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    commands = {"version", "check", "assign-rsid", "rsid-to-chrpos", "to-format"}
    has_command = any((a in commands) for a in argv if not a.startswith("-"))
    if not has_command:
        argv = ["check"] + list(argv)

    parser = argparse.ArgumentParser(prog="gwaslab", add_help=True)
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_version = subparsers.add_parser("version")
    p_version.set_defaults(func=_cmd_version)

    p_check = subparsers.add_parser("check")
    p_check.add_argument("--input", required=True)
    p_check.add_argument("--fmt", default="auto")
    p_check.add_argument("--nrows", type=int)
    p_check.add_argument("--n-cores", type=int, default=1)
    p_check.add_argument("--remove", action="store_true")
    p_check.add_argument("--remove-dup", action="store_true")
    p_check.add_argument("--normalize", action="store_true")
    p_check.add_argument("--out")
    p_check.add_argument("--quiet", action="store_true")
    p_check.set_defaults(func=_cmd_check)

    p_assign = subparsers.add_parser("assign-rsid")
    p_assign.add_argument("--input", required=True)
    p_assign.add_argument("--fmt", default="auto")
    p_assign.add_argument("--nrows", type=int)
    p_assign.add_argument("--ref-tsv")
    p_assign.add_argument("--ref-vcf")
    p_assign.add_argument("--overwrite", choices=["all", "invalid", "empty"], default="empty")
    p_assign.add_argument("--n-cores", type=int, default=1)
    p_assign.add_argument("--out")
    p_assign.add_argument("--quiet", action="store_true")
    p_assign.set_defaults(func=_cmd_assign_rsid)

    p_rtc = subparsers.add_parser("rsid-to-chrpos")
    p_rtc.add_argument("--input", required=True)
    p_rtc.add_argument("--fmt", default="auto")
    p_rtc.add_argument("--nrows", type=int)
    p_rtc.add_argument("--ref-tsv")
    p_rtc.add_argument("--build", default="19")
    p_rtc.add_argument("--overwrite", action="store_true")
    p_rtc.add_argument("--remove", action="store_true")
    p_rtc.add_argument("--chunksize", type=int, default=5000000)
    p_rtc.add_argument("--out")
    p_rtc.add_argument("--quiet", action="store_true")
    p_rtc.set_defaults(func=_cmd_rsid_to_chrpos)

    p_to = subparsers.add_parser("to-format")
    p_to.add_argument("--input", required=True)
    p_to.add_argument("--fmt", default="auto")
    p_to.add_argument("--nrows", type=int)
    p_to.add_argument("--out", required=True)
    p_to.add_argument("--out-fmt", default="gwaslab")
    p_to.add_argument("--tab-fmt", default="tsv")
    p_to.add_argument("--no-gzip", action="store_true")
    p_to.add_argument("--bgzip", action="store_true")
    p_to.add_argument("--tabix", action="store_true")
    p_to.add_argument("--hapmap3", action="store_true")
    p_to.add_argument("--exclude-hla", action="store_true")
    p_to.add_argument("--hla-lower", type=int, default=25)
    p_to.add_argument("--hla-upper", type=int, default=34)
    p_to.add_argument("--n", type=float)
    p_to.add_argument("--chr-prefix", default="")
    p_to.add_argument("--xymt-number", action="store_true")
    p_to.add_argument("--quiet", action="store_true")
    p_to.set_defaults(func=_cmd_to_format)

    args = parser.parse_args(argv)
    args.func(args)

if __name__ == "__main__":
    main()
