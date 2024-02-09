from gwaslab.ldsc_sumstats import estimate_h2
from gwaslab.ldsc_sumstats import estimate_rg
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.qc_fix_sumstats import skipped
from gwaslab.io_read_ldsc import parse_ldsc_summary

class ARGS():
    def __init__(self, **args):
        
        self.out = "ldsc"
        self.bfile = None
        self.l2 = False
        self.extract = None
        self.keep = None
        self.ld_wind_snps = None
        self.ld_wind_kb = None
        self.ld_wind_cm = None
        self.print_snps = None
        self.annot =None
        self.thin_annot = False
        self.cts_bin = None
        self.cts_breaks = None
        self.cts_names = None
        self.per_allele = False
        self.pq_exp =None
        self.no_print_annot = False

        if "h2" in args.keys():
            self.h2 = args["h2"]
        else:
            self.h2 = None

        self.h2_cts = None

        if "rg" in args.keys():
            self.rg = args["rg"]
        else:
            self.rg = None

        if "ref_ld" in args.keys():
            self.ref_ld = args["ref_ld"]
        else:
            self.ref_ld = None

        if "ref_ld_chr" in args.keys():
            self.ref_ld_chr = args["ref_ld_chr"]
        else:
            self.ref_ld_chr = None

        if "w_ld" in args.keys():
            self.w_ld = args["w_ld"]
        else:
            self.w_ld = None
        
        if "w_ld_chr" in args.keys():
            self.w_ld_chr = args["w_ld_chr"]
        else:
            self.w_ld_chr = None     

        self.overlap_annot = False
        self.print_coefficients = "ldsc"
        self.frqfile = None
        self.frqfile_chr = None
        self.no_intercept = None
        self.intercept_h2 = None
        self.intercept_gencov = None
        self.M = None
        self.two_step = None
        self.chisq_max = None
        self.ref_ld_chr_cts = None
        self.print_cov = None
        self.print_delete_vals = False
        self.chunk_size = 50
        self.pickle = False
        self.yes_really = False
        self.invert_anyway = False
        self.n_blocks = 200
        self.not_M_5_50 = False
        self.no_check_alleles = False
        self.return_silly_things = False
        self.samp_prev =None
        self.pop_prev =None

def _estimate_h2_by_ldsc(insumstats, log, verbose=True, **args):
    sumstats = insumstats.copy()
    ##start function with col checking##########################################################
    _start_line = "run LD score regression"
    _end_line = "running LD score regression"
    _start_cols =[]
    _start_function = ".estimate_h2_by_ldsc()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None
    ############################################################################################
    log.write(" -Run single variate LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)
    
    for key, value in args.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**args)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    log.write(" -LDSC log:", verbose=verbose)
    summary = estimate_h2(sumstats, default_args, log)
    
    log.write(" -Results have been stored in .ldsc_h2", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return parse_ldsc_summary(summary)

def _estimate_rg_by_ldsc(insumstats, other_traits ,log, verbose=True, **args):
    sumstats = insumstats.copy()
    ##start function with col checking##########################################################
    _start_line = "run LD score regression for genetic correlation"
    _end_line = "running LD score regression for genetic correlation"
    _start_cols =[]
    _start_function = ".estimate_rg_by_ldsc()"
    _must_args ={}

    is_enough_info = start_to(sumstats=sumstats,
                            log=log,
                            verbose=verbose,
                            start_line=_start_line,
                            end_line=_end_line,
                            start_cols=_start_cols,
                            start_function=_start_function,
                            **_must_args)
    if is_enough_info == False: return None
    ############################################################################################
    log.write(" -Run cross-trait LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)
    
    for key, value in args.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**args)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    other_traits_to_use = []
    
    for each_other_sumstats in other_traits:
        log.write("Processing {}".format(each_other_sumstats.meta["gwaslab"]["study_name"]))
        if "rsID" not in each_other_sumstats.data.columns:
            to_append = each_other_sumstats.filter_hapmap3(verbose=False).data.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
        else:        
            to_append = each_other_sumstats.data.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
        
        if "Z" not in to_append.columns:
            to_append["Z"] = to_append["BETA"]/to_append["SE"]

        other_traits_to_use.append(to_append[["SNP","A1","A2","Z","N"]])    

    log.write(" -LDSC log:", verbose=verbose)
    summary = estimate_rg(sumstats[["SNP","A1","A2","Z","N"]], other_traits_to_use, default_args, log)[1]
    
    log.write(" -Results have been stored in .ldsc_rg", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return summary