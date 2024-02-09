from gwaslab.ldsc_sumstats import estimate_h2
from gwaslab.g_Log import Log

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

def _estimate_h2_by_ldsc(insumstats, log, **args):
    
    default_args = ARGS(**args)

    sumstats = insumstats.copy()

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})

    summary = estimate_h2(sumstats, default_args, log)
    

    return summary