from gwaslab.ldsc_sumstats import estimate_h2
from gwaslab.ldsc_sumstats import estimate_rg
from gwaslab.ldsc_sumstats import cell_type_specific
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.qc_fix_sumstats import skipped
from gwaslab.io_read_ldsc import parse_ldsc_summary
from gwaslab.io_read_ldsc import parse_partitioned_ldsc_summary
class ARGS():
    def __init__(self, **args):
        
        self.out = "ldsc"

        if "bfile" in args.keys():
            self.bfile = args["bfile"]
        else:
            self.bfile = None 
        
        if "l2" in args.keys():
            self.l2 = args["l2"]
        else:
            self.l2 = None 

        if "extract" in args.keys():
            self.extract = args["extract"]
        else:
            self.extract = None 

        if "keep" in args.keys():
            self.keep = args["keep"]
        else:
            self.keep = None 

        if "ld_wind_snps" in args.keys():
            self.ld_wind_snps = args["ld_wind_snps"]
        else:
            self.ld_wind_snps = None 

        if "ld_wind_kb" in args.keys():
            self.ld_wind_kb = args["ld_wind_kb"]
        else:
            self.ld_wind_kb = None 

        if "ld_wind_cm" in args.keys():
            self.ld_wind_cm = args["ld_wind_cm"]
        else:
            self.ld_wind_cm = None 

        if "print_snps" in args.keys():
            self.print_snps = args["print_snps"]
        else:
            self.print_snps = None 

        if "annot" in args.keys():
            self.annot = args["annot"]
        else:
            self.annot = None 

        if "thin_annot" in args.keys():
            self.thin_annot = args["thin_annot"]
        else:
            self.thin_annot = None 

        if "cts_bin" in args.keys():
            self.cts_bin = args["cts_bin"]
        else:
            self.cts_bin = None 

        if "cts_breaks" in args.keys():
            self.cts_breaks = args["cts_breaks"]
        else:
            self.cts_breaks = None 

        if "cts_names" in args.keys():
            self.cts_names = args["cts_names"]
        else:
            self.cts_names = None 

        if "per_allele" in args.keys():
            self.per_allele = args["per_allele"]
        else:
            self.per_allele = None 

        if "pq_exp" in args.keys():
            self.pq_exp = args["pq_exp"]
        else:
            self.pq_exp = None 

        if "no_print_annot" in args.keys():
            self.no_print_annot = args["no_print_annot"]
        else:
            self.no_print_annot = None 

        if "h2" in args.keys():
            self.h2 = args["h2"]
        else:
            self.h2 = None

        if "h2_cts" in args.keys():
            self.h2_cts = args["h2_cts"]
        else:
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

        if "overlap_annot" in args.keys():
            self.overlap_annot = args["overlap_annot"]
        else:
            self.overlap_annot = None 

        if "print_coefficients" in args.keys():
            self.print_coefficients = args["print_coefficients"]
        else:
            self.print_coefficients = "ldsc" 

        if "frqfile" in args.keys():
            self.frqfile = args["frqfile"]
        else:
            self.frqfile = None 

        if "frqfile_chr" in args.keys():
            self.frqfile_chr = args["frqfile_chr"]
        else:
            self.frqfile_chr = None 

        if "no_intercept" in args.keys():
            self.no_intercept = args["no_intercept"]
        else:
            self.no_intercept = None 

        if "intercept_h2" in args.keys():
            self.intercept_h2 = args["intercept_h2"]
        else:
            self.intercept_h2 = None 

        if "intercept_gencov" in args.keys():
            self.intercept_gencov = args["intercept_gencov"]
        else:
            self.intercept_gencov = None 

        if "M" in args.keys():
            self.M = args["M"]
        else:
            self.M = None 

        if "two_step" in args.keys():
            self.two_step = args["two_step"]
        else:
            self.two_step = None 

        if "chisq_max" in args.keys():
            self.chisq_max = args["chisq_max"]
        else:
            self.chisq_max= None 

        if "ref_ld_chr_cts" in args.keys():
            self.ref_ld_chr_cts = args["ref_ld_chr_cts"]
        else:
            self.ref_ld_chr_cts = None

        if "print_all_cts" in args.keys():
            self.print_all_cts = args["print_all_cts"]
        else:
            self.print_all_cts = False 

        if "print_cov" in args.keys():
            self.print_cov = args["print_cov"]
        else:
            self.print_cov = None 

        self.print_delete_vals = False
        if "print_delete_vals" in args.keys():
            self.print_delete_vals = args["print_delete_vals"]
        else:
            self.print_delete_vals = False 

        if "chunk_size" in args.keys():
            self.chunk_size = args["chunk_size"]
        else:
            self.chunk_size = 50 

        if "pickle" in args.keys():
            self.pickle = args["pickle"]
        else:
            self.pickle = False 

        if "yes_really" in args.keys():
            self.yes_really = args["yes_really"]
        else:
            self.yes_really = False 

        if "invert_anyway" in args.keys():
            self.invert_anyway = args["invert_anyway"]
        else:
            self.invert_anyway = False 

        if "n_blocks" in args.keys():
            self.n_blocks = args["n_blocks"]
        else:
            self.n_blocks = 200 

        if "not_M_5_50" in args.keys():
            self.not_M_5_50 = args["not_M_5_50"]
        else:
            self.not_M_5_50 = False 

        if "no_check_alleles" in args.keys():
            self.no_check_alleles = args["no_check_alleles"]
        else:
            self.no_check_alleles = False 
        
        if "return_silly_things" in args.keys():
            self.return_silly_things = args["return_silly_things"]
        else:
            self.return_silly_things = False 

        if "samp_prev" in args.keys():
            self.samp_prev = args["samp_prev"]
        else:
            self.samp_prev = None 
        
        if "pop_prev" in args.keys():
            self.pop_prev = args["pop_prev"]
        else:
            self.pop_prev = None 


####################################################################################################################


def _estimate_h2_by_ldsc(insumstats, log, verbose=True, **args):
    sumstats = insumstats.copy()
    ##start function with col checking##########################################################
    _start_line = "run LD score regression"
    _end_line = "running LD score regression"
    _start_cols =["CHR","POS","EA","NEA"]
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


####################################################################################################################

def _estimate_partitioned_h2_by_ldsc(insumstats, log, verbose=True, **args):
    sumstats = insumstats.copy()
    ##start function with col checking##########################################################
    _start_line = "run LD score regression"
    _end_line = "running LD score regression"
    _start_cols =["CHR","POS","EA","NEA"]
    _start_function = ".estimate_partitioned_h2_by_ldsc()"
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
    log.write(" -Run partitioned LD score regression:", verbose=verbose)
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
    summary,results = estimate_h2(sumstats, default_args, log)
    
    log.write(" -Results have been stored in .ldsc_h2", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return parse_partitioned_ldsc_summary(summary), results


####################################################################################################################



def _estimate_rg_by_ldsc(insumstats, other_traits ,log, verbose=True, **args):
    sumstats = insumstats.copy()
    ##start function with col checking##########################################################
    _start_line = "run LD score regression for genetic correlation"
    _end_line = "running LD score regression for genetic correlation"
    _start_cols =["CHR","POS","EA","NEA"]
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
    log.write("  -Please cite LDSC: Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)
    
    for key, value in args.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**args)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    other_traits_to_use = []
    alias = default_args.rg.split(",")[1:]

    for index, each_other_sumstats in enumerate(other_traits):
        log.write(" -Processing sumstats with alias {} ({})".format(alias[index], each_other_sumstats.meta["gwaslab"]["study_name"]))
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


####################################################################################################################


def _estimate_h2_cts_by_ldsc(insumstats, log, verbose=True, **args):
    sumstats = insumstats.copy()
    ##start function with col checking##########################################################
    _start_line = "run LD score regression"
    _end_line = "running LD score regression"
    _start_cols =["CHR","POS","EA","NEA"]
    _start_function = ".estimate_h2_cts_by_ldsc()"
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
    log.write(" -Run cell type specific LD score regression:", verbose=verbose)
    log.write("  -Adopted from LDSC source code: https://github.com/bulik/ldsc", verbose=verbose)
    log.write("  -Please cite LDSC: Finucane, H. K., Reshef, Y. A., Anttila, V., Slowikowski, K., Gusev, A., Byrnes, A., ... & Price, A. L. (2018). Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types. Nature genetics, 50(4), 621-629.", verbose=verbose)
    log.write(" -Arguments:", verbose=verbose)
    
    for key, value in args.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**args)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    log.write(" -LDSC log:", verbose=verbose)
    summary = cell_type_specific(sumstats, default_args, log)
    
    log.write(" -Results have been stored in .ldsc_partitioned_h2", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return summary