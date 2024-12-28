from gwaslab.ldsc_sumstats import estimate_h2
from gwaslab.ldsc_sumstats import estimate_rg
from gwaslab.ldsc_sumstats import cell_type_specific
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.qc_fix_sumstats import skipped
from gwaslab.io_read_ldsc import parse_ldsc_summary
from gwaslab.io_read_ldsc import parse_partitioned_ldsc_summary
from gwaslab.util_in_filter_value import filtervalues
from gwaslab.util_in_filter_value import _filter_palindromic
from gwaslab.util_in_filter_value import _exclude_hla
from gwaslab.util_in_filter_value import _exclude_sexchr

class ARGS():
    def __init__(self, **kwargs):
        
        self.out = "ldsc"

        if "bfile" in kwargs.keys():
            self.bfile = kwargs["bfile"]
        else:
            self.bfile = None 
        
        if "l2" in kwargs.keys():
            self.l2 = kwargs["l2"]
        else:
            self.l2 = None 

        if "extract" in kwargs.keys():
            self.extract = kwargs["extract"]
        else:
            self.extract = None 

        if "keep" in kwargs.keys():
            self.keep = kwargs["keep"]
        else:
            self.keep = None 

        if "ld_wind_snps" in kwargs.keys():
            self.ld_wind_snps = kwargs["ld_wind_snps"]
        else:
            self.ld_wind_snps = None 

        if "ld_wind_kb" in kwargs.keys():
            self.ld_wind_kb = kwargs["ld_wind_kb"]
        else:
            self.ld_wind_kb = None 

        if "ld_wind_cm" in kwargs.keys():
            self.ld_wind_cm = kwargs["ld_wind_cm"]
        else:
            self.ld_wind_cm = None 

        if "print_snps" in kwargs.keys():
            self.print_snps = kwargs["print_snps"]
        else:
            self.print_snps = None 

        if "annot" in kwargs.keys():
            self.annot = kwargs["annot"]
        else:
            self.annot = None 

        if "thin_annot" in kwargs.keys():
            self.thin_annot = kwargs["thin_annot"]
        else:
            self.thin_annot = None 

        if "cts_bin" in kwargs.keys():
            self.cts_bin = kwargs["cts_bin"]
        else:
            self.cts_bin = None 

        if "cts_breaks" in kwargs.keys():
            self.cts_breaks = kwargs["cts_breaks"]
        else:
            self.cts_breaks = None 

        if "cts_names" in kwargs.keys():
            self.cts_names = kwargs["cts_names"]
        else:
            self.cts_names = None 

        if "per_allele" in kwargs.keys():
            self.per_allele = kwargs["per_allele"]
        else:
            self.per_allele = None 

        if "pq_exp" in kwargs.keys():
            self.pq_exp = kwargs["pq_exp"]
        else:
            self.pq_exp = None 

        if "no_print_annot" in kwargs.keys():
            self.no_print_annot = kwargs["no_print_annot"]
        else:
            self.no_print_annot = None 

        if "h2" in kwargs.keys():
            self.h2 = kwargs["h2"]
        else:
            self.h2 = None

        if "h2_cts" in kwargs.keys():
            self.h2_cts = kwargs["h2_cts"]
        else:
            self.h2_cts = None

        if "rg" in kwargs.keys():
            self.rg = kwargs["rg"]
        else:
            self.rg = None

        if "ref_ld" in kwargs.keys():
            self.ref_ld = kwargs["ref_ld"]
        else:
            self.ref_ld = None

        if "ref_ld_chr" in kwargs.keys():
            self.ref_ld_chr = kwargs["ref_ld_chr"]
        else:
            self.ref_ld_chr = None

        if "w_ld" in kwargs.keys():
            self.w_ld = kwargs["w_ld"]
        else:
            self.w_ld = None
        
        if "w_ld_chr" in kwargs.keys():
            self.w_ld_chr = kwargs["w_ld_chr"]
        else:
            self.w_ld_chr = None     

        if "overlap_annot" in kwargs.keys():
            self.overlap_annot = kwargs["overlap_annot"]
        else:
            self.overlap_annot = None 

        if "print_coefficients" in kwargs.keys():
            self.print_coefficients = kwargs["print_coefficients"]
        else:
            self.print_coefficients = "ldsc" 

        if "frqfile" in kwargs.keys():
            self.frqfile = kwargs["frqfile"]
        else:
            self.frqfile = None 

        if "frqfile_chr" in kwargs.keys():
            self.frqfile_chr = kwargs["frqfile_chr"]
        else:
            self.frqfile_chr = None 

        if "no_intercept" in kwargs.keys():
            self.no_intercept = kwargs["no_intercept"]
        else:
            self.no_intercept = None 

        if "intercept_h2" in kwargs.keys():
            self.intercept_h2 = kwargs["intercept_h2"]
        else:
            self.intercept_h2 = None 

        if "intercept_gencov" in kwargs.keys():
            self.intercept_gencov = kwargs["intercept_gencov"]
        else:
            self.intercept_gencov = None 

        if "M" in kwargs.keys():
            self.M = kwargs["M"]
        else:
            self.M = None 

        if "two_step" in kwargs.keys():
            self.two_step = kwargs["two_step"]
        else:
            self.two_step = None 

        if "chisq_max" in kwargs.keys():
            self.chisq_max = kwargs["chisq_max"]
        else:
            self.chisq_max= None 

        if "ref_ld_chr_cts" in kwargs.keys():
            self.ref_ld_chr_cts = kwargs["ref_ld_chr_cts"]
        else:
            self.ref_ld_chr_cts = None

        if "print_all_cts" in kwargs.keys():
            self.print_all_cts = kwargs["print_all_cts"]
        else:
            self.print_all_cts = False 

        if "print_cov" in kwargs.keys():
            self.print_cov = kwargs["print_cov"]
        else:
            self.print_cov = None 

        self.print_delete_vals = False
        if "print_delete_vals" in kwargs.keys():
            self.print_delete_vals = kwargs["print_delete_vals"]
        else:
            self.print_delete_vals = False 

        if "chunk_size" in kwargs.keys():
            self.chunk_size = kwargs["chunk_size"]
        else:
            self.chunk_size = 50 

        if "pickle" in kwargs.keys():
            self.pickle = kwargs["pickle"]
        else:
            self.pickle = False 

        if "yes_really" in kwargs.keys():
            self.yes_really = kwargs["yes_really"]
        else:
            self.yes_really = False 

        if "invert_anyway" in kwargs.keys():
            self.invert_anyway = kwargs["invert_anyway"]
        else:
            self.invert_anyway = False 

        if "n_blocks" in kwargs.keys():
            self.n_blocks = kwargs["n_blocks"]
        else:
            self.n_blocks = 200 

        if "not_M_5_50" in kwargs.keys():
            self.not_M_5_50 = kwargs["not_M_5_50"]
        else:
            self.not_M_5_50 = False 

        if "no_check_alleles" in kwargs.keys():
            self.no_check_alleles = kwargs["no_check_alleles"]
        else:
            self.no_check_alleles = False 
        
        if "return_silly_things" in kwargs.keys():
            self.return_silly_things = kwargs["return_silly_things"]
        else:
            self.return_silly_things = False 

        if "samp_prev" in kwargs.keys():
            self.samp_prev = kwargs["samp_prev"]
        else:
            self.samp_prev = None 
        
        if "pop_prev" in kwargs.keys():
            self.pop_prev = kwargs["pop_prev"]
        else:
            self.pop_prev = None 


####################################################################################################################


def _estimate_h2_by_ldsc(insumstats, log, verbose=True, munge=False, munge_args=None, **kwargs):
    sumstats = insumstats.copy()
    
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].astype("int64")

    if munge:
        if munge_args is None:
            munge_args={}
        log.write("Start to munge sumstats.")
        sumstats = _munge_sumstats(sumstats, log=log, verbose=verbose,**munge_args)
        log.write("Finished munging sumstats.")

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
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    default_args = ARGS(**kwargs)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    log.write(" -LDSC log:", verbose=verbose)
    summary = estimate_h2(sumstats, default_args, log)
    
    results_table = None
    if type(summary) is tuple:
        results_table = summary[1]
        summary = summary[0]
        log.write(" -Coefficient results have been stored in .ldsc_h2_results", verbose=verbose)
        

    log.write(" -Results have been stored in .ldsc_h2", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return parse_ldsc_summary(summary), results_table


####################################################################################################################

def _estimate_partitioned_h2_by_ldsc(insumstats, log, verbose=True, **kwargs):
    sumstats = insumstats.copy()
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].astype("int64")
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
    
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**kwargs)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    log.write(" -LDSC log:", verbose=verbose)
    summary,results = estimate_h2(sumstats, default_args, log)
    
    log.write(" -Results have been stored in .ldsc_h2", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return parse_partitioned_ldsc_summary(summary), results


####################################################################################################################



def _estimate_rg_by_ldsc(insumstats, other_traits ,log, verbose=True, **kwargs):
    sumstats = insumstats.copy()
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].astype("int64")
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
    
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**kwargs)

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


def _estimate_h2_cts_by_ldsc(insumstats, log, verbose=True, **kwargs):
    sumstats = insumstats.copy()
    if "N" in sumstats.columns:
        sumstats["N"] = sumstats["N"].astype("int64")
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
    
    for key, value in kwargs.items():
        log.write("  -{}:{}".format(key, value), verbose=verbose)
    
    default_args = ARGS(**kwargs)

    if "Z" not in sumstats.columns:
        sumstats["Z"] = sumstats["BETA"]/sumstats["SE"]

    sumstats = sumstats.rename(columns={"EA":"A1","NEA":"A2","rsID":"SNP"})
    
    log.write(" -LDSC log:", verbose=verbose)
    summary = cell_type_specific(sumstats, default_args, log)
    
    log.write(" -Results have been stored in .ldsc_partitioned_h2", verbose=verbose)
    finished(log=log,verbose=verbose,end_line=_end_line)
    return summary



def _munge_sumstats(sumstats, log, 
                    info=0.9, maf=0.01, 
                    n=None, nopalindromic=True,
                    exclude_hla=True, exclude_sexchr=True,  
                    verbose=True, **kwargs):
    if "CHR" in sumstats.columns and "POS" in sumstats.columns:
        if exclude_hla == True:
            sumstats = _exclude_hla(sumstats, verbose=verbose, log=log)
    
    if "CHR" in sumstats.columns:
        if exclude_sexchr == True:
            sumstats = _exclude_sexchr(sumstats, verbose=verbose, log=log)
    
    # filter_info
    if "INFO" in sumstats.columns:
        sumstats = filtervalues(sumstats, 'INFO >={}'.format(info) ,verbose=verbose, log=log)
    
    # frequency
    if "EAF" in sumstats.columns:
        sumstats = filtervalues(sumstats,'EAF>={} and EAF<={}'.format(maf, 1-maf),verbose=verbose, log=log)
    
    # N
    if "N" in sumstats.columns:
        if n is None:
            min_n = sumstats.N.quantile(0.9) / 1.5
        else:
            min_n = n
        sumstats = filtervalues(sumstats,'N>={}'.format(min_n),verbose=verbose, log=log)
    
    # remove strand-unambiguous SNPs
    if "EA" in sumstats.columns and "NEA" in sumstats.columns:
        if nopalindromic==True:
            sumstats = _filter_palindromic(sumstats, mode="out", verbose=verbose, log=log)
    
    return sumstats
