import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.g_version import _checking_r_version
from gwaslab.g_version import _check_susie_version
from gwaslab.qc_fix_sumstats import start_to
from gwaslab.qc_fix_sumstats import finished
from gwaslab.util_ex_calculate_ldmatrix import _extract_variants_in_locus
from gwaslab.util_in_get_sig import getsig

def _run_hyprcoloc(  sumstats_multi, 
                     r="Rscript",
                     study="Group1",
                     traits=None,
                     types=None, 
                     loci=None,
                     nstudy=2,
                     windowsizekb=1000,
                     build="99",
                     log=Log(), 
                     verbose=True):
    
    log.write(" Start to run hyprcoloc from command line:", verbose=verbose)
    if traits is None:
         traits_to_form_string = [ '"trait_{}"'.format(i+1) for i in range(nstudy)]
    else:
         traits_to_form_string = ['"{}"'.format(i) for i in traits]
    
    hyprcoloc_res_combined = pd.DataFrame()

    if loci is None:
        log.write(" -Loci were not provided. All significant loci will be automatically extracted...",verbose=verbose)
        sig_df = getsig(sumstats_multi,id="SNPID",chrom="CHR",pos="POS",p="P_MIN", build=build)
    else:
        sig_df = sumstats_multi.loc[sumstats_multi["SNPID"].isin(loci),:]

    for index,row  in  sig_df.iterrows():
    #row = sig_df.iloc[0,:]
        
        # extract locus
        locus = row["SNPID"]
        log.write(" -Running hyprcoloc for locus : {}...".format(locus),verbose=verbose)

        # prepare input files sumstats_multi
        output_beta_cols = []
        output_se_cols = []

        for i in range(nstudy):
            output_beta_cols.append("BETA_{}".format(i+1))
            output_se_cols.append("SE_{}".format(i+1))

        matched_sumstats = _extract_variants_in_locus(sumstats_multi, windowsizekb, locus = (row["CHR"],row["POS"]))
        
        to_export = matched_sumstats[["SNPID"] + output_se_cols + output_beta_cols].dropna()
        
        if len(to_export)>0:
            log.write(" -Number of shared variants in locus {} : {}...".format(locus, len(to_export)),verbose=verbose)
            to_export[["SNPID"] + output_beta_cols].to_csv("{}_{}_beta_cols.tsv.gz".format(study,locus), index=None,sep="\t")
            to_export[["SNPID"] + output_se_cols].to_csv("{}_{}_se_cols.tsv.gz".format(study,locus), index=None,sep="\t")
        else:
            log.write(" -No shared variants in locus {}...skipping".format(locus),verbose=verbose)
            continue

        r_log=""
        log = _checking_r_version(r, log)

        rscript='''
library(hyprcoloc)

betas<-read.csv("{study}_{locus}_beta_cols.tsv.gz",row.names = 1,sep="\\t")
ses  <-read.csv("{study}_{locus}_se_cols.tsv.gz",row.names = 1,sep="\\t")

betas <- as.matrix(betas)
ses <- as.matrix(ses)

traits <- c({traits_string})
snpid <- rownames(betas)

res <- hyprcoloc(betas, 
          ses, 
          trait.names = traits, 
          snp.id = snpid,
          snpscore=TRUE)

write.csv(res[[1]], "{study}_{locus}_{nstudy}studies.res",row.names = FALSE) 
        '''.format(
            study=study,
            locus=locus,
            nstudy=nstudy,
            traits_string = ','.join(traits_to_form_string)
        )

        output_path = "{study}_{locus}_{nstudy}studies.res".format(study=study, locus=locus, nstudy=nstudy)
        output_prefix = "{study}_{locus}_{nstudy}studies".format(study=study, locus=locus, nstudy=nstudy)
        
        with open("_{}_{}_gwaslab_hyprcoloc_temp.R".format(study,locus),"w") as file:
                file.write(rscript)

        script_run_r = "{} _{}_{}_gwaslab_hyprcoloc_temp.R".format(r, study,locus)
        
        try:
            log.write(" Running hyprcoloc from command line...", verbose=verbose)
            output = subprocess.check_output(script_run_r, stderr=subprocess.STDOUT, shell=True,text=True)
            r_log+= output + "\n"
            #os.remove("_{}_{}_gwaslab_hyprcoloc_temp.R".format(study,locus))
            hyprcoloc_res = pd.read_csv(output_path)    
            hyprcoloc_res["PREFIX"] = output_prefix 
            hyprcoloc_res_combined = pd.concat([hyprcoloc_res_combined, hyprcoloc_res],ignore_index=True)
        except subprocess.CalledProcessError as e:
            log.write(e.output)
            #os.remove("_{}_{}_gwaslab_hyprcoloc_temp.R".format(study,locus))
        log.write(" -Finishing hyprcoloc for locus : {}...".format(locus),verbose=verbose)
    log.write("Finished clocalization using hyprcoloc.", verbose=verbose)
    return hyprcoloc_res_combined