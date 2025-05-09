import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
from gwaslab.util_in_filter_value import _exclude_hla

def _run_magma(sumstats, 
                magma="magma",
                study="Study1",
                exclude_hla=True,
                window="35,10",
                id_to_use="rsID",
                ref=None,
                ncbi=None,
                set_annot=None,
                out="./",
                delete=True, 
                ncol="N",
                build="19",
                log=Log(), 
                verbose=True):
    
    log.write(" Start to run magma from command line:", verbose=verbose)

    if exclude_hla==True:
        sumstats = _exclude_hla(sumstats, build =build)

    snploc="{}{}.rsid.chr.pos.tsv".format(out,study)
    pval="{}{}.rsid.p.n.tsv".format(out, study)

    log.write(f" -writing temp file for --snp-loc:{snploc}", verbose=verbose)
    sumstats.dropna()[[id_to_use,"CHR","POS"]].rename(columns={id_to_use:"SNP"}).to_csv("{}{}.rsid.chr.pos.tsv".format(out,study),index=None, sep="\t")
    
    log.write(f" -writing temp file for --pval:{pval}", verbose=verbose)
    sumstats.dropna()[[id_to_use,"P","N"]].rename(columns={id_to_use:"SNP"}).to_csv("{}{}.rsid.p.n.tsv".format(out,study),index=None, sep="\t")
    
    log.write(f" --annotate window: {window}", verbose=verbose)
    log.write(f" --gene-loc: {ncbi}", verbose=verbose)
    log.write(f" --bfile: {ref}", verbose=verbose)
    log.write(f" Output prefix: {out}", verbose=verbose)

    bash_script=f'''#!/bin/bash

{magma} --annotate window={window}  --snp-loc {snploc} --gene-loc {ncbi} --out {study}

{magma} --bfile {ref} --pval {pval} ncol={ncol} --gene-annot {study}.genes.annot --out {study}
'''
    
    if set_annot is not None:
        bash_script+=f'''
{magma} --gene-results {study}.genes.raw --set-annot {set_annot} --out {study}
'''
    log.write(f"Script: {bash_script}")

    try:
        log.write(" Running magma from command line...", verbose=verbose)
        output = subprocess.check_output(bash_script, stderr=subprocess.STDOUT, shell=True,text=True)
        output =  output + "\n"
        
        if delete == True:
            os.remove(snploc)
            os.remove(pval)

    except subprocess.CalledProcessError as e:
        log.warning("ERROR!")
        log.write(e.output)

    log.write("Finished running magma.", verbose=verbose)

