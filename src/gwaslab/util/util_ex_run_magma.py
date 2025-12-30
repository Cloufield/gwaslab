from typing import TYPE_CHECKING, Optional
import subprocess
import os
import gc
import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_filter_value import _exclude_hla

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

def _run_magma(
    gls: 'Sumstats',
    magma: str = "magma",
    study: str = "Study1",
    exclude_hla: bool = True,
    window: str = "35,10",
    id_to_use: str = "rsID",
    ref: Optional[str] = None,
    ncbi: Optional[str] = None,
    set_annot: Optional[str] = None,
    out: Optional[str] = None,
    delete: bool = True,
    ncol: str = "N",
    build: str = "19",
    log: Log = Log(),
    verbose: bool = True
) -> None:
    
    log.write(" Start to run magma from command line:", verbose=verbose)

    sumstats = gls.data
    gls.offload() 
    if exclude_hla==True:
        sumstats = _exclude_hla(sumstats, build =build)
    
    if out is None:
        out = os.path.join("./", study)
    else:
        out = os.path.join(out, study)

    snploc="{}.rsid.chr.pos.tsv".format(out)
    pval="{}.rsid.p.n.tsv".format(out)

    log.write(f" -writing temp file for --snp-loc:{snploc}", verbose=verbose)
    sumstats.dropna()[[id_to_use,"CHR","POS"]].rename(columns={id_to_use:"SNP"}).to_csv("{}.rsid.chr.pos.tsv".format(out),index=None, sep="\t")
    
    log.write(f" -writing temp file for --pval:{pval}", verbose=verbose)
    sumstats.dropna()[[id_to_use,"P","N"]].rename(columns={id_to_use:"SNP"}).to_csv("{}.rsid.p.n.tsv".format(out),index=None, sep="\t")
    
    log.write(f" --annotate window: {window}", verbose=verbose)
    log.write(f" --gene-loc: {ncbi}", verbose=verbose)
    log.write(f" --bfile: {ref}", verbose=verbose)
    log.write(f" Output prefix: {out}", verbose=verbose)
    
    
    
    bash_script=f'''#!/bin/bash

{magma} --annotate window={window}  --snp-loc {snploc} --gene-loc {ncbi} --out {out}

{magma} --bfile {ref} --pval {pval} ncol={ncol} --gene-annot {out}.genes.annot --out {out}
'''
    
    if set_annot is not None:
        bash_script+=f'''
{magma} --gene-results {out}.genes.raw --set-annot {set_annot} --out {out}
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

    gls.reload() 
    log.write("Finished running magma.", verbose=verbose)

