import pandas as pd
import numpy as np
import subprocess
from gwaslab.Log import Log
import os

def _process_vcf_and_bfile(chrlist, bfile, vcf, n_cores, plink_log, log, overwrite, load_bim=False):
    # return bfile_path and plink_log
    ref_bims = []
    if bfile is None:
        if vcf is not None:
            log.write(" -Processing VCF : {}...".format(vcf))
            for i in chrlist:
                log.write("  -Processing VCF for CHR {}...".format(i))
                if "@" in vcf:
                    vcf_to_load = vcf.replace("@",str(i))
                else:
                    vcf_to_load = vcf

                bfile_root = vcf.replace(".vcf.gz","")
                bfile_prefix = bfile_root + ".{}".format(i)
                
                if (not os.path.exists(bfile_prefix+".bed")) or overwrite:
                    script_vcf_to_bfile = """
                    plink2 \
                        --vcf {} \
                        --chr {} \
                        --make-bed \
                        --rm-dup force-first \
                        --threads {}\
                        --out {}
                    """.format(vcf_to_load, i, n_cores, bfile_prefix)
                    
                    try:
                        log.write("  -Converting VCF to bed: {}.bim/bed/fam...".format(bfile_prefix))
                        output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                        plink_log+=output + "\n"
                    except subprocess.CalledProcessError as e:
                        log.write(e.output)
                        
                else:
                    log.write("  -Plink bfile for CHR {} exists. Skipping...".format(i))
                
                if load_bim == True:
                    ref_bims = _load_bim(bfile_prefix, ref_bims, log)

            return bfile_root+".@", plink_log, ref_bims
        else:
            log.write("  -Please provide PLINK bfile or VCF as reference!")
    else:
        log.write(" -PLINK bfile as LD reference panel: {}".format(bfile))   
        bfile_prefix = bfile
        if load_bim==True:  
            ref_bims = _load_bim(bfile_prefix, ref_bims, log)
        return bfile, plink_log, ref_bims

def _load_bim(bfile_prefix, ref_bims, log):
    bim_path =bfile_prefix+".bim"
    single_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
    log.write(" -#variants in ref file: {}".format(len(single_bim))) 
    ref_bims.append(single_bim)
    return ref_bims






###########################################################################################################
#def _process_vcf_and_bfile(bfile, vcf, row, n_cores, log):
#    plink_log =""
#    if bfile is None:
#        if vcf is not None:
#            log.write(" -Processing VCF : {}...".format(vcf))
#            bfile_gwaslab = vcf.replace(".vcf.gz","") + ".{}".format(row["CHR"])
#            bfile_prefix= vcf.split("/")[-1].replace(".vcf.gz","") + ".{}".format(row["CHR"])
#            log.write("  -Processing VCF for CHR {}...".format(row["CHR"]))
#            
#            if not os.path.exists(bfile_gwaslab+".bed"):
#                script_vcf_to_bfile = """
#                plink2 \
#                    --vcf {} \
#                    --chr {} \
#                    --make-bed \
#                    --rm-dup force-first \
#                    --threads {}\
#                    --out {}
#                """.format(vcf, row["CHR"], n_cores, bfile_gwaslab)
#                
#                try:
#                    output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
#                    plink_log+=output + "\n"
#                except subprocess.CalledProcessError as e:
#                    log.write(e.output)
#                    
#            else:
#                log.write("  -Plink bfile for CHR {} exists. Skipping...".format(row["CHR"]))
#            ## load bim file
#            bim_path =bfile_gwaslab+".bim"
#            ref_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
#            log.write("#variants in ref file: {}".format(len(ref_bim)))
#        else:
#            log.write("  -Please provide PLINK bfile or VCF as reference!")
#    else:
#        log.write(" -PLINK bfile as LD reference panel: {}".format(bfile))  
#        ## load bim file
#        bim_path =bfile_gwaslab+".bim"
#        ref_bim = pd.read_csv(bim_path,sep="\s+",usecols=[1,3,4,5],header=None,dtype={1:"string",3:"int",4:"string",5:"string"}).rename(columns={1:"SNPID",4:"NEA_bim",5:"EA_bim"})
#        log.write(" -#variants in ref file: {}".format(len(ref_bim)))
#    
#    return ref_bim , plink_log, bfile_gwaslab
#
#def _process_reference(chrlist, bfile, vcf, n_cores, plink_log, log, overwrite):
#    # return bfile_path and plink_log
#    if bfile is None:
#        if vcf is not None:
#            log.write(" -Processing VCF : {}...".format(vcf))
#            for i in chrlist:
#                log.write("  -Processing VCF for CHR {}...".format(i))
#                if "@" in vcf:
#                    vcf_to_load = vcf.replace("@",str(i))
#                else:
#                    vcf_to_load = vcf
#
#                bfile_root = vcf.replace(".vcf.gz","")
#                bfile_prefix = bfile_root + ".{}".format(i)
#                
#                if (not os.path.exists(bfile_prefix+".bed")) or overwrite:
#                    script_vcf_to_bfile = """
#                    plink2 \
#                        --vcf {} \
#                        --chr {} \
#                        --make-bed \
#                        --rm-dup force-first \
#                        --threads {}\
#                        --out {}
#                    """.format(vcf_to_load, i, n_cores, bfile_prefix)
#                    
#                    try:
#                        log.write("  -Converting VCF to bed: {}.bim/bed/fam...".format(bfile_prefix))
#                        output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
#                        plink_log+=output + "\n"
#                    except subprocess.CalledProcessError as e:
#                        log.write(e.output)
#                        
#                else:
#                    log.write("  -Plink bfile for CHR {} exists. Skipping...".format(i))
#                
#            return bfile_root+".@", plink_log
#        else:
#            log.write("  -Please provide PLINK bfile or VCF as reference!")
#    else:
#        log.write(" -PLINK bfile as LD reference panel: {}".format(bfile))      
#        return bfile, plink_log   