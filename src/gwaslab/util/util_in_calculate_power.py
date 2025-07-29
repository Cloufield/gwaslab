import pandas as pd
import numpy as np
import scipy.stats as ss
from gwaslab.g_Log import Log
import scipy as sp

def get_power(
              mode="b",
              genotype_rr=None ,
              genotype_or=None ,
              beta=0.3,
              eaf=0.1,
              n=10000,
              ncase= 2000,
              ncontrol= 15000,
              prevalence= 0.15,
              or_to_rr=False,
              daf = 0.2,
              sig_level= 5e-8,
              vary=1,
              log=Log(),
              verbose=True
             ):
    log.write(" Start to calculate statistical power...", verbose=verbose)
    if mode=="b":
        log.write(" -Input settings (b mode):", verbose=verbose)
        log.write("  -Number of cases:{}".format(ncase), verbose=verbose)
        log.write("  -Number of controls:{}".format(ncontrol), verbose=verbose)
        if genotype_rr is not None:
            log.write("  -Risk allele RR:{:.3f}".format(genotype_rr), verbose=verbose)
        elif genotype_or is not None:
            log.write("  -Risk allele OR:{:.3f}".format(genotype_or), verbose=verbose)
        elif beta is not None:
            log.write("  -Risk allele beta:{:.3f}".format(beta), verbose=verbose)
        else:
            genotype_rr = 0.1
            log.write("  -Risk allele RR:{:.3f}".format(genotype_rr), verbose=verbose)
        log.write("  -Disease prevalence:{:.3f}".format(prevalence), verbose=verbose)
        log.write("  -Risk allele frequency: {:.3f}".format(daf), verbose=verbose)
        log.write("  -Significance level: {:.3e}".format(sig_level), verbose=verbose)
        # Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). Joint analysis is more efficient than replication-based analysis for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
        aaf = daf**2
        abf = 2 * (daf) * (1 - daf)
        bbf = (1- daf)**2

        if genotype_rr is None:
            if or_to_rr == False:
                if genotype_or is None:
                    genotype_or = np.exp(beta)
                genotype_rr = genotype_or
            else:
                if genotype_or is None:
                    genotype_or = np.exp(beta)
                genotype_rr = genotype_or/ ((1-prevalence)+(genotype_or*prevalence))
                # https://jamanetwork.com/journals/jama/fullarticle/188182
        
            if or_to_rr ==False:
                log.write(" -Alogorithm: Skol, Andrew D., et al. Nature genetics 38.2 (2006): 209-213....", verbose=verbose)
                log.write(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....", verbose=verbose)
            else:
                log.write(" -OR is converted to GRR using base prevalence: {}".format(prevalence), verbose=verbose)
                log.write(" -Alogorithm: Zhang, J., & Kai, F. Y. (1998). What's the relative risk?: A method of correcting the odds ratio in cohort studies of common outcomes. Jama, 280(19), 1690-1691.....", verbose=verbose)
            
        # additive
        x = [ 2*genotype_rr-1, genotype_rr, 1 ] 

        aap= x[0] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
        abp= x[1] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
        bbp= x[2] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
        log.write("Probability of disease :", verbose=verbose)
        log.write(" - Individuals with AA genotype: {:.3f}".format(aap), verbose=verbose)
        log.write(" - Individuals with AB genotype: {:.3f}".format(abp), verbose=verbose)
        log.write(" - Individuals with BB genotype: {:.3f}".format(bbp), verbose=verbose)
        
        pcase= (aap * aaf + abp * abf*0.5) / prevalence
        pcontrol=((1-aap )* aaf + (1-abp )* abf*0.5) / (1 - prevalence)

        vcase = pcase *(1-pcase)
        vcontrol =pcontrol *(1-pcontrol)
        log.write("Expected risk allele frequency:", verbose=verbose)
        log.write(" - In cases: {:.3f}".format(pcase), verbose=verbose)
        log.write(" - In controls: {:.3f}".format(pcontrol), verbose=verbose)

        num= (pcase - pcontrol)
        den= np.sqrt( (vcase/ncase +  vcontrol/ncontrol)*0.5 )
        u = num / den

        c = ss.norm.isf(sig_level/2)
        power = 1 - ss.norm.cdf(c-u) + ss.norm.cdf(-c-u)
        log.write("Expected power: {:.3f}".format(power), verbose=verbose)

    elif mode=="q":
        if beta is None:
            beta = 0.1

        log.write(" -Input settings (q mode):", verbose=verbose)
        log.write("  -Significance level: {}".format(sig_level), verbose=verbose)
        log.write("  -EAF: {}".format(eaf), verbose=verbose)
        log.write("  -BETA: {}".format(beta), verbose=verbose)
        log.write("  -N: {}".format(n), verbose=verbose)
        log.write("  -SNPR2: {}".format(2*eaf*(1-eaf)*(beta**2)), verbose=verbose)
        c = ss.chi2.isf(sig_level,df=1)
        NCP = n * 2*eaf*(1-eaf)*(beta**2)/vary
        power = 1 - ss.ncx2.cdf(c, df=1, nc=NCP)
    log.write("Finished calculating statistical power.", verbose=verbose)
    return power

def get_beta(
              mode="q",
              eaf_range=(0.0001,0.5),
              beta_range=(0.0001,10),
              t=0,
              n=None,
              sig_level= 5e-8,
              vary=1,
              log=Log(),
              verbose=True,
              n_matrix=500
             ):
    if mode=="q":
        if t >0:
            def calculate_power_single(
                                beta, 
                                eaf, 
                                n, 
                                sig_level=5e-8,
                                vary=1):
                
                c = ss.chi2.isf(sig_level,df=1)
                h2 = 2*eaf*(1-eaf)*(beta**2)
                NCP = n * h2/vary
                power = 1 - ss.ncx2.cdf(c,df=1,nc=NCP)
                return power
            
            eaf_beta_matrix = np.zeros((n_matrix,n_matrix),dtype=float)
            eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
            betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
            
            log.write(" -Updating eaf-beta matrix...", verbose=verbose)
            for i in range(n_matrix):
                    eaf_beta_matrix[i,] = calculate_power_single(beta=betas,eaf=eafs[i],n=n,sig_level=sig_level,vary=vary)
            
            log.write(" -Extracting eaf-beta combinations with power = {}...".format(t), verbose=verbose)
            i,j=1,1
            eaf_beta = []
            while i<n_matrix-1 and j<n_matrix-1:
                if eaf_beta_matrix[i,j] < t:
                    j+=1
                else:
                    i+=1
                    eaf_beta.append((eafs[i],betas[j]))
        return pd.DataFrame(eaf_beta)

def get_beta_binary(
              prevalence=None,
              ncase=None,
              ncontrol=None,
              eaf_range=(0.0001,0.5),
              beta_range=(0.0001,10),
              t=0,
              sig_level= 5e-8,
              vary=1,
              log=Log(),
              verbose=True,
              n_matrix=500,
              or_to_rr=False
             ):
    if t >0:
        def calculate_power_single(
                            beta, 
                            daf, 
                            prevalence,
                            ncase, 
                            ncontrol, 
                            sig_level=5e-8,
                            or_to_rr=False):
                
                aaf = daf**2
                abf = 2 * (daf) * (1 - daf)
                bbf = (1- daf)**2

                if or_to_rr == False:
                    genotype_or = np.exp(beta)
                    genotype_rr = genotype_or
                else:
                    genotype_or = np.exp(beta)
                    genotype_rr = genotype_or/ ((1-prevalence)+(genotype_or*prevalence))
                    # https://jamanetwork.com/journals/jama/fullarticle/188182
                # additive
                x = [ 2*genotype_rr-1, genotype_rr, 1 ] 
                aap= x[0] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
                abp= x[1] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
                bbp= x[2] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
                pcase= (aap * aaf + abp * abf*0.5) / prevalence
                pcontrol=((1-aap )* aaf + (1-abp )* abf*0.5) / (1 - prevalence)
                vcase = pcase *(1-pcase)
                vcontrol =pcontrol *(1-pcontrol)
                num= (pcase - pcontrol)
                den= np.sqrt( (vcase/ncase +  vcontrol/ncontrol)*0.5 )
                u = num / den
                c = ss.norm.isf(sig_level/2)
                power = 1 - ss.norm.cdf(c-u) + ss.norm.cdf(-c-u)
                return power
        
        eaf_beta_matrix = np.zeros((n_matrix,n_matrix),dtype=float)
        eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
        betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
        
        log.write(" -Updating eaf-beta matrix...", verbose=verbose)
        if or_to_rr ==False:
            log.write(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....", verbose=verbose)
        else:
            log.write(" -OR is converted to GRR using base prevalence: {}".format(prevalence), verbose=verbose)
        
        for i in range(n_matrix):
                eaf_beta_matrix[i,] = calculate_power_single(beta=betas,
                                                                daf=eafs[i],
                                                                ncase=ncase,
                                                                ncontrol=ncontrol,
                                                                prevalence=prevalence,
                                                                sig_level=sig_level,
                                                                or_to_rr=or_to_rr)
        
        log.write(" -Extracting eaf-beta combinations with power = {}...".format(t), verbose=verbose)
        i,j=1,1
        eaf_beta = []
        while i<n_matrix-1 and j<n_matrix-1:
            if eaf_beta_matrix[i,j] < t:
                j+=1
            else:
                i+=1
                eaf_beta.append((eafs[i],betas[j]))
    return pd.DataFrame(eaf_beta)
    
    