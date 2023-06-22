import pandas as pd
import numpy as np
import scipy.stats as ss
from gwaslab.Log import Log
import scipy as sp

def get_power(
              mode="b",
              t=0,
              genotype_or=1.3 ,
              beta=0.3,
              eaf=0.1,
              n=10000,
              scase= 2000,
              scontrol= 15000,
              prevalence= 0.15,
              prevalence_base=None,
              daf = 0.2,
              sig_level= 5e-8,
              vary=1,
              log=Log(),
              verbose=True
             ):
    if mode=="b":
        print("Input settings:{}".format(daf))
        print(" -Number of cases:{}".format(scase))
        print(" -Number of controls:{}".format(scontrol))
        print(" -Risk allele OR:{:.3f}".format(genotype_or))
        print(" -Disease prevalence:{:.3f}".format(prevalence))
        print(" -Risk allele frequency: {:.3f}".format(daf))
        print(" -Significance level: {:.3e}".format(sig_level))
        # Skol, A. D., Scott, L. J., Abecasis, G. R., & Boehnke, M. (2006). Joint analysis is more efficient than replication-based analysis for two-stage genome-wide association studies. Nature genetics, 38(2), 209-213.
        aaf = daf**2
        abf = 2 * (daf) * (1 - daf)
        bbf = (1- daf)**2

        if prevalence_base is None:
            genotype_or = np.exp(beta)
            genotype_rr = genotype_or
        else:
            genotype_or = np.exp(beta)
            genotype_rr = genotype_or/ ((1-prevalence_base)+(genotype_or*prevalence_base))
            # https://jamanetwork.com/journals/jama/fullarticle/188182
        
        # additive
        x = [ 2*genotype_rr-1, genotype_rr, 1 ] 

        aap= x[0] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
        abp= x[1] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
        bbp= x[2] * prevalence / (x[0]*aaf + x[1]*abf + x[2]*bbf)
        print("Probability of disease :")
        print(" - Individuals with AA genotype: {:.3f}".format(aap))
        print(" - Individuals with AB genotype: {:.3f}".format(abp))
        print(" - Individuals with BB genotype: {:.3f}".format(bbp))
        
        pcase= (aap * aaf + abp * abf*0.5) / prevalence
        pcontrol=((1-aap )* aaf + (1-abp )* abf*0.5) / (1 - prevalence)

        vcase = pcase *(1-pcase)
        vcontrol =pcontrol *(1-pcontrol)
        print("Expected risk allele frequency:")
        print(" - In cases: {:.3f}".format(pcase))
        print(" - In controls: {:.3f}".format(pcontrol))

        num= (pcase - pcontrol)
        den= np.sqrt( (vcase/scase +  vcontrol/scontrol)*0.5 )
        u = num / den

        c = ss.norm.isf(sig_level/2)
        power = 1 - ss.norm.cdf(c-u) + ss.norm.cdf(-c-u)
        print("Expected power: {:.3f}".format(power))

    elif mode=="q":
        if verbose:
            log.write("Significance level: {}".format(sig_level))
            log.write("EAF: {}".format(eaf))
            log.write("BETA: {}".format(beta))
            log.write("N: {}".format(n))
            log.write("H2: {}".format(2*eaf*(1-eaf)*(beta**2)))
        c = ss.chi2.isf(sig_level/2,df=1)
        NCP = n * 2*eaf*(1-eaf)*(beta**2)/vary
        power = 1 - ss.ncx2.cdf(c,df=1,nc=NCP)

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
                
                c = ss.chi2.isf(sig_level/2,df=1)
                h2 = 2*eaf*(1-eaf)*(beta**2)
                NCP = n * h2/vary
                power = 1 - ss.ncx2.cdf(c,df=1,nc=NCP)
                return power
            
            eaf_beta_matrix = np.zeros((n_matrix,n_matrix),dtype=float)
            eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
            betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
            
            if verbose: log.write(" -Updating eaf-beta matrix...")
            for i in range(n_matrix):
                    eaf_beta_matrix[i,] = calculate_power_single(beta=betas,eaf=eafs[i],n=n,sig_level=sig_level,vary=vary)
            
            if verbose: log.write(" -Extracting eaf-beta combinations with power = {}...".format(t))
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
              scase=None,
              scontrol=None,
              eaf_range=(0.0001,0.5),
              beta_range=(0.0001,10),
              t=0,
              sig_level= 5e-8,
              vary=1,
              log=Log(),
              verbose=True,
              n_matrix=500,
              prevalence_base=None
             ):
    if t >0:
        def calculate_power_single(
                            beta, 
                            daf, 
                            prevalence,
                            scase, 
                            scontrol, 
                            sig_level=5e-8,
                            prevalence_base=None):
                
                aaf = daf**2
                abf = 2 * (daf) * (1 - daf)
                bbf = (1- daf)**2

                if prevalence_base is None:
                    genotype_or = np.exp(beta)
                    genotype_rr = genotype_or
                else:
                    genotype_or = np.exp(beta)
                    genotype_rr = genotype_or/ ((1-prevalence_base)+(genotype_or*prevalence_base))
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
                den= np.sqrt( (vcase/scase +  vcontrol/scontrol)*0.5 )
                u = num / den
                c = ss.norm.isf(sig_level/2)
                power = 1 - ss.norm.cdf(c-u) + ss.norm.cdf(-c-u)
                return power
        
        eaf_beta_matrix = np.zeros((n_matrix,n_matrix),dtype=float)
        eafs = np.linspace(eaf_range[1],eaf_range[0],n_matrix)
        betas =  np.linspace(beta_range[0],beta_range[1],n_matrix)
        
        if verbose: log.write(" -Updating eaf-beta matrix...")
        if prevalence_base is None:
            if verbose: log.write(" -Alogorithm: Skol, Andrew D., et al. Nature genetics 38.2 (2006): 209-213....")
            if verbose: log.write(" -GRR is approximated using OR. For prevalence < 10%, GRR is very similar to OR....")
        else:
            if verbose: log.write(" -OR is converted to GRR using base prevalence: {}".format(prevalence_base))

        for i in range(n_matrix):
                eaf_beta_matrix[i,] = calculate_power_single(beta=betas,
                                                                daf=eafs[i],
                                                                scase=scase,
                                                                scontrol=scontrol,
                                                                prevalence=prevalence,
                                                                sig_level=sig_level,
                                                                prevalence_base=prevalence_base)
        
        if verbose: log.write(" -Extracting eaf-beta combinations with power = {}...".format(t))
        i,j=1,1
        eaf_beta = []
        while i<n_matrix-1 and j<n_matrix-1:
            if eaf_beta_matrix[i,j] < t:
                j+=1
            else:
                i+=1
                eaf_beta.append((eafs[i],betas[j]))
    return pd.DataFrame(eaf_beta)
    
    