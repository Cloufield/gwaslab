def get_power(genotype_or=1.3 ,
              scase= 2000,
              scontrol= 15000,
              prevalence= 0.15,
              daf = 0.2,
              sig_level= 5e-8
             ):
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

    # additive
    x = [ 2*genotype_or-1, genotype_or, 1 ] 

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
    return power