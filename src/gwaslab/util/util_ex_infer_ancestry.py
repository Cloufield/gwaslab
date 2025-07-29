
import pandas as pd
from gwaslab.g_Log import Log

def _infer_ancestry(sumstats, 
                    ancestry_af=None,
                    build="19",
                    log=Log(),
                    verbose=True):
    log.write("Start to infer ancestry based on Fst...", verbose=verbose)
    ref_af = pd.read_csv(ancestry_af, sep="\t")
    
    data_af = pd.merge(sumstats[["CHR","POS","EA","NEA","EAF"]] ,ref_af,on=["CHR","POS"],how="inner") 

    log.write(f"  -Estimating Fst using {len(data_af)} variants...", verbose=verbose)

    is_filp = data_af["EA"] == data_af["ALT"]
    data_af.loc[is_filp, ["EA","NEA"]] = data_af.loc[is_filp, ["NEA","EA"]]
    data_af.loc[is_filp, "EAF"] = 1 - data_af.loc[is_filp, "EAF"]
    
    headers = []
    for i in ['GBR', 'FIN', 'CHS', 'PUR', 'CDX',
        'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'ACB', 'GWD', 'ESN', 'BEB', 'MSL',
        'STU', 'ITU', 'CEU', 'YRI', 'CHB', 'JPT', 'LWK', 'ASW', 'MXL', 'TSI',
        'GIH', 'EUR', 'EAS', 'AMR', 'SAS', 'AFR']:
        headers.append(f"FST_{i}")
        data_af[f"FST_{i}"] = data_af.apply(lambda x: calculate_fst(x["EAF"], x[i]), axis=1)
    
    for i,value in data_af[headers].mean().sort_values().items():
        log.write( f"  -{i} : {value}", verbose=verbose)
        
        closest_ancestry = data_af[headers].mean().sort_values().idxmin()

    log.write(f"  -Closest Ancestry: {closest_ancestry.split('_')[1]}", verbose=verbose)
    log.write("Finished inferring ancestry.", verbose=verbose)
    return closest_ancestry.split("_")[1]

def calculate_fst(p_1, p_2):
    # https://bios1140.github.io/understanding-fst-the-fixation-index.html
    # calculate q1 and q2
    q_1 = 1 - p_1
    q_2 = 1 - p_2

    # calculate total allele frequency
    p_t = (p_1 + p_2)/2
    q_t = 1 - p_t

    # calculate expected heterozygosity
    # first calculate expected heterozygosity for the two populations
    # pop1
    hs_1 = 2*p_1*q_1
    # pop2
    hs_2 = 2*p_2*q_2
    # then take the mean of this
    hs = (hs_1 + hs_2)/2

    # next calculate expected heterozygosity for the metapopulations
    ht = 2*p_t*q_t

    # calculate fst
    fst = (ht - hs)/ht

    # return output
    return fst
    