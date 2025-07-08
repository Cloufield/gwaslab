import pandas as pd
from gwaslab.g_Log import Log
import pandas as pd
import numpy as np

def _extract_associations(sumstats, rsid="rsID", log = Log(), verbose=True):

    assoc, traits, studies = get_associations_from_gwascatalog(sumstats, rsid=rsid, log=log, verbose=verbose)
    
    assoc = _fix_beta(assoc)

    traits_agg = traits.groupby("associationId")[["trait","shortForm"]].agg(lambda x: ",".join(x)).reset_index()
    
    assoc_traits_agg= pd.merge(assoc, traits_agg, on ="associationId",how="left")
    
    assoc_traits_agg= pd.merge(assoc_traits_agg, studies, on ="associationId",how="left")
    
    summary_columns=['associationId', 'rsID', 'RA', 'riskFrequency','betaNum', 'pvalue','trait','cohort','initialSampleSize','publicationInfo.pubmedId']
    
    assoc_traits_agg_summary = assoc_traits_agg[summary_columns]
    
    return assoc_traits_agg, assoc_traits_agg_summary

def get_associations_from_gwascatalog(sumstats, rsid="rsID", log=Log(), verbose=True):
    from pandasgwas import get_associations
    from pandasgwas import get_traits
    from pandasgwas import get_studies

    association = pd.DataFrame()
    strongest_risk_alleles=pd.DataFrame()
    
    unique_sumstats = sumstats.dropna(subset=[rsid]).drop_duplicates(subset=[rsid])
    
    for index,row in unique_sumstats.iterrows():
        if index==0:
            log.write(f"Getting associations from GWAS Catalog for {row[rsid]}...", end="",verbose=verbose)
        elif index < len(unique_sumstats)-1:
            log.write(f",{row[rsid]}...",end="", show_time=False, verbose=verbose)
        else:
            log.write(f",{row[rsid]}...", show_time=False, verbose=verbose)

        df = get_associations(variant_id = row[rsid])
        
        if len(df.associations)>0:
            df.associations[rsid] = row[rsid]
            association = pd.concat([association, df.associations],ignore_index=True)
        
            df.strongest_risk_alleles[rsid] = row[rsid]
            strongest_risk_alleles = pd.concat([strongest_risk_alleles, df.strongest_risk_alleles],ignore_index=True)
    
    if len(strongest_risk_alleles)>0:
        strongest_risk_alleles["RA"] = strongest_risk_alleles["riskAlleleName"].str.split("-").str[-1]
    
    if len(association)>0:
        association = pd.merge(association, strongest_risk_alleles[["associationId","RA"]],on="associationId",how="left")

    n_assoc = len( association.drop_duplicates(subset=["associationId"]))
    traits = pd.DataFrame()
    for index,row in association.drop_duplicates(subset=["associationId"]).iterrows():
        if index==0:
            log.write(f"Getting traits from GWAS Catalog for {row["associationId"]}...", end="",verbose=verbose)
        elif index < n_assoc-1:
            log.write(f",{row["associationId"]}...",end="", show_time=False, verbose=verbose)
        else:
            log.write(f",{row["associationId"]}...", show_time=False, verbose=verbose)

        df = get_traits(association_id = row["associationId"])
        df.efo_traits["associationId"] = row["associationId"]
        traits = pd.concat([traits, df.efo_traits],ignore_index=True)
    
    studies = pd.DataFrame()
    for index,row in association.drop_duplicates(subset=["associationId"]).iterrows():
        if index==0:
            log.write(f"Getting studies from GWAS Catalog for {row["associationId"]}...", end="",verbose=verbose)
        elif index < n_assoc-1:
            log.write(f",{row["associationId"]}...",end="", show_time=False, verbose=verbose)
        else:
            log.write(f",{row["associationId"]}...", show_time=False, verbose=verbose)

        df = get_studies(association_id = row["associationId"])
        df.studies["associationId"] = row["associationId"]
        
        studies = pd.concat([studies, df.studies],ignore_index=True)

    return association, traits, studies

def _fix_beta(association):
    

    is_or_available = (association["betaNum"].isna()) & (~association["orPerCopyNum"].isna())
    is_range_available = (association["betaNum"].isna()) & (association["orPerCopyNum"].isna()) & (~association["range"].isna())

    association.loc[is_or_available ,"betaNum"] = np.log(association.loc[is_or_available,"orPerCopyNum"]) 
    association.loc[is_range_available ,"betaNum"] = association.loc[is_range_available,"range"].apply(lambda x: parse_range(x))
    return association

def parse_range(x):
    range_list = x.strip("[|]").split("-")
    high = np.log(range_list[1])
    low = np.log(range_list[0])
    beta = (high + low)/2
    return beta