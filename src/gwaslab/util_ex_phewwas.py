import pandas as pd
from gwaslab.g_Log import Log
import pandas as pd
import numpy as np

def _extract_associations(sumstats, rsid="rsID", log = Log(), verbose=True):

    assoc, traits, studies, variants = get_associations_from_gwascatalog(sumstats, rsid=rsid, log=log, verbose=verbose)
    
    assoc = _fix_beta(assoc)

    traits_agg = traits.groupby("associationId")[["trait","shortForm"]].agg(lambda x: ",".join(x)).reset_index()

    assoc_traits_agg= pd.merge(assoc, traits_agg, on ="associationId",how="left")
    
    assoc_traits_agg= pd.merge(assoc_traits_agg, studies, on ="associationId", how="left")

    assoc_traits_agg= pd.merge(assoc_traits_agg, variants, on ="associationId",how="left")
    
    assoc_traits_agg = assoc_traits_agg.rename(columns={"trait":"GWASCATALOG_TRAIT",
                                                        "riskFrequency":"RAF",
                                                        "betaNum":"Beta",
                                                        "pvalue":"P-value"
                                                        })
    
    summary_columns=['GWASCATALOG_TRAIT','associationId', 'rsID', "geneName", 
                     'RA', 'RAF','Beta', 'P-value','cohort','initialSampleSize','publicationInfo.pubmedId',
                     "functionalClass","gene.geneName"]
    
    assoc_traits_agg_summary = assoc_traits_agg[summary_columns]

    return assoc_traits_agg, assoc_traits_agg_summary

def get_associations_from_gwascatalog(sumstats, rsid="rsID", log=Log(), verbose=True):
    from pandasgwas import get_associations
    from pandasgwas import get_traits
    from pandasgwas import get_studies
    from pandasgwas import get_variants

    association = pd.DataFrame()
    strongest_risk_alleles=pd.DataFrame()
    author_reported_genes = pd.DataFrame()
    unique_sumstats = sumstats.dropna(subset=[rsid]).drop_duplicates(subset=[rsid])
    
    for index,row in unique_sumstats.iterrows():
        log.write(f"Getting associations from GWAS Catalog for {row[rsid]}...",verbose=verbose)

        df = get_associations(variant_id = row[rsid])
        
        empty=[]
        if len(df.associations)>0:
            df.associations[rsid] = row[rsid]
            association = pd.concat([association, df.associations],ignore_index=True)
        
            df.strongest_risk_alleles[rsid] = row[rsid]
            strongest_risk_alleles = pd.concat([strongest_risk_alleles, df.strongest_risk_alleles],ignore_index=True)
            
            try:
                author_reported_genes = pd.concat([author_reported_genes, df.author_reported_genes],ignore_index=True)
            except:
                pass
            log.write("", show_time=False, verbose=verbose)
        else:
            empty.append(row[rsid])

    log.write(f"No associations: {empty}", verbose=verbose)

    if len(strongest_risk_alleles)>0:
        strongest_risk_alleles["RA"] = strongest_risk_alleles["riskAlleleName"].str.split("-").str[-1]
    
    if len(association)>0:
        association = pd.merge(association, strongest_risk_alleles[["associationId","RA"]],on="associationId",how="left")

        author_reported_genes = author_reported_genes.groupby("associationId")["geneName"].agg(lambda x: ",".join(x))
        association = pd.merge(association, author_reported_genes,on="associationId",how="left")

    log.write(f"Retrieved {len(association)} associations from GWAS Catalog...", verbose=verbose)

    traits = pd.DataFrame()
    studies = pd.DataFrame()
    variants = pd.DataFrame()

    for index,row in association.drop_duplicates(subset=["associationId"]).iterrows():
        log.write(f"Getting traits/studies/variants from GWAS Catalog for associationId: {row["associationId"]}...",verbose=verbose)

        df = get_traits(association_id = row["associationId"])
        df.efo_traits["associationId"] = row["associationId"]
        traits = pd.concat([traits, df.efo_traits],ignore_index=True)
        
        df = get_studies(association_id = row["associationId"])
        df.studies["associationId"] = row["associationId"]
        studies = pd.concat([studies, df.studies],ignore_index=True)

        df = get_variants(association_id = row["associationId"])
        df.variants["associationId"] = row["associationId"]
        min_distance = df.genomic_contexts["distance"].min()
        df.genomic_contexts = df.genomic_contexts.loc[df.genomic_contexts["distance"]==min_distance,:].drop_duplicates("gene.geneName").groupby("rsId")["gene.geneName"].agg(lambda x: ",".join(x))
        df.variants = pd.merge(df.variants[["rsId","functionalClass","associationId"]],df.genomic_contexts, on="rsId")
        variants = pd.concat([variants, df.variants[["associationId","functionalClass","gene.geneName"]]],ignore_index=True)

    return association, traits, studies, variants

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