from typing import TYPE_CHECKING, Dict, Any, Optional
import pandas as pd
import time
import numpy as np

if TYPE_CHECKING:
    from gwaslab.g_Sumstats import Sumstats

status_dic_12={
"13":"CHM13",
"19":"hg19",
"38":"hg38",
"97":"Unmapped",
"98":"Unknown",
"99":"Unchecked"
    }  
status_dic_3={
"0":"rsid valid & SNPID valid",
"1":"rsid valid & SNPID invalid",
"2":"rsid invalid & SNPID valid",
"3":"rsid invalid & SNPID invalid",
"4":"rsid valid & SNPID valid",
"5":"rsid valid & SNPID unknown",
"6":"rsid unknown & SNPID valid",
"7":"rsid invalid & SNPID unknown",
"8":"rsid unknown & SNPID invalid",
"9":"Unchecked"
    }  
status_dic_4={
"0":"CHR valid & POS valid",
"2":"CHR invalid & POS invalid",
"3":"CHR invalid & POS valid",
"4":"CHR valid & POS invalid",
"5":"CHR valid & POS unknown",
"6":"CHR unknown & POS valid",
"7":"CHR invalid & POS unknown",
"8":"CHR unknown & POS invalid",
"9":"Unchecked"
    }  
status_dic_5={
"0":"standardized SNP",
"1":"standardized & normalized insertion",
"2":"standardized & normalized deletion",
"3":"standardized & normalized indel",
"4":"standardized indel",
"5":"indistinguishable or not normalized allele",
"6":"invalid allele notation",
"7":"Unknown",
"9":"Unchecked"
    }  
status_dic_6={
"0":"Match: NEA=REF",
"1":"Flipped_fixed",
"2":"Reverse_complementary_fixed",
"3":"Flipped",
"4":"Reverse_complementary",
"5":"Reverse_complementary+Flipped",
"6":"Both_alleles_on_ref+indistinguishable",
"8":"Not_on_reference_genome",
"9":"Unchecked"
    }  
status_dic_7={
"0":"Not_palindromic_SNPs",
"1":"Palindromic+strand",
"2":"Palindromic-strand_fixed",
"3":"Indel_match",
"4":"Indel_flipped_fixed",
"5":"Palindromic-strand",
"6":"Indel_flipped",
"7":"Indistinguishable",
"8":"No_matching_or_no_info",
"9":"Unchecked"
    }   

def summarize(
    self: 'Sumstats',
    chrom: str = "CHR",
    snpid: str = "SNPID",
    rsid: str = "rsID",
    eaf: str = "EAF",
    p: str = "P",
    beta: str = "BETA",
    n: str = "N",
    status: str = "STATUS",
) -> Dict[str, Any]:
    """
    Generate a structured quality-control summary for GWAS summary statistics.

    This function computes descriptive metrics across several categories:
        - META: Basic dataset information
        - CHR: Chromosome counts and distribution
        - MISSING: Missing value counts across relevant fields
        - MAF: Minor allele frequency distribution
        - P: P-value significance levels
        - BETA: Effect-size magnitude among genome-wide significant variants
        - STATUS: Variant QC status summary (see `lookup_status`)

    Returns
    -------
    dict
        A hierarchical dictionary stored under `self.meta["gwaslab"]["summary"]`,
        organized by categories listed above.
    """

    # -------------------------------------------------------------------------
    # Prepare data and initialize
    # -------------------------------------------------------------------------
    insumstats = self.data
    self.meta["gwaslab"]["summary"] = {}
    output = {}

    selected_cols = [c for c in [snpid, rsid, eaf, p, beta, n, status] if c in insumstats.columns]
    sumstats = insumstats[selected_cols].copy()
    total_variants = len(sumstats)

    # -------------------------------------------------------------------------
    # META information
    # -------------------------------------------------------------------------
    meta_info = {
        "Row_num": str(insumstats.shape[0]),
        "Column_num": str(insumstats.shape[1]),
        "Column_names": ",".join(insumstats.columns.astype(str)),
        "Last_checked_time": str(time.ctime(time.time())),
        "QC and Harmonization": self.check_sumstats_qc_status()
    }

    self.meta["gwaslab"]["summary"]["overview"] = meta_info
    output["META"] = meta_info

    # -------------------------------------------------------------------------
    # Chromosome distribution
    # -------------------------------------------------------------------------
    chr_info = {
        "Chromosomes_notations": sorted(insumstats[chrom].unique()),
        "Chromosomes_numbers": len(insumstats[chrom].unique()),
    }

    chr_counts = insumstats.groupby(chrom)[chrom].count().to_dict()
    chr_info.update({f"chr{key}": val for key, val in chr_counts.items()})

    self.meta["gwaslab"]["summary"]["chromosomes"] = chr_info
    output["CHR"] = chr_info

    # -------------------------------------------------------------------------
    # Missing values
    # -------------------------------------------------------------------------
    missing_info = {}
    missing_mask = sumstats.isna()

    if missing_mask.any(axis=None):
        missing_info["Missing_total"] = missing_mask.any(axis=1).sum()
        for col in missing_mask.columns:
            count = missing_mask[col].sum()
            if count > 0:
                missing_info[f"Missing_{col}"] = count
    else:
        missing_info["Missing_total"] = 0

    self.meta["gwaslab"]["summary"]["missing_values"] = missing_info
    output["MISSING"] = missing_info

    # -------------------------------------------------------------------------
    # MAF categories
    # -------------------------------------------------------------------------
    if eaf in sumstats.columns:
        maf = pd.to_numeric(sumstats[eaf], errors="coerce")
        maf = maf.where(maf <= 0.5, 1 - maf)

        maf_info = {
            "Common (MAF>=0.05)": (maf >= 0.05).sum(),
            "Low_frequency (0.01<MAF<=0.05)": ((maf > 0.01) & (maf <= 0.05)).sum(),
            "Rare (0.001<MAF<=0.01)": ((maf > 0.001) & (maf <= 0.01)).sum(),
            "Ultra Rare (MAF<=0.001)": (maf <= 0.001).sum(),
        }

        self.meta["gwaslab"]["summary"]["MAF"] = maf_info
        output["MAF"] = maf_info

    # -------------------------------------------------------------------------
    # P-value significance thresholds
    # -------------------------------------------------------------------------
    if p in sumstats.columns and len(sumstats) > 0:
        p_values = sumstats[p].dropna()
        if len(p_values) > 0:
            p_info = {
                "Minimum": float(p_values.min()),
                "P<5e-8": (p_values < 5e-8).sum(),
                "P<5e-6": (p_values < 5e-6).sum(),
            }
            self.meta["gwaslab"]["summary"]["p_values"] = p_info
            output["P"] = p_info

    # Beta magnitude among significant variants
    if p in sumstats.columns and beta in sumstats.columns and len(sumstats) > 0:
        sig = sumstats.loc[sumstats[p] < 5e-8, beta].abs()
        if len(sig) > 0:
            beta_info = {
                f"P<5e-8 with ABS(BETA)>{thr}": (sig > thr).sum()
                for thr in [10, 3, 1, 0.5, 0.3, 0.2, 0.1, 0.05]
            }

            self.meta["gwaslab"]["summary"]["beta_for_significant_variants"] = beta_info
            output["BETA"] = beta_info

    # -------------------------------------------------------------------------
    # STATUS categories
    # -------------------------------------------------------------------------
    if status in sumstats.columns:
        tmp_id = "__tmp_idx__"
        sumstats[tmp_id] = np.arange(total_variants)

        status_summary = sum_status(tmp_id, sumstats[[tmp_id, status]])
        sumstats.drop(columns=tmp_id, inplace=True)

        status_info = {str(idx): {"count": row.iloc[0], "explanation":_explain_status(str(idx))} for idx, row in status_summary.iterrows()}

        self.meta["gwaslab"]["summary"]["variant_status"] = status_info
        output["STATUS"] = status_info

    # -------------------------------------------------------------------------
    # Add percentages
    # -------------------------------------------------------------------------
    for section, stats in self.meta["gwaslab"]["summary"].items():
        for key, value in list(stats.items()):
            if key != "Minimum" and isinstance(value, (int, float, np.number)):
                if total_variants > 0:
                    stats[f"{key} percentage"] = value / total_variants
                else:
                    stats[f"{key} percentage"] = 0.0

    # -------------------------------------------------------------------------
    # Additional variant / sample metadata
    # -------------------------------------------------------------------------
    self.meta["gwaslab"]["summary"].setdefault("variants", {})
    self.meta["gwaslab"]["summary"]["variants"]["variant_number"] = total_variants

    if p in sumstats.columns and len(sumstats) > 0:
        p_values = sumstats[p].dropna()
        if len(p_values) > 0:
            self.meta["gwaslab"]["summary"]["variants"]["min_P"] = float(np.nanmin(p_values))

    if eaf in sumstats.columns and len(sumstats) > 0:
        eaf_values = sumstats[eaf].dropna()
        if len(eaf_values) > 0:
            self.meta["gwaslab"]["summary"]["variants"]["min_minor_allele_freq"] = float(min(
                np.nanmin(eaf_values), 1 - np.nanmax(eaf_values)
            ))

    if n in sumstats.columns and len(sumstats) > 0:
        n_values = sumstats[n].dropna()
        if len(n_values) > 0:
            self.meta["gwaslab"]["summary"].setdefault("samples", {})
            self.meta["gwaslab"]["summary"]["samples"]["sample_size"] = int(n_values.max())
            self.meta["gwaslab"]["summary"]["samples"]["sample_size_median"] = float(n_values.median())
            self.meta["gwaslab"]["summary"]["samples"]["sample_size_min"] = int(n_values.min())



    return to_python(self.meta["gwaslab"]["summary"])



from typing import Any, Union, List, Dict
import numpy as np
import pandas as pd

def to_python(obj: Any) -> Any:
    # Dict
    if isinstance(obj, dict):
        return {k: to_python(v) for k, v in obj.items()}

    # List / Tuple
    elif isinstance(obj, (list, tuple)):
        return [to_python(v) for v in obj]

    # Pandas DataFrame -> list of dict
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient="records")

    # Pandas Series -> dict
    elif isinstance(obj, pd.Series):
        return obj.to_dict()

    # NumPy array -> Python list
    elif isinstance(obj, np.ndarray):
        return obj.tolist()

    # NumPy scalar -> Python scalar
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)

    # Everything else unchanged
    else:
        return obj

def sum_status(id_to_use: str, sumstats: pd.DataFrame) -> pd.DataFrame:
        results = sumstats.groupby("STATUS",observed=True).count()
        results = results.loc[results[id_to_use]>0,:].sort_values(id_to_use,ascending=False)
        return results

def _explain_status(i: str) -> Dict[str, str]:
    title = ["Genome_Build","rsID&SNPID","CHR&POS","Stadardize&Normalize","Align","Panlidromic_SNP&Indel"]
    value = [status_dic_12[i[:2]] , status_dic_3[i[2]] , status_dic_4[i[3]], status_dic_5[i[4]], status_dic_6[i[5]], status_dic_7[i[6]]]
    explanation = {i:j for i,j in zip(title, value)}
    return explanation

def lookupstatus(status: pd.Series) -> pd.DataFrame:
    """
    Decode and analyze variant status codes.
    
    Processes status codes that encode multiple layers of variant validation information 
    using a string-based format. Returns a structured DataFrame with detailed status 
    breakdown and frequency statistics.
    
    Status Code Structure:
    Each status code string contains 7+ digits encoding:
    - 1st 2 digits: Genome build mapping (CHM13/hg19/hg38)
    - 3rd digit: rsID and SNPID validation status
    - 4th digit: Chromosome and position validation
    - 5th digit: Standardization and normalization status
    - 6th digit: Alignment to reference genome
    - 7th digit: Palindromic SNP/indel status
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with columns:
        - Genome_Build: Reference genome build information
        - rsID&SNPID: Validation status of ID fields
        - CHR&POS: Chromosome/position validation status
        - Stadardize&Normalize: Standardization status
        - Align: Reference alignment status
        - Panlidromic_SNP&Indel: Variant type characteristics
        - Count: Absolute count of each status code
        - Percentage(%): Relative percentage of total variants
    """
    uniq_status = status.unique()
 
    status_dic={}
    for i in uniq_status:
        count = sum(status==i)
        percentage =  np.around(count*100 / len(status),2)
        # Convert integer status to string for digit extraction
        i_str = str(int(i)).zfill(7)  # Ensure 7 digits
        description = [status_dic_12[i_str[:2]] , status_dic_3[i_str[2]] , status_dic_4[i_str[3]], status_dic_5[i_str[4]], status_dic_6[i_str[5]], status_dic_7[i_str[6]],count,percentage]
        status_dic[i] = description
    df = pd.DataFrame.from_dict({i: status_dic[i] 
                           for i in status_dic.keys()},
                           orient='index',columns=["Genome_Build","rsID&SNPID","CHR&POS","Stadardize&Normalize","Align","Panlidromic_SNP&Indel","Count","Percentage(%)"], dtype="string")
    df = df.sort_index()
    return df
