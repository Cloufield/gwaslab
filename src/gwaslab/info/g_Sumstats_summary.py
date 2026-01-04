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
        A hierarchical dictionary containing comprehensive quality-control summary
        statistics. The dictionary is organized into the following sections:
        
        **META** (key: "overview")
            Basic dataset metadata including:
            - "Row_num": Total number of variants (rows) in the dataset
            - "Column_num": Total number of columns in the dataset
            - "Column_names": Comma-separated list of all column names
            - "Last_checked_time": Timestamp of when the summary was generated
            - "QC and Harmonization": Status of QC and harmonization checks
        
        **CHR** (key: "chromosomes")
            Chromosome distribution information:
            - "Chromosomes_notations": Sorted list of unique chromosome identifiers
            - "Chromosomes_numbers": Count of unique chromosomes
            - "chr{chr_id}": Count of variants per chromosome (one key per chromosome)
        
        **MISSING** (key: "missing_values")
            Missing value statistics:
            - "Missing_total": Total number of variants with at least one missing value
            - "Missing_{col_name}": Count of missing values per column (only included if > 0)
            - Each count also includes a corresponding "{key} percentage" entry
        
        **MAF** (key: "MAF", only if EAF column exists)
            Minor allele frequency distribution:
            - "Common (MAF>=0.05)": Count of common variants
            - "Low_frequency (0.01<MAF<=0.05)": Count of low-frequency variants
            - "Rare (0.001<MAF<=0.01)": Count of rare variants
            - "Ultra Rare (MAF<=0.001)": Count of ultra-rare variants
            - Each count also includes a corresponding "{key} percentage" entry
        
        **P** (key: "p_values", only if P column exists)
            P-value significance statistics:
            - "Minimum": Minimum p-value in the dataset
            - "P<5e-8": Count of genome-wide significant variants (p < 5e-8)
            - "P<5e-6": Count of suggestive significant variants (p < 5e-6)
            - Each count also includes a corresponding "{key} percentage" entry
        
        **BETA** (key: "beta_for_significant_variants", only if both P and BETA columns exist)
            Effect size magnitude for genome-wide significant variants (p < 5e-8):
            - "P<5e-8 with ABS(BETA)>{threshold}": Count of significant variants with
              absolute beta exceeding threshold, for thresholds: 10, 3, 1, 0.5, 0.3, 0.2, 0.1, 0.05
            - Each count also includes a corresponding "{key} percentage" entry
        
        **STATUS** (key: "variant_status", only if STATUS column exists)
            Variant QC status summary, organized by status code:
            - Each status code (as string key) contains:
              - "count": Number of variants with this status code
              - "explanation": Dictionary explaining each digit of the status code:
                - "Genome_Build": Reference genome build (CHM13/hg19/hg38/Unmapped/Unknown/Unchecked)
                - "rsID&SNPID": Validation status of rsID and SNPID fields
                - "CHR&POS": Validation status of chromosome and position fields
                - "Stadardize&Normalize": Standardization and normalization status
                - "Align": Alignment status to reference genome
                - "Panlidromic_SNP&Indel": Palindromic SNP and indel status
            - Each count also includes a corresponding "{key} percentage" entry
        
        **variants** (key: "variants")
            Additional variant-level metadata:
            - "variant_number": Total number of variants
            - "min_P": Minimum p-value (only if P column exists and has values)
            - "min_minor_allele_freq": Minimum minor allele frequency (only if EAF column exists)
        
        **samples** (key: "samples", only if N column exists)
            Sample size statistics:
            - "sample_size": Maximum sample size
            - "sample_size_median": Median sample size
            - "sample_size_min": Minimum sample size
        
        All numeric values are converted to native Python types (int, float, str, list, dict)
        via the `to_python` function. Percentages are calculated as fractions of the total
        variant count and are included for all count-based statistics.
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
        A DataFrame with one row per unique status code, sorted by status code index.
        All columns are of string dtype. The DataFrame contains the following columns:
        
        **Index**: 
            Status codes (as strings), sorted in ascending order. Each status code
            is a 7-digit string where each digit position encodes different validation
            information.
        
        **Genome_Build** (str):
            Reference genome build mapping status. Possible values:
            - "CHM13": Mapped to CHM13 reference
            - "hg19": Mapped to hg19/GRCh37 reference
            - "hg38": Mapped to hg38/GRCh38 reference
            - "Unmapped": Variant could not be mapped to any reference
            - "Unknown": Genome build status unknown
            - "Unchecked": Genome build not yet checked
            Encoded in the 1st-2nd digits of the status code.
        
        **rsID&SNPID** (str):
            Validation status of rsID and SNPID identifier fields. Possible values:
            - "rsid valid & SNPID valid": Both identifiers are valid
            - "rsid valid & SNPID invalid": rsID valid, SNPID invalid
            - "rsid invalid & SNPID valid": rsID invalid, SNPID valid
            - "rsid invalid & SNPID invalid": Both identifiers invalid
            - "rsid valid & SNPID valid": Both valid (alternative encoding)
            - "rsid valid & SNPID unknown": rsID valid, SNPID status unknown
            - "rsid unknown & SNPID valid": rsID status unknown, SNPID valid
            - "rsid invalid & SNPID unknown": rsID invalid, SNPID status unknown
            - "rsid unknown & SNPID invalid": rsID status unknown, SNPID invalid
            - "Unchecked": Validation not yet performed
            Encoded in the 3rd digit of the status code.
        
        **CHR&POS** (str):
            Validation status of chromosome and position fields. Possible values:
            - "CHR valid & POS valid": Both chromosome and position are valid
            - "CHR invalid & POS invalid": Both chromosome and position invalid
            - "CHR invalid & POS valid": Chromosome invalid, position valid
            - "CHR valid & POS invalid": Chromosome valid, position invalid
            - "CHR valid & POS unknown": Chromosome valid, position status unknown
            - "CHR unknown & POS valid": Chromosome status unknown, position valid
            - "CHR invalid & POS unknown": Chromosome invalid, position status unknown
            - "CHR unknown & POS invalid": Chromosome status unknown, position invalid
            - "Unchecked": Validation not yet performed
            Encoded in the 4th digit of the status code.
        
        **Stadardize&Normalize** (str):
            Standardization and normalization status of variant alleles. Possible values:
            - "standardized SNP": Variant is a standardized SNP
            - "standardized & normalized insertion": Standardized insertion variant
            - "standardized & normalized deletion": Standardized deletion variant
            - "standardized & normalized indel": Standardized indel variant
            - "standardized indel": Standardized but not normalized indel
            - "indistinguishable or not normalized allele": Allele cannot be distinguished or normalized
            - "invalid allele notation": Allele notation is invalid
            - "Unknown": Standardization status unknown
            - "Unchecked": Standardization not yet performed
            Encoded in the 5th digit of the status code.
        
        **Align** (str):
            Alignment status to reference genome. Possible values:
            - "Match: NEA=REF": Non-effect allele matches reference
            - "Flipped_fixed": Alleles were flipped and fixed
            - "Reverse_complementary_fixed": Reverse complement applied and fixed
            - "Flipped": Alleles need to be flipped
            - "Reverse_complementary": Reverse complement needed
            - "Reverse_complementary+Flipped": Both reverse complement and flip needed
            - "Both_alleles_on_ref+indistinguishable": Both alleles on reference, indistinguishable
            - "Not_on_reference_genome": Variant not found on reference genome
            - "Unchecked": Alignment not yet checked
            Encoded in the 6th digit of the status code.
        
        **Panlidromic_SNP&Indel** (str):
            Palindromic SNP and indel status. Possible values:
            - "Not_palindromic_SNPs": Variant is not a palindromic SNP
            - "Palindromic+strand": Palindromic SNP on positive strand
            - "Palindromic-strand_fixed": Palindromic SNP on negative strand, fixed
            - "Indel_match": Indel matches reference
            - "Indel_flipped_fixed": Indel flipped and fixed
            - "Palindromic-strand": Palindromic SNP on negative strand
            - "Indel_flipped": Indel needs to be flipped
            - "Indistinguishable": Variant type indistinguishable
            - "No_matching_or_no_info": No matching information available
            - "Unchecked": Status not yet checked
            Encoded in the 7th digit of the status code.
        
        **Count** (str):
            Absolute count of variants with this status code. Represented as a string
            but contains numeric count values.
        
        **Percentage(%)** (str):
            Relative percentage of variants with this status code, calculated as
            (count / total_variants) * 100. Rounded to 2 decimal places. Represented
            as a string but contains numeric percentage values.
        
        The DataFrame is sorted by the status code index in ascending order. Only
        status codes that appear in the input Series are included in the output.
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
