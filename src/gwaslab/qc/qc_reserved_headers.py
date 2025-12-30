researved_header="""

Here is a list of all headers (gl.Sumstats.data) used in GWASLab functions. These headers are fixed to avoid ambiguity within GWASLab framework.

| Header              | Keyword of args in preformat | Datatype   | Description                                    | Note                                         |
|---------------------|------------------------------|------------|------------------------------------------------|----------------------------------------------|
| `SNPID`             | snpid                        | `string`   | variant ID (CHR:POS:NEA:EA)                    | GWASLab assumes EA=ALT and NEA=REF           |
| `rsID`              | rsid                         | `string`   | dbSNP rsID                                     | -                                            |
| `CHR`               | chrom                        | `Int64`    | chromosome number (X 23, Y 24, MT 25)          | -                                            |
| `POS`               | pos                          | `Int64`    | base pair position                             | -                                            |
| `EA`                | ea                           | `category` | effect allele                                  | -                                            |
| `NEA`               | nea                          | `category` | non-effect allele                              | -                                            |
| `STATUS`            | status                       | `category` | variant standardization & harmonization status | lookup_status() to check                                            |
| `REF`               | ref                          | `category` | reference allele in reference genome           | -                                            |
| `ALT`               | alt                          | `category` | alternative allele                             | -                                            |
| `EAF`               | eaf                          | `float64`  | effect allele frequency                        | -                                            |
| `NEAF`              | neaf                         | `float64`  | non-effect allele frequency                    | -                                            |
| `MAF`               | maf                          | `float64`  | minor allele frequency                         | if EAF, use fill_data to get MAF                                            |
| `INFO`              | info                         | `float32`  | imputation INFO/RSQ                            | -                                            |
| `BETA`              | beta                         | `float64`  | effect size beta                               | -                                            |
| `SE`                | se                           | `float64`  | standard error of beta                         | -                                            |
| `BETA_95U`          | beta_95U                     | `float64`  | upper bound of beta 95% condidence interval    | -                                            |
| `BETA_95L`          | beta_95L                     | `float64`  | lower bound of beta 95% condidence interval    | -                                            |
| `OR`                | OR                           | `float64`  | odds ratio                                     | -                                            |
| `OR_95U`            | OR_95U                       | `float64`  | upper bound of OR 95% condidence interval      | -                                            |
| `OR_95L`            | OR_95L                       | `float64`  | lower bound of OR 95% condidence interval      | -                                            |
| `HR`                | HR                           | `float64`  | hazard ratio                                   | -                                            |
| `HR_95U`            | HR_95U                       | `float64`  | upper bound of HR 95% condidence interval      | -                                            |
| `HR_95L`            | HR_95L                       | `float64`  | lower bound of HR 95% condidence interval      | -                                            |
| `CHISQ`             | chisq                        | `float64`  | chi square                                     | -                                            |
| `Z`                 | z                            | `float64`  | z score                                        | -                                            |
| `T`                 | t                            | `float64`  | t statistics                                   | -                                            |
| `F`                 | f                            | `float64`  | F statistics                                   | -                                            |
| `P`                 | p                            | `float64`  | P value                                        | -                                            |
| `P_MANTISSA`        | p_mantissa                   | `float64`  | P mantissa                                     | -                                            |
| `P_EXPONENT`        | p_exponent                   | `float64`  | P exponent                                     | -                                            |
| `MLOG10P`           | mlog10p                      | `float64`  | $-log_{10}(P)$                                 | -                                            |
| `SNPR2`             | snpr2                        | `float64`  | per variant R2                                 | -                                            |
| `DOF`               | dof                          | `Int64`    | degree of freedom                              | -                                            |
| `P_HET`             | phet                         | `float64`  | heterogeneity test P value                     | -                                            |
| `I2_HET`            | i2                           | `float64`  | heterogeneity I2                               | -                                            |
| `I2`                | i2                           | `float64`  | heterogeneity I2 (alias for I2_HET)            | -                                            |
| `DIRECTION`         | direction                    | `string`   | direction of effect (+, -, ?)                  | -                                            |
| `DENSITY`           | density                      | `Int64`    | signal density                                 | -                                            |
| `N`                 | n                            | `Int64`    | total sample size                              | -                                            |
| `N_CASE`            | ncase                        | `Int64`    | number of cases                                | -                                            |
| `N_CONTROL`         | ncontrol                     | `Int64`    | number of controls                             | -                                            |
| `GENENAME`          |                       | `string`   | nearest gene symbol                            | -                                            |
| `CIS/TRANS`         |                      | `string`   | whether the variant is in cis or trans region  | `Cis`,`Trans`,`NoReference`                  |
| `DISTANCE_TO_KNOWN` |              | `Int64`    | distance to nearest known variants             | -                                            |
| `LOCATION_OF_KNOWN` |              | `string`   | relative location to nearest known variants    | `Same`,`Upstream`,`Downstream`,`NoReference` |
| `KNOWN_ID`          |                       | `string`   | nearest known variant ID                       | -                                            |
| `KNOWN_PUBMED_ID`   |                | `string`   | pubmed ID of the known variant                 | -                                            |
| `KNOWN_AUTHOR`      |                   | `string`   | author of the study                            | -                                            |
| `KNOWN_SET_VARIANT` |              | `string`   | known set and overlapping variant              | -                                            |
| `KNOWN_VARIANT`     |                  | `string`   | known variant overlapping with the variant     | -                                            |
| `KNOWN_SET`         |                      | `string`   | variant set of the known variant               | -                                            |
| `NOVEL`             |                          | `string`   | if the identified variants are novel           | -                                            |

"""

from typing import Dict, List, Optional, Tuple, Any, Union
import json
import os

"""
Reserved Headers JSON File Structure (qc_researved_header_json)

This JSON file serves as the single source of truth for all reserved headers in GWASLab.
Each header entry contains the following fields:

Fields:
- Keyword (str): The keyword used in preformat arguments to specify this column
- Datatype (str): The preferred/standard datatype for this column (e.g., "string", "Int64", "float64", "category")
- Category (str): Classification of the header into one of three categories:
    - "info": Core variant identification information (SNPID, rsID, CHR, POS, EA, NEA, STATUS)
    - "stats": Statistical and effect size columns (BETA, SE, P, MLOG10P, N, OR, HR, etc.)
    - "others": Additional metadata and annotation columns (REF, ALT, GENENAME, known variant info, etc.)
- ValidDtypes (list): List of acceptable data types for validation (e.g., ["string", "object"], ["Int64", "int64", "int32"])
- ValidDtypesPolars (list): List of Polars type identifiers for validation (e.g., ["String"], ["Int64"], ["Float64"])
- ValidRange (list or null): Valid value range for sanity checking [min, max]. 
    - For numeric values: [min, max] where min and max are numbers or "Inf" for infinity
    - For headers without range checking: null
    - Used in _sanity_check_stats() function for range validation
- Description (str): Human-readable description of what the column represents
- Note (str): Additional notes or important information about the column

Example entries:
{
  "SNPID": {
    "Keyword": "snpid",
    "Datatype": "string",
    "Category": "info",
    "ValidDtypes": ["string", "object"],
    "ValidDtypesPolars": ["String"],
    "ValidRange": null,
    "Description": "variant ID (CHR:POS:NEA:EA)",
    "Note": "GWASLab assumes EA=ALT and NEA=REF"
  },
  "BETA": {
    "Keyword": "beta",
    "Datatype": "float64",
    "Category": "stats",
    "ValidDtypes": ["float64"],
    "ValidDtypesPolars": ["Float64"],
    "ValidRange": [-100, 100],
    "Description": "effect size beta",
    "Note": "-"
  }
}

This JSON file is used to:
- Generate dtype_dict for data type validation (pandas)
- Generate dtype_dict_polars for Polars data type validation
- Generate dtype_dic and description_dic dictionaries
- Categorize headers into info/stats/others via _get_headers()
- Define DEFAULT_COLUMN_ORDER
- Provide ValidRange for sanity checking in _sanity_check_stats()
- Build the reserved_header markdown documentation
"""

# Load reserved headers JSON from separate file
_json_file_path = os.path.join(os.path.dirname(__file__), "qc_researved_header.json")
with open(_json_file_path, 'r') as f:
    researved_header_json = f.read()
def _build_reserved_header(j: str) -> str:
    d = json.loads(j)
    rows = []
    rows.append("Here is a list of all headers (gl.Sumstats.data) used in GWASLab functions. These headers are fixed to avoid ambiguity within GWASLab framework.")
    rows.append("")
    rows.append("| Header              | Keyword of args in preformat | Datatype   | Description                                    | Note                                         |")
    rows.append("|---------------------|------------------------------|------------|------------------------------------------------|----------------------------------------------|")
    for k, v in d.items():
        keyword = v.get("Keyword","")
        datatype = v.get("Datatype","")
        description = v.get("Description","")
        note = v.get("Note","")
        if keyword == "":
            kw = " "
        else:
            kw = keyword
        if datatype == "":
            dt = " "
        else:
            dt = f"`{datatype}`"
        rows.append(f"| `{k}` | {kw} | {dt} | {description} | {note} |")
    rows.append("")
    return "\n".join(rows)

# Extract default column order from reserved headers JSON
_reserved_headers_dict = json.loads(researved_header_json)
DEFAULT_COLUMN_ORDER = list(_reserved_headers_dict.keys())

researved_header = _build_reserved_header(researved_header_json)

# Generate dtype_dic and description_dic from researved_header_json (single source of truth)
dtype_dic = {k: v.get("Datatype", "") for k, v in _reserved_headers_dict.items()}
description_dic = {k: v.get("Description", "").strip() for k, v in _reserved_headers_dict.items()}

# Generate dtype_dict from researved_header_json (single source of truth)
dtype_dict = {k: v.get("ValidDtypes", []) for k, v in _reserved_headers_dict.items() if v.get("ValidDtypes")}

def get_default_sanity_ranges() -> Dict[str, Optional[Tuple[float, float]]]:
    """
    Get default sanity check ranges from JSON file.
    
    Returns
    -------
    dict
        Dictionary mapping header names to (min, max) tuples for range validation.
        Values are converted from JSON format (where "Inf" is a string) to Python format
        (where float("Inf") is used). Returns None for headers without range checking.
    """
    ranges = {}
    for header, info in _reserved_headers_dict.items():
        valid_range = info.get("ValidRange")
        if valid_range is None:
            ranges[header] = None
        elif isinstance(valid_range, list) and len(valid_range) == 2:
            # Convert "Inf" strings to float("Inf")
            min_val = float("Inf") if valid_range[0] == "Inf" else valid_range[0]
            max_val = float("Inf") if valid_range[1] == "Inf" else valid_range[1]
            ranges[header] = (min_val, max_val)
    return ranges

# Generate dtype_dict_polars from researved_header_json (single source of truth)
def get_dtype_dict_polars() -> Dict[str, List[Any]]:
    """
    Generate Polars dtype dictionary from JSON file.
    Converts string identifiers from ValidDtypesPolars to actual Polars type objects.
    
    Returns
    -------
    dict
        Dictionary mapping header names to lists of Polars type objects.
    """
    import polars as pl
    type_mapping = {
        "String": pl.String(),
        "Int64": pl.Int64(),
        "Float64": pl.Float64(),
    }
    
    dtype_dict_polars = {}
    for k, v in _reserved_headers_dict.items():
        polars_types = v.get("ValidDtypesPolars", [])
        if polars_types:
            dtype_dict_polars[k] = [type_mapping.get(t, pl.String()) for t in polars_types]
    return dtype_dict_polars

def _get_headers(mode: str = "all") -> List[str]:
    if mode=="info":
        return [k for k, v in _reserved_headers_dict.items() if v.get("Category") == "info"]
    elif mode=="stats":
        return [k for k, v in _reserved_headers_dict.items() if v.get("Category") == "stats"]
    else:
        return list(_reserved_headers_dict.keys())

def _check_overlap_with_reserved_keys(other):
    overlapped=[]
    for i in other:
        if i in _get_headers():
            overlapped.append(i)
    return overlapped
