"""
Format detection module for GWAS summary statistics files.

This module provides functionality to automatically detect the format of sumstats files
by comparing their headers against a formatbook of known formats. The detection algorithm
uses a robust scoring system that handles:
- Extra or missing columns
- Renamed/aliased headers
- Different delimiters (tab, comma, space)
- Incomplete headers

The main function is `detect_sumstats_format()`, which analyzes a sumstats file and returns
the most likely format along with confidence scores and alternative candidates.
"""

import csv
import gzip
import io
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from os import path

from gwaslab.bd.bd_config import options
from gwaslab.bd.bd_download import update_formatbook


def _open_text_maybe_gz(file_path: Union[str, Path], encoding: str = "utf-8") -> io.TextIOWrapper:
    """
    Open a text file, handling gzip compression if the path ends with .gz.
    
    Parameters:
    -----------
    file_path : str or Path
        Path to the file to open
    encoding : str, default "utf-8"
        Text encoding to use
        
    Returns:
    --------
    io.TextIOWrapper
        Text file object
    """
    file_path = str(file_path)
    if file_path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(file_path, "rb"), encoding=encoding, newline="")
    return open(file_path, "r", encoding=encoding, newline="")


_norm_re = re.compile(r"[^a-z0-9#]+")


def norm_header(s: str) -> str:
    """
    Normalize header string: keep # (for VCF/#CHROM), remove punctuation/underscores/spaces, lowercase.
    
    Parameters:
    -----------
    s : str
        Header string to normalize
        
    Returns:
    --------
    str
        Normalized header string
    """
    return _norm_re.sub("", s.strip().lower())


def sniff_delimiter(sample: str) -> str:
    """
    Detect the delimiter used in a text sample.
    
    Tries csv.Sniffer first, then falls back to a heuristic that picks the delimiter
    which yields the most columns in the first non-empty line.
    
    Parameters:
    -----------
    sample : str
        Sample text to analyze
        
    Returns:
    --------
    str
        Detected delimiter character
    """
    # try csv.Sniffer first; fallback to common delimiters
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", " "])
        return dialect.delimiter
    except Exception:
        # heuristic: pick the delimiter which yields most columns in the first non-empty line
        best = "\t"
        best_n = 0
        for d in [",", "\t", " "]:
            for line in sample.splitlines():
                if line.strip() and not line.startswith("##"):
                    n = len(line.rstrip("\n").split(d))
                    if n > best_n:
                        best_n, best = n, d
                    break
        return best


def read_header_and_rows(file_path: Union[str, Path], max_lines: int = 2000) -> Tuple[Optional[List[str]], List[List[str]], Dict]:
    """
    Read header and sample rows from a sumstats file.
    
    Parameters:
    -----------
    file_path : str or Path
        Path to the sumstats file
    max_lines : int, default 2000
        Maximum number of data rows to read
        
    Returns:
    --------
    tuple:
        - header: list[str] or None - Column headers if found
        - rows: list[list[str]] - Sample data rows
        - meta: dict - Metadata about the file (delimiter, is_vcf, etc.)
    """
    with _open_text_maybe_gz(file_path) as f:
        buf = []
        for _ in range(200):  # enough for sniffer
            line = f.readline()
            if not line:
                break
            buf.append(line)
        sample = "".join(buf)
        delim = sniff_delimiter(sample)

    header = None
    rows = []
    meta = {"delimiter": delim, "is_vcf": False}

    with _open_text_maybe_gz(file_path) as f:
        for line in f:
            if len(rows) >= max_lines:
                break
            s = line.rstrip("\n")
            if not s:
                continue

            # VCF meta/header
            if s.startswith("##"):
                if s.startswith("##fileformat=VCF"):
                    meta["is_vcf"] = True
                continue

            if s.startswith("#") and not s.startswith("#CHROM"):
                # comment line (non-VCF)
                continue

            parts = s.split(delim) if delim != " " else s.split()

            if s.startswith("#CHROM"):
                meta["is_vcf"] = True
                header = parts
                continue

            if header is None:
                # Decide whether this line is header or data:
                # Heuristic: if many tokens are non-numeric-ish, treat as header
                nonnum = 0
                for p in parts:
                    p2 = p.strip()
                    if not p2:
                        continue
                    if re.fullmatch(r"[-+]?\d+(\.\d+)?([eE][-+]?\d+)?", p2):
                        continue
                    nonnum += 1
                # header usually has lots of non-numeric tokens
                if nonnum >= max(2, int(0.5 * len(parts))):
                    header = parts
                    continue
                else:
                    # no header
                    header = None
                    rows.append(parts)
                    continue
            else:
                rows.append(parts)

    return header, rows, meta


def build_idf_weights(formatbook: Dict, candidate_formats: List[str]) -> Tuple[Dict[str, float], Dict[str, set]]:
    """
    Compute IDF (Inverse Document Frequency) weights for format keys.
    
    Headers that appear in fewer formats are more discriminative and get higher weights.
    
    Parameters:
    -----------
    formatbook : dict
        Format book dictionary loaded from JSON
    candidate_formats : list[str]
        List of format names to consider
        
    Returns:
    --------
    tuple:
        - idf: dict mapping normalized keys to IDF weights
        - fmt_keys: dict mapping format names to sets of normalized keys
    """
    # compute df for each normalized original header key across formats
    df = Counter()
    fmt_keys = {}
    for fmt in candidate_formats:
        keys = list(formatbook[fmt].get("format_dict", {}).keys())
        nkeys = set(norm_header(k) for k in keys)
        fmt_keys[fmt] = nkeys
        for k in nkeys:
            df[k] += 1

    N = len(candidate_formats)
    idf = {}
    for k, d in df.items():
        # smooth idf
        idf[k] = math.log((N + 1) / (d + 1)) + 1.0
    return idf, fmt_keys


def detect_sumstats_format(
    sumstats_path: Union[str, Path],
    formatbook_path: Optional[Union[str, Path]] = None,
    topk: int = 5
) -> Dict:
    """
    Detect the format of a sumstats file by comparing its headers against formatbook definitions.
    
    Algorithm Overview
    -----------------
    The detection uses a multi-stage scoring system that is robust to:
    - Extra or missing columns
    - Renamed/aliased column headers
    - Different delimiters (tab, comma, space)
    - Incomplete headers
    
    Detection Process:
    
    1. **File Parsing & Delimiter Detection**
       - Automatically detects delimiter (tab, comma, or space) using csv.Sniffer
       - Falls back to heuristic: picks delimiter that yields most columns
       - Handles gzip-compressed files
       - Identifies header line (skips comment lines, handles VCF metadata)
       - Reads up to 2000 sample rows for analysis
    
    2. **VCF Fast-Path Detection**
       - If file starts with "##fileformat=VCF" or has "#CHROM" header, immediately
         returns "vcf" format with high confidence (0.99)
    
    3. **Candidate Format Filtering**
       - Excludes all formats starting with "auto*" (e.g., "auto", "auto_detect")
       - Only considers formats that have P, MLOG10P, or LOG10P in their canonical fields
       - This filters out non-sumstats formats
    
    4. **Header Normalization**
       - Normalizes file headers: lowercase, remove punctuation/underscores/spaces
       - Keeps "#" character for VCF-style headers
       - Maps file headers to canonical fields via "auto" format dictionary
       - Creates two sets: normalized original headers and canonical field mappings
    
    5. **IDF Weighting (Inverse Document Frequency)**
       - Computes weights for each header key across all candidate formats
       - Headers that appear in fewer formats get higher weights (more discriminative)
       - Formula: IDF = log((N + 1) / (df + 1)) + 1.0, where N = number of formats, df = document frequency
       - Example: "P_BOLT_LMM" is highly discriminative (only in BOLT-LMM format) → high weight
       - Example: "P" appears in many formats → lower weight
    
    6. **Scoring for Each Candidate Format**
       For each candidate format, computes a score:
       
       a) **Direct Match Score**: Weighted sum of normalized header matches
          - Finds intersection of file's normalized headers with format's expected headers
          - Each match is weighted by its IDF value
          - More discriminative headers contribute more to the score
       
       b) **Soft Match Score**: Canonical field overlap (0.5 × number of matches)
          - Maps file headers through "auto" format → canonical fields
          - Compares canonical fields in file vs. format's expected canonical fields
          - Helps when headers are aliased (e.g., "chr" → "CHR", "position" → "POS")
       
       c) **Penalties**: Subtracted for missing required fields
          - Missing CHR+POS: -3.0 points
          - Missing P (or LOG10P/MLOG10P): -2.0 points
          - Missing BETA or OR: -1.0 points
       
       Final Score = Direct Match Score + Soft Match Score - Penalties
    
    7. **Confidence Calculation**
       - Confidence based on gap between top-1 and top-2 scores
       - Uses sigmoid function: confidence = 1 / (1 + exp(-gap))
       - If absolute score < 1.0, confidence is reduced by 60%
       - If confidence < 0.6, returns "unknown" format
    
    8. **Result Filtering**
       - Removes any "auto*" formats from results (safety check)
       - Returns top-k formats sorted by score
    
    Parameters:
    -----------
    sumstats_path : str or Path
        Path to the sumstats file to analyze
    formatbook_path : str or Path, optional
        Path to formatbook.json. If None, uses options.paths["formatbook"]
    topk : int, default 5
        Number of top candidate formats to return
        
    Returns:
    --------
    dict:
        - best_format: str - Detected format name, or "unknown" if confidence is low
        - confidence: float - Confidence score (0-1)
        - top: list[tuple[str, float]] - Top-k formats with scores
        - best_detail: dict - Detailed information about the best match
        - meta: dict - File metadata (delimiter, is_vcf, etc.)
        - notes: list[str] - Notes about the detection
    
    Examples:
    --------
    >>> result = detect_sumstats_format("sumstats.txt")
    >>> print(result["best_format"])
    'plink2'
    >>> print(result["confidence"])
    0.85
    >>> print(result["top"])
    [('plink2', 12.345), ('bolt_lmm', 8.234), ('saige', 6.123)]
    """
    # Load formatbook
    if formatbook_path is None:
        formatbook_path = options.paths["formatbook"]
    
    if not path.exists(formatbook_path):
        update_formatbook()
    
    formatbook = json.load(open(formatbook_path, "r", encoding="utf-8"))

    # Read header and sample rows
    header, rows, meta = read_header_and_rows(sumstats_path)

    # Candidate formats: skip "auto*" formats; keep those that mention P or LOG10P/MLOG10P in canonical targets
    candidates = []
    for fmt, spec in formatbook.items():
        # Exclude any format that starts with "auto"
        if fmt.startswith("auto"):
            continue
        fdict = spec.get("format_dict", {})
        canon_targets = set(fdict.values())
        if ("P" in canon_targets) or ("MLOG10P" in canon_targets) or ("LOG10P" in fdict):
            candidates.append(fmt)

    # VCF fast path
    if meta.get("is_vcf") or (header and norm_header(header[0]) == "#chrom"):
        if "vcf" in formatbook:
            return {
                "best_format": "vcf",
                "confidence": 0.99,
                "top": [("vcf", 999.0)],
                "meta": meta,
                "notes": ["Detected VCF-like header/meta (#CHROM / ##fileformat)."]
            }

    # Prepare header sets
    file_cols = header if header else []
    file_norm = set(norm_header(c) for c in file_cols)

    # auto mapping for soft canonical coverage
    auto_dict = formatbook.get("auto", {}).get("format_dict", {})
    auto_norm_map = {norm_header(k): v for k, v in auto_dict.items()}
    file_canon = set()
    for c in file_cols:
        k = norm_header(c)
        if k in auto_norm_map:
            file_canon.add(auto_norm_map[k])

    # weights for direct matches
    idf, fmt_norm_keys = build_idf_weights(formatbook, candidates)

    # Define required canonical fields for "sumstats-like" files
    required_base = {"CHR", "POS"}
    required_stats = {"P"}  # allow LOG10P/MLOG10P as substitute later

    results = []
    for fmt in candidates:
        fdict = formatbook[fmt].get("format_dict", {})
        expected_keys_norm = fmt_norm_keys[fmt]
        expected_canon = set(fdict.values())

        # Direct score: weighted overlap on original keys
        direct = file_norm.intersection(expected_keys_norm)
        direct_score = sum(idf.get(k, 1.0) for k in direct)

        # Soft score: canonical overlap (helps when headers are aliased)
        canon_overlap = file_canon.intersection(expected_canon)
        soft_score = 0.5 * len(canon_overlap)

        # Substitute for P if file has LOG10P/MLOG10P mapped
        has_p = ("P" in file_canon)
        has_logp = ("MLOG10P" in file_canon)
        p_ok = has_p or has_logp

        # Check minimum plausibility
        base_ok = required_base.issubset(file_canon) or required_base.issubset(expected_canon.intersection(file_canon))
        effect_ok = ("BETA" in file_canon) or ("OR" in file_canon) or ("BETA" in expected_canon and "BETA" in canon_overlap) or ("OR" in expected_canon and "OR" in canon_overlap)

        penalty = 0.0
        if not base_ok:
            penalty += 3.0
        if not p_ok:
            penalty += 2.0
        if not effect_ok:
            penalty += 1.0

        score = direct_score + soft_score - penalty

        results.append((fmt, score, {
            "direct_matches": sorted(direct),
            "canon_overlap": sorted(canon_overlap),
            "missing_base": sorted(required_base - file_canon),
            "has_p_or_logp": p_ok,
            "meta": formatbook[fmt].get("meta_data", {})
        }))

    results.sort(key=lambda x: x[1], reverse=True)
    # Filter out any auto* formats from results (safety check)
    results = [(fmt, score, info) for fmt, score, info in results if not fmt.startswith("auto")]
    top = results[:topk]

    # confidence: gap between top-1 and top-2, and absolute score
    if not top:
        return {"best_format": "unknown", "confidence": 0.0, "top": [], "meta": meta}

    best_fmt, best_score, best_info = top[0]
    second_score = top[1][1] if len(top) > 1 else best_score - 1e-6
    gap = best_score - second_score

    # simple confidence heuristic
    confidence = 1 / (1 + math.exp(-(gap)))  # sigmoid(gap)
    if best_score < 1.0:
        confidence *= 0.4  # low absolute evidence

    return {
        "best_format": best_fmt if confidence >= 0.6 else "unknown",
        "confidence": round(float(confidence), 3),
        "top": [(fmt, round(float(score), 3)) for fmt, score, _ in top],
        "best_detail": best_info,
        "meta": meta,
        "notes": [
            "If best_format=unknown, use format='auto' mapping and proceed with canonicalized columns."
        ]
    }
