import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
import subprocess
import shutil
import os
import tempfile
from gwaslab.io.io_vcf import auto_check_vcf_chr_dict, is_vcf_file

def _assign_rsid(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    lookup_path: str | None = None,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    n_cores: int = 6,
    rsid: str = "rsID",
    chrom: str = "CHR",
    pos: str = "POS",
    ea: str = "EA",
    nea: str = "NEA",
    overwrite: str = "empty",
    chr_dict: dict | None = None,
    log: "Log" = Log(),
    verbose: bool = True,
):
    """
    Assign rsIDs to GWAS summary statistics using reference data with allele matching and STATUS filtering.

    This function assigns rsIDs to a GWAS summary statistics DataFrame by matching variants against a reference
    VCF or TSV file. It performs allele-aware matching and applies STATUS-based filtering to determine which
    variants should be assigned rsIDs. The function handles various reference formats and allows control over
    overwrite behavior for existing rsID values.

    Parameters
    ----------
    path : str or None, optional
        Path to reference file (VCF/BCF or TSV). Overrides `tsv_path`.
    vcf_path : str or None, optional
        Path to VCF/BCF file. Overrides `path` and `tsv_path`.
    tsv_path : str or None, optional
        Path to precomputed lookup TSV file. If not provided, generated from VCF.
    reuse_lookup : bool, optional
        If True, reuse existing lookup TSV if available.
    convert_to_bcf : bool, optional
        If True, convert VCF to BCF before processing.
    strip_info : bool, optional
        If True, strip INFO fields when converting VCF to BCF.
    n_cores : int, optional
        Number of threads for bcftools operations.
    overwrite : str, optional
        Overwrite mode: "all", "invalid", or "empty". Determines which existing rsID values to overwrite.
        Default is "empty".
    chr_dict : dict or None, optional
        Dictionary mapping sumstats chromosome names to reference chromosome names.
    log : gwaslab.g_Log.Log, optional
        Log object for recording progress. Default is a new Log instance.
    verbose : bool, optional
        If True, log detailed progress messages. Default is True.

    Returns
    -------
    pd.DataFrame
        The input `sumstats` DataFrame with rsID column updated.

    Raises
    ------
    ValueError
        If required columns are missing in `sumstats` or invalid `overwrite` value is provided.
    FileNotFoundError
        If specified reference file is not found.

    Notes
    -----
    - The function first checks for required columns in `sumstats`.
    - STATUS filtering uses a regex pattern to identify variants eligible for rsID assignment.
    - Overwrite modes:
      * "all": overwrite all rsIDs for eligible variants
      * "invalid": overwrite only non-rsID formatted values (e.g., not matching "rs[0-9]+")
      * "empty": only fill missing rsID values
    - If `ref_mode` is "auto", the function determines whether to use VCF or TSV based on file extension.
    """
    import os
    import re


    log.write(" -Starting to assign rsID from reference...", verbose=verbose)

    # ---------------------------
    # Required column checks
    # ---------------------------
    for col in [chrom, pos, ea, nea, "STATUS"]:
        if col not in sumstats.columns:
            raise ValueError(f"sumstats missing required column: {col}")

    # Ensure rsID exists
    if rsid not in sumstats.columns:
        sumstats[rsid] = pd.Series(pd.NA, dtype="string")

    total_before = len(sumstats)
    had_rsid_before = (~sumstats[rsid].isna()).sum()

    # ---------------------------
    # Determine lookup TSV
    # ---------------------------

    if vcf_path is not None:
        if is_vcf_file(vcf_path):
            ref_mode = "vcf/bcf"
            path_to_use = vcf_path
    elif path is not None:
        if is_vcf_file(path):
            ref_mode = "vcf/bcf" 
            path_to_use = path
        else:
            ref_mode = "tsv"
            path_to_use = path
    else:
        ref_mode = "tsv"
        path_to_use = tsv_path

    log.write(" -Determining reference mode: {}...".format(ref_mode), verbose=verbose)

    if ref_mode == "tsv":
        # path_to_use is tsv
        if not os.path.exists(path_to_use):
            raise FileNotFoundError(f"Lookup TSV not found: {path_to_use}")
        log.write(" -Using TSV directly for lookup: {}...".format(path_to_use), verbose=verbose)

    else:  
        # path_to_use is vcf/bcf

        if convert_to_bcf:
            if path_to_use[-3]!="bcf":
                log.write(" -Converting VCF to BCF (strip_info={})...".format(strip_info), verbose=verbose)
                path_to_use = _convert_vcf_to_bcf(path, threads=n_cores, strip=strip_info)
            else:
                log.write(" -Already bcf")

        log.write(" -Extracting new lookup TSV from: {}...".format(path_to_use), verbose=verbose)

        lookup_tsv, rm_tmp_lookup = _extract_lookup_table_from_vcf_bcf(
            vcf_path=path_to_use,
            sumstats=sumstats[[chrom, pos]].rename(columns={chrom: "CHR", pos: "POS"}),
            chr_dict=chr_dict,
            assign_cols=["rsID"],  # ensures the lookup carries ID/rsID usable by assigner
            out_lookup=lookup_path,
            threads=n_cores,
            verbose=verbose,
            log=log,
        )

    # ---------------------------
    # STATUS filter (your original logic)
    # ---------------------------
    standardized_normalized = sumstats["STATUS"].str.match(
        r"\w\w\w[0][01234]\w\w", case=False, flags=0, na=False
    )
    to_assign_mask = standardized_normalized.copy()

    # ---------------------------
    # Overwrite modes (your logic)
    # ---------------------------
    if overwrite == "all":
        sumstats.loc[to_assign_mask, rsid] = pd.NA
    elif overwrite == "invalid":
        invalid_mask = ~sumstats[rsid].astype("string").str.match(r"^rs[0-9]+$", na=False)
        sumstats.loc[to_assign_mask & invalid_mask, rsid] = pd.NA
    elif overwrite == "empty":
        # only fill missing → do nothing here
        pass
    else:
        raise ValueError("overwrite must be 'all', 'invalid', or 'empty'")

    before_missing = sumstats[rsid].isna().sum()

    # ---------------------------
    # Apply new allele-aware assignment
    # ---------------------------
    log.write(" -Running _assign_from_lookup(rsID)...", verbose=verbose)

    sumstats = _assign_from_lookup(
        sumstats=sumstats,
        lookup_table=lookup_tsv,
        assign_cols=("rsID",),   # function will map ID → rsID automatically
        chrom=chrom,
        pos=pos,
        ea=ea,
        nea=nea,
        verbose=verbose,
        log=log,
        rm_tmp_lookup=rm_tmp_lookup
    )

    after_missing = sumstats[rsid].isna().sum()
    filled = before_missing - after_missing

    had_rsid_after = (~sumstats[rsid].isna()).sum()

    log.write(" -Filled {} rsIDs...".format(filled), verbose=verbose)
    log.write(" -rsID count: {} → {} / {}...".format(had_rsid_before, had_rsid_after, total_before), verbose=verbose)
    log.write(" -Finished assigning rsID from reference.", verbose=verbose)

    return sumstats


def _annotate_sumstats(
    sumstats: pd.DataFrame,
    path: str | None = None,
    vcf_path: str | None = None,
    tsv_path: str | None = None,
    lookup_path: str | None = None,
    assign_cols=("rsID",),
    chr_dict=None,
    threads=6,
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    reuse_lookup=True,
    convert_to_bcf=False,
    strip_info=True,
    verbose=True,
    log=Log()
):
    """
    Annotate GWAS summary statistics by assigning fields (e.g., rsID, AF)
    from a lookup table extracted from a VCF/BCF.

    Two modes:
      (1) If tsv_path exists and reuse_lookup=True → skip extraction.
      (2) Otherwise extract from VCF → create tsv_path → annotate.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Must contain CHR, POS, EA, NEA.
    vcf_path : str or None
        Required if lookup table needs to be generated.
    tsv_path : str
        Lookup table file (tsv or tsv.gz).
    assign_cols : tuple[str]
        Columns to assign (e.g., ("rsID","AF")).
    chr_dict : dict, optional
        Map sumstats CHR → VCF CHR naming.
    threads : int
        bcftools threads.
    chrom, pos, ea, nea : str
        Column names in sumstats.
    reuse_lookup : bool
        If True, reuse lookup if exists.

    Returns
    -------
    sumstats : pd.DataFrame
    tsv_path : str
    """
    import os

    # ----------------------------------------------
    # Step 1 — reuse or extract lookup table
    # ----------------------------------------------
    if vcf_path is not None:
        if is_vcf_file(vcf_path):
            ref_mode = "vcf/bcf"
            path_to_use = vcf_path
    elif path is not None:
        if is_vcf_file(path):
            ref_mode = "vcf/bcf" 
            path_to_use = path
        else:
            ref_mode = "tsv"
            path_to_use = path
    else:
        ref_mode = "tsv"
        path_to_use = tsv_path

    log.write(" -Determining reference mode: {}...".format(ref_mode), verbose=verbose)
    if ref_mode == "tsv":
        # path_to_use is tsv
        if not os.path.exists(path_to_use):
            raise FileNotFoundError(f"Lookup TSV not found: {path_to_use}")
        log.write(" -Using TSV directly for lookup: {}...".format(path_to_use), verbose=verbose)

    else:  
        # path_to_use is vcf/bcf
        if convert_to_bcf:
            if path_to_use[-3]!="bcf":
                log.write(" -Converting VCF to BCF (strip_info={})...".format(strip_info), verbose=verbose)
                path_to_use = _convert_vcf_to_bcf(path, threads=threads, strip=strip_info)
            else:
                log.write(" -Already bcf")

        log.write(" -Creating lookup table from VCF: {}...".format(vcf_path), verbose=verbose)
        lookup_tsv, rm_tmp_lookup = _extract_lookup_table_from_vcf_bcf(
            vcf_path   = path_to_use,
            sumstats   = sumstats,
            chr_dict   = chr_dict,
            assign_cols= assign_cols,
            out_lookup = lookup_path,
            threads    = threads,
            log        =log,
            verbose    = verbose
        )

    # ----------------------------------------------
    # Step 2 — assign annotation fields
    # ----------------------------------------------
    sumstats = _assign_from_lookup(
        sumstats     = sumstats,
        lookup_table = lookup_tsv,
        assign_cols  = assign_cols,
        chrom        = chrom,
        pos          = pos,
        ea           = ea,
        nea          = nea,
        log          = log,
        verbose      =verbose,
        rm_tmp_lookup=rm_tmp_lookup
    )

    return sumstats


def _extract_lookup_table_from_vcf_bcf(
    vcf_path,
    sumstats,
    chr_dict = None,
    assign_cols   = None,
    out_lookup=None,
    threads=6,
    verbose=True,
    rm_out_lookup=False,
    log=Log()
    ):

    if chr_dict is None:
        chr_dict = auto_check_vcf_chr_dict(vcf_path, None, verbose, log)
        log.write(" -Auto-determined chr_dict: {}...".format(chr_dict), verbose=verbose)

    if assign_cols is None:
        assign_cols = []
    if not shutil.which("bcftools"):
        raise RuntimeError("bcftools not found in PATH")

    if out_lookup is None:
        out_lookup_tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".lookup.txt")
        out_lookup = out_lookup_tmp.name


    tmp_targets = tempfile.NamedTemporaryFile(delete=False, suffix=".targets.txt")
    _targets_path = tmp_targets.name
    tmp_targets.close()

    # ---------------- NEW: Save original CHR for later restore ----------------
    if chr_dict is not None:
        inv_chr_dict = {v: k for k, v in chr_dict.items()}  # inverse map
        log.write(" -Converting chromosome notation using chr_dict to reference notation...", verbose=verbose)
        sumstats_for_targets = sumstats[["CHR","POS"]].copy()
        sumstats_for_targets["CHR"] = sumstats_for_targets["CHR"].map(chr_dict)
    else:
        sumstats_for_targets = sumstats[["CHR","POS"]].copy()

    # Create targets
    sumstats_for_targets.dropna(subset=["CHR"]) \
        .drop_duplicates() \
        .sort_values(["CHR","POS"]) \
        .to_csv(_targets_path, sep="\t", header=False, index=False)

    log.write(" -Created target list: {}...".format(_targets_path), verbose=verbose)

    # ---- extract from VCF/BCF ----
    _tmp_filtered_bcf = tempfile.NamedTemporaryFile(delete=False, suffix=".filtered.bcf").name
    cmd_filter = [
        "bcftools", "view",
        "-T", _targets_path,
        "-Ob", "-o", _tmp_filtered_bcf,
        "--threads", str(threads),
        vcf_path
    ]
    log.write(" -Extracting target sites...", verbose=verbose)
    log.write(" -Calling: {}".format(" ".join(cmd_filter)), verbose=verbose)
    
    subprocess.check_call(cmd_filter)
    subprocess.check_call(["bcftools", "index", "-f", _tmp_filtered_bcf])

    info_tags = [col for col in assign_cols if col not in ("ID","rsID")]
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID"
    for tag in info_tags:
        fmt += f"\t%INFO/{tag}"
    fmt += "\n"

    id_col = "rsID" if ("rsID" in assign_cols) else "ID"
    header_cols = ["CHR", "POS", "REF", "ALT", id_col] + info_tags
    header_line = "\t".join(header_cols) + "\n"
    cmd_query = ["bcftools", "query", "-f", fmt, _tmp_filtered_bcf]

    log.write(" -Writing lookup to {}...".format(out_lookup), verbose=verbose)
    if out_lookup.endswith(".gz"):
        import gzip
        out_handle = gzip.open(out_lookup, "wt")
    else:
        out_handle = open(out_lookup, "w")

    with out_handle as out_f:
        out_f.write(header_line)
        log.write(" -Calling: {}".format(" ".join(cmd_query)), verbose=verbose)
        p = subprocess.Popen(cmd_query, stdout=subprocess.PIPE, text=True)

        shutil.copyfileobj(p.stdout, out_f)
        p.wait()

    # ---- Post-process output ----
    df = pd.read_csv(out_lookup, sep="\t", compression="infer")

    # Convert ID → rsID if needed
    if id_col == "rsID":
        df = df.rename(columns={"ID": "rsID"})

    # ---------------- NEW: Convert CHR back using inverse dictionary ----------------
    if chr_dict is not None:
        log.write(" -Converting CHR back to original sumstats notation...", verbose=verbose)
        df["CHR"] = df["CHR"].astype(str).map(inv_chr_dict).astype("category")

    df.to_csv(out_lookup, sep="\t", index=False, compression="infer")
    
    to_clean_up = [_targets_path,
             _tmp_filtered_bcf,
             _tmp_filtered_bcf + ".csi", 
             _tmp_filtered_bcf + ".tbi"]
    
    # ---- Cleanup ----

    for f in to_clean_up:
        if isinstance(f, str) and os.path.exists(f):
            try:
                os.remove(f)
                log.write(" -Cleaninig up : {}...".format(f), verbose=verbose)
            except:
                pass

    log.write(" -Lookup table created: {}...".format(out_lookup), verbose=verbose)
    return out_lookup, rm_out_lookup



import pandas as pd
import numpy as np

def _assign_from_lookup(
    sumstats,
    lookup_table,
    assign_cols=("rsID",),
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    verbose=True,
    log=Log(),
    rm_tmp_lookup=False
):
    import pandas as pd
    import numpy as np

    chunksize = 5_000_000

    log.write(" -Starting to assign annotation from lookup: {}...".format(lookup_table), verbose=verbose)

    lookup_header = pd.read_csv(lookup_table, sep="\t", nrows=0).columns.tolist()

    # -------------------------------
    # Detect allele column mode
    # -------------------------------
    if ("REF" in lookup_header) and ("ALT" in lookup_header):
        lookup_nea_col = "REF"
        lookup_ea_col  = "ALT"
        mode = "REF_ALT"
    elif ("NEA" in lookup_header) and ("EA" in lookup_header):
        lookup_nea_col = "NEA"
        lookup_ea_col  = "EA"
        mode = "EA_NEA"
    else:
        raise ValueError(
            f"Lookup must contain either (REF,ALT) or (NEA,EA). Found: {lookup_header}"
        )

    log.write(" -Detected allele mode: {} using {}(EA/ALT) / {}(NEA/REF)...".format(mode, lookup_ea_col, lookup_nea_col), verbose=verbose)
    original_assign_cols = assign_cols
    # -------------------------------
    # Normalize ID vs rsID behavior
    # -------------------------------
    # Case 1: lookup has both → prefer rsID
    if "rsID" in lookup_header and "ID" in lookup_header:
        id_col = "rsID"
    # Case 2: lookup has only one
    elif "rsID" in lookup_header:
        id_col = "rsID"
    elif "ID" in lookup_header:
        id_col = "ID"
    else:
        id_col = None  # no rsID or ID in lookup

    # Replace requested assign_cols accordingly
    normalized_assign_cols = []
    for col in assign_cols:
        if col in ("rsID", "ID"):
            if id_col is not None:
                normalized_assign_cols.append(id_col)
        elif col in lookup_header:
            normalized_assign_cols.append(col)

    # remove duplicates but keep order
    assign_cols = tuple(dict.fromkeys(normalized_assign_cols))

    # Report dropped annotation fields
    dropped = set(original_assign_cols) - set(assign_cols)

    if dropped:
        log.warning("Annotation columns not available in lookup, skipped: {}...".format(dropped), verbose=verbose)

    # Initialize missing annotation columns in sumstats
    for col in assign_cols:
        if col not in sumstats.columns:
            sumstats[col] = pd.NA

    if "ALLELE_FLIPPED" not in sumstats.columns:
        sumstats["ALLELE_FLIPPED"] = False

    # Required columns
    usecols = [chrom, pos, lookup_ea_col, lookup_nea_col] + list(assign_cols)

    dtype = {
        chrom: "int64",
        pos:   "int64",
        lookup_ea_col:  "category",
        lookup_nea_col: "category"
    }

    processed_variants = 0
    flipped_count = 0

    # -------------------------------
    # Stream lookup data
    # -------------------------------
    for chunk in pd.read_csv(
        lookup_table,
        sep="\t",
        usecols=usecols,
        dtype=dtype,
        chunksize=chunksize,
    ):
        log.write(" -Loaded {:,} lookup rows...".format(len(chunk)), verbose=verbose)

        ss_sub = sumstats[sumstats[chrom].isin(chunk[chrom].unique())]
        if ss_sub.empty:
            log.write(" -No matching chromosomes in this chunk...", verbose=verbose)
            continue

        ss_sub = ss_sub.copy()

        # unify alleles
        alleles = pd.CategoricalDtype(
            categories=list(
                set(ss_sub[ea].unique())
                | set(ss_sub[nea].unique())
                | set(chunk[lookup_ea_col].unique())
                | set(chunk[lookup_nea_col].unique())
            ),
            ordered=False
        )

        ss_sub[ea]  = ss_sub[ea].astype(alleles)
        ss_sub[nea] = ss_sub[nea].astype(alleles)
        chunk[lookup_ea_col]  = chunk[lookup_ea_col].astype(alleles)
        chunk[lookup_nea_col] = chunk[lookup_nea_col].astype(alleles)

        lookup = (
            chunk.drop_duplicates([chrom, pos, lookup_nea_col, lookup_ea_col])
                 .set_index([chrom, pos, lookup_nea_col, lookup_ea_col])[list(assign_cols)]
        )

        key_fwd = pd.MultiIndex.from_arrays([ss_sub[chrom], ss_sub[pos], ss_sub[nea], ss_sub[ea]])
        key_rev = pd.MultiIndex.from_arrays([ss_sub[chrom], ss_sub[pos], ss_sub[ea], ss_sub[nea]])

        vals_fwd = lookup.reindex(key_fwd).to_numpy()
        vals_rev = lookup.reindex(key_rev).to_numpy()

        flipped = pd.isna(vals_fwd).all(axis=1) & pd.notna(vals_rev).any(axis=1)
        assigned = np.where(pd.notna(vals_fwd), vals_fwd, vals_rev)

        mask = sumstats.loc[ss_sub.index, assign_cols].isna().any(axis=1).to_numpy()
        sumstats.loc[ss_sub.index[mask], assign_cols] = assigned[mask]
        sumstats.loc[ss_sub.index[mask], "ALLELE_FLIPPED"] |= flipped[mask]

        processed_variants += mask.sum()
        flipped_count += flipped[mask].sum()

        log.write(" -Assigned: {:,} variants | Flipped: {:,}...".format(mask.sum(), flipped[mask].sum()), verbose=verbose)

    log.write(" -Total assigned: {:,}...".format(processed_variants), verbose=verbose)
    log.write(" -Total flipped alleles: {:,}...".format(flipped_count), verbose=verbose)
    log.write(" -Finished annotation.", verbose=verbose)
    if rm_tmp_lookup is True:
        for f in [lookup_table]:
            if isinstance(f, str) and os.path.exists(f):
                try:
                    os.remove(f)
                    log.write(" -Cleaninig up : {}...".format(f), verbose=verbose)
                except:
                    pass
    return sumstats

import subprocess
from pathlib import Path

def _convert_vcf_to_bcf(reference, 
                        threads=6, 
                        strip=True, 
                        ref_fa=None, 
                        log=Log(), 
                        verbose=True):
    """
    Normalize a reference VCF (multi-allelic splitting) with optional INFO/FORMAT stripping
    and optional left-normalization using a reference FASTA.

    Parameters
    ----------
    reference : str or Path
        Path to the reference VCF (bgzipped). If missing '.gz', it will be appended.
    threads : int
        Number of threads for bcftools processing.
    strip : bool, default=True
        If True, remove all INFO and FORMAT fields and name output as <reference>.strip.bcf.
        If False, keep INFO/FORMAT and name output <reference>.bcf.
    ref_fa : str or Path, optional
        Reference FASTA for left-normalization (`bcftools norm -f ref.fa`). If None, skip.
        FASTA must be indexed (.fai).
    log : Log
    verbose : bool

    Returns
    -------
    str
        Path to the generated .bcf file.
    """

    reference = str(reference)
    if not reference.endswith(".gz"):
        reference_vcf = reference + ".gz"
    else:
        reference_vcf = reference

    # Output name
    out_bcf = reference + (".strip.bcf" if strip else ".bcf")

    # ---- Build normalization command ----
    # Always split multi-allelics
    cmd_norm = f"bcftools norm -m - {reference_vcf} -Ou"

    # Add FASTA for left normalization if provided
    if ref_fa is not None:
        cmd_norm += f" -f {ref_fa}"

    # ---- Strip or keep INFO/FORMAT ----
    if strip:
        cmd_annot = f" | bcftools annotate -x INFO,FORMAT -Ob --threads {threads} -o {out_bcf}"
    else:
        cmd_annot = f" | bcftools view -Ob --threads {threads} -o {out_bcf}"

    cmd1 = ["bash", "-c", cmd_norm + cmd_annot]
    if verbose:
        log.write(" -Running: {}...".format(cmd1[2]), verbose=verbose)
    subprocess.check_call(cmd1)

    # ---- Index ----
    cmd2 = ["bcftools", "index", "-f", out_bcf]
    if verbose:
        log.write(" -Running: {}...".format(" ".join(cmd2)), verbose=verbose)
    subprocess.check_call(cmd2)

    if verbose:
        log.write(" -Done. Output: {} and index...".format(out_bcf), verbose=verbose)
    return out_bcf
