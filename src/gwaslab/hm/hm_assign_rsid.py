import pandas as pd
import numpy as np
from gwaslab.g_Log import Log
import subprocess
import shutil
import os
import tempfile

def _assign_rsid(
    sumstats: pd.DataFrame,
    path: str,
    ref_mode: str = "auto",
    tsv_path: str | None = None,
    reuse_lookup: bool = True,
    convert_to_bcf: bool = False,
    strip_info: bool = True,
    threads: int = 6,
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
    Assign rsID using new lookup + allele-matching + STATUS filtering.
    """
    import os
    import re

    def _looks_like_vcf(p: str) -> bool:
        p = p.lower()
        return p.endswith(".vcf") or p.endswith(".vcf.gz") or p.endswith(".bcf")

    log.write("[assign_rsid_from_ref] Start", verbose=verbose)

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
    if ref_mode == "auto":
        ref_mode = "vcf" if _looks_like_vcf(path) else "tsv"

    log.write(f"[assign] reference mode = {ref_mode}", verbose=verbose)

    if ref_mode == "tsv":
        lookup_tsv = path
        if not os.path.exists(lookup_tsv):
            raise FileNotFoundError(f"Lookup TSV not found: {lookup_tsv}")
        log.write(f"[lookup] Using TSV directly: {lookup_tsv}", verbose=verbose)

    else:  # ref_mode == "vcf"
        if tsv_path:
            tsv_exists = os.path.exists(tsv_path) or os.path.exists(tsv_path + ".gz")
        else:
            tsv_exists = False

        if reuse_lookup and tsv_exists:
            lookup_tsv = tsv_path if os.path.exists(tsv_path) else tsv_path + ".gz"
            log.write(f"[lookup] Reusing existing lookup TSV: {lookup_tsv}", verbose=verbose)
        else:
            vcf_input = path
            if convert_to_bcf:
                log.write(f"[vcf→bcf] Converting (strip_info={strip_info})", verbose=verbose)
                vcf_input = _convert_vcf_to_bcf(path, threads=threads, strip=strip_info)

            lookup_tsv = tsv_path or "gwaslab_sumstats_lookup_table.tsv.gz"
            log.write(f"[lookup] Extracting new lookup TSV from: {vcf_input}", verbose=verbose)

            lookup_tsv = _extract_lookup_table_from_vcf_bcf(
                vcf_path=vcf_input,
                sumstats=sumstats[[chrom, pos]].rename(columns={chrom: "CHR", pos: "POS"}),
                chr_dict=chr_dict,
                assign_cols=["rsID"],  # ensures the lookup carries ID/rsID usable by assigner
                out_lookup=lookup_tsv,
                threads=threads,
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
    log.write("[assign] Running _assign_from_lookup(rsID)...", verbose=verbose)

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
    )

    after_missing = sumstats[rsid].isna().sum()
    filled = before_missing - after_missing

    had_rsid_after = (~sumstats[rsid].isna()).sum()

    log.write(f"[assign] Filled {filled} rsIDs.", verbose=verbose)
    log.write(f"[total] rsID count: {had_rsid_before} → {had_rsid_after} / {total_before}", verbose=verbose)
    log.write("[assign_rsid_from_ref] Done", verbose=verbose)

    return sumstats


def _annotate_sumstats(
    sumstats: pd.DataFrame,
    vcf_path=None,
    tsv_path="gwaslab_sumstats_lookup_table.txt.gz",
    assign_cols=("rsID",),
    chr_dict=None,
    threads=6,
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    overwrite = False,
    is_id_rsid = True,
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
    if tsv_path is not None:
        tsv_exists = os.path.exists(tsv_path) or os.path.exists(tsv_path + ".gz")
    else:
        tsv_exists = False

    if overwrite==True:
        tsv_exists = False

    if reuse_lookup and tsv_exists:
        # Normalize to actual existing path
        tsv_path = tsv_path if os.path.exists(tsv_path) else tsv_path + ".gz"
        log.write(f"[info] Reusing existing lookup table: {tsv_path}")

    else:
        if vcf_path is None:
            raise ValueError(
                "vcf_path must be provided when lookup table does not exist or reuse_lookup=False."
            )

        if convert_to_bcf:
            log.write(f"[info] Converting VCF → BCF (strip_info={strip_info})")
            vcf_path = _convert_vcf_to_bcf(vcf_path, threads=threads, strip=strip_info)

        log.write(f"[info] Creating lookup table from VCF: {vcf_path}")
        tsv_path = _extract_lookup_table_from_vcf_bcf(
            vcf_path   = vcf_path,
            sumstats   = sumstats,
            chr_dict   = chr_dict,
            assign_cols     = assign_cols,
            out_lookup = tsv_path,
            threads    = threads,
            log=log,
            verbose=verbose
        )

    # ----------------------------------------------
    # Step 2 — assign annotation fields
    # ----------------------------------------------
    sumstats = _assign_from_lookup(
        sumstats     = sumstats,
        lookup_table = tsv_path,
        assign_cols  = assign_cols,
        is_id_rsid=is_id_rsid,
        chrom        = chrom,
        pos          = pos,
        ea           = ea,
        nea          = nea,
                log=log,
            verbose=verbose
    )

    return sumstats


def _extract_lookup_table_from_vcf_bcf(
    vcf_path,
    sumstats,
    chr_dict = None,
    assign_cols   = None,
    out_lookup="gwaslab_sumstats_lookup_table.txt.gz",
    threads=6,
    verbose=True,
    log=Log()
    ):

    if assign_cols is None:
        assign_cols = []
    if not shutil.which("bcftools"):
        raise RuntimeError("bcftools not found in PATH")

    tmp_targets = tempfile.NamedTemporaryFile(delete=False, suffix=".targets.txt")
    targets_path = tmp_targets.name
    tmp_targets.close()

    # ---------------- NEW: Save original CHR for later restore ----------------
    if chr_dict is not None:
        inv_chr_dict = {v: k for k, v in chr_dict.items()}  # inverse map
        log.write("Converting chromosome notation using chr_dict → reference notation...")
        sumstats_for_targets = sumstats[["CHR","POS"]].copy()
        sumstats_for_targets["CHR"] = (
            sumstats_for_targets["CHR"].astype(str).map(chr_dict)
        )
    else:
        sumstats_for_targets = sumstats[["CHR","POS"]].copy()

    # Create targets
    sumstats_for_targets.dropna(subset=["CHR"]) \
        .drop_duplicates() \
        .sort_values(["CHR","POS"]) \
        .to_csv(targets_path, sep="\t", header=False, index=False)

    log.write(f"Created target list → {targets_path}")

    # ---- extract from VCF/BCF ----
    tmp_filtered_bcf = tempfile.NamedTemporaryFile(delete=False, suffix=".filtered.bcf").name
    cmd_filter = [
        "bcftools", "view",
        "-T", targets_path,
        "-Ob", "-o", tmp_filtered_bcf,
        "--threads", str(threads),
        vcf_path
    ]
    log.write("Extracting target sites: " + " ".join(cmd_filter))
    subprocess.check_call(cmd_filter)
    subprocess.check_call(["bcftools", "index", "-f", tmp_filtered_bcf])

    info_tags = [col for col in assign_cols if col not in ("ID","rsID")]
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID"
    for tag in info_tags:
        fmt += f"\t%INFO/{tag}"
    fmt += "\n"

    id_col = "rsID" if ("rsID" in assign_cols) else "ID"
    header_cols = ["CHR", "POS", "REF", "ALT", id_col] + info_tags
    header_line = "\t".join(header_cols) + "\n"
    cmd_query = ["bcftools", "query", "-f", fmt, tmp_filtered_bcf]

    log.write(f"Writing lookup → {out_lookup}")
    if out_lookup.endswith(".gz"):
        import gzip
        out_handle = gzip.open(out_lookup, "wt")
    else:
        out_handle = open(out_lookup, "w")

    with out_handle as out_f:
        out_f.write(header_line)
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
        log.write("Converting CHR back to original sumstats notation...")
        df["CHR"] = df["CHR"].astype(str).map(inv_chr_dict).astype("category")

    df.to_csv(out_lookup, sep="\t", index=False, compression="infer")

    log.write("lookup_table created:", out_lookup)
    return out_lookup

    # ---- Cleanup ----
    #for f in [tmp_targets, tmp_filtered_bcf, tmp_filtered_bcf + ".csi", tmp_filtered_bcf + ".tbi"]:
    #    if isinstance(f, str) and os.path.exists(f):
    #        try:
    #            os.remove(f)
    #        except:
    #            pass

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
    is_id_rsid=True,
    log=Log()
):
    import pandas as pd
    import numpy as np

    chunksize = 5_000_000

    log.write(f"[START] Assigning annotation from lookup: {lookup_table}")

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

    log.write(f"[DETECT] Allele mode detected: {mode} → using {lookup_ea_col}(EA/ALT) / {lookup_nea_col}(NEA/REF)")
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
        log.write(f"[WARN] annotation columns not available in lookup → skipped: {dropped}")

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
        log.write(f"[CHUNK] Loaded {len(chunk):,} lookup rows.")

        ss_sub = sumstats[sumstats[chrom].isin(chunk[chrom].unique())]
        if ss_sub.empty:
            log.write("[SKIP] No matching chromosomes in this chunk.")
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

        log.write(f"[UPDATE] Assigned: {mask.sum():,} variants | Flipped: {flipped[mask].sum():,}")

    log.write(f"[DONE] Total assigned: {processed_variants:,}")
    log.write(f"[DONE] Total flipped alleles: {flipped_count:,}")
    log.write("[END] Annotation complete.")

    return sumstats

import subprocess
from pathlib import Path

def _convert_vcf_to_bcf(reference, threads=6, strip=True, ref_fa=None, log=Log(), verbose=True):
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
        log.write("Running:", cmd1[2])
    subprocess.check_call(cmd1)

    # ---- Index ----
    cmd2 = ["bcftools", "index", "-f", out_bcf]
    if verbose:
        log.write("Running:", " ".join(cmd2))
    subprocess.check_call(cmd2)

    if verbose:
        log.write(f"Done. Output: {out_bcf} and index.")
    return out_bcf