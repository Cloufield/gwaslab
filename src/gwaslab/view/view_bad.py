import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.bd.bd_common_data import get_chr_list
from gwaslab.bd.bd_get_hapmap3 import gethapmap3
from gwaslab.qc.qc_fix_sumstats import fixchr, fixpos, fixallele

@with_logging(
        start_to_msg="view bad variants",
        finished_msg="viewing bad variants",
        start_function=".view_bad()",
        start_cols=["STATUS"]
)
def view_bad(sumstats,
             status="STATUS",
             chrom="CHR",
             pos="POS",
             ea="EA",
             nea="NEA",
             snpid="SNPID",
             rsid="rsID",
             select=("chr","pos","allele","align","palindromic"),
             auto_qc_detect=True,
             limit=100,
             log=Log(),
             verbose=True):
    """
    Return small tables of problematic variants grouped by category.

    When `STATUS` is usable, masks are derived from STATUS digits:
    - digit 4 (index 3): CHR/POS issues
    - digit 5 (index 4): allele normalization/validity issues
    - digit 6 (index 5): alignment to reference (8 = not on reference)
    - digit 7 (index 6): palindromic/indel indistinguishable/no info

    If `STATUS` is missing or all 9 and `auto_qc_detect=True`, minimal QC
    checks are run on a copy to compute masks without modifying input.

    Parameters
    - select: tuple of categories to include ("chr","pos","allele","align","palindromic")
    - auto_qc_detect: derive masks via quick QC when STATUS is not informative
    - limit: number of rows to return per category
    """
    cols = [c for c in [snpid, rsid, chrom, pos, nea, ea, status] if c in sumstats.columns]
    results = {}

    has_status = status in sumstats.columns
    status_str = sumstats[status].astype("string") if has_status else pd.Series([], dtype="string")
    initial_all_nine = False
    if has_status and len(sumstats) > 0:
        initial_all_nine = status_str.str.match(r"^..99999$", na=False).all()

    status_for_masks = status_str
    used_qc = False
    # Run minimal QC on a copy to obtain STATUS-based masks when unavailable
    if auto_qc_detect and (not has_status or initial_all_nine):
        qc_df = sumstats.copy()
        if status not in qc_df.columns:
            qc_df[status] = pd.Series(["9999999"] * len(qc_df), dtype="string")
        else:
            qc_df[status] = qc_df[status].astype("string").fillna("9999999")
        try:
            qc_df = fixchr(qc_df, chrom=chrom, status=status, verbose=False, log=log)
        except Exception:
            pass
        try:
            qc_df = fixpos(qc_df, pos=pos, status=status, verbose=False, log=log)
        except Exception:
            pass
        try:
            qc_df = fixallele(qc_df, ea=ea, nea=nea, status=status, verbose=False, log=log)
        except Exception:
            pass
        if status in qc_df.columns:
            status_for_masks = qc_df[status].astype("string")
            used_qc = True

    def _derive_chr_bad():
        s_chr = sumstats[chrom].astype("string") if chrom in sumstats.columns else pd.Series([], dtype="string")
        is_numeric = s_chr.str.isnumeric().fillna(False)
        chr_list = get_chr_list()
        invalid = (~is_numeric) | (~s_chr.isin(chr_list))
        return invalid

    def _derive_pos_bad():
        s_pos = sumstats[pos]
        s_pos = pd.to_numeric(s_pos, errors="coerce")
        lower, upper = 0, 250000000
        invalid = s_pos.isna() | (s_pos <= lower) | (s_pos >= upper)
        return invalid

    def _derive_allele_bad():
        s_ea = sumstats[ea].astype("string") if ea in sumstats.columns else pd.Series([], dtype="string")
        s_nea = sumstats[nea].astype("string") if nea in sumstats.columns else pd.Series([], dtype="string")
        bad_ea = s_ea.str.contains("[^ACTG]", na=True, regex=True) | s_ea.isna()
        bad_nea = s_nea.str.contains("[^ACTG]", na=True, regex=True) | s_nea.isna()
        same = s_ea == s_nea
        invalid = bad_ea | bad_nea | same
        return invalid

    if ("chr" in select) or ("pos" in select):
        # CHR/POS bad (digit 4: any of 2,3,4,5,6,7,8)
        if (has_status and not initial_all_nine) or used_qc:
            chrpos_bad_mask = status_for_masks.str[3].str.match(r"[2345678]", na=False)
        else:
            if ("chr" in select) and ("pos" in select):
                chrpos_bad_mask = _derive_chr_bad() | _derive_pos_bad()
            elif "chr" in select:
                chrpos_bad_mask = _derive_chr_bad()
            else:
                chrpos_bad_mask = _derive_pos_bad()
        results["CHR_POS_BAD"] = sumstats.loc[chrpos_bad_mask, cols].head(limit).copy()
        log.write(" -CHR/POS bad: {}".format(int(chrpos_bad_mask.sum() if len(chrpos_bad_mask)>0 else 0)), verbose=verbose)

    if "allele" in select:
        # Allele bad (digit 5: 5/6/7)
        if (has_status and not initial_all_nine) or used_qc:
            allele_bad_mask = status_for_masks.str[4].str.match(r"[567]", na=False)
        else:
            allele_bad_mask = _derive_allele_bad()
        results["ALLELE_BAD"] = sumstats.loc[allele_bad_mask, cols].head(limit).copy()
        log.write(" -ALLELE bad: {}".format(int(allele_bad_mask.sum() if len(allele_bad_mask)>0 else 0)), verbose=verbose)

    # Alignment masks (digit 6): 8 = not on reference; 3/4/5/6 need review
    if "align" in select and ((has_status and not initial_all_nine) or used_qc):
        align_bad_mask = status_for_masks.str[5].str.match(r"8", na=False)
        align_review_mask = status_for_masks.str[5].str.match(r"[3456]", na=False)
        results["ALIGN_BAD"] = sumstats.loc[align_bad_mask, cols].head(limit).copy()
        results["ALIGN_REVIEW"] = sumstats.loc[align_review_mask, cols].head(limit).copy()
        log.write(" -ALIGN bad: {}".format(int(align_bad_mask.sum() if len(align_bad_mask)>0 else 0)), verbose=verbose)
        log.write(" -ALIGN review: {}".format(int(align_review_mask.sum() if len(align_review_mask)>0 else 0)), verbose=verbose)

    # Palindromic/indel problematic (digit 7: 7/8)
    if "palindromic" in select and ((has_status and not initial_all_nine) or used_qc):
        pal_bad_mask = status_for_masks.str[6].str.match(r"[78]", na=False)
        results["PALINDROMIC_INDEL_BAD"] = sumstats.loc[pal_bad_mask, cols].head(limit).copy()
        log.write(" -PALINDROMIC/INDEL bad: {}".format(int(pal_bad_mask.sum() if len(pal_bad_mask)>0 else 0)), verbose=verbose)

    return results

@with_logging(
        start_to_msg="view HLA region variants",
        finished_msg="viewing HLA region variants",
        start_function=".view_hla()",
        start_cols=["CHR","POS"]
)
def view_hla(sumstats,
             chrom="CHR",
             pos="POS",
             mode="xmhc",
             build=None,
             limit=100,
             log=Log(),
             verbose=True):
    """
    Return variants within HLA/xMHC region on chr6.

    Parameters
    - mode: "xmhc" for extended MHC, "hla"/"mhc" for classical HLA bounds
    - build: optional genome build to select bounds; defaults use general ranges
    - limit: number of rows to return
    """
    lower, upper = None, None
    if build is None:
        if mode == "xmhc":
            lower, upper = 25000000, 34000000
        else:
            lower, upper = 29500000, 33500000
    else:
        if mode == "xmhc":
            lower, upper = 25000000, 34000000
        else:
            lower, upper = 29500000, 33500000
    if str(sumstats[chrom].dtype) == "string":
        mask = (sumstats[chrom].astype("string") == "6") & (sumstats[pos] > lower) & (sumstats[pos] < upper)
    else:
        mask = (sumstats[chrom] == 6) & (sumstats[pos] > lower) & (sumstats[pos] < upper)
    return sumstats.loc[mask, :].head(limit).copy()

@with_logging(
        start_to_msg="view HapMap3 variants",
        finished_msg="viewing HapMap3 variants",
        start_function=".view_hapmap3()",
        start_cols=None
)
def view_hapmap3(sumstats,
                 build=None,
                 limit=100,
                 log=Log(),
                 verbose=True):
    """
    Return HapMap3 subset of the input summary statistics.

    Parameters
    - build: optional genome build hint forwarded to `gethapmap3`
    - limit: number of rows to return
    """
    out = gethapmap3(sumstats, build=build, log=log)
    return out.head(limit).copy()

@with_logging(
        start_to_msg="view variants on specific chromosome",
        finished_msg="viewing variants on specific chromosome",
        start_function=".view_chr()",
        start_cols=["CHR"]
)
def view_chr(sumstats,
             chrom_value=1,
             chrom="CHR",
             limit=100,
             log=Log(),
             verbose=True):
    """
    Return variants on a specific chromosome.

    Parameters
    - chrom_value: target chromosome number (int)
    - limit: number of rows to return
    """
    if str(sumstats[chrom].dtype) == "string":
        mask = sumstats[chrom].astype("string") == str(chrom_value)
    else:
        mask = sumstats[chrom] == chrom_value
    return sumstats.loc[mask, :].head(limit).copy()

@with_logging(
        start_to_msg="view variants on chr1",
        finished_msg="viewing variants on chr1",
        start_function=".view_chr1()",
        start_cols=["CHR"]
)
def view_chr1(sumstats,
              chrom="CHR",
              limit=100,
              log=Log(),
              verbose=True):
    """
    Return variants on chromosome 1.
    """
    return view_chr(sumstats=sumstats, chrom_value=1, chrom=chrom, limit=limit, log=log, verbose=verbose)

@with_logging(
        start_to_msg="view bad variants and region sets",
        finished_msg="viewing bad variants and region sets",
        start_function=".view_all()",
        start_cols=None
)
def view_all(sumstats,
             status="STATUS",
             chrom="CHR",
             pos="POS",
             ea="EA",
             nea="NEA",
             snpid="SNPID",
             rsid="rsID",
             select=("chr","pos","allele","align","palindromic"),
             auto_qc_detect=True,
             include_hla=False,
             hla_mode="xmhc",
             include_hapmap3=False,
             include_chr=None,
             build=None,
             limit=100,
             log=Log(),
             verbose=True):
    """
    Return a dictionary aggregating bad-variant views and optional region sets.

    Keys may include: CHR_POS_BAD, ALLELE_BAD, ALIGN_BAD, ALIGN_REVIEW,
    PALINDROMIC_INDEL_BAD, HLA, HAPMAP3, and CHR_{n} depending on inputs.

    Parameters
    - select: categories to include for bad variants
    - auto_qc_detect: compute masks when STATUS is not informative
    - include_hla/include_hapmap3/include_chr: add region-based subsets
    - limit: number of rows per returned table
    """
    output = {}
    output.update(view_bad(sumstats=sumstats,
                           status=status,
                           chrom=chrom,
                           pos=pos,
                           ea=ea,
                           nea=nea,
                           snpid=snpid,
                           rsid=rsid,
                           select=select,
                           auto_qc_detect=auto_qc_detect,
                           limit=limit,
                           log=log,
                           verbose=verbose))
    if include_hla:
        output["HLA"] = view_hla(sumstats=sumstats, chrom=chrom, pos=pos, mode=hla_mode, build=build, limit=limit, log=log, verbose=verbose)
    if include_hapmap3:
        output["HAPMAP3"] = view_hapmap3(sumstats=sumstats, build=build, limit=limit, log=log, verbose=verbose)
    if include_chr is not None:
        output[f"CHR_{include_chr}"] = view_chr(sumstats=sumstats, chrom_value=include_chr, chrom=chrom, limit=limit, log=log, verbose=verbose)
    return output

@with_logging(
        start_to_msg="view variants by simple key",
        finished_msg="viewing variants by simple key",
        start_function=".view()",
        start_cols=None
)
def view(sumstats,
         key,
         status="STATUS",
         chrom="CHR",
         pos="POS",
         ea="EA",
         nea="NEA",
         snpid="SNPID",
         rsid="rsID",
         limit=100,
         log=Log(),
         verbose=True):
    """
    Convenience wrapper to view by simple key returning a DataFrame.

    Supported keys:
    - "CHR": returns CHR/POS-related bad variants (auto QC masks if needed)
    - "CHR_n": returns variants on chromosome n
    - "HLA": returns HLA/xMHC subset
    - "HAPMAP3": returns HapMap3 subset
    - "BAD"/"STATUS": returns concatenated bad-variant categories
    - Combined tokens like "BAD CHR" or "BAD_CHR" return concatenated views
      with a column `VIEW` indicating source category.
    """
    k = str(key).upper()
    underscore_join = "_" in k
    import re as _re
    tokens = [t for t in _re.split(r"[\s_]+", k) if t]
    if len(tokens) >= 1:
        parts = []
        has_bad = ("BAD" in tokens) or ("STATUS" in tokens)
        has_chr_token = ("CHR" in tokens)
        chr_specific = None
        for t in tokens:
            if t.startswith("CHR_"):
                try:
                    chr_specific = int(t.split("_")[1])
                except Exception:
                    chr_specific = None

        if has_bad:
            sel = ("chr","pos") if (has_chr_token or underscore_join) else ("chr","pos","allele","align","palindromic")
            bad_dict = view_bad(sumstats=sumstats,
                                status=status,
                                chrom=chrom,
                                pos=pos,
                                ea=ea,
                                nea=nea,
                                snpid=snpid,
                                rsid=rsid,
                                select=sel,
                                auto_qc_detect=True,
                                limit=limit,
                                log=log,
                                verbose=verbose)
            for kname, df in bad_dict.items():
                if isinstance(df, pd.DataFrame) and len(df)>0:
                    df = df.copy()
                    df["VIEW"] = kname
                    parts.append(df)

        if chr_specific is not None:
            df_chr = view_chr(sumstats=sumstats, chrom_value=chr_specific, chrom=chrom, limit=limit, log=log, verbose=verbose)
            if isinstance(df_chr, pd.DataFrame) and len(df_chr)>0:
                df_chr = df_chr.copy()
                df_chr["VIEW"] = f"CHR_{chr_specific}"
                parts.append(df_chr)

        if (not has_bad) and has_chr_token and chr_specific is None:
            bad_chr_only = view_bad(sumstats=sumstats,
                                    status=status,
                                    chrom=chrom,
                                    pos=pos,
                                    ea=ea,
                                    nea=nea,
                                    snpid=snpid,
                                    rsid=rsid,
                                    select=("chr","pos"),
                                    auto_qc_detect=True,
                                    limit=limit,
                                    log=log,
                                    verbose=verbose)
            df = bad_chr_only.get("CHR_POS_BAD", sumstats.iloc[0:0,:].copy())
            if len(df)>0:
                df = df.copy()
                df["VIEW"] = "CHR_POS_BAD"
                parts.append(df)

        if "HLA" in tokens:
            df_hla = view_hla(sumstats=sumstats, chrom=chrom, pos=pos, mode="xmhc", limit=limit, log=log, verbose=verbose)
            if isinstance(df_hla, pd.DataFrame) and len(df_hla)>0:
                df_hla = df_hla.copy()
                df_hla["VIEW"] = "HLA"
                parts.append(df_hla)

        if "HAPMAP3" in tokens:
            df_hm3 = view_hapmap3(sumstats=sumstats, limit=limit, log=log, verbose=verbose)
            if isinstance(df_hm3, pd.DataFrame) and len(df_hm3)>0:
                df_hm3 = df_hm3.copy()
                df_hm3["VIEW"] = "HAPMAP3"
                parts.append(df_hm3)

        if len(parts)==0:
            return sumstats.iloc[0:0,:].copy()
        return pd.concat(parts, ignore_index=True)
    if k == "CHR" or k.startswith("CHR_"):
        if k == "CHR":
            bad_chr_only = view_bad(sumstats=sumstats,
                                    status=status,
                                    chrom=chrom,
                                    pos=pos,
                                    ea=ea,
                                    nea=nea,
                                    snpid=snpid,
                                    rsid=rsid,
                                    select=("chr","pos"),
                                    auto_qc_detect=True,
                                    limit=limit,
                                    log=log,
                                    verbose=verbose)
            df = bad_chr_only.get("CHR_POS_BAD", sumstats.iloc[0:0,:].copy())
            return df
        else:
            try:
                chrom_value = int(k.split("_")[1])
            except Exception:
                chrom_value = 1
            return view_chr(sumstats=sumstats, chrom_value=chrom_value, chrom=chrom, limit=limit, log=log, verbose=verbose)
    if k == "HLA":
        return view_hla(sumstats=sumstats, chrom=chrom, pos=pos, mode="xmhc", limit=limit, log=log, verbose=verbose)
    if k == "HAPMAP3":
        return view_hapmap3(sumstats=sumstats, limit=limit, log=log, verbose=verbose)
    if k == "BAD" or k == "STATUS":
        out_dict = view_bad(sumstats=sumstats,
                            status=status,
                            chrom=chrom,
                            pos=pos,
                            ea=ea,
                            nea=nea,
                            snpid=snpid,
                            rsid=rsid,
                            select=("chr","pos","allele","align","palindromic"),
                            auto_qc_detect=True,
                            limit=limit,
                            log=log,
                            verbose=verbose)
        parts=[]
        for kname, df in out_dict.items():
            if isinstance(df, pd.DataFrame) and len(df)>0:
                df = df.copy()
                df["VIEW"] = kname
                parts.append(df)
        if len(parts)==0:
            return sumstats.iloc[0:0,:].copy()
        return pd.concat(parts, ignore_index=True)
    return sumstats.iloc[0:0,:].copy()
