from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gwaslab.g_Sumstats import Sumstats
from gwaslab.info.g_Log import Log
from gwaslab.util.util_in_get_sig import _get_sig
from gwaslab.qc.qc_check_datatype import categorical_safe_str
from gwaslab.viz.viz_aux_save_figure import save_figure
from gwaslab.viz.viz_aux_style_options import set_plot_style


def plot_lead_overlap(
    objects: Sequence[Sumstats],
    titles: Optional[Sequence[str]] = None,
    mode: str = "auto",
    sig_level: float = 5e-8,
    windowsizekb: int = 500,
    windowsizekb_for_overlap: int = 1000,
    use_p: bool = False,
    get_lead_kwargs: Optional[Dict[str, Any]] = None,
    anno: bool = True,
    build: Optional[Union[str, Sequence[str]]] = None,
    source: str = "ensembl",
    gtf_path: Optional[str] = None,
    wc_correction: bool = False,
    show_counts: bool = True,
    show_genes: bool = True,
    max_gene_labels: int = 30,
    sort_by: str = "count",
    venn_kwargs: Optional[Dict[str, Any]] = None,
    upset_kwargs: Optional[Dict[str, Any]] = None,
    fig_kwargs: Optional[Dict[str, Any]] = None,
    save: Optional[Union[bool, str]] = False,
    save_kwargs: Optional[Dict[str, Any]] = None,
    font_kwargs: Optional[Dict[str, Any]] = None,
    legend_kwargs: Optional[Dict[str, Any]] = None,
    title: Optional[str] = None,
    log: Log = Log(),
    verbose: bool = True,
    **kwargs: Any,
) -> List[Any]:
    """Plot overlap of annotated lead loci across multiple Sumstats objects.

    Parameters
    ----------
    objects : sequence of Sumstats
        At least two GWASLab Sumstats objects to compare.
    titles : sequence of str, optional
        Display names for each study. Defaults to ``meta["gwaslab"]["study_name"]``
        or ``Sumstats_N``.
    mode : str
        ``"auto"`` (Venn for 2–3 studies, UpSet for 4+), ``"venn"``, or ``"upset"``.
    sig_level : float
        Genome-wide significance threshold for lead extraction.
    windowsizekb : int
        Sliding window (kb) for lead extraction within each study.
    windowsizekb_for_overlap : int
        Window (kb) for merging leads from different studies into shared loci.
    use_p : bool
        If True, rank leads by P instead of MLOG10P.
    get_lead_kwargs : dict, optional
        Extra kwargs forwarded to lead extraction (``_get_sig`` / ``get_lead``).
    anno : bool
        Annotate leads with nearest gene names.
    build : str or sequence of str, optional
        Genome build for annotation per object.
    source : str
        Gene annotation backend (``"ensembl"`` or ``"refseq"``).
    gtf_path : str, optional
        Custom GTF path for gene annotation.
    wc_correction : bool
        Apply Winner's Curse correction during lead extraction.
    show_counts : bool
        Show intersection counts on the plot.
    show_genes : bool
        Show gene labels on UpSet intersections.
    max_gene_labels : int
        Maximum gene labels on UpSet plot.
    sort_by : str
        UpSet row sort key (e.g. ``"count"``).
    venn_kwargs, upset_kwargs : dict, optional
        Extra kwargs for Venn or UpSet drawing.
    fig_kwargs, save_kwargs, font_kwargs, legend_kwargs : dict, optional
        Figure styling and save options.
    save : bool or str
        Save figure to default or custom path.
    title : str, optional
        Figure title.
    log : Log
        GWASLab log object.
    verbose : bool
        Print progress messages.

    Returns
    -------
    list
        ``[overlap_df, fig, log]``. ``overlap_df`` has one row per locus group with
        columns such as ``MEMBERSHIP_KEY``, ``N_STUDIES``, ``GENE``, and ``IN_<study>``.
        For UpSet plots, ``overlap_df.attrs["set_list"]`` maps set IDs to memberships.

    See Also
    --------
    Sumstats.get_lead : Single-study lead extraction.
    """
    log.write("Start to plot lead overlap across sumstats objects...", verbose=verbose)
    labels = _normalize_labels(objects, titles)
    selected_mode = _resolve_mode(mode, len(objects))

    # Reuse the central GWASLab plotting style machinery so this behaves like
    # the other top-level visualization APIs.
    style = set_plot_style(
        plot="plot_lead_overlap",
        fig_kwargs=fig_kwargs,
        save_kwargs=save_kwargs,
        font_kwargs=font_kwargs,
        legend_kwargs=legend_kwargs,
        save=save,
        verbose=verbose,
        log=log,
    )
    fig_kwargs = style.get("fig_kwargs", {})
    save_kwargs = style.get("save_kwargs", {})
    font_kwargs = style.get("font_kwargs", {})
    legend_kwargs = style.get("legend_kwargs", {})

    # Step 1: extract one lead-variant table per study using the same
    # significance/window logic as Sumstats.get_lead().
    leads = _extract_leads_for_objects(
        objects=objects,
        labels=labels,
        sig_level=sig_level,
        windowsizekb=windowsizekb,
        use_p=use_p,
        get_lead_kwargs=get_lead_kwargs,
        anno=anno,
        build=build,
        source=source,
        gtf_path=gtf_path,
        wc_correction=wc_correction,
        log=log,
        verbose=verbose,
    )
    # Step 2: determine overlap by CHR/POS proximity. Leads on the same
    # chromosome are assigned to the same locus group when the next lead falls
    # within windowsizekb_for_overlap kb of the current group's end coordinate.
    overlap_df = _cluster_lead_loci(
        leads=leads,
        labels=labels,
        windowsizekb_for_overlap=windowsizekb_for_overlap,
        show_genes=show_genes,
    )
    log.write(" -Identified {} overlapping lead-locus groups.".format(len(overlap_df)), verbose=verbose)

    if selected_mode == "venn":
        fig = _draw_venn_overlap(
            overlap_df=overlap_df,
            labels=labels,
            show_counts=show_counts,
            fig_kwargs=fig_kwargs,
            font_kwargs=font_kwargs,
            title=title,
            venn_kwargs=venn_kwargs,
        )
    else:
        # UpSet plots use compact Set1/Set2/... labels on the data side only;
        # the visual matrix intentionally has no x tick labels.
        overlap_df, set_list = _assign_upset_set_ids(
            overlap_df=overlap_df,
            labels=labels,
            sort_by=sort_by,
        )
        fig = _draw_upset_overlap(
            overlap_df=overlap_df,
            set_list=set_list,
            labels=labels,
            show_counts=show_counts,
            show_genes=show_genes,
            max_gene_labels=max_gene_labels,
            sort_by=sort_by,
            fig_kwargs=fig_kwargs,
            font_kwargs=font_kwargs,
            legend_kwargs=legend_kwargs,
            title=title,
            upset_kwargs=upset_kwargs,
        )

    save_figure(fig=fig, save=save, keyword="lead_overlap", save_kwargs=save_kwargs, log=log, verbose=verbose)
    log.write("Finished plotting lead overlap.", verbose=verbose)
    return [overlap_df, fig, log]


def _resolve_mode(mode: str, n_objects: int) -> str:
    if n_objects < 2:
        raise ValueError("Please provide at least two GWASLab Sumstats objects.")
    if mode not in {"auto", "venn", "upset"}:
        raise ValueError('Please select mode from "auto", "venn", or "upset".')
    if mode == "auto":
        return "venn" if n_objects in {2, 3} else "upset"
    if mode == "venn" and n_objects not in {2, 3}:
        raise ValueError('mode="venn" supports exactly two or three Sumstats objects.')
    return mode


def _normalize_labels(objects: Sequence[Sumstats], titles: Optional[Sequence[str]]) -> List[str]:
    if not isinstance(objects, Sequence) or isinstance(objects, (str, bytes)) or len(objects) < 2:
        raise ValueError("Please provide at least two GWASLab Sumstats objects in `objects`.")
    for index, obj in enumerate(objects):
        if not isinstance(obj, Sumstats):
            raise ValueError("Please provide GWASLab Sumstats objects. Object #{} is not a Sumstats.".format(index + 1))
    if titles is not None:
        if len(titles) != len(objects):
            raise ValueError("Please provide one title for each Sumstats object.")
        return [str(title) for title in titles]

    labels = []
    seen = set()
    for index, obj in enumerate(objects):
        label = None
        meta = getattr(obj, "meta", {})
        if isinstance(meta, dict):
            label = meta.get("gwaslab", {}).get("study_name")
        if not label:
            label = "Sumstats_{}".format(index + 1)
        label = str(label)
        if label in seen:
            label = "{}_{}".format(label, index + 1)
        seen.add(label)
        labels.append(label)
    return labels


def _get_build_for_object(build: Optional[Union[str, Sequence[str]]], obj: Sumstats, index: int) -> str:
    if isinstance(build, Sequence) and not isinstance(build, str):
        return str(build[index])
    if build is not None:
        return str(build)
    obj_build = getattr(obj, "build", None)
    if obj_build is not None:
        return str(obj_build)
    meta = getattr(obj, "meta", {})
    if isinstance(meta, dict):
        meta_build = meta.get("gwaslab", {}).get("genome_build")
        if meta_build is not None:
            return str(meta_build)
    return "19"


def _extract_leads_for_objects(
    objects: Sequence[Sumstats],
    labels: Sequence[str],
    sig_level: float,
    windowsizekb: int,
    use_p: bool,
    get_lead_kwargs: Optional[Dict[str, Any]],
    anno: bool,
    build: Optional[Union[str, Sequence[str]]],
    source: str,
    gtf_path: Optional[str],
    wc_correction: bool,
    log: Log,
    verbose: bool,
) -> pd.DataFrame:
    lead_frames = []
    # Keep get_lead-compatible option names here so callers can pass familiar
    # kwargs either directly or through get_lead_kwargs.
    base_kwargs = {
        "windowsizekb": windowsizekb,
        "sig_level": sig_level,
        "use_p": use_p,
        "anno": anno,
        "source": source,
        "gtf_path": gtf_path,
        "wc_correction": wc_correction,
    }
    if get_lead_kwargs is not None:
        base_kwargs.update(get_lead_kwargs)
    build_override = base_kwargs.pop("build", build)

    for index, obj in enumerate(objects):
        label = labels[index]
        build_for_object = _get_build_for_object(build_override, obj, index)
        lead_kwargs = dict(base_kwargs)
        log.write(" -Extracting lead variants for {}...".format(label), verbose=verbose)
        leads = _get_sig(obj, build=build_for_object, log=log, verbose=verbose, **lead_kwargs)
        if leads is None or len(leads) == 0:
            leads = pd.DataFrame(columns=list(obj.data.columns))
        else:
            leads = leads.copy()
        leads["STUDY"] = label
        leads["STUDY_INDEX"] = index
        lead_frames.append(leads)
        log.write("  -{} lead variants for {}.".format(len(leads), label), verbose=verbose)

    if not lead_frames:
        return pd.DataFrame()
    return pd.concat(lead_frames, ignore_index=True, sort=False)


def _cluster_lead_loci(
    leads: pd.DataFrame,
    labels: Sequence[str],
    windowsizekb_for_overlap: int,
    show_genes: bool,
) -> pd.DataFrame:
    membership_cols = ["IN_{}".format(label) for label in labels]
    if leads.empty:
        return pd.DataFrame(columns=["LOCUS_ID", "CHR", "START", "END", "POS", "STUDIES", "N_STUDIES", "GENE", "LEAD_SNPS", "MEMBERSHIP_KEY"] + membership_cols)

    required = {"CHR", "POS", "STUDY", "STUDY_INDEX"}
    missing = required.difference(leads.columns)
    if missing:
        raise ValueError("Lead extraction output is missing required columns: {}".format(", ".join(sorted(missing))))

    working = leads.copy()
    working["CHR"] = pd.to_numeric(working["CHR"], errors="coerce")
    working["POS"] = pd.to_numeric(working["POS"], errors="coerce")
    working = working.dropna(subset=["CHR", "POS"]).sort_values(["CHR", "POS"])
    window_bp = int(windowsizekb_for_overlap) * 1000

    rows = []
    locus_index = 0
    for chrom, chrom_df in working.groupby("CHR", sort=True):
        current_rows = []
        current_start = None
        current_end = None
        for _, row in chrom_df.iterrows():
            pos = int(row["POS"])
            if current_start is None:
                # Start the first locus group on this chromosome.
                current_rows = [row]
                current_start = pos
                current_end = pos
                continue
            if pos <= current_end + window_bp:
                # Same locus: the lead is close enough to the current group.
                # Extending current_end lets chains of nearby leads merge.
                current_rows.append(row)
                current_end = max(current_end, pos)
            else:
                # New locus: the next lead is outside the overlap window.
                rows.append(_summarize_locus(locus_index, chrom, current_start, current_end, current_rows, labels, membership_cols, show_genes))
                locus_index += 1
                current_rows = [row]
                current_start = pos
                current_end = pos
        if current_start is not None:
            rows.append(_summarize_locus(locus_index, chrom, current_start, current_end, current_rows, labels, membership_cols, show_genes))
            locus_index += 1

    return pd.DataFrame(rows)


def _summarize_locus(
    locus_index: int,
    chrom: Any,
    start: int,
    end: int,
    rows: Sequence[pd.Series],
    labels: Sequence[str],
    membership_cols: Sequence[str],
    show_genes: bool,
) -> Dict[str, Any]:
    locus_df = pd.DataFrame(rows)
    study_set = set(categorical_safe_str(locus_df["STUDY"]))
    membership = [label in study_set for label in labels]
    genes = []
    if show_genes and "GENE" in locus_df.columns:
        for value in categorical_safe_str(locus_df["GENE"].dropna()):
            if value and value != "Unknown":
                genes.extend([gene.strip() for gene in value.split(",") if gene.strip()])
    gene_label = ",".join(dict.fromkeys(genes)) if genes else pd.NA

    snp_col = "SNPID" if "SNPID" in locus_df.columns else None
    lead_snps = {}
    for label in labels:
        study_rows = locus_df.loc[locus_df["STUDY"] == label]
        if snp_col is not None and len(study_rows) > 0:
            lead_snps[label] = ",".join(categorical_safe_str(study_rows[snp_col].dropna()).tolist())
        else:
            lead_snps[label] = ""

    out = {
        "LOCUS_ID": "Locus_{}".format(locus_index + 1),
        "CHR": int(chrom) if pd.notna(chrom) and float(chrom).is_integer() else chrom,
        "START": int(start),
        "END": int(end),
        "POS": int(round(float(locus_df["POS"].median()))),
        "STUDIES": ",".join([label for label, present in zip(labels, membership) if present]),
        "N_STUDIES": int(sum(membership)),
        "GENE": gene_label,
        "LEAD_SNPS": lead_snps,
        # One bit per study in label order; e.g. 1|0|1 means present in
        # studies 1 and 3 but absent from study 2.
        "MEMBERSHIP_KEY": "|".join(["1" if present else "0" for present in membership]),
    }
    for col, present in zip(membership_cols, membership):
        out[col] = present
    return out


def _membership_sets(overlap_df: pd.DataFrame, labels: Sequence[str]) -> Dict[str, set]:
    sets = {}
    for label in labels:
        col = "IN_{}".format(label)
        if col in overlap_df.columns:
            sets[label] = set(overlap_df.loc[overlap_df[col], "LOCUS_ID"])
        else:
            sets[label] = set()
    return sets


def _assign_upset_set_ids(
    overlap_df: pd.DataFrame,
    labels: Sequence[str],
    sort_by: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Convert each unique membership pattern to Set1, Set2, ... so the UpSet
    # plot can stay compact while the mapping remains available to users.
    set_list = _build_upset_set_list(overlap_df=overlap_df, sort_by=sort_by)
    overlap_df = overlap_df.copy()
    if "MEMBERSHIP_KEY" in overlap_df.columns and not set_list.empty:
        key_to_set_id = dict(zip(set_list["MEMBERSHIP_KEY"], set_list["SET_ID"]))
        overlap_df["SET_ID"] = overlap_df["MEMBERSHIP_KEY"].map(key_to_set_id)
    else:
        overlap_df["SET_ID"] = pd.Series(dtype="string")
    overlap_df.attrs["set_list"] = set_list
    return overlap_df, set_list


def _build_upset_set_list(overlap_df: pd.DataFrame, sort_by: str) -> pd.DataFrame:
    if overlap_df.empty or "MEMBERSHIP_KEY" not in overlap_df.columns:
        return pd.DataFrame(columns=["SET_ID", "MEMBERSHIP_KEY", "STUDIES", "COUNT", "GENE", "LOCUS_IDS"])

    rows = []
    for membership_key, group in overlap_df.groupby("MEMBERSHIP_KEY", sort=False):
        # Each row in set_list summarizes all locus groups with the same
        # membership pattern.
        genes = []
        if "GENE" in group.columns:
            for value in categorical_safe_str(group["GENE"].dropna()):
                if value and value != "Unknown":
                    genes.extend([gene.strip() for gene in value.split(",") if gene.strip()])
        rows.append(
            {
                "MEMBERSHIP_KEY": membership_key,
                "STUDIES": group["STUDIES"].iloc[0] if "STUDIES" in group.columns else "",
                "COUNT": len(group),
                "GENE": ",".join(dict.fromkeys(genes)) if genes else pd.NA,
                "LOCUS_IDS": ",".join(categorical_safe_str(group["LOCUS_ID"]).tolist()) if "LOCUS_ID" in group.columns else "",
                "N_STUDIES": str(membership_key).count("1"),
            }
        )

    set_list = pd.DataFrame(rows)
    if sort_by == "degree":
        set_list = set_list.sort_values(["N_STUDIES", "COUNT"], ascending=[False, False])
    elif sort_by == "key":
        set_list = set_list.sort_values("MEMBERSHIP_KEY")
    else:
        set_list = set_list.sort_values(["COUNT", "N_STUDIES"], ascending=[False, False])
    set_list = set_list.reset_index(drop=True)
    set_list.insert(0, "SET_ID", ["Set{}".format(index + 1) for index in range(len(set_list))])
    return set_list.drop(columns=["N_STUDIES"])


def _split_membership_key(key: Any) -> List[str]:
    key = str(key)
    if "|" in key:
        return key.split("|")
    return list(key)


def _draw_venn_overlap(
    overlap_df: pd.DataFrame,
    labels: Sequence[str],
    show_counts: bool,
    fig_kwargs: Dict[str, Any],
    font_kwargs: Dict[str, Any],
    title: Optional[str],
    venn_kwargs: Optional[Dict[str, Any]],
):
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError as exc:
        raise ImportError("Please install matplotlib-venn to use mode='venn'.") from exc

    fig, ax = plt.subplots(**fig_kwargs)
    sets = _membership_sets(overlap_df, labels)
    venn_kwargs = venn_kwargs or {}
    if len(labels) == 2:
        venn2([sets[labels[0]], sets[labels[1]]], set_labels=tuple(labels), ax=ax, **venn_kwargs)
    else:
        venn3([sets[labels[0]], sets[labels[1]], sets[labels[2]]], set_labels=tuple(labels), ax=ax, **venn_kwargs)
    if not show_counts:
        for text in ax.texts:
            if text.get_text().isdigit():
                text.set_text("")
    if title is not None:
        ax.set_title(title, fontdict=font_kwargs)
    return fig


def _draw_upset_overlap(
    overlap_df: pd.DataFrame,
    set_list: pd.DataFrame,
    labels: Sequence[str],
    show_counts: bool,
    show_genes: bool,
    max_gene_labels: int,
    sort_by: str,
    fig_kwargs: Dict[str, Any],
    font_kwargs: Dict[str, Any],
    legend_kwargs: Dict[str, Any],
    title: Optional[str],
    upset_kwargs: Optional[Dict[str, Any]],
):
    del legend_kwargs
    upset_kwargs = upset_kwargs or {}
    del overlap_df, show_genes, max_gene_labels, sort_by
    combo_df = set_list.reset_index(drop=True)

    fig = plt.figure(**fig_kwargs)
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1.4], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(len(combo_df))
    color = upset_kwargs.get("bar_color", "#597FBD")
    connector_gap = upset_kwargs.get("connector_gap", 0.13)
    ax_bar.bar(x, combo_df["COUNT"] if len(combo_df) else [], color=color)
    ax_bar.set_ylabel("Lead loci", fontdict=font_kwargs)
    if show_counts:
        for xpos, count in zip(x, combo_df["COUNT"]):
            ax_bar.text(xpos, count, str(int(count)), ha="center", va="bottom", fontsize=font_kwargs.get("fontsize", 9))
    if title is not None:
        ax_bar.set_title(title, fontdict=font_kwargs)

    for xpos, key in zip(x, combo_df["MEMBERSHIP_KEY"]):
        bits = _split_membership_key(key)
        active = [len(labels) - i - 1 for i, bit in enumerate(bits) if bit == "1"]
        if len(active) > 1:
            # Draw the connector only between the edges of active dots, not
            # through their centers, to avoid a line extending past the dot.
            y_min = min(active) + connector_gap
            y_max = max(active) - connector_gap
            if y_min > y_max:
                y_min, y_max = min(active), max(active)
            ax_matrix.plot(
                [xpos, xpos],
                [y_min, y_max],
                color=color,
                lw=1,
                zorder=2,
                solid_capstyle="butt",
            )

    for study_index, label in enumerate(labels):
        y = len(labels) - study_index - 1
        ax_matrix.text(-0.6, y, label, ha="right", va="center", fontsize=font_kwargs.get("fontsize", 9))
        for xpos, key in zip(x, combo_df["MEMBERSHIP_KEY"]):
            bits = _split_membership_key(key)
            present = study_index < len(bits) and bits[study_index] == "1"
            ax_matrix.scatter(
                xpos,
                y,
                s=60 if present else 25,
                color=color if present else "#cccccc",
                zorder=3,
            )

    ax_matrix.set_yticks([])
    ax_matrix.set_xticks([])
    ax_matrix.set_ylim(-0.75, len(labels) - 0.25)
    ax_matrix.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax_bar.spines[["top", "right"]].set_visible(False)
    plt.setp(ax_bar.get_xticklabels(), visible=False)

    return fig
