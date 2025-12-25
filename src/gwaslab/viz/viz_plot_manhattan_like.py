import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from gwaslab.info.g_Log import Log
from gwaslab.viz.viz_aux_quickfix import _quick_fix_p_value, _quick_fix_mlog10p
from gwaslab.util.util_in_get_sig import _get_sig, _get_top, _anno_gene


def _configure_fig_save_kwargs(mode="m",
                                save=None,
                                fig_kwargs=None,
                                scatter_kwargs=None,
                                qq_scatter_kwargs=None,
                                save_kwargs=None,
                                log=Log(),
                                verbose=True):
    if fig_kwargs is None:
        fig_kwargs = dict()
    if scatter_kwargs is None:
        scatter_kwargs = dict()
    if qq_scatter_kwargs is None:
        qq_scatter_kwargs = dict()
    if save_kwargs is None:
        save_kwargs = dict()

    if save is not None:
        if type(save) is not bool:
            if len(save) > 3:
                if save[-3:] == "pdf" or save[-3:] == "svg":
                    fig_kwargs["dpi"] = 72
                    if mode != "r":
                        scatter_kwargs["rasterized"] = True
                        qq_scatter_kwargs["rasterized"] = True
                        qq_scatter_kwargs["antialiased"] = False,
                        log.write("Saving as pdf/svg: scatter plot will be rasterized for mqq...", verbose=verbose)
                    else:
                        scatter_kwargs["rasterized"] = False
                        qq_scatter_kwargs["rasterized"] = False
                else:
                    fig_kwargs["dpi"] = save_kwargs["dpi"]
    return fig_kwargs, scatter_kwargs, qq_scatter_kwargs, save_kwargs


def _add_pad_to_x_axis(ax1, xpad, xpadl, xpadr, sumstats, pos, chrpad, xtight, log, verbose):
    if xtight == True:
        log.write(" -Adjusting X padding on both side : tight mode", verbose=verbose)
        xmax = sumstats["i"].max()
        xmin = sumstats["i"].min()
        ax1.set_xlim([xmin, xmax])
    else:
        chrpad_to_remove = sumstats[pos].max() * chrpad
        if ax1 is not None:
            xmax = sumstats["i"].max()
            xmin = sumstats["i"].min()
            length = xmax - xmin
            if xpad is not None:
                log.write(" -Adjusting X padding on both side: {}".format(xpad), verbose=verbose)
                pad = xpad * length
                ax1.set_xlim([xmin - pad + chrpad_to_remove, xmax + pad - chrpad_to_remove])
            if xpad is None and xpadl is not None:
                log.write(" -Adjusting X padding on left side: {}".format(xpadl), verbose=verbose)
                xmax = ax1.get_xlim()[1]
                pad = xpadl * length
                ax1.set_xlim([xmin - pad + chrpad_to_remove, xmax])
            if xpad is None and xpadr is not None:
                log.write(" -Adjusting X padding on right side: {}".format(xpadr), verbose=verbose)
                xmin = ax1.get_xlim()[0]
                pad = xpadr * length
                ax1.set_xlim([xmin, xmax + pad - chrpad_to_remove])
    return ax1


def _configure_cols_to_use(insumstats, snpid, chrom, pos, ea, nea, eaf, p, mlog10p, scaled, mode, stratified, anno, anno_set, anno_alias, _chrom_df_for_i, highlight, pinpoint, density_color):
    usecols = []
    if mlog10p in insumstats.columns and scaled is True:
        usecols.append(mlog10p)
    elif p in insumstats.columns:
        usecols.append(p)
    elif "b" in mode:
        pass
    else:
        raise ValueError("Please make sure " + p + " column is in input sumstats.")

    if (chrom is not None) and (pos is not None) and (("qq" in mode) or ("m" in mode) or ("r" in mode)):
        if chrom in insumstats.columns:
            usecols.append(chrom)
        else:
            raise ValueError("Please make sure " + chrom + " column is in input sumstats.")
        if pos in insumstats.columns:
            usecols.append(pos)
        else:
            raise ValueError("Please make sure " + pos + " column is in input sumstats.")

    if ("r" in mode) and (ea in insumstats.columns) and (nea in insumstats.columns):
        try:
            usecols.append(ea)
            usecols.append(nea)
        except:
            raise ValueError("Please make sure {} and {} columns are in input sumstats.".format(nea, ea))

    if len(highlight) > 0 or len(pinpoint) > 0 or (snpid is not None):
        if snpid in insumstats.columns:
            usecols.append(snpid)
        else:
            raise ValueError("Please make sure " + snpid + " column is in input sumstats.")

    if (stratified is True) and (eaf is not None):
        if eaf in insumstats.columns:
            usecols.append(eaf)
        else:
            raise ValueError("Please make sure " + eaf + " column is in input sumstats.")

    if _chrom_df_for_i is not None:
        usecols.append("i")

    if (anno_set is not None and len(anno_set) > 0) or (anno_alias is not None and len(anno_alias) > 0):
        if anno is None:
            anno = True

    if isinstance(anno,str):
        if anno == "GENENAME":
            pass
        elif (anno in insumstats.columns):
            if (anno not in usecols):
                usecols.append(anno)
        else:
            raise ValueError("Please make sure " + anno + " column is in input sumstats.")

    if (density_color == True) and ("b" in mode and "DENSITY" in insumstats.columns):
        usecols.append("DENSITY")
    return usecols


def _sanity_check(sumstats, mode, chrom, pos, stratified, _if_quick_qc, log, verbose):
    if ("m" in mode or "r" in mode) and _if_quick_qc:
        pre_number = len(sumstats)
        sumstats = sumstats.dropna(subset=[chrom, pos])
        after_number = len(sumstats)
        removed_count = pre_number - after_number
        log.write(" -Removed " + str(removed_count) + " variants with nan in CHR or POS column ...", verbose=verbose if removed_count > 0 else False)
        out_of_range_chr = sumstats[chrom] <= 0
        removed_chr_count = sum(out_of_range_chr)
        log.write(" -Removed {} variants with CHR <=0...".format(removed_chr_count), verbose=verbose if removed_chr_count > 0 else False)
        sumstats = sumstats.loc[~out_of_range_chr, :]

    if stratified is True and _if_quick_qc:
        pre_number = len(sumstats)
        sumstats = sumstats.dropna(subset=["MAF"])
        after_number = len(sumstats)
        log.write(" -Removed " + str(pre_number - after_number) + " variants with nan in EAF column ...", verbose=verbose)

    if "b" not in mode and _if_quick_qc:
        pre_number = len(sumstats)
        sumstats = sumstats.dropna(subset=["raw_P"])
        after_number = len(sumstats)
        removed_count = pre_number - after_number
        log.write(" -Removed " + str(removed_count) + " variants with nan in P column ...", verbose=verbose if removed_count > 0 else False)
    return sumstats


def _process_p_value(sumstats, mode, p, mlog10p, scaled, log, verbose):
    if "b" in mode:
        sumstats["scaled_P"] = sumstats["DENSITY"].copy()
        sumstats["raw_P"] = -np.log10(sumstats["DENSITY"].copy() + 2)
    elif scaled is True:
        log.write(" -P values are already converted to -log10(P)!", verbose=verbose)
        sumstats["scaled_P"] = sumstats["raw_P"].copy()
        sumstats["raw_P"] = np.power(10, -sumstats["scaled_P"].astype("float64"))
    else:
        if not scaled:
            sumstats = _quick_fix_p_value(sumstats, p=p, mlog10p=mlog10p, verbose=verbose, log=log)
        sumstats = _quick_fix_mlog10p(sumstats, p=p, mlog10p=mlog10p, scaled=scaled, verbose=verbose, log=log)
    return sumstats


def _process_highlight(sumstats, highlight, highlight_chrpos, highlight_windowkb, highlight_lim, highlight_lim_mode, snpid, chrom, pos):
    if pd.api.types.is_list_like(highlight[0]):
        if highlight_chrpos == False:
            for i, highlight_set in enumerate(highlight):
                to_highlight = sumstats.loc[sumstats[snpid].isin(highlight_set), :]
                for j, (index, row) in enumerate(to_highlight.iterrows()):
                    target_chr = int(row[chrom])
                    target_pos = int(row[pos])
                    right_chr = sumstats[chrom] == target_chr
                    if (highlight_lim is not None and i < len(highlight_lim) and j < len(highlight_lim[i]) and highlight_lim[i][j] is not None):
                        if highlight_lim_mode == "absolute":
                            start_pos, end_pos = highlight_lim[i][j]
                            up_pos = sumstats[pos] >= start_pos
                            low_pos = sumstats[pos] <= end_pos
                        else:
                            lower_kb, upper_kb = highlight_lim[i][j]
                            up_pos = sumstats[pos] > target_pos + lower_kb * 1000
                            low_pos = sumstats[pos] < target_pos + upper_kb * 1000
                    else:
                        up_pos = sumstats[pos] > target_pos - highlight_windowkb * 1000
                        low_pos = sumstats[pos] < target_pos + highlight_windowkb * 1000
                    sumstats.loc[right_chr & up_pos & low_pos, "HUE"] = i
        else:
            for i, highlight_chrpos_tuple in enumerate(highlight):
                if len(highlight_chrpos_tuple) == 2:
                    target_chr = int(highlight_chrpos_tuple[0])
                    target_pos = int(highlight_chrpos_tuple[1])
                    right_chr = sumstats[chrom] == target_chr
                    up_pos = sumstats[pos] > target_pos - highlight_windowkb * 1000
                    low_pos = sumstats[pos] < target_pos + highlight_windowkb * 1000
                    sumstats.loc[right_chr & up_pos & low_pos, "HUE"] = 0
                elif len(highlight_chrpos_tuple) == 3:
                    target_chr = int(highlight_chrpos_tuple[0])
                    target_pos_low = int(highlight_chrpos_tuple[1])
                    target_pos_up = int(highlight_chrpos_tuple[2])
                    right_chr = sumstats[chrom] == target_chr
                    up_pos = sumstats[pos] > target_pos_low
                    low_pos = sumstats[pos] < target_pos_up
                    sumstats.loc[right_chr & up_pos & low_pos, "HUE"] = 0
    else:
        to_highlight = sumstats.loc[sumstats[snpid].isin(highlight), :]
        for j, (index, row) in enumerate(to_highlight.iterrows()):
            target_chr = int(row[chrom])
            target_pos = int(row[pos])
            right_chr = sumstats[chrom] == target_chr
            if (highlight_lim is not None and j < len(highlight_lim) and highlight_lim[j] is not None):
                if highlight_lim_mode == "absolute":
                    start_pos, end_pos = highlight_lim[j]
                    up_pos = sumstats[pos] >= start_pos
                    low_pos = sumstats[pos] <= end_pos
                else:
                    lower_kb, upper_kb = highlight_lim[j]
                    up_pos = sumstats[pos] > target_pos + lower_kb * 1000
                    low_pos = sumstats[pos] < target_pos + upper_kb * 1000
            else:
                up_pos = sumstats[pos] > target_pos - highlight_windowkb * 1000
                low_pos = sumstats[pos] < target_pos + highlight_windowkb * 1000
            sumstats.loc[right_chr & up_pos & low_pos, "HUE"] = 0
    return sumstats


def _process_line(ax1, sig_line, suggestive_sig_line, additional_line, lines_to_plot, sc_linewidth, sig_line_color, suggestive_sig_line_color, additional_line_color, mode, bmean, bmedian, log=Log(), verbose=True):
    log.write(" -Processing lines...", verbose=False)
    if sig_line is True:
        ax1.axhline(y=lines_to_plot[0], linewidth=sc_linewidth, linestyle="--", color=sig_line_color, zorder=1)
    if suggestive_sig_line is True:
        ax1.axhline(y=lines_to_plot[1], linewidth=sc_linewidth, linestyle="--", color=suggestive_sig_line_color, zorder=1)
    if additional_line is not None:
        for index, level in enumerate(lines_to_plot[2:2 + len(additional_line)].values):
            ax1.axhline(y=level, linewidth=sc_linewidth, linestyle="--", color=additional_line_color[index % len(additional_line_color)], zorder=1)
    if "b" in mode:
        bmean = lines_to_plot.iat[-2]
        bmedian = lines_to_plot.iat[-1]
        log.write(" -Plotting horizontal line (  mean DENISTY): y = {}".format(bmean), verbose=verbose)
        ax1.axhline(y=bmean, linewidth=sc_linewidth, linestyle="-", color=sig_line_color, zorder=1000)
        log.write(" -Plotting horizontal line ( median DENISTY): y = {}".format(bmedian), verbose=verbose)
        ax1.axhline(y=bmedian, linewidth=sc_linewidth, linestyle="--", color=sig_line_color, zorder=1000)
    return ax1


def _process_cbar(cbar, cbar_fontsize, cbar_font_family, cbar_title, log=Log(), verbose=True):
    log.write(" -Processing color bar...", verbose=False)
    cbar_yticklabels = cbar.get_yticklabels()
    cbar.set_yticklabels(cbar_yticklabels, fontsize=cbar_fontsize, family=cbar_font_family)
    cbar_xticklabels = cbar.get_xticklabels()
    cbar.set_xticklabels(cbar_xticklabels, fontsize=cbar_fontsize, family=cbar_font_family)
    cbar.set_title(cbar_title, fontsize=cbar_fontsize, family=cbar_font_family, loc="center", y=1.00)
    return cbar


def _process_xtick(ax1, mode, chrom_df, xtick_chr_dict, fontsize, font_family="Arial", ax3=None, log=Log(), verbose=True):
    log.write(" -Processing X ticks...", verbose=False)
    if mode != "r":
        ax1.set_xticks(chrom_df.astype("float64"))
        ax1.set_xticklabels(chrom_df.index.astype("Int64").map(xtick_chr_dict), fontsize=fontsize, family=font_family)
    if ax3 is not None:
        ax3.tick_params(axis='x', labelsize=fontsize, labelfontfamily=font_family)
    return ax1, ax3


def _process_ytick(ax1, fontsize, font_family, ax4, log=Log(), verbose=True):
    log.write(" -Processing Y labels...", verbose=False)
    ax1.tick_params(axis='y', labelsize=fontsize, labelfontfamily=font_family)
    if ax4 is not None:
        ax4.tick_params(axis='y', labelsize=fontsize, labelfontfamily=font_family)
    return ax1, ax4


def _process_xlabel(region, xlabel, ax1, gtf_path, mode, fontsize, font_family="Arial", ax3=None, log=Log(), verbose=True):
    log.write(" -Processing X labels...", verbose=False)
    if region is not None:
        if xlabel is None:
            xlabel = "Chromosome " + str(region[0]) + " (MB)"
        if (gtf_path is not None) and ("r" in mode):
            ax3.set_xlabel(xlabel, fontsize=fontsize, family=font_family)
        else:
            ax1.set_xlabel(xlabel, fontsize=fontsize, family=font_family)
    else:
        if xlabel is None:
            xlabel = "Chromosome"
        ax1.set_xlabel(xlabel, fontsize=fontsize, family=font_family)
    return ax1, ax3


def _process_ylabel(ylabel, ax1, mode, bwindowsizekb, fontsize, font_family, math_fontfamily, ax4=None, log=Log(), verbose=True):
    log.write(" -Processing Y labels...", verbose=False)
    if "b" in mode:
        if ylabel is None:
            ylabel = "Density of GWAS \n SNPs within " + str(bwindowsizekb) + " kb"
        ax1.set_ylabel(ylabel, ha="center", va="bottom", fontsize=fontsize, family=font_family, math_fontfamily=math_fontfamily)
    else:
        if ylabel is None:
            ylabel = "$\\mathregular{-log_{10}(P)}$"
        ax1.set_ylabel(ylabel, fontsize=fontsize, family=font_family, math_fontfamily=math_fontfamily)
    if ax4 is not None:
        ax4_ylabel = ax4.get_ylabel()
        ax4.set_ylabel(ax4_ylabel, fontsize=fontsize, family=font_family, math_fontfamily=math_fontfamily)
    return ax1, ax4


def _extract_to_annotate(sumstats,
                         snpid,
                         chrom,
                         pos,
                         scaled_threhosld,
                         anno_sig_level,
                         windowsizekb,
                         scaled,
                         anno,
                         anno_set,
                         anno_source,
                         anno_gtf_path,
                         build,
                         mode,
                         log=Log(),
                         verbose=True):
    log.write("Start to extract variants for annotation...", verbose=verbose)
    # Initialize empty DataFrame to store variants to annotate
    to_annotate = pd.DataFrame()
    
    # Check if annotation is enabled (anno is a string/column name) or if specific annotation set is provided
    if (anno and anno != True) or (len(anno_set) > 0):
        # Case 1: User provided a specific set of variants to annotate
        if len(anno_set) > 0:
            # Filter sumstats to only include variants in the specified annotation set
            to_annotate = sumstats.loc[sumstats[snpid].isin(anno_set), :]
            # Log the number of specified variants found
            if (to_annotate.empty is not True):
                log.write(" -Found " + str(len(to_annotate)) + " specified variants to annotate...", verbose=verbose)
        else:
            # Case 2: Annotation is enabled but no specific set provided - auto-detect variants
            # For bubble plot mode: get top variants by density
            if "b" in mode:
                to_annotate = _get_top(
                    sumstats,
                    id="i",
                    chrom=chrom,
                    pos=pos,
                    by="DENSITY",
                    threshold=scaled_threhosld,
                    windowsizekb=windowsizekb,
                    verbose=False,
                )
            else:
                # For regular Manhattan plot mode: get significant variants using sliding window
                # Filter to variants at or above threshold, then find lead variants
                # Use >= instead of > to ensure anno_sig_level threshold is properly applied
                to_annotate = _get_sig(
                    sumstats.loc[sumstats["scaled_P"] >= scaled_threhosld, :],
                    snpid,
                    chrom,
                    pos,
                    "raw_P",
                    sig_level=anno_sig_level,
                    windowsizekb=windowsizekb,
                    scaled=scaled,
                    mlog10p="scaled_P",
                    verbose=False,
                )
            # Handle case where _get_top or _get_sig returns None
            if to_annotate is None:
                to_annotate = pd.DataFrame()
            # Log the number of significant variants found (only for non-bubble mode)
            if ("b" not in mode) and (to_annotate.empty is not True):
                log.write(" -Found " + str(len(to_annotate)) + " significant variants with a sliding window size of " + str(windowsizekb) + " kb...", verbose=verbose)
    else:
        # Case 3: Annotation not explicitly enabled, but still extract variants for potential annotation
        # For bubble plot mode: get top variants by density
        if "b" in mode:
            to_annotate = _get_top(
                sumstats,
                variant_id="i",
                chrom=chrom,
                pos=pos,
                by="DENSITY",
                threshold=scaled_threhosld,
                windowsizekb=windowsizekb,
                verbose=False,
            )
        else:
            # For regular Manhattan plot mode: get significant variants using sliding window
            # Use "i" (index column) instead of snpid for variant identification
            # Use >= instead of > to ensure anno_sig_level threshold is properly applied
            to_annotate = _get_sig(
                sumstats.loc[sumstats["scaled_P"] >= scaled_threhosld, :],
                "i",
                chrom,
                pos,
                "raw_P",
                windowsizekb=windowsizekb,
                scaled=scaled,
                verbose=False,
                mlog10p="scaled_P",
                sig_level=anno_sig_level,
            )
        if to_annotate is None:
            to_annotate = pd.DataFrame()
        if ("b" not in mode) and (to_annotate.empty is not True):
            log.write(" -Found " + str(len(to_annotate)) + " significant variants with a sliding window size of " + str(windowsizekb) + " kb...", verbose=verbose)
    if (to_annotate is not None) and (to_annotate.empty is not True) and anno == "GENENAME":
        to_annotate = _anno_gene(to_annotate,
                               id=snpid,
                               chrom=chrom,
                               pos=pos,
                               log=log,
                               build=build,
                               source=anno_source,
                               gtf_path=anno_gtf_path,
                               verbose=verbose).rename(columns={"GENE": "Annotation"})
        if "b" in mode and (to_annotate.empty is not True):
            from gwaslab.viz.viz_plot_density_mode import b_annotation_log
            b_annotation_log(to_annotate, snpid, log=log, verbose=verbose)
    log.write("Finished extracting variants for annotation...", verbose=verbose)
    return to_annotate


def _process_spine(ax1, mode):
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(True)
    if mode == "r":
        ax1.spines["top"].set_visible(True)
        ax1.spines["top"].set_zorder(1)
        ax1.spines["right"].set_visible(True)
    return ax1


def _process_layout(mode, figax, fig_kwargs, mqqratio, region_hspace):
    explicit = {"gridspec_kw"}
    fig_kwargs = {k: v for k, v in fig_kwargs.items() if k not in explicit}
    if mode == "qqm":
        fig, (ax2, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, mqqratio]}, **fig_kwargs)
        ax3 = None
    elif mode == "mqq":
        if figax is not None:
            fig = figax[0]
            ax1 = figax[1]
            ax2 = figax[2]
        else:
            fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [mqqratio, 1]}, **fig_kwargs)
        ax3 = None
    elif mode == "m":
        if figax is not None:
            fig = figax[0]
            ax1 = figax[1]
        else:
            fig, ax1 = plt.subplots(1, 1, **fig_kwargs)
        ax2 = None
        ax3 = None
    elif mode == "qq":
        fig, ax2 = plt.subplots(1, 1, **fig_kwargs)
        ax1 = None
        ax3 = None
    elif mode == "r":
        if figax is not None:
            fig = figax[0]
            ax1 = figax[1]
            ax3 = figax[2]
            ax2 = None
        else:
            fig_kwargs["figsize"] = (15, 10)
            fig, (ax1, ax3) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [mqqratio, 1]}, **fig_kwargs)
            ax2 = None
            plt.subplots_adjust(hspace=region_hspace)
    elif mode == "b":
        if figax is not None:
            fig = figax[0]
            ax1 = figax[1]
            ax3 = None
            ax2 = None
        else:
            fig_kwargs["figsize"] = (15, 5)
            fig, ax1 = plt.subplots(1, 1, **fig_kwargs)
            ax2 = None
            ax3 = None
    else:
        raise ValueError("Please select one from the 5 modes: mqq/qqm/m/qq/r/b")
    ax4 = None
    cbar = None
    return fig, ax1, ax2, ax3, ax4, cbar
