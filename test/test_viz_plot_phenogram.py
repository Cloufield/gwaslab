import gzip
import os
import sys
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.viz.viz_plot_phenogram import _plot_phenogram
from gwaslab.viz.viz_aux_phenogram import (
    _AnnotationLayoutBlock,
    _build_lead_text_kwargs,
    _block_bottom,
    _block_bottom_from_center,
    _block_top,
    _block_top_from_center,
    _centromere_cap_path,
    _centromere_fill_patches,
    _chr_ideogram_fill_patches,
    _chr_ideogram_outline_path,
    _chr_ideogram_top_y,
    _chr_number_label_position,
    _compute_chr_ideogram_geometry,
    _compute_group_marker_positions,
    _compute_phenogram_xlim,
    _cap_phenogram_arrow_shaft_pt,
    _effective_marker_radius_pt,
    _format_phenogram_label,
    _label_eligible_mask,
    _label_width_data_dx,
    _marker_legend_row_handles,
    _make_marker_color_legend_handles,
    _make_marker_shape_legend_handles,
    _marker_center_step_pt,
    _marker_anchor_from_block_center,
    _marker_group_block_extents_pt,
    _marker_label_x_and_ha,
    _marker_label_position,
    _marker_radius_pt,
    _marker_row_x_bounds,
    _plot_chr_ideogram_outline,
    _plot_phenogram_text_connector,
    _points_to_data_delta,
    _phenogram_provisional_xmax,
    _phenogram_text_fontproperties,
    _prepare_marker_group_display_label,
    _prepare_phenogram_display_label,
    _resolve_lead_connector_style,
    _resolve_marker_group_label,
    _resolve_phenogram_chr_width,
    _resolve_phenogram_label,
    _ANNOTATION_PATH_PAD_PT,
    _ANNOTATION_Y_MARGIN_DATA,
    _chr_panel_bottom_from_row_top,
    _compute_phenogram_row_bands,
    _unit_bottom_data_y,
    _phenogram_inter_unit_gap_pt,
    _phenogram_panel_ceiling_y,
    _phenogram_panel_floor_y,
    _place_marker_unit_at_center,
    _pt_to_data_height,
    _solve_bounded_stack_layout,
    _spread_group_blocks_display,
    _TELOMERE_LENGTH,
    _spread_phenogram_text_labels,
    _text_bbox_attach_point,
    _wrap_phenogram_label,
    _wrap_phenogram_label_chars,
    _DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
    _DEFAULT_MARKER_MAX_PER_ROW,
    _LEGEND_FONTSIZE,
    _LEGEND_MARKER_SIZE,
)


class TestPhenogramAnno(unittest.TestCase):
    def _make_leads(self):
        return pd.DataFrame(
            {
                "SNPID": ["rs1", "rs2", "rs3"],
                "CHR": [1, 1, 2],
                "POS": [1000, 2000, 3000],
                "P": [1e-10, 1e-9, 1e-8],
                "MLOG10P": [10.0, 9.0, 8.0],
                "GENE": ["AAA", "BBB", "CCC"],
            }
        )

    def test_resolve_label_none(self):
        lead = {"chrom": 1, "pos": 1000, "snpid": "rs1"}
        self.assertIsNone(_resolve_phenogram_label(lead, None))

    def test_resolve_label_true(self):
        lead = {"chrom": 19, "pos": 46214297, "snpid": "rs1"}
        self.assertEqual(
            _resolve_phenogram_label(lead, True),
            "Chr19:46214297",
        )

    def test_resolve_label_column(self):
        lead = {"chrom": 1, "pos": 1000, "snpid": "rs1", "annotation": "FTO"}
        self.assertEqual(_resolve_phenogram_label(lead, "GENE"), "FTO")

    def test_resolve_label_alias(self):
        lead = {"chrom": 1, "pos": 1000, "snpid": "rs1", "annotation": "FTO"}
        self.assertEqual(
            _resolve_phenogram_label(lead, "GENE", {"rs1": "Custom"}),
            "Custom",
        )

    def test_anno_set_filters_labels(self):
        leads = self._make_leads()
        mask = _label_eligible_mask(
            leads,
            snpid="SNPID",
            anno_set=["rs1", "rs3"],
            anno_max_rows=40,
            p="P",
            mlog10p="MLOG10P",
            verbose=False,
        )
        self.assertTrue(mask.loc[0])
        self.assertFalse(mask.loc[1])
        self.assertTrue(mask.loc[2])

    def test_anno_set_zero_match_all_false(self):
        leads = self._make_leads()
        mask = _label_eligible_mask(
            leads,
            snpid="SNPID",
            anno_set=["missing_id"],
            anno_max_rows=40,
            p="P",
            mlog10p="MLOG10P",
            verbose=False,
        )
        self.assertFalse(mask.any())

    def test_format_phenogram_label_truncates(self):
        long_text = "A" * 30
        self.assertEqual(_format_phenogram_label(long_text), "A" * 17 + "...")

    def test_compute_phenogram_xlim_keeps_chr_fraction(self):
        xmin = _phenogram_provisional_xmax(0.0, 0.35, 0.18) * -0.1
        xmin = -0.09
        wide_label_xmax = 5.0
        xmax = _compute_phenogram_xlim(xmin, 0.0, 0.35, wide_label_xmax, anno_x_pad=0.18)
        self.assertGreaterEqual(0.35 / (xmax - xmin), 0.25)

    def test_resolve_phenogram_chr_width_scales_with_ncols(self):
        self.assertGreater(
            _resolve_phenogram_chr_width(0.35, 6),
            _resolve_phenogram_chr_width(0.35, 12),
        )

    def test_text_bbox_attach_left_of_anchor(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        text = ax.text(0.5, 0.5, "GENE1", ha="left", va="center")
        fig.canvas.draw()
        attach_x, attach_y = _text_bbox_attach_point(
            ax, text, fig.canvas.get_renderer(), side="left"
        )
        self.assertLessEqual(attach_x, 0.5)
        self.assertAlmostEqual(attach_y, 0.5, delta=0.05)
        plt.close(fig)

    def test_anno_max_rows_limits_labels(self):
        leads = self._make_leads()
        mask = _label_eligible_mask(
            leads,
            snpid="SNPID",
            anno_set=[],
            anno_max_rows=2,
            p="P",
            mlog10p="MLOG10P",
            verbose=False,
        )
        self.assertEqual(int(mask.sum()), 2)
        self.assertTrue(mask.loc[0])
        self.assertTrue(mask.loc[1])
        self.assertFalse(mask.loc[2])

    def test_marker_center_step_fixed(self):
        step = _marker_center_step_pt(42, 16)
        expected = 2 * _effective_marker_radius_pt(42, 0.6) + 16
        self.assertAlmostEqual(step, expected, places=4)

    def test_effective_marker_radius_includes_linewidth(self):
        base = _marker_radius_pt(42)
        effective = _effective_marker_radius_pt(42, 0.6)
        self.assertAlmostEqual(effective, base + 0.3, places=6)

    def test_marker_display_spacing_matches_step(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        marker_size = 42.0
        marker_gap_pt = 16.0
        marker_linewidth = 0.6
        step_pt = _marker_center_step_pt(marker_size, marker_gap_pt, marker_linewidth)
        items = [{"lead": {"pos": 1}}, {"lead": {"pos": 2}}]
        positions, _ = _compute_group_marker_positions(
            ax,
            0.5,
            0.3,
            items,
            marker_size,
            marker_gap_pt,
            marker_max_per_row=6,
            marker_row_gap_pt=12.0,
            marker_linewidth=marker_linewidth,
        )
        self.assertEqual(len(positions), 2)
        x1, y1, _ = positions[0]
        x2, y2, _ = positions[1]
        fig.canvas.draw()
        disp1 = ax.transData.transform((x1, y1))
        disp2 = ax.transData.transform((x2, y2))
        dist_px = abs(disp2[0] - disp1[0])
        expected_px = step_pt * fig.dpi / 72.0
        self.assertAlmostEqual(dist_px, expected_px, delta=1.0)
        plt.close(fig)

    def test_two_groups_same_marker_pixel_gap(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        marker_size = 42.0
        marker_gap_pt = 16.0
        marker_linewidth = 0.6
        items = [{"lead": {"pos": i}} for i in range(3)]

        pos_a, _ = _compute_group_marker_positions(
            ax, 0.5, 0.2, items, marker_size, marker_gap_pt, 6, 12.0, marker_linewidth
        )
        pos_b, _ = _compute_group_marker_positions(
            ax, 0.5, 0.8, items, marker_size, marker_gap_pt, 6, 12.0, marker_linewidth
        )
        fig.canvas.draw()

        def pair_gap(positions):
            x1, y1, _ = positions[0]
            x2, y2, _ = positions[1]
            d1 = ax.transData.transform((x1, y1))
            d2 = ax.transData.transform((x2, y2))
            return abs(d2[0] - d1[0])

        self.assertAlmostEqual(pair_gap(pos_a), pair_gap(pos_b), delta=0.5)
        plt.close(fig)

    def test_text_connector_fixed_horizontal_arm(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        chr_right_x = 0.35
        locus_y = 0.4
        text_x = 0.55
        text_y = 0.7
        arrow_shaft_pt = 18.0

        _plot_phenogram_text_connector(
            ax,
            chr_right_x,
            locus_y,
            text_x,
            text_y,
            arrow_shaft_pt=arrow_shaft_pt,
        )
        fig.canvas.draw()

        shaft_pt = _cap_phenogram_arrow_shaft_pt(
            ax, chr_right_x, locus_y, text_x, arrow_shaft_pt
        )
        dx_arm, _ = _points_to_data_delta(
            ax, chr_right_x, locus_y, dx_points=shaft_pt
        )
        elbow_x = chr_right_x + dx_arm
        elbow_disp = ax.transData.transform((elbow_x, locus_y))
        locus_disp = ax.transData.transform((chr_right_x, locus_y))
        arm_px = abs(elbow_disp[0] - locus_disp[0])
        expected_px = shaft_pt * fig.dpi / 72.0
        self.assertAlmostEqual(arm_px, expected_px, delta=1.5)
        self.assertLess(elbow_x, text_x)

        text_disp = ax.transData.transform((text_x, text_y))
        self.assertGreater(abs(text_disp[1] - elbow_disp[1]), 1.0)
        plt.close(fig)

    def test_arrow_shaft_capped_before_text(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        chr_right_x = 0.35
        text_x = 0.45
        capped = _cap_phenogram_arrow_shaft_pt(
            ax, chr_right_x, 0.4, text_x, arrow_shaft_pt=100.0
        )
        gap_pt = (
            abs(
                ax.transData.transform((text_x, 0.4))[0]
                - ax.transData.transform((chr_right_x, 0.4))[0]
            )
            * 72.0
            / fig.dpi
        )
        self.assertLess(capped, gap_pt)
        self.assertLessEqual(capped, 100.0)
        plt.close(fig)

    def test_anno_kwargs_single_overrides_text(self):
        base = {"fontsize": 8, "color": "black"}
        lead = {"snpid": "rs1"}
        single = {"rs1": {"fontsize": 12, "color": "red", "arrow_shaft": 22}}
        merged = _build_lead_text_kwargs(base, single, lead)
        self.assertEqual(merged["fontsize"], 12)
        self.assertEqual(merged["color"], "red")
        self.assertNotIn("arrow_shaft", merged)

    def test_arrow_kwargs_and_single_merge(self):
        anno = {"arrow_shaft": 18, "arrow_pad": 10}
        single = {"rs1": {"arrow_shaft": 24, "color": "blue"}}
        shaft, pad, shrink, arrowprops, line_kwargs = _resolve_lead_connector_style(
            anno,
            {"lw": 1.2, "color": "gray"},
            single,
            "rs1",
        )
        self.assertEqual(shaft, 24)
        self.assertEqual(arrowprops["lw"], 1.2)
        self.assertEqual(line_kwargs["color"], "gray")

    def test_marker_group_block_grows_with_rows(self):
        one_row = _marker_group_block_extents_pt(42, 14, 10, 1.5, n_marker_rows=1)
        two_rows = _marker_group_block_extents_pt(42, 14, 10, 1.5, n_marker_rows=2)
        self.assertGreater(two_rows[1], one_row[1])

    def test_chr_ideogram_top_includes_telomere_cap(self):
        self.assertAlmostEqual(_chr_ideogram_top_y(0.6, 0.02), 0.62)

    def test_chr_number_label_above_ideogram(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.5)
        ax.set_ylim(-0.2, 1.2)
        fig.canvas.draw()

        body_top_y = _chr_ideogram_top_y(0.5, 0.02)
        label_y, label_va = _chr_number_label_position(
            ax, body_top_y, ref_x=0.35, fontsize=10, pad_pt=4
        )
        self.assertEqual(label_va, "bottom")
        body_disp = ax.transData.transform((0.35, body_top_y))
        label_disp = ax.transData.transform((0.35, label_y))
        self.assertGreater(label_disp[1], body_disp[1])

        t = ax.text(
            0.35, label_y, "22", ha="center", va=label_va, fontsize=10
        )
        fig.canvas.draw()
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        self.assertGreater(bb.y0, body_disp[1])
        plt.close(fig)

    def test_chr_side_edges_drawn(self):
        fig, ax = plt.subplots(figsize=(4, 6))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 2)
        geom = _compute_chr_ideogram_geometry(
            chr_x=0.1,
            chr_width=0.3,
            chr_size=100_000_000,
            chr_centromere_u=44_000_000,
            chr_centromere_l=46_000_000,
            max_chr_size=250_000_000,
            offset=0.2,
        )
        _plot_chr_ideogram_outline(ax, geom)
        outline = ax.patches[-1]
        self.assertAlmostEqual(outline.get_linewidth(), 0.8)
        plt.close(fig)

    def test_chr_geometry_band_width_matches_body(self):
        geom = _compute_chr_ideogram_geometry(
            chr_x=0.0,
            chr_width=0.35,
            chr_size=100_000_000,
            chr_centromere_u=44_000_000,
            chr_centromere_l=46_000_000,
            max_chr_size=250_000_000,
            offset=0.0,
        )
        self.assertAlmostEqual(geom.right_x - geom.left_x, geom.width)
        self.assertAlmostEqual(geom.center_x, geom.left_x + geom.width / 2.0)

    def test_marker_group_label_uses_anno(self):
        group_items = [
            {
                "lead": {
                    "annotation": "FTO",
                    "anno_group": "LOCUS1",
                    "snpid": "rs1",
                }
            }
        ]
        self.assertEqual(
            _resolve_marker_group_label("LOCUS1", group_items, anno="GENE"),
            "FTO",
        )
        self.assertEqual(
            _resolve_marker_group_label("LOCUS1", group_items, anno=None),
            "FTO",
        )

    def test_marker_label_at_fixed_anno_column(self):
        row_left, row_center, row_right = 0.5, 0.7, 0.9
        anno_col_x = 0.53
        label_x, ha = _marker_label_x_and_ha(
            row_left, row_center, row_right, "center", anno_col_x
        )
        self.assertEqual(label_x, anno_col_x)
        self.assertEqual(ha, "left")

    def test_marker_label_align_center_without_column(self):
        row_left, row_center, row_right = 0.5, 0.7, 0.9
        label_x, ha = _marker_label_x_and_ha(
            row_left, row_center, row_right, "center"
        )
        self.assertEqual(label_x, row_center)
        self.assertEqual(ha, "center")

    def test_marker_connector_reaches_row(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()
        items = [{"lead": {"pos": 1}}, {"lead": {"pos": 2}}]
        positions, _ = _compute_group_marker_positions(
            ax, 0.55, 0.3, items, 42, 16, 6, 12.0, 0.6
        )
        _, row_center, _ = _marker_row_x_bounds(ax, positions, 42, 0.6)
        self.assertGreaterEqual(row_center, min(x for x, _, _ in positions))
        plt.close(fig)

    def test_marker_multiline_label_sits_lower(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        _, y_one = _marker_label_position(
            ax, 0.55, 0.35, 81, 2.75, 0.6, n_text_lines=1, marker_fontsize=11
        )
        _, y_two = _marker_label_position(
            ax, 0.55, 0.35, 81, 2.75, 0.6, n_text_lines=2, marker_fontsize=11
        )
        self.assertLess(y_two, y_one)
        plt.close(fig)

    def test_marker_label_clear_of_marker_row(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        marker_size = 81
        anno_col_x = 0.55
        bottom_row_y = 0.35
        group_items = [
            {"lead": {"annotation": "GENE"}, "marker_shape": "o", "marker_color": "red"},
        ]
        positions, _ = _compute_group_marker_positions(
            ax, anno_col_x, bottom_row_y, group_items,
            marker_size, 4.0, 4, 2.0, marker_linewidth=0.6,
        )
        _, label_y = _marker_label_position(
            ax, anno_col_x, bottom_row_y, marker_size, 2.75, 0.6,
            n_text_lines=2, marker_fontsize=11,
        )
        text = ax.text(
            anno_col_x, label_y, "Line one\nLine two",
            ha="left", va="bottom", fontsize=11,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        text_bb = text.get_window_extent(renderer)
        radius_pt = _effective_marker_radius_pt(marker_size, 0.6)
        scale = fig.dpi / 72.0
        _, marker_bottom_disp = ax.transData.transform(
            (anno_col_x, bottom_row_y)
        )
        marker_bottom_disp -= radius_pt * scale
        self.assertLessEqual(
            text_bb.y1, marker_bottom_disp - 2.75 * scale + 1.0
        )
        plt.close(fig)

    def test_place_marker_unit_label_clear_of_markers(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 1.2)
        fig.canvas.draw()

        ref_x = 0.55
        group_items = [
            {"lead": {"annotation": "A"}, "marker_shape": "o", "marker_color": "red"},
            {"lead": {"annotation": "B"}, "marker_shape": "^", "marker_color": "blue"},
            {"lead": {"annotation": "C"}, "marker_shape": "s", "marker_color": "green"},
            {"lead": {"annotation": "D"}, "marker_shape": "D", "marker_color": "black"},
            {"lead": {"annotation": "E"}, "marker_shape": "v", "marker_color": "grey"},
        ]
        placed = _place_marker_unit_at_center(
            ax, ref_x, 0.5, group_items,
            marker_size=81, marker_gap_pt=4.0, marker_label_gap_pt=2.75,
            marker_max_per_row=4, marker_row_gap_pt=2.0, marker_linewidth=0.6,
            marker_fontsize=11, n_text_lines=2, group_label_box_pad_pt=1.5,
            marker_label_align="left",
        )
        scale = fig.dpi / 72.0
        radius_pt = _effective_marker_radius_pt(81, 0.6)
        min_marker_bottom = min(
            ax.transData.transform((x, y))[1] - radius_pt * scale
            for x, y, _ in placed["marker_positions"]
        )
        text = ax.text(
            placed["label_x"], placed["label_y"], "GENE\nLINE2",
            ha=placed["label_ha"], va="bottom", fontsize=11,
        )
        fig.canvas.draw()
        text_bb = text.get_window_extent(fig.canvas.get_renderer())
        self.assertLessEqual(text_bb.y1, min_marker_bottom - 2.75 * scale + 1.0)
        plt.close(fig)

    def test_marker_label_below_row(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()

        marker_y = 0.35
        anno_col_x = 0.55
        marker_size = 81
        _, label_y = _marker_label_position(
            ax, anno_col_x, marker_y, marker_size, 2.75, 0.6
        )
        marker_disp_y = ax.transData.transform((anno_col_x, marker_y))[1]
        text = ax.text(
            anno_col_x,
            label_y,
            "GROUP",
            ha="left",
            va="bottom",
            fontsize=11,
        )
        fig.canvas.draw()
        bbox = text.get_window_extent(fig.canvas.get_renderer())
        self.assertLess(bbox.y0, marker_disp_y)
        plt.close(fig)

    def test_panel_floor_caps_spread_below_stacked_chr(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 2.0)
        fig.canvas.draw()

        ref_x = 0.55
        below_pt = 40.0
        panel_floor_y = 1.2
        y_vals = [0.35, 0.38, 0.41]
        spread = _spread_group_blocks_display(
            ax,
            y_vals,
            ref_x=ref_x,
            above_pt=10.0,
            below_pt=below_pt,
            gap_pt=4.0,
            max_iter=50,
            panel_floor_y=panel_floor_y,
        )
        scale = fig.dpi / 72.0
        for y in spread:
            bottom_y = _block_bottom_from_center(
                ax, ref_x, y, 10.0, below_pt
            )
            bottom_disp = ax.transData.transform((ref_x, bottom_y))[1]
            floor_disp = ax.transData.transform((ref_x, panel_floor_y))[1]
            self.assertGreaterEqual(bottom_disp, floor_disp - 1e-6)
        plt.close(fig)

    def test_panel_floor_cascade_keeps_distinct_label_rows(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 2.0)
        fig.canvas.draw()

        y_vals = list(np.linspace(0.15, 0.45, 8))
        spread = _spread_phenogram_text_labels(
            ax,
            y_vals,
            ref_x=0.55,
            anno_fontsize=9,
            gap_pt=4.0,
            group_label_box_pad_pt=1.5,
            max_iter=100,
            panel_floor_y=0.69,
        )
        self.assertGreater(len(np.unique(np.round(spread, 4))), 1)
        plt.close(fig)

    def test_phenogram_panel_floor_y_uses_row_band_bottom(self):
        row_band_bottom = [2.0, 1.0, 0.0]
        floor = _phenogram_panel_floor_y(
            visual_row=0,
            row_band_bottom=row_band_bottom,
            has_chr_below=True,
        )
        self.assertIsNotNone(floor)
        self.assertGreater(floor, row_band_bottom[0])
        self.assertIsNone(
            _phenogram_panel_floor_y(
                visual_row=2,
                row_band_bottom=row_band_bottom,
                has_chr_below=False,
            )
        )

    def test_phenogram_panel_ceiling_y_uses_row_above(self):
        row_band_bottom = [2.0, 1.0, 0.0]
        self.assertIsNone(
            _phenogram_panel_ceiling_y(
                visual_row=0,
                row_band_bottom=row_band_bottom,
                has_chr_above=False,
            )
        )
        ceiling = _phenogram_panel_ceiling_y(
            visual_row=1,
            row_band_bottom=row_band_bottom,
            has_chr_above=True,
        )
        self.assertIsNotNone(ceiling)
        self.assertLess(ceiling, row_band_bottom[0])

    def test_chr1_row_is_visual_top(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 3.0)
        fig.canvas.draw()

        chr_sizes = {f"chr{i}": float(i * 10_000_000) for i in range(1, 23)}
        chr_list = [f"chr{i}" for i in range(1, 23)]
        max_chr_size = max(chr_sizes.values())
        n_grid_rows, row_ideogram_top, row_band_bottom, _, _, _, _ = (
            _compute_phenogram_row_bands(
                chr_list,
                chr_sizes,
                max_chr_size,
                ncols=11,
                ax=ax,
                ref_x=0.35,
            )
        )
        self.assertEqual(n_grid_rows, 2)
        self.assertGreater(row_ideogram_top[0], row_ideogram_top[1])
        chr1_offset = _chr_panel_bottom_from_row_top(
            row_ideogram_top[0],
            chr_sizes["chr1"] / max_chr_size,
            _TELOMERE_LENGTH,
        )
        chr12_offset = _chr_panel_bottom_from_row_top(
            row_ideogram_top[1],
            chr_sizes["chr12"] / max_chr_size,
            _TELOMERE_LENGTH,
        )
        self.assertGreater(chr1_offset, chr12_offset)
        plt.close(fig)

    def test_row_top_aligns_ideogram_tops(self):
        max_chr_size = 250_000_000.0
        row_top = 2.0
        tel = _TELOMERE_LENGTH
        geom_chr1 = _compute_chr_ideogram_geometry(
            0.0, 0.35, 250_000_000, 110_000_000, 114_000_000,
            max_chr_size,
            _chr_panel_bottom_from_row_top(row_top, 1.0, tel),
            telomere_h=tel,
        )
        geom_chr21 = _compute_chr_ideogram_geometry(
            0.0, 0.35, 48_000_000, 21_000_000, 22_000_000,
            max_chr_size,
            _chr_panel_bottom_from_row_top(row_top, 48_000_000 / max_chr_size, tel),
            telomere_h=tel,
        )
        self.assertAlmostEqual(geom_chr1.ideogram_top, row_top)
        self.assertAlmostEqual(geom_chr21.ideogram_top, row_top)
        self.assertGreater(geom_chr21.y_arm2_bottom, geom_chr1.y_arm2_bottom)

    def test_panel_ceiling_blocks_spread_into_row_above(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 2.0)
        fig.canvas.draw()

        ref_x = 0.55
        above_pt = 12.0
        panel_ceiling_y = 0.72
        y_vals = list(np.linspace(0.78, 0.95, 6))
        spread = _spread_group_blocks_display(
            ax,
            y_vals,
            ref_x=ref_x,
            above_pt=above_pt,
            below_pt=40.0,
            gap_pt=4.0,
            max_iter=50,
            panel_ceiling_y=panel_ceiling_y,
        )
        for y in spread:
            top_y = _block_top_from_center(ax, ref_x, y, above_pt, 40.0)
            top_disp = ax.transData.transform((ref_x, top_y))[1]
            ceiling_disp = ax.transData.transform((ref_x, panel_ceiling_y))[1]
            self.assertLessEqual(top_disp, ceiling_disp + 1e-6)
        plt.close(fig)

    def test_text_spread_uses_per_label_line_count(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 1.0)
        fig.canvas.draw()

        y_vals = [0.3, 0.32]
        uniform = _spread_phenogram_text_labels(
            ax,
            y_vals,
            ref_x=0.55,
            anno_fontsize=9,
            gap_pt=4.0,
            group_label_box_pad_pt=1.5,
            max_iter=50,
            n_lines=2,
        )
        per_label = _spread_phenogram_text_labels(
            ax,
            y_vals,
            ref_x=0.55,
            anno_fontsize=9,
            gap_pt=4.0,
            group_label_box_pad_pt=1.5,
            max_iter=50,
            n_lines_list=[1, 2],
        )
        gap_uniform = abs(uniform[1] - uniform[0])
        gap_per_label = abs(per_label[1] - per_label[0])
        self.assertLess(gap_per_label, gap_uniform)
        plt.close(fig)

    def test_phenogram_inter_unit_gap_uses_marker_label_gap(self):
        gap = _phenogram_inter_unit_gap_pt(
            marker_mode=True,
            marker_label_gap_pt=2.75,
            anno_fontsize=11,
            group_min_vertical_gap_pt=0,
            group_marker_to_marker_gap_pt=0,
            repel_force=1.0,
        )
        self.assertAlmostEqual(gap, 2.75)
        floor_gap = _phenogram_inter_unit_gap_pt(
            marker_mode=True,
            marker_label_gap_pt=2.0,
            anno_fontsize=11,
            group_min_vertical_gap_pt=5.0,
            group_marker_to_marker_gap_pt=0,
            repel_force=1.0,
        )
        self.assertAlmostEqual(floor_gap, 5.0)

    def test_inter_unit_spread_gap_matches_marker_label_gap(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 1.2)
        fig.canvas.draw()

        ref_x = 0.55
        gap_pt = 8.0
        above_pt = 12.0
        below_pt = 28.0
        spread = _spread_group_blocks_display(
            ax,
            [0.35, 0.37],
            ref_x=ref_x,
            above_pt=above_pt,
            below_pt=below_pt,
            gap_pt=gap_pt,
            max_iter=80,
        )
        scale = fig.dpi / 72.0
        order = np.argsort([ax.transData.transform((ref_x, y))[1] for y in spread])
        rank_lo = int(order[0])
        rank_hi = int(order[1])
        d_lo = ax.transData.transform((ref_x, spread[rank_lo]))[1]
        d_hi = ax.transData.transform((ref_x, spread[rank_hi]))[1]
        edge_gap = (d_hi - below_pt * scale) - (d_lo + above_pt * scale)
        self.assertAlmostEqual(edge_gap, gap_pt * scale, delta=scale * 0.5)
        plt.close(fig)

    def test_marker_unit_bottom_extent_for_layout(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 1.2)
        fig.canvas.draw()

        anchor_y = 0.45
        _, below_pt = _marker_group_block_extents_pt(
            81, 2.75, 11, 1.5, n_text_lines=2
        )
        unit_bottom = _unit_bottom_data_y(
            ax,
            0.55,
            anchor_y,
            below_pt + _ANNOTATION_PATH_PAD_PT,
        )
        layout_ymax = unit_bottom + _ANNOTATION_Y_MARGIN_DATA
        self.assertLess(unit_bottom, anchor_y)
        ax.set_ylim(0, layout_ymax)
        fig.canvas.draw()
        ylim_lo, ylim_hi = ax.get_ylim()
        self.assertLessEqual(unit_bottom, max(ylim_lo, ylim_hi) + 1e-6)
        plt.close(fig)

    def test_spread_units_fit_within_ylim(self):
        fig, ax = plt.subplots(figsize=(6, 8))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 1.5)
        fig.canvas.draw()

        ref_x = 0.55
        below_pt = 35.0
        spread = _spread_group_blocks_display(
            ax,
            list(np.linspace(0.85, 1.05, 4)),
            ref_x=ref_x,
            above_pt=12.0,
            below_pt=below_pt,
            gap_pt=2.75,
            max_iter=100,
        )
        ymax = 1.5
        above_pt = 12.0
        for y in spread:
            unit_bottom = _block_bottom_from_center(
                ax, ref_x, y, above_pt, below_pt
            )
            self.assertLessEqual(
                unit_bottom + _ANNOTATION_Y_MARGIN_DATA,
                ymax + 1e-6,
            )
        plt.close(fig)

    def test_default_marker_max_per_row(self):
        self.assertEqual(_DEFAULT_MARKER_MAX_PER_ROW, 4)

    def test_marker_legend_row_handles(self):
        items = _make_marker_shape_legend_handles({"A": "o"})
        row = _marker_legend_row_handles("Shape", items)
        self.assertEqual(len(row), 2)
        self.assertEqual(row[0].get_label(), "Shape")
        self.assertEqual(row[1].get_label(), "A")

    def test_marker_legend_handles_split(self):
        shape_handles = _make_marker_shape_legend_handles({"A": "o", "B": "^"})
        color_handles = _make_marker_color_legend_handles({"X": "red", "Y": "blue"})
        self.assertEqual(len(shape_handles), 2)
        self.assertEqual(len(color_handles), 2)
        self.assertEqual(shape_handles[0].get_marker(), "o")
        self.assertEqual(color_handles[0].get_markerfacecolor(), "red")
        self.assertEqual(shape_handles[0].get_markersize(), _LEGEND_MARKER_SIZE)
        self.assertEqual(
            shape_handles[0].get_markersize(),
            color_handles[0].get_markersize(),
        )

    def test_marker_group_display_label_wraps(self):
        group_items = [{"lead": {"annotation": "Very Long Gene Name Example"}}]
        display, n_lines = _prepare_marker_group_display_label(
            "LOCUS1",
            group_items,
            anno="GENE",
            anno_wrap_chars_per_line=9,
        )
        self.assertIn("\n", display)
        self.assertGreater(n_lines, 1)
        for line in display.split("\n"):
            self.assertLessEqual(len(line), 9)

    def test_wrap_long_label_multiline(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        fp = _phenogram_text_fontproperties({"fontsize": 9})
        wrapped = _wrap_phenogram_label(
            "Very Long Gene Name Example",
            max_width_pt=40.0,
            fontproperties=fp,
            renderer=renderer,
        )
        self.assertIn("\n", wrapped)
        display, n_lines, _ = _prepare_phenogram_display_label(
            "Very Long Gene Name Example",
            anno_wrap=True,
            anno_wrap_width_pt=40.0,
            anno_wrap_chars_per_line=None,
            anno_max_len=None,
            fontproperties=fp,
            renderer=renderer,
        )
        self.assertGreater(n_lines, 1)
        self.assertNotIn("...", display)
        plt.close(fig)

    def test_wrap_label_chars_per_line(self):
        text = "Very Long Gene Name Example"
        wrapped = _wrap_phenogram_label_chars(text, 9)
        for line in wrapped.split("\n"):
            self.assertLessEqual(len(line), 9)
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        fp = _phenogram_text_fontproperties({"fontsize": 9})
        display, n_lines, _ = _prepare_phenogram_display_label(
            text,
            anno_wrap=True,
            anno_wrap_width_pt=None,
            anno_wrap_chars_per_line=_DEFAULT_ANNO_WRAP_CHARS_PER_LINE,
            anno_max_len=None,
            fontproperties=fp,
            renderer=renderer,
        )
        self.assertEqual(display, wrapped)
        self.assertGreater(n_lines, 1)
        for line in display.split("\n"):
            self.assertLessEqual(len(line), 9)
        plt.close(fig)

    def test_xlim_expands_for_wrapped_label(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        fp = _phenogram_text_fontproperties({"fontsize": 9})
        display, _, width_pt = _prepare_phenogram_display_label(
            "Very Long Gene Name Example",
            anno_wrap=True,
            anno_wrap_width_pt=80.0,
            anno_wrap_chars_per_line=None,
            anno_max_len=None,
            fontproperties=fp,
            renderer=renderer,
        )
        _, right_x = _label_width_data_dx(ax, 0.55, 0.4, width_pt, ha="left")
        self.assertGreater(right_x, 0.55)
        plt.close(fig)

    def test_centromere_telomere_style_caps(self):
        geom = _compute_chr_ideogram_geometry(
            0.0, 0.35, 250_000_000, 110_000_000, 114_000_000,
            250_000_000, 0.5, telomere_h=_TELOMERE_LENGTH,
        )
        band = _centromere_fill_patches(geom)
        self.assertEqual(len(band), 1)
        self.assertAlmostEqual(band[0].get_y(), geom.y_cent_bottom)
        self.assertAlmostEqual(band[0].get_height(), geom.centromere_h)
        self.assertGreater(geom.y_cent_mid, geom.y_cent_bottom)
        self.assertLess(geom.y_cent_mid, geom.y_cent_top)

        top_cap = _centromere_cap_path(geom, at_cent_top=True)
        bot_cap = _centromere_cap_path(geom, at_cent_top=False)
        top_ys = [v[1] for v in top_cap.vertices if v[1] != 0.0]
        bot_ys = [v[1] for v in bot_cap.vertices if v[1] != 0.0]
        self.assertAlmostEqual(
            min(top_ys), geom.y_cent_top - geom.telomere_h / 2.0, places=4
        )
        self.assertAlmostEqual(
            max(bot_ys), geom.y_cent_bottom + geom.telomere_h / 2.0, places=4
        )

        outline = _chr_ideogram_outline_path(geom)
        ys = [v[1] for v in outline.vertices]
        self.assertTrue(
            any(abs(y - geom.y_cent_top) < 1e-5 for y in ys),
            "outline should pass through y_cent_top",
        )
        self.assertTrue(
            any(abs(y - geom.y_cent_bottom) < 1e-5 for y in ys),
            "outline should pass through y_cent_bottom",
        )

    def test_bounded_stack_no_overlap(self):
        blocks = [
            _AnnotationLayoutBlock(
                block_id=0, target_y=0.2, height=0.15, y=0.2,
                min_gap=0.02, panel_min=0.0, panel_max=2.0,
            ),
            _AnnotationLayoutBlock(
                block_id=1, target_y=0.25, height=0.15, y=0.25,
                min_gap=0.02, panel_min=0.0, panel_max=2.0,
            ),
            _AnnotationLayoutBlock(
                block_id=2, target_y=0.28, height=0.15, y=0.28,
                min_gap=0.02, panel_min=0.0, panel_max=2.0,
            ),
        ]
        solved, fits = _solve_bounded_stack_layout(blocks, max_iterations=100)
        self.assertTrue(fits)
        solved.sort(key=lambda b: b.y)
        for j in range(1, len(solved)):
            self.assertGreaterEqual(
                _block_bottom(solved[j]),
                _block_top(solved[j - 1]) + solved[j].min_gap - 1e-9,
            )

    def test_bounded_stack_respects_panel(self):
        blocks = [
            _AnnotationLayoutBlock(
                block_id=i,
                target_y=0.6 + i * 0.02,
                height=0.08,
                y=0.6 + i * 0.02,
                min_gap=0.01,
                panel_min=0.55,
                panel_max=0.95,
            )
            for i in range(4)
        ]
        solved, fits = _solve_bounded_stack_layout(blocks, max_iterations=200)
        self.assertTrue(fits)
        for block in solved:
            self.assertGreaterEqual(_block_bottom(block), 0.55 - 1e-9)
            self.assertLessEqual(_block_top(block), 0.95 + 1e-9)

    def test_bounded_stack_minimizes_displacement(self):
        panel_min, panel_max = 0.0, 2.0
        gap = 0.05
        heights = [0.12, 0.12, 0.12]
        targets = [0.4, 0.42, 0.44]

        def _naive_forward(targets_in, heights_in):
            order = np.argsort(targets_in)
            ys = [float(targets_in[i]) for i in order]
            hs = [heights_in[i] for i in order]
            for j in range(1, len(ys)):
                overlap = (
                    ys[j - 1] + hs[j - 1] / 2.0 + gap - (ys[j] - hs[j] / 2.0)
                )
                if overlap > 0:
                    ys[j] += overlap
            return ys

        naive_ys = _naive_forward(targets, heights)
        naive_disp = sum(abs(y - t) for y, t in zip(naive_ys, sorted(targets)))

        blocks = [
            _AnnotationLayoutBlock(
                block_id=i,
                target_y=targets[i],
                height=heights[i],
                y=targets[i],
                min_gap=gap,
                panel_min=panel_min,
                panel_max=panel_max,
            )
            for i in range(3)
        ]
        solved, fits = _solve_bounded_stack_layout(blocks, max_iterations=200)
        self.assertTrue(fits)
        solved_disp = sum(abs(b.y - b.target_y) for b in solved)
        self.assertLessEqual(solved_disp, naive_disp + 1e-9)

    def test_marker_unit_rigid_placement(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(0, 1.2)
        fig.canvas.draw()

        ref_x = 0.55
        above_pt, below_pt = _marker_group_block_extents_pt(
            81, 2.75, 11, 1.5, n_text_lines=1
        )
        height = _pt_to_data_height(ax, ref_x, 0.5, above_pt + below_pt)
        gap_data = _pt_to_data_height(ax, ref_x, 0.5, 8.0)
        blocks = [
            _AnnotationLayoutBlock(
                block_id=0,
                target_y=0.35,
                height=height,
                y=0.35,
                min_gap=gap_data,
                panel_min=0.0,
                panel_max=1.2,
                is_marker=True,
            ),
            _AnnotationLayoutBlock(
                block_id=1,
                target_y=0.37,
                height=height,
                y=0.37,
                min_gap=gap_data,
                panel_min=0.0,
                panel_max=1.2,
                is_marker=True,
            ),
        ]
        solved, fits = _solve_bounded_stack_layout(blocks, max_iterations=100)
        self.assertTrue(fits)
        solved.sort(key=lambda b: b.y)

        group_items = [
            {"lead": {"annotation": "A"}, "marker_shape": "o", "marker_color": "red"},
        ]
        for block in solved:
            placed = _place_marker_unit_at_center(
                ax,
                ref_x,
                block.y,
                group_items,
                marker_size=81,
                marker_gap_pt=4.0,
                marker_label_gap_pt=2.75,
                marker_max_per_row=4,
                marker_row_gap_pt=2.0,
                marker_linewidth=1.0,
                marker_fontsize=11,
                n_text_lines=1,
                group_label_box_pad_pt=1.5,
                marker_label_align="left",
            )
            anchor = _marker_anchor_from_block_center(
                ax,
                ref_x,
                block.y,
                placed["above_pt"],
                placed["below_pt"],
            )
            self.assertAlmostEqual(placed["anchor_y"], anchor)
            unit_top = _block_top_from_center(
                ax, ref_x, block.y, placed["above_pt"], placed["below_pt"]
            )
            unit_bottom = _block_bottom_from_center(
                ax, ref_x, block.y, placed["above_pt"], placed["below_pt"]
            )
            self.assertAlmostEqual(
                unit_top - unit_bottom, block.height, delta=0.01
            )

        for j in range(1, len(solved)):
            self.assertGreaterEqual(
                _block_bottom(solved[j]),
                _block_top(solved[j - 1]) + solved[j].min_gap - 1e-6,
            )
        plt.close(fig)

    def test_connector_uses_bbox_for_wrapped(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.7)
        ax.set_ylim(-0.1, 1.2)
        fig.canvas.draw()
        text = ax.text(
            0.55,
            0.4,
            "Line one\nLine two",
            ha="left",
            va="center",
            fontsize=9,
        )
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        attach_x, attach_y = _text_bbox_attach_point(
            ax, text, renderer, side="left"
        )
        anchor_disp = ax.transData.transform((0.55, 0.4))
        attach_disp = ax.transData.transform((attach_x, attach_y))
        self.assertAlmostEqual(attach_disp[1], anchor_disp[1], delta=3.0)
        plt.close(fig)

    def test_plot_phenogram_smoke(self):
        rows = []
        for chrom in range(1, 23):
            chr_name = "chr{}".format(chrom)
            rows.append("{}\t0\t250000000\tq\tgpos".format(chr_name))
            rows.append("{}\t110000000\t114000000\tacen\tacen".format(chr_name))
        with tempfile.NamedTemporaryFile(suffix=".txt.gz", delete=False) as tmp:
            cytoband_path = tmp.name
        try:
            with gzip.open(cytoband_path, "wt") as gz:
                gz.write("\n".join(rows))

            sumstats = pd.DataFrame(
                {
                    "SNPID": ["rs1"],
                    "CHR": [1],
                    "POS": [1_000_000],
                    "P": [1e-10],
                    "MLOG10P": [10.0],
                }
            )
            fig = _plot_phenogram(
                sumstats,
                cytoband_path=cytoband_path,
                use_lead_extraction=False,
                anno=True,
                ncols=6,
                figsize=(10, 12),
                dpi=100,
                verbose=False,
            )
            self.assertIsNotNone(fig)
            self.assertGreater(len(fig.axes), 0)
            plt.close(fig)
        finally:
            os.unlink(cytoband_path)


if __name__ == "__main__":
    unittest.main()
