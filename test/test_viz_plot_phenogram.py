import os
import sys
import unittest

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

import pandas as pd

from gwaslab.viz.viz_plot_phenogram import (
    _build_lead_text_kwargs,
    _chr_ideogram_top_y,
    _chr_number_label_position,
    _compute_group_marker_positions,
    _cap_phenogram_arrow_shaft_pt,
    _effective_marker_radius_pt,
    _label_eligible_mask,
    _marker_center_step_pt,
    _marker_group_block_extents_pt,
    _marker_radius_pt,
    _plot_phenogram_text_connector,
    _points_to_data_delta,
    _resolve_lead_connector_style,
    _resolve_phenogram_label,
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
        ax.invert_yaxis()
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
        ax.invert_yaxis()
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
        ax.invert_yaxis()
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
        ax.invert_yaxis()
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
        self.assertAlmostEqual(_chr_ideogram_top_y(0.6, 0.02), 0.59)

    def test_chr_number_label_above_ideogram_inverted_y(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0, 1.5)
        ax.set_ylim(-0.2, 1.2)
        ax.invert_yaxis()
        fig.canvas.draw()

        body_top_y = _chr_ideogram_top_y(0.5, 0.02)
        label_y, label_va = _chr_number_label_position(
            ax, body_top_y, ref_x=0.35, fontsize=10, pad_pt=4
        )
        self.assertEqual(label_va, "top")
        body_disp = ax.transData.transform((0.35, body_top_y))
        label_disp = ax.transData.transform((0.35, label_y))
        self.assertLess(label_disp[1], body_disp[1])

        t = ax.text(
            0.35, label_y, "22", ha="center", va=label_va, fontsize=10
        )
        fig.canvas.draw()
        bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
        self.assertLess(bb.y1, body_disp[1])
        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
