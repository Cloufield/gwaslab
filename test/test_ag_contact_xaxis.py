"""X-axis alignment tests for ag_contact stacked panels."""

import os
import sys
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
FIXTURES_DIR = os.path.join(ROOT, "test", "fixtures")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
if FIXTURES_DIR not in sys.path:
    sys.path.insert(0, FIXTURES_DIR)

import gwaslab as gl
from gwaslab.viz.viz_plot_alphagenome import plot_ag_contact
from panel_demo_fixtures import setup_demo


class TestAgContactXaxis(unittest.TestCase):
    def setUp(self):
        self.demo = setup_demo()
        self.region = self.demo["region"]

    def tearDown(self):
        plt.close("all")

    def test_plot_ag_contact_limits_match_region(self):
        fig, ax = plt.subplots()
        plot_ag_contact(
            self.demo["ag_contact"], [ax], region=self.region, verbose=False,
        )
        self.assertAlmostEqual(ax.get_xlim()[0], self.region[1], places=3)
        self.assertAlmostEqual(ax.get_xlim()[1], self.region[2], places=3)
        self.assertEqual(ax.get_xlim(), ax.get_ylim())

    def test_ag_contact_aligned_in_plot_panels(self):
        panels = [
            gl.Panel(
                "track",
                track_path=self.demo["gtf"],
                region=self.region,
                feature_type="gene",
                verbose=False,
            ),
            gl.Panel(
                "ag_contact",
                bundle=self.demo["ag_contact"],
                region=self.region,
                verbose=False,
            ),
        ]
        fig, _ = gl.plot_panels(
            panels,
            region=self.region,
            fig_kwargs={"figsize": (10, 6), "dpi": 80},
            verbose=False,
        )
        track_ax, contact_ax = fig.axes[0], fig.axes[1]
        self.assertAlmostEqual(track_ax.get_xlim()[0], contact_ax.get_xlim()[0], places=3)
        self.assertAlmostEqual(track_ax.get_xlim()[1], contact_ax.get_xlim()[1], places=3)
        self.assertEqual(contact_ax.get_xlim(), contact_ax.get_ylim())


if __name__ == "__main__":
    unittest.main()
