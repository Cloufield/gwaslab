import unittest

import numpy as np

from gwaslab.algorithm.core.ess import metal_effective_sample_size
from gwaslab.algorithm.core.simulation import norm_sf, p_from_z, z_to_mlog10p
from gwaslab.algorithm.population.fst import fst_from_allele_frequencies


class TestAlgorithmPhase4(unittest.TestCase):
    def test_metal_effective_sample_size(self):
        n_eff = metal_effective_sample_size(np.array([1000.0]), np.array([1000.0]))
        self.assertAlmostEqual(float(n_eff[0]), 2000.0)

    def test_norm_sf_zero(self):
        self.assertAlmostEqual(float(norm_sf(np.array([0.0]))[0]), 0.5)

    def test_p_from_z(self):
        p = p_from_z(np.array([1.96]))
        self.assertLess(p[0], 0.06)
        self.assertGreater(p[0], 0.04)

    def test_z_to_mlog10p(self):
        m = z_to_mlog10p(np.array([1.96]))
        self.assertGreater(float(m[0]), 0.0)

    def test_fst_bounds(self):
        fst = fst_from_allele_frequencies(0.8, 0.2)
        self.assertGreaterEqual(fst, 0.0)
        self.assertLessEqual(fst, 1.0)


if __name__ == "__main__":
    unittest.main()
