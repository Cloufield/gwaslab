import unittest
import numpy as np

from gwaslab.algorithm.core.conversions import betase_to_z, p_to_mlog10p
from gwaslab.algorithm.core.genomic_control import lambda_gc_from_z


class TestAlgorithmCore(unittest.TestCase):
    def test_betase_to_z(self):
        z = betase_to_z(np.array([0.2]), np.array([0.1]))
        self.assertAlmostEqual(float(z[0]), 2.0)

    def test_p_to_mlog10p(self):
        m = p_to_mlog10p(np.array([0.01]))
        self.assertAlmostEqual(float(m[0]), 2.0)

    def test_lambda_gc_from_z_null(self):
        rng = np.random.default_rng(0)
        z = rng.normal(0, 1, 1000)
        lam = lambda_gc_from_z(z)
        self.assertGreater(lam, 0.8)
        self.assertLess(lam, 1.2)


if __name__ == "__main__":
    unittest.main()
