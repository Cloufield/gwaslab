import os
import sys
import unittest

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from gwaslab.extension.ldsc.ldsc_regressions import Hsq, _as_python_float
import gwaslab.extension.ldsc.ldsc_jackknife as jk


def _numpy_at_least_2_4():
    return np.lib.NumpyVersion(np.__version__) >= "2.4.0"


class TestLdscNumpyScalar(unittest.TestCase):
    def test_float_on_1x1_array_is_the_numpy_2_4_failure_mode(self):
        arr = np.ones((1, 1))
        if _numpy_at_least_2_4():
            with self.assertRaises(TypeError):
                float(arr)
        self.assertIsInstance(_as_python_float(arr), float)
        self.assertEqual(_as_python_float(arr), 1.0)

    def test_as_python_float_shapes(self):
        cases = [
            0.5,
            np.float64(0.5),
            np.array(0.5),
            np.array([0.5]),
            np.array([[0.5]]),
        ]
        for value in cases:
            with self.subTest(shape=np.asarray(value).shape, dtype=type(value)):
                out = _as_python_float(value)
                self.assertIsInstance(out, float)
                self.assertAlmostEqual(out, 0.5)

    def test_ratio_jackknife_scalars(self):
        rg_ratio = np.array([[0.5]])
        numer = np.array([[0.40], [0.55], [0.70]])
        denom = np.array([[0.90], [1.00], [0.80]])
        rg = jk.RatioJackknife(rg_ratio, numer, denom)
        self.assertEqual(rg.jknife_est.shape[0], 1)
        jknife = _as_python_float(rg.jknife_est)
        se = _as_python_float(rg.jknife_se)
        ratio = _as_python_float(rg_ratio)
        self.assertIsInstance(jknife, float)
        self.assertIsInstance(se, float)
        self.assertIsInstance(ratio, float)
        self.assertGreater(se, 0.0)
        self.assertAlmostEqual(ratio, 0.5)

    def test_hsq_weights_accepts_1x1_M(self):
        n_snp = 4
        ld = np.ones((n_snp, 1))
        w_ld = np.ones((n_snp, 1))
        N = np.full((n_snp, 1), 1000.0)
        M = np.array([[1e4]])
        w = Hsq.weights(ld, w_ld, N, M, hsq=0.2)
        self.assertEqual(w.shape, (n_snp, 1))
        self.assertTrue(np.all(np.isfinite(w)))


if __name__ == "__main__":
    unittest.main()
