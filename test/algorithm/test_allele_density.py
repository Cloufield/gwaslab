import unittest

import numpy as np

from gwaslab.algorithm.allele.complement import reverse_complement
from gwaslab.algorithm.density.signal import density_all_variants


class TestAlgorithmAlleleDensity(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AC"), "GT")

    def test_density_all_variants(self):
        pos = np.array([1000, 1100, 5000], dtype=np.int64)
        counts = density_all_variants(pos, window_bp=200)
        self.assertEqual(list(counts), [1, 1, 0])


if __name__ == "__main__":
    unittest.main()
