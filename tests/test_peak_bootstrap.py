"""Checks of the joint-temperature bootstrap used for peak uncertainty."""
import unittest
import numpy as np
from cu3au.statistics import bootstrap_peak


class PeakBootstrapTests(unittest.TestCase):
    def test_known_peak_and_energy_offset_invariance(self):
        rng = np.random.default_rng(17)
        x = rng.normal(size=(2, 240, 1))
        energies = x*np.array([1., 4., 1.])[None, None, :]
        t = np.array([100., 200., 300.])
        peaks = bootstrap_peak(energies, t, 108, 12, 80)
        np.testing.assert_array_equal(peaks, np.full(80, 200.))
        np.testing.assert_array_equal(
            bootstrap_peak(energies+10000, t, 108, 12, 80), peaks)

    def test_grid_support_reproducibility_and_invalid_parameters(self):
        rng = np.random.default_rng(18)
        e = rng.normal(size=(2, 150, 3))
        t = np.array([900., 1000., 1100.])
        a = bootstrap_peak(e, t, 864, 10, 50)
        np.testing.assert_array_equal(a, bootstrap_peak(e, t, 864, 10, 50))
        self.assertTrue(np.isin(a, t).all())
        with self.assertRaises(ValueError):
            bootstrap_peak(e, t, 864, 151, 50)
        with self.assertRaises(ValueError):
            bootstrap_peak(e, t[::-1], 864, 10, 50)
