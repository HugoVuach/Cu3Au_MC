import unittest
import numpy as np
from cu3au.lattice import build_lattice, ordered, random_configuration
from cu3au.model import ordre
from cu3au.correlation import (local_order, order_fft, spectra, radial_correlation,
                                second_moment_length, SIGNS, within_variant_spectrum)


class CorrelationTests(unittest.TestCase):
    def test_four_variants_and_order_normalization(self):
        lat = build_lattice(3)
        for variant in range(4):
            s = ordered(lat, variant)
            psi = local_order(s, lat.L)
            np.testing.assert_allclose(psi, np.broadcast_to(SIGNS[:, variant], psi.shape))
            raw, connected, _ = spectra(np.array([s, s]), lat.L)
            self.assertAlmostEqual(raw[0, 0, 0], lat.L**3)
            self.assertLess(np.max(connected), 1e-12)

    def test_fft_matches_direct_sum_and_real_space_pairs(self):
        lat = build_lattice(3)
        rng = np.random.default_rng(7)
        states = np.array([random_configuration(lat, rng) for _ in range(12)])
        field = local_order(states, lat.L)
        transform = order_fft(states, lat.L)
        xyz = np.stack(np.meshgrid(*([np.arange(lat.L)]*3), indexing='ij'), axis=-1)
        direct = np.sum(field*np.exp(-2j*np.pi*xyz[..., 0]/lat.L)[None, ..., None], axis=(1, 2, 3))
        np.testing.assert_allclose(transform[:, 1, 0, 0], direct, atol=1e-12)
        raw, connected, _ = spectra(states, lat.L)
        inverse = np.fft.ifftn(connected).real
        centered = field-field.mean(axis=0)
        pair = np.mean(np.sum(centered*np.roll(centered, -1, axis=1), axis=-1))/3
        self.assertAlmostEqual(inverse[1, 0, 0], pair)
        q2 = np.array([ordre(s, lat.sublattice)[1]**2 for s in states]).mean()
        self.assertAlmostEqual(raw[0, 0, 0], lat.L**3*q2)

    def test_length_known_spectrum_and_undefined_cases(self):
        L, xi = 6, 1.7
        s0 = 10.
        s1 = s0/(1+(2*np.sin(np.pi/L)*xi)**2)
        self.assertAlmostEqual(second_moment_length(s0, s1, L), xi)
        self.assertTrue(np.isnan(second_moment_length(0, 0, L)))
        self.assertTrue(np.isnan(second_moment_length(1, 2, L)))
        radius, values = radial_correlation(np.ones((L, L, L)))
        self.assertEqual(radius[0], 0)
        self.assertAlmostEqual(values[0], 1)
        np.testing.assert_allclose(values[1:], 0)

    def test_variant_mixture_is_not_intra_variant_fluctuation(self):
        lat = build_lattice(3)
        states = np.array([ordered(lat, v) for v in range(4) for _ in range(10)])
        raw, connected, _ = spectra(states, lat.L)
        self.assertAlmostEqual(connected[0, 0, 0], lat.L**3)
        conditional, coverage = within_variant_spectrum(states, lat.L)
        self.assertEqual(coverage, 1)
        np.testing.assert_allclose(conditional, 0, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
