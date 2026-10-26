"""Tests physiques indépendants du résultat de la campagne."""
import itertools
import unittest
import numpy as np
from cu3au.lattice import build_lattice, ordered, random_configuration
from cu3au.model import (energy, energy_int, delta_reference, cowley, ordre,
                         random_energy_limit, KB, SCALE)
from cu3au.kernel import delta_fast, run_segment
from cu3au.meanfield import coexistence


class PhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lat = build_lattice(6)

    def test_geometry_and_periodic_distances(self):
        lat = self.lat
        self.assertEqual(lat.n, 864)
        self.assertTrue(np.all(lat.xyz.sum(axis=1) % 2 == 0))
        for i, row in enumerate(lat.neighbors):
            self.assertEqual(len(set(row)), 18)
            d = (lat.xyz[row]-lat.xyz[i]+lat.L) % (2*lat.L)-lat.L
            np.testing.assert_array_equal((d[:12]**2).sum(axis=1), np.full(12, 2))
            np.testing.assert_array_equal((d[12:]**2).sum(axis=1), np.full(6, 4))
            for j in row[:12]:
                self.assertIn(i, lat.neighbors[j, :12])
            for j in row[12:]:
                self.assertIn(i, lat.neighbors[j, 12:])

    def test_ordered_exact_energy_and_cowley(self):
        s = ordered(self.lat)
        self.assertEqual(int(s.sum()), 216)
        self.assertEqual(energy_int(s, self.lat.neighbors), -3294432000)
        np.testing.assert_allclose(cowley(s, self.lat.neighbors), [-1/3, 1], atol=1e-14)
        eta, q, p = ordre(s, self.lat.sublattice)
        np.testing.assert_allclose([eta, q], [1, 1])

    def test_four_equivalent_variants(self):
        for v in range(4):
            s = ordered(self.lat, v)
            self.assertAlmostEqual(energy(s, self.lat.neighbors), -3294.432)
            eta, q, _ = ordre(s, self.lat.sublattice)
            self.assertAlmostEqual(q, 1)
            self.assertAlmostEqual(eta, 1 if v == 0 else -1/3)

    def test_nearest_neighbor_defect(self):
        s = ordered(self.lat)
        i = int(np.flatnonzero(s)[0])
        j = int(self.lat.neighbors[i, 0])
        self.assertEqual(delta_reference(s, i, j, self.lat.neighbors), 558600)
        self.assertEqual(delta_fast(s, i, j, self.lat.neighbors), 558600)

    def test_random_swaps_against_full_energy(self):
        rng = np.random.default_rng(99)
        for L in [3, 6]:
            lat = build_lattice(L)
            s = random_configuration(lat, rng)
            for k in range(250):
                # Inclut explicitement des voisins 1, 2 et des sites quelconques.
                i = int(rng.integers(lat.n))
                j = int(lat.neighbors[i, k % 18] if k % 3 else rng.integers(lat.n))
                before = energy_int(s, lat.neighbors)
                fast = int(delta_fast(s, i, j, lat.neighbors))
                self.assertEqual(fast, delta_reference(s, i, j, lat.neighbors))
                s[i], s[j] = s[j], s[i]
                self.assertEqual(energy_int(s, lat.neighbors)-before, fast)
                self.assertEqual(delta_fast(s, i, j, lat.neighbors), -fast)
                self.assertEqual(s.sum(), lat.n//4)

    def test_ground_state_bound_identity(self):
        rng = np.random.default_rng(52)
        for _ in range(25):
            s = random_configuration(self.lat, rng)
            aa = int(np.sum(s[:, None]*s[self.lat.neighbors[:, :12]]))//2
            ac = int(np.sum((s[:, None] == 1) &
                            (s[self.lat.neighbors[:, 12:]] == 0)))
            self.assertEqual(energy_int(s, self.lat.neighbors)+3294432000,
                             161000*aa+6300*ac)

    def test_exact_restricted_canonical_detailed_balance(self):
        # Toutes les configurations d'un ensemble restreint : 1 Au mobile
        # parmi quatre sites, tous les autres sites figés. Aucun MC ici.
        lat = build_lattice(3)
        base = ordered(lat)
        active = [0, int(lat.neighbors[0, 0]), int(lat.neighbors[0, 12]), 13]
        active = list(dict.fromkeys(active))
        self.assertEqual(len(active), 4)
        base[active] = 0
        states = []
        for site in active:
            s = base.copy()
            s[site] = 1
            states.append(s)
        e = np.array([energy(s, lat.neighbors) for s in states])
        T = 1200
        p = np.exp(-(e-e.min())/(KB*T)); p /= p.sum()
        P = np.zeros((4, 4))
        for i, j in itertools.permutations(range(4), 2):
            P[i, j] = min(1., np.exp(-(e[j]-e[i])/(KB*T)))/3
        np.fill_diagonal(P, 1-P.sum(axis=1))
        np.testing.assert_allclose(p@P, p, atol=1e-13)
        np.testing.assert_allclose(p[:, None]*P, (p[:, None]*P).T, atol=1e-13)

    def test_replica_exchange_sign_and_balance(self):
        for E1, E2 in [(-100., -90.), (-90., -100.), (-90., -90.)]:
            b1, b2 = 1/(KB*900), 1/(KB*1600)
            loga = (b1-b2)*(E1-E2)
            forward, reverse = min(0., loga), min(0., -loga)
            self.assertAlmostEqual(forward-reverse, loga)
            self.assertAlmostEqual(loga, (-b1*E2-b2*E1)-(-b1*E1-b2*E2))

    def test_short_actual_kernel_and_reproducibility(self):
        lat = build_lattice(3)
        def run():
            states = np.array([ordered(lat, v) for v in range(3)])
            energies = np.array([energy_int(s, lat.neighbors) for s in states])
            labels = np.arange(3, dtype=np.int32)
            res = run_segment(states, energies, np.array([500., 1400., 4000.]),
                              lat.neighbors, lat.sublattice, 50, 5, 5, 123, labels)
            return states, energies, res
        a, b = run(), run()
        np.testing.assert_array_equal(a[0], b[0])
        np.testing.assert_array_equal(a[2][0], b[2][0])
        for s, e in zip(a[0], a[1]):
            self.assertEqual(energy_int(s, lat.neighbors), e)
            self.assertEqual(s.sum(), lat.n//4)

    def test_high_temperature_limit(self):
        lat = build_lattice(3)
        rng = np.random.default_rng(18)
        s = np.array([random_configuration(lat, rng)])
        e = np.array([energy_int(s[0], lat.neighbors)])
        records = run_segment(s, e, np.array([1e12]), lat.neighbors, lat.sublattice,
                              3000, 0, 1, 1234, np.array([0], dtype=np.int32))[0][:, 0]
        self.assertLess(abs(records[:, 0].mean()/lat.n-random_energy_limit(lat.n)), .001)
        np.testing.assert_allclose(records[:, 3:5].mean(axis=0), -1/(lat.n-1), atol=.012)
        self.assertLess(abs(records[:, 1].mean()), .015)

    def test_meanfield_reference(self):
        r = coexistence()
        self.assertAlmostEqual(r['T_MF_K'], 1712.17574575, places=5)
        self.assertGreater(r['T_MF_K'], r['T_spinodal_K'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
