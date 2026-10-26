import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from cu3au.statistics import tau_int, split_rhat, moments
from cu3au.simulation import load_config
from cu3au.model import KB


class StatisticsTests(unittest.TestCase):
    def test_autocorrelation_known_ar1(self):
        rng = np.random.default_rng(33)
        x = np.zeros(100000)
        noise = rng.normal(size=len(x))
        for i in range(1, len(x)):
            x[i] = .8*x[i-1]+noise[i]
        # Tau exact de AR(1), convention 1/2+sum rho = 4.5.
        self.assertLess(abs(tau_int(x)-4.5), .8)
        self.assertEqual(tau_int(np.ones(100)), .5)

    def test_rhat_detects_separated_chains(self):
        rng = np.random.default_rng(32)
        x = rng.normal(size=(3, 4000))
        self.assertLess(split_rhat(x), 1.01)
        x[0] += 4
        self.assertGreater(split_rhat(x), 1.2)

    def test_total_energy_heat_capacity_units(self):
        a = np.zeros((100, 9))
        a[:, 0] = np.tile([-1., 1.], 50)
        a[:, 2] = .5
        T, n = 1000, 864
        expected = 1/(n*(KB*T)**2)
        self.assertAlmostEqual(moments(a, T, n)[5], expected)

    def test_invalid_temperature_grid_rejected(self):
        base = Path(__file__).resolve().parents[1]
        cfg = json.loads((base/'configs/quick.json').read_text())
        cfg['temperatures'] = [1000, 900, 1200]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'bad.json'
            path.write_text(json.dumps(cfg))
            with self.assertRaises(ValueError):
                load_config(path)


if __name__ == '__main__':
    unittest.main(verbosity=2)
