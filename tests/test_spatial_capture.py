"""Les états supplémentaires correspondent aux mesures, même aux bornes de segment."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from cu3au.simulation import campaign, load_config
from cu3au.model import observables
from cu3au.lattice import build_lattice


class SpatialCaptureTests(unittest.TestCase):
    def test_captures_match_measured_observables(self):
        cfg = dict(L=3, temperatures=[300, 1100, 4000], seeds=[61, 62],
                   burn_sweeps=10, sample_sweeps=30, sample_every=2,
                   exchange_every=2, chunk_sweeps=10, store_states_every=6, bootstrap=12)
        progress = []
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            campaign(cfg, tmp, progress=lambda *v: progress.append(v))
            lat = build_lattice(3)
            for path in sorted(Path(tmp).glob('chain_*.npz')):
                with np.load(path) as a:
                    np.testing.assert_array_equal(a['spatial_sample_sweeps'], [6, 12, 18, 24, 30])
                    self.assertEqual(a['spatial_states'].shape, (5, 3, 108))
                    for i, sweep in enumerate(a['spatial_sample_sweeps']):
                        for t in range(3):
                            reference = observables(a['spatial_states'][i, t], lat)
                            np.testing.assert_allclose(reference, a['samples'][sweep//2-1, t], atol=1e-11)
                    np.testing.assert_array_equal(a['spatial_states'][-1], a['final_states'])
        self.assertEqual(progress[-1][2:], (30, 30, 5))

    def test_capture_config_rejects_incompatible_intervals(self):
        cfg = dict(L=3, temperatures=[300, 1100, 4000], seeds=[61, 62],
                   burn_sweeps=10, sample_sweeps=30, sample_every=2,
                   exchange_every=2, chunk_sweeps=10)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'config.json'
            for invalid in [-1, True, 3, 8]:
                path.write_text(json.dumps(dict(cfg, store_states_every=invalid)))
                with self.assertRaises(ValueError):
                    load_config(path)
            path.write_text(json.dumps(dict(cfg, store_states_every=6)))
            self.assertEqual(load_config(path)['store_states_every'], 6)
