"""Physical continuity, complete measurements, and exact plateau restart."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from cu3au.thermal_cycle import (run_cycle, schedule, stage_path,
                                 load_cycle_config, analyze_cycle, _save_stage, write_status)


def small_config():
    return {'L': 3, 'temperatures': [300, 1400, 4000], 'seeds': [71, 72],
            'burn_sweeps': 10, 'sample_sweeps': 20, 'sample_every': 2,
            'exchange_every': 0, 'chunk_sweeps': 10, 'bootstrap': 12,
            'cycle': 'heating_then_cooling'}


class ThermalCycleTests(unittest.TestCase):
    def test_a_locked_monitoring_file_does_not_abort_the_cycle(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch.object(Path, 'replace', side_effect=PermissionError('locked')), \
                    patch('cu3au.thermal_cycle.time.sleep'):
                self.assertFalse(write_status(temp, 'simulating'))
            self.assertTrue(write_status(temp, 'simulating'))
            data = json.loads((Path(temp)/'execution.json').read_text(encoding='utf8'))
            self.assertEqual(data['state'], 'simulating')

    def test_schedule_and_invalid_replica_exchange(self):
        cfg = small_config()
        self.assertEqual([v[2] for v in schedule(cfg)], [300, 1400, 4000, 4000, 1400, 300])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'cfg.json'
            path.write_text(json.dumps(cfg), encoding='utf8')
            self.assertEqual(load_cycle_config(path)['L'], 3)
            cfg['exchange_every'] = 2
            path.write_text(json.dumps(cfg), encoding='utf8')
            with self.assertRaises(ValueError):
                load_cycle_config(path)

    def test_complete_cycle_continuity_and_observable_archives(self):
        with tempfile.TemporaryDirectory() as temp:
            result = run_cycle(small_config(), temp)
            self.assertEqual(result['plateaus'], 12)
            for c in range(2):
                previous = None
                for step in range(6):
                    with np.load(stage_path(temp, c, step)) as a:
                        self.assertEqual(a['samples'].shape, (10, 9))
                        self.assertEqual(a['burn_samples'].shape, (5, 9))
                        self.assertEqual(int(a['final_state'].sum()), 27)
                        if previous is not None:
                            np.testing.assert_array_equal(a['initial_state'], previous)
                        previous = a['final_state'].copy()
            self.assertTrue((Path(temp)/'paired_differences.csv').exists())

    def test_restart_replays_the_identical_trajectory(self):
        cfg = small_config()
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            run_cycle(cfg, first)
            def interrupt_after_two(c, step, path):
                if c == 0 and step == 1:
                    raise RuntimeError('Intentional interruption after a saved plateau.')
            with self.assertRaises(RuntimeError):
                run_cycle(cfg, second, progress=interrupt_after_two)
            run_cycle(cfg, second, resume=True)
            for c in range(2):
                for step in range(6):
                    with np.load(stage_path(first, c, step)) as a, np.load(stage_path(second, c, step)) as b:
                        np.testing.assert_array_equal(a['samples'], b['samples'])
                        np.testing.assert_array_equal(a['final_state'], b['final_state'])
            cfg['sample_sweeps'] = 40
            with self.assertRaises(ValueError):
                run_cycle(cfg, second, resume=True)

    def test_validation_rejects_a_broken_handoff(self):
        with tempfile.TemporaryDirectory() as temp:
            run_cycle(small_config(), temp)
            path = stage_path(temp, 0, 3)
            with np.load(path) as a:
                arrays = {key: a[key] for key in a.files}
            arrays['initial_state'] = 1-arrays['initial_state']
            _save_stage(path, **arrays)
            with self.assertRaises(AssertionError):
                analyze_cycle(temp)
