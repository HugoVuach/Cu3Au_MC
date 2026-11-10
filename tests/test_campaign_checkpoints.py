import contextlib
import io
from pathlib import Path
import tempfile
import unittest
import numpy as np
from cu3au.simulation import campaign


def config():
    return dict(L=3, temperatures=[300, 1100, 4000], seeds=[71, 72],
                burn_sweeps=12, sample_sweeps=24, sample_every=2,
                exchange_every=2, chunk_sweeps=6, store_states_every=4, bootstrap=12)


class CampaignCheckpointTests(unittest.TestCase):
    def test_restart_matches_uninterrupted_in_both_phases_and_second_chain(self):
        cfg = config()
        for interrupted_chain, interrupted_phase, interrupted_sweep in [
                (1, 'equilibration', 6), (1, 'mesure', 8),
                (1, 'mesure', 24), (2, 'mesure', 8)]:
            with self.subTest(chain=interrupted_chain, phase=interrupted_phase), \
                    tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
                reference, resumed = Path(temp)/'reference', Path(temp)/'resumed'
                expected = campaign(cfg, reference)

                def interrupt(chain, phase, done, total, captures):
                    if (chain, phase, done) == (interrupted_chain, interrupted_phase, interrupted_sweep):
                        raise InterruptedError('simulated external stop')

                with self.assertRaises(InterruptedError):
                    campaign(cfg, resumed, checkpoint=True, progress=interrupt)
                # Un fichier temporaire laissé par une écriture interrompue n'est pas publié.
                folder = resumed/'checkpoints'/f'chain_{interrupted_chain:02d}'
                next_index = len(list(folder.glob('segment_*.npz')))
                (folder/f'segment_{next_index:06d}.npz.tmp').write_bytes(b'incomplete')
                actual = campaign(cfg, resumed, checkpoint=True, resume=True)
                self.assertEqual(expected['chains'], actual['chains'])
                self.assertEqual(expected['T_cv_peak_grid_K'], actual['T_cv_peak_grid_K'])
                for chain in [1, 2]:
                    with np.load(reference/f'chain_{chain:02d}.npz') as left, \
                            np.load(resumed/f'chain_{chain:02d}.npz') as right:
                        for key in left.files:
                            np.testing.assert_array_equal(left[key], right[key], err_msg=key)
                self.assertEqual((reference/'observables.csv').read_bytes(),
                                 (resumed/'observables.csv').read_bytes())

    def test_mismatch_gap_and_implicit_overwrite_are_rejected(self):
        cfg = config()
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            out = Path(temp)

            def interrupt(chain, phase, done, total, captures):
                if phase == 'mesure' and done == 8:
                    raise InterruptedError()

            with self.assertRaises(InterruptedError):
                campaign(cfg, out, checkpoint=True, progress=interrupt)
            with self.assertRaises(FileExistsError):
                campaign(cfg, out, checkpoint=True)
            with self.assertRaisesRegex(ValueError, 'configuration mismatch'):
                campaign(dict(cfg, sample_sweeps=48), out, checkpoint=True, resume=True)
            # Une sauvegarde publiée manquante n'est pas ignorée.
            (out/'checkpoints/chain_01/segment_000001.npz').rename(
                out/'checkpoints/chain_01/segment_000001.removed')
            with self.assertRaisesRegex(ValueError, 'missing'):
                campaign(cfg, out, checkpoint=True, resume=True)

    def test_completed_chain_is_reused_after_process_stops(self):
        cfg = config()
        with tempfile.TemporaryDirectory() as temp, contextlib.redirect_stdout(io.StringIO()):
            out = Path(temp)
            first = campaign(cfg, out, checkpoint=True)
            archives = {p.name: p.read_bytes() for p in out.glob('chain_*.npz')}
            second = campaign(cfg, out, checkpoint=True, resume=True)
            self.assertEqual(first['chains'], second['chains'])
            self.assertEqual(archives, {p.name: p.read_bytes() for p in out.glob('chain_*.npz')})
