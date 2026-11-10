"""Checkpoints incrémentaux atomiques aux frontières des segments Monte-Carlo.

Chaque fichier contient l'état de reprise et les nouvelles mesures du segment.
Les fichiers .tmp incomplets sont ignorés ; un checkpoint publié corrompu ou
manquant provoque une erreur explicite, jamais une reprise silencieuse au hasard.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import numpy as np
import numba

SCHEMA = 1


def trajectory_signature():
    digest = hashlib.sha256()
    for name in ['kernel.py', 'model.py', 'lattice.py', 'simulation.py', 'checkpoints.py']:
        digest.update((Path(__file__).parent/name).read_bytes())
    digest.update(f'{platform.python_version()}|{np.__version__}|{numba.__version__}'.encode())
    return digest.hexdigest()


class ChainCheckpoints:
    def __init__(self, output, chain, cfg, resume=False):
        self.folder = Path(output)/'checkpoints'/f'chain_{chain+1:02d}'
        self.config = json.dumps(cfg, sort_keys=True)
        self.cfg = cfg
        self.signature = trajectory_signature()
        self.files = sorted(self.folder.glob('segment_*.npz'))
        if self.files and not resume:
            raise FileExistsError('Checkpoints exist: explicit resume=True is required.')
        self.folder.mkdir(parents=True, exist_ok=True)
        self.next_index = len(self.files)

    def validate(self, archive):
        if int(archive['checkpoint_schema']) != SCHEMA:
            raise ValueError('Unsupported checkpoint schema.')
        if str(archive['config_json']) != self.config:
            raise ValueError('Checkpoint configuration mismatch.')
        if str(archive['trajectory_signature']) != self.signature:
            raise ValueError('Simulation code/environment differs from the checkpoint.')

    def restore(self):
        if not self.files:
            return None
        records, burn, traces, captures, capture_sweeps = [], [], [], [], []
        last = None
        previous_mode, previous_done = 0, 0
        for index, path in enumerate(self.files):
            if path.name != f'segment_{index:06d}.npz':
                raise ValueError('A published checkpoint is missing from the sequence.')
            with np.load(path, allow_pickle=False) as a:
                self.validate(a)
                mode, done = int(a['mode_index']), int(a['sweeps_done'])
                limits = [self.cfg['burn_sweeps'], self.cfg['sample_sweeps']]
                if mode not in (0, 1) or not 0 <= done <= limits[mode]:
                    raise ValueError('Checkpoint progress outside configured phase.')
                if index:
                    if mode < previous_mode or mode > previous_mode+1:
                        raise ValueError('Invalid checkpoint phase sequence.')
                    old_done = previous_done if mode == previous_mode else 0
                    if mode != previous_mode and previous_done != limits[previous_mode]:
                        raise ValueError('Checkpoint switched phase before equilibration ended.')
                    expected = (done-old_done)//self.cfg['sample_every']
                    if done <= old_done or len(a['observations']) != expected:
                        raise ValueError('Checkpoint progress and recorded samples disagree.')
                    if len(a['record_labels']) != expected:
                        raise ValueError('Checkpoint label history is incomplete.')
                elif done != 0 or len(a['observations']):
                    raise ValueError('Missing initial checkpoint.')
                if len(a['observations']):
                    if mode == 0:
                        burn.append(a['observations'].copy())
                    else:
                        records.append(a['observations'].copy())
                        traces.append(a['record_labels'].copy())
                if len(a['captured_states']):
                    captures.extend(a['captured_states'].copy())
                    capture_sweeps.extend(a['captured_sweeps'].tolist())
                last = {key: a[key].copy() for key in [
                    'states', 'energies', 'replica_labels', 'exchange_phase',
                    'rng_state_json', 'mode_index', 'sweeps_done',
                    'count_ok', 'count_n', 'local_ok', 'elapsed_seconds']}
                previous_mode, previous_done = mode, done
        last.update(records=records, burn_records=burn, lab_records=traces,
                    spatial_states=captures, spatial_sweeps=capture_sweeps)
        return last

    def save(self, *, states, energies, labels, phase, rng, mode, done,
             count_ok, count_n, local_ok, observations, record_labels,
             captures, capture_sweeps, elapsed):
        path = self.folder/f'segment_{self.next_index:06d}.npz'
        if path.exists():
            raise FileExistsError(f'Checkpoint overwrite refused: {path}')
        temporary = path.with_suffix('.npz.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream, checkpoint_schema=SCHEMA, config_json=self.config,
                trajectory_signature=self.signature, states=states, energies=energies,
                replica_labels=labels, exchange_phase=phase,
                rng_state_json=json.dumps(rng.bit_generator.state), mode_index=mode,
                sweeps_done=done, count_ok=count_ok, count_n=count_n, local_ok=local_ok,
                observations=observations, record_labels=record_labels,
                captured_states=captures, captured_sweeps=capture_sweeps,
                elapsed_seconds=elapsed)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
        self.files.append(path)
        self.next_index += 1
        return path
