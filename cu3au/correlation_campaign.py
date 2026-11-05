"""Campagne canonique dédiée aux corrélations avec conservation des états.

python -m cu3au.correlation_campaign
Réutilise seulement les chaînes complètes dont la configuration est identique.
"""
import json
import time
from pathlib import Path
import numpy as np
from .simulation import load_config
from .lattice import build_lattice, ordered, random_configuration
from .model import energy_int, SCALE
from .kernel import run_segment


def campaign(config='configs/correlation_L6.json', output='results/correlation_L6'):
    cfg = load_config(config)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'config.json').exists():
        prior = json.loads((out/'config.json').read_text(encoding='utf8'))
        if cfg != prior:
            raise ValueError('Different configuration in existing output directory.')
    (out/'config.json').write_text(json.dumps(cfg, indent=2), encoding='utf8')
    lat = build_lattice(cfg['L'])
    temps = np.asarray(cfg['temperatures'], dtype=float)
    clock = time.perf_counter()
    for c, seed in enumerate(cfg['seeds']):
        dest = out/f'chain_{c+1:02d}.npz'
        if dest.exists():
            with np.load(dest) as archive:
                if str(archive['config_json']) != json.dumps(cfg, sort_keys=True):
                    raise ValueError('Archive configuration mismatch.')
            print(f'Chain {c+1} reused.', flush=True)
            continue
        rng = np.random.default_rng(seed)
        # Réutiliser les extrémités d'une campagne validée pour limiter le transitoire,
        # puis effectuer une nouvelle équilibration avec des graines indépendantes.
        source = Path(f'results/production/chain_{c+1:02d}.npz')
        if source.exists():
            with np.load(source) as old:
                if not np.array_equal(old['xyz'], lat.xyz):
                    raise ValueError('Source geometry mismatch.')
                indices = [int(np.argmin(abs(old['temperatures']-t))) for t in temps]
                states = old['final_states'][indices].copy()
            initialization = str(source)
        else:
            states = np.array([ordered(lat, 2*c%4) if t < 1100 else
                               random_configuration(lat, rng) for t in temps])
            initialization = 'ordered cold / random hot'
        initial = states.copy()
        energies = np.array([energy_int(s, lat.neighbors) for s in states], dtype=np.int64)
        labels = np.arange(len(temps), dtype=np.int32)
        phase = 0
        saved, records, traces = [], [], []
        swap_ok = np.zeros(len(temps)-1, dtype=np.int64)
        swap_attempts = swap_ok.copy()
        for mode, length in [('burn', cfg['burn_sweeps']), ('production', cfg['sample_sweeps'])]:
            stride = cfg['chunk_sweeps'] if mode == 'burn' else cfg['sample_every']
            for done in range(stride, length+1, stride):
                rec, lab, _, ok, attempted, phase = run_segment(
                    states, energies, temps, lat.neighbors, lat.sublattice, stride,
                    cfg['exchange_every'], cfg['sample_every'],
                    int(rng.integers(1, 2**31-1)), labels, phase)
                check = np.array([energy_int(s, lat.neighbors) for s in states])
                if (not np.array_equal(check, energies) or
                        not np.all(states.sum(axis=1) == lat.n//4) or
                        np.any(check < round(-3.813*lat.n*SCALE))):
                    raise RuntimeError('Exact energy/composition check failed.')
                if mode == 'production':
                    saved.append(states.copy())
                    records.append(rec[-1])
                    traces.append(lab[-1])
                    swap_ok += ok
                    swap_attempts += attempted
                if done % 1000 == 0:
                    print(f'Correlation chain {c+1}/{len(cfg["seeds"])} {mode} '
                          f'{done}/{length}; elapsed {time.perf_counter()-clock:.1f}s', flush=True)
        temp = dest.with_suffix('.npz.tmp')
        with temp.open('wb') as stream:
            np.savez_compressed(stream, states=np.array(saved), samples=np.array(records),
                labels=np.array(traces), temperatures=temps, xyz=lat.xyz,
                sublattice=lat.sublattice, initial_states=initial, final_states=states,
                swap_ok=swap_ok, swap_attempts=swap_attempts, seed=seed,
                rng_state_json=json.dumps(rng.bit_generator.state), initialization=initialization,
                config_json=json.dumps(cfg, sort_keys=True))
        temp.replace(dest)
    print('Correlation configuration archives complete.', flush=True)


if __name__ == '__main__':
    campaign()
