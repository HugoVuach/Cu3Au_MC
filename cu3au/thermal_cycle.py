"""Sequential Metropolis heating/cooling, with a restart after every plateau.

There is no replica exchange: the last state at T becomes the initial state
at the next T. Cooling begins from the actual final hot configuration.
"""
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import time
import numpy as np
from .kernel import run_segment
from .lattice import build_lattice, ordered
from .model import energy_int, SCALE
from .simulation import load_config, save_csv
from .statistics import summarize


def load_cycle_config(path):
    cfg = load_config(path)
    if cfg['exchange_every'] != 0:
        raise ValueError('A sequential thermal cycle requires exchange_every=0.')
    if cfg['sample_sweeps']//cfg['sample_every'] < 8:
        raise ValueError('At least eight production measurements per plateau.')
    if cfg.get('cycle', 'heating_then_cooling') != 'heating_then_cooling':
        raise ValueError('Unsupported cycle order.')
    return cfg


def schedule(cfg):
    t = cfg['temperatures']
    # The upper temperature is sampled again on the return branch, without
    # reinitialization: it is a second hot dwell, not a duplicated sample.
    return [('heating', r, float(T)) for r, T in enumerate(t)]+[
        ('cooling', r, float(t[r])) for r in range(len(t)-1, -1, -1)]


def config_fingerprint(cfg):
    return hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()


def write_status(output, state, **details):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    data = {'state': state, 'pid': os.getpid(),
            'updated_at': datetime.now().astimezone().isoformat(), **details}
    temporary = out/f'execution_{os.getpid()}.tmp'
    # Windows readers/antivirus can temporarily deny replacing an open file.
    # Monitoring is best effort: it must never stop the physical trajectory.
    for attempt in range(10):
        try:
            temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding='utf8')
            temporary.replace(out/'execution.json')
            return True
        except PermissionError:
            if attempt < 9:
                time.sleep(.2)
    return False


def stage_path(output, chain, step):
    return Path(output)/f'chain_{chain+1:02d}'/f'step_{step:03d}.npz'


def _save_stage(path, **arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    with temporary.open('wb') as f:
        np.savez_compressed(f, **arrays)
    temporary.replace(path)


def run_cycle(cfg, output, resume=False, progress=None):
    """Save every completed plateau atomically; partial plateaus are replayed.

    The random-generator state is saved at each plateau, so a resumed run
    follows the same trajectory as an uninterrupted run on the same runtime.
    """
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    fingerprint = config_fingerprint(cfg)
    existing = sorted(out.glob('chain_*/step_*.npz'))
    if existing and not resume:
        raise FileExistsError('Existing plateaus; use --resume or a new output directory.')
    config_path = out/'config.json'
    if config_path.exists():
        previous = json.loads(config_path.read_text(encoding='utf8'))
        if config_fingerprint(previous) != fingerprint:
            raise ValueError('Configuration differs from saved cycle.')
    config_path.write_text(json.dumps(cfg, indent=2), encoding='utf8')
    lat = build_lattice(cfg['L'])
    steps = schedule(cfg)
    expected = {stage_path(out, c, j) for c in range(len(cfg['seeds']))
                for j in range(len(steps))}
    if any(p not in expected for p in existing):
        raise ValueError('Unexpected plateau archives in output directory.')
    start = time.perf_counter()
    completed = 0
    for c, seed in enumerate(cfg['seeds']):
        rng = np.random.default_rng(seed)
        states = ordered(lat, (2*c) % 4)[None, :]
        energies = np.array([energy_int(states[0], lat.neighbors)], dtype=np.int64)
        found_gap = False
        for step, (branch, index, temperature) in enumerate(steps):
            path = stage_path(out, c, step)
            if path.exists():
                if found_gap:
                    raise ValueError('Missing plateau before an existing archive.')
                with np.load(path) as a:
                    if str(a['config_fingerprint']) != fingerprint:
                        raise ValueError('Plateau configuration fingerprint mismatch.')
                    np.testing.assert_array_equal(a['initial_state'], states[0])
                    if int(a['temperature_index']) != index or float(a['temperature_K']) != temperature:
                        raise ValueError('Plateau schedule mismatch.')
                    states[0] = a['final_state']
                    energies[0] = int(a['final_energy_int'])
                    rng.bit_generator.state = json.loads(str(a['rng_state_json']))
                completed += 1
                continue
            found_gap = True
            initial_state = states[0].copy()
            initial_energy = energies[0]
            records = {'equilibration': [], 'production': []}
            accepted = {'equilibration': 0, 'production': 0}
            labels = np.zeros(1, dtype=np.int32)
            for mode, length in [('equilibration', cfg['burn_sweeps']),
                                 ('production', cfg['sample_sweeps'])]:
                done = 0
                while done < length:
                    chunk = min(cfg['chunk_sweeps'], length-done)
                    rec, _, acc, _, _, _ = run_segment(
                        states, energies, np.array([temperature]), lat.neighbors,
                        lat.sublattice, chunk, 0, cfg['sample_every'],
                        int(rng.integers(1, 2**31-1)), labels)
                    check = energy_int(states[0], lat.neighbors)
                    if (check != energies[0] or states[0].sum() != lat.n//4 or
                            check < round(-3.813*lat.n*SCALE)):
                        raise RuntimeError('Energy, composition, or lower-bound check failed.')
                    records[mode].append(rec[:, 0, :])
                    accepted[mode] += int(acc[0])
                    done += chunk
                    write_status(out, 'simulating', L=cfg['L'], n_atoms=lat.n,
                        chain=c+1, chains=len(cfg['seeds']), branch=branch,
                        temperature_K=temperature, step=step+1, steps_per_chain=len(steps),
                        phase=mode, sweeps_done=done, sweeps_target=length,
                        completed_plateaus=completed,
                        total_plateaus=len(steps)*len(cfg['seeds']))
            _save_stage(path, samples=np.concatenate(records['production']),
                burn_samples=np.concatenate(records['equilibration']),
                initial_state=initial_state, final_state=states[0],
                initial_energy_int=initial_energy, final_energy_int=energies[0],
                temperature_K=temperature, temperature_index=index, branch=branch,
                chain=c, step=step, seed=seed,
                accepted_production=accepted['production'],
                accepted_equilibration=accepted['equilibration'],
                rng_state_json=json.dumps(rng.bit_generator.state),
                config_fingerprint=fingerprint, xyz=lat.xyz, sublattice=lat.sublattice)
            completed += 1
            print(f'cycle {c+1}/{len(cfg["seeds"])} | {branch} | {temperature:g} K | '
                  f'palier {step+1}/{len(steps)} sauvegardé | '
                  f'{time.perf_counter()-start:.1f} s', flush=True)
            if progress is not None:
                progress(c, step, path)
    write_status(out, 'analyzing', L=cfg['L'], completed_plateaus=completed,
                 total_plateaus=len(steps)*len(cfg['seeds']))
    return analyze_cycle(out)


def load_cycle_archives(output):
    out = Path(output)
    cfg = json.loads((out/'config.json').read_text(encoding='utf8'))
    steps = schedule(cfg)
    data = {}
    for branch in ('heating', 'cooling'):
        per_chain = []
        for c in range(len(cfg['seeds'])):
            records = [None]*len(cfg['temperatures'])
            for step, (direction, index, _) in enumerate(steps):
                if direction == branch:
                    with np.load(stage_path(out, c, step)) as a:
                        records[index] = {key: a[key] for key in a.files}
            per_chain.append(records)
        data[branch] = per_chain
    return cfg, data


def analyze_cycle(output):
    out = Path(output)
    cfg, data = load_cycle_archives(out)
    lat = build_lattice(cfg['L'])
    steps = schedule(cfg)
    # Validate time continuity across *all* heating and cooling plateaus.
    for c in range(len(cfg['seeds'])):
        previous = ordered(lat, (2*c) % 4)
        for branch, index, temperature in steps:
            a = data[branch][c][index]
            np.testing.assert_array_equal(previous, a['initial_state'])
            if energy_int(a['initial_state'], lat.neighbors) != int(a['initial_energy_int']):
                raise RuntimeError('Initial plateau energy mismatch.')
            if energy_int(a['final_state'], lat.neighbors) != int(a['final_energy_int']):
                raise RuntimeError('Final plateau energy mismatch.')
            if a['final_state'].sum() != lat.n//4:
                raise RuntimeError('Composition mismatch.')
            for key in ('samples', 'burn_samples'):
                s = a[key]
                if not np.isfinite(s).all():
                    raise RuntimeError('Non-finite observation.')
                p = s[:, 5:9]
                np.testing.assert_allclose(p.sum(axis=1), 1, atol=2e-13)
                np.testing.assert_allclose(s[:, 1], (p[:, 0]-.25)/.75, atol=2e-13)
                np.testing.assert_allclose(s[:, 2], np.sqrt((4/3)*((p-.25)**2).sum(axis=1)), atol=2e-13)
                b1 = (1-s[:, 3])*(lat.n/4)*12*.75
                b2 = (1-s[:, 4])*(lat.n/4)*6*.75
                expected = -3.813*lat.n+.161*(12*(lat.n/4)-b1)/2+.0063*b2
                np.testing.assert_allclose(s[:, 0], expected, atol=2e-9, rtol=0)
            previous = a['final_state']
    temperatures = np.asarray(cfg['temperatures'], dtype=float)
    all_rows = []
    branch_rows = {}
    for branch, chains in data.items():
        samples = np.array([[a['samples'] for a in chain] for chain in chains]).transpose(0, 2, 1, 3)
        rows = summarize(samples, temperatures, lat.n, bootstrap=cfg.get('bootstrap', 250))
        for index, row in enumerate(rows):
            row['branch'] = branch
            row['local_acceptance'] = float(np.mean([
                int(chain[index]['accepted_production'])/(lat.n*cfg['sample_sweeps']) for chain in chains]))
            for mu in range(4):
                row[f'p{mu}'] = float(samples[:, :, index, 5+mu].mean())
        save_csv(out/f'observables_{branch}.csv', rows)
        all_rows.extend(rows)
        branch_rows[branch] = rows
    save_csv(out/'observables.csv', all_rows)
    # Heating/cooling belong to the same cycles: differences are paired.
    # Publish cycle-to-cycle dispersion, not an independence-based error bar.
    from .statistics import moments, NAMES
    differences = []
    for k, T in enumerate(temperatures):
        delta = np.array([moments(data['cooling'][c][k]['samples'], T, lat.n)-
                          moments(data['heating'][c][k]['samples'], T, lat.n)
                          for c in range(len(cfg['seeds']))])
        row = {'T_K': float(T)}
        for j, name in enumerate(NAMES):
            row[f'delta_{name}'] = float(delta[:, j].mean())
            row[f'delta_{name}_between_cycles_se'] = float(delta[:, j].std(ddof=1)/np.sqrt(len(delta)))
        differences.append(row)
    save_csv(out/'paired_differences.csv', differences)
    peaks = {branch: float(temperatures[np.argmax([r['cv_kB'] for r in rows])])
             for branch, rows in branch_rows.items()}
    result = {'config': cfg, 'n_atoms': lat.n, 'plateaus': len(steps)*len(cfg['seeds']),
        'peaks_cv_grid_K': peaks,
        'peak_difference_cooling_minus_heating_K': peaks['cooling']-peaks['heating'],
        'diagnostic_failed_temperatures_K': {
            branch: [r['T_K'] for r in rows if not r['diagnostic_ok']]
            for branch, rows in branch_rows.items()},
        'validation': 'energy, composition, all observables, and plateau continuity passed',
        'notes': [
            'No replica exchange; temperature changes only between completed plateaus.',
            'Cooling begins with a second independent dwell at the maximum T, without a state reset.',
            'Branch variances describe fluctuations; equilibrium interpretation requires diagnostics.',
            'Finite-time hysteresis is protocol dependent, not a proof of thermodynamic first order.',
            'Paired-difference SE describes between-cycle dispersion only, with few independent cycles.',
            'Archives store all sampled burn-in/production measurements and both endpoint states.',
            'Resume starts at the last completed plateau; an interrupted partial plateau is replayed.']}
    (out/'summary.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf8')
    return result
