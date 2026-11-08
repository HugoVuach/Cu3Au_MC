"""Campagnes reproductibles et contrôles d'intégrité de chaque segment."""
import csv
import json
import platform
import time
from pathlib import Path
import numpy as np
from .lattice import build_lattice, ordered, random_configuration
from .model import energy_int, observables, SCALE, random_energy_limit
from .kernel import run_segment
from .statistics import summarize


def load_config(path):
    cfg = json.loads(Path(path).read_text(encoding='utf8'))
    required = ['L', 'temperatures', 'seeds', 'burn_sweeps', 'sample_sweeps',
                'sample_every', 'exchange_every', 'chunk_sweeps']
    for k in required:
        if k not in cfg:
            raise ValueError(f'Paramètre manquant : {k}')
    t = np.asarray(cfg['temperatures'], dtype=float)
    if len(t) < 3 or np.any(t <= 0) or np.any(np.diff(t) <= 0) or not np.all(np.isfinite(t)):
        raise ValueError('Au moins trois températures strictement croissantes et positives.')
    if len(set(cfg['seeds'])) != len(cfg['seeds']) or len(cfg['seeds']) < 2:
        raise ValueError('Au moins deux graines distinctes pour les diagnostics.')
    for k in ['burn_sweeps', 'sample_sweeps', 'sample_every', 'chunk_sweeps']:
        if not isinstance(cfg[k], int) or cfg[k] <= 0:
            raise ValueError(f'{k} doit être un entier strictement positif.')
    if not isinstance(cfg['exchange_every'], int) or cfg['exchange_every'] < 0:
        raise ValueError('exchange_every doit être un entier positif ou nul.')
    capture = cfg.get('store_states_every', 0)
    if (not isinstance(capture, int) or isinstance(capture, bool) or capture < 0 or
            (capture and (capture % cfg['sample_every'] or cfg['sample_sweeps'] % capture))):
        raise ValueError('store_states_every doit être nul ou un multiple de sample_every divisant sample_sweeps.')
    for k in ['burn_sweeps', 'sample_sweeps', 'chunk_sweeps']:
        if cfg[k] % cfg['sample_every']:
            raise ValueError(f'{k} doit être multiple de sample_every.')
        if cfg['exchange_every'] and cfg[k] % cfg['exchange_every']:
            raise ValueError(f'{k} doit être multiple de exchange_every.')
    build_lattice(cfg['L'])
    return cfg


def save_csv(path, rows):
    with Path(path).open('w', newline='', encoding='utf8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def campaign(cfg, output, progress=None):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    lat = build_lattice(cfg['L'])
    temperatures = np.array(cfg['temperatures'], dtype=np.float64)
    nr = len(temperatures)
    n = lat.n
    all_samples, all_labels, diagnostics = [], [], []
    start = time.perf_counter()
    for chain, seed in enumerate(cfg['seeds']):
        rng = np.random.default_rng(seed)
        # Deux variantes différentes et, pour la seconde chaîne, départs chauds
        # désordonnés. Aux très basses T, un départ aléatoire serait piégé.
        states = np.array([ordered(lat, (2*chain) % 4) if
                           (chain == 0 or t < cfg.get('random_start_above_K', 1300))
                           else random_configuration(lat, rng)
                           for t in temperatures])
        energies = np.array([energy_int(s, lat.neighbors) for s in states], dtype=np.int64)
        labels = np.arange(nr, dtype=np.int32)
        phase = 0
        serial = 0
        records, lab_records = [], []
        spatial_states, spatial_sweeps = [], []
        burn_records = []
        count_ok = np.zeros(nr-1, dtype=np.int64)
        count_n = np.zeros(nr-1, dtype=np.int64)
        local_ok = np.zeros(nr, dtype=np.int64)
        for mode, length in [('equilibration', cfg['burn_sweeps']),
                             ('mesure', cfg['sample_sweeps'])]:
            done = 0
            while done < length:
                chunk = min(cfg['chunk_sweeps'], length-done)
                capture = cfg.get('store_states_every', 0) if mode == 'mesure' else 0
                if capture:
                    chunk = min(chunk, capture-done % capture)
                serial += 1
                res = run_segment(states, energies, temperatures, lat.neighbors,
                                  lat.sublattice, chunk, cfg['exchange_every'],
                                  cfg['sample_every'], int(rng.integers(1, 2**31-1)),
                                  labels, phase)
                rec, lab, acc, ok, attempted, phase = res
                # Vérification indépendante, exacte en arithmétique entière.
                check = np.array([energy_int(s, lat.neighbors) for s in states])
                if not np.array_equal(check, energies):
                    raise RuntimeError('Dérive de l’énergie détectée.')
                if not np.all(states.sum(axis=1) == n//4):
                    raise RuntimeError('La composition a changé.')
                if np.any(energies < round(-3.813*n*SCALE)):
                    raise RuntimeError('Énergie inférieure à la borne analytique.')
                if mode == 'mesure':
                    records.append(rec)
                    lab_records.append(lab)
                    count_ok += ok
                    count_n += attempted
                    local_ok += acc
                else:
                    burn_records.append(rec)
                done += chunk
                if capture and done % capture == 0:
                    spatial_states.append(states.copy())
                    spatial_sweeps.append(done)
                if progress is not None:
                    progress(chain+1, mode, done, length, len(spatial_states))
                print(f'chaîne {chain+1}/{len(cfg["seeds"])} | {mode} '
                      f'{done}/{length} balayages | {time.perf_counter()-start:.1f} s', flush=True)
        samples = np.concatenate(records)
        trace = np.concatenate(lab_records)
        all_samples.append(samples)
        all_labels.append(trace)
        extra = {}
        if cfg.get('store_states_every', 0):
            extra = dict(spatial_states=np.array(spatial_states),
                         spatial_sample_sweeps=np.array(spatial_sweeps),
                         spatial_capture_every=cfg['store_states_every'])
        destination = out/f'chain_{chain+1:02d}.npz'
        temporary = destination.with_suffix('.npz.tmp')
        with temporary.open('wb') as stream:
            np.savez_compressed(stream, samples=samples,
                                burn_samples=np.concatenate(burn_records), labels=trace,
                                final_states=states, temperatures=temperatures,
                                xyz=lat.xyz, sublattice=lat.sublattice, **extra)
        temporary.replace(destination)
        # Nombre de trajets froid -> chaud -> froid, par identité de réplique.
        roundtrips = np.zeros(nr, dtype=int)
        stage = np.zeros(nr, dtype=int)
        for lab in trace:
            cold, hot = int(lab[0]), int(lab[-1])
            if stage[cold] == 2:
                roundtrips[cold] += 1
            stage[cold] = 1
            if stage[hot] == 1:
                stage[hot] = 2
        diagnostics.append({
            'seed': seed,
            'swap_acceptance': np.divide(count_ok, count_n,
                                         out=np.zeros(nr-1), where=count_n > 0).tolist(),
            'swap_attempts': count_n.tolist(),
            'local_acceptance': (local_ok/(cfg['sample_sweeps']*n)).tolist(),
            'roundtrips': roundtrips.tolist(),
            'energy_checks': 'exact integer agreement after every chunk',
            'composition_checks': 'N_Au=N/4 after every chunk'})
    samples = np.array(all_samples)
    rows = summarize(samples, temperatures, n, bootstrap=cfg.get('bootstrap', 200))
    save_csv(out/'observables.csv', rows)
    peak = int(np.argmax([r['cv_kB'] for r in rows]))
    left = temperatures[max(0, peak-1)]
    right = temperatures[min(nr-1, peak+1)]
    critical_indices = list(range(max(0, peak-2), min(nr, peak+3)))
    summary = {
        'model': 'Cu3Au, CFC rigide, interactions 1er et 2e voisins',
        'n_atoms': n, 'config': cfg,
        'T_cv_peak_grid_K': float(temperatures[peak]),
        'T_cv_neighbor_bracket_K': [float(left), float(right)],
        'T_peak_interpretation': 'Maximum on sampled grid; neighbor bracket is not a confidence interval.',
        'critical_diagnostics_ok': all(rows[i]['diagnostic_ok'] for i in critical_indices),
        'critical_window_rows': critical_indices,
        'ordered_energy_eV_atom': -3.813,
        'random_limit_eV_atom_finite_N': random_energy_limit(n),
        'elapsed_seconds': time.perf_counter()-start,
        'runtime': {'python': platform.python_version(), 'platform': platform.platform(),
                    'numpy': np.__version__},
        'chains': diagnostics,
        'notes': [
            'Les deux chaînes ne sont pas indépendantes entre températures ; elles sont indépendantes entre graines.',
            'R-hat est calculé sur énergie et amplitude q, pas sur eta signé qui dépend de la variante.',
            'Un bon diagnostic ne prouve pas l’exploration de tous les états métastables.',
            'Aucune extrapolation à la limite thermodynamique avec une seule taille.']}
    (out/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding='utf8')
    return summary
