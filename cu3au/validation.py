"""Validation de cohérence des résultats livrés ; pas de courbes synthétiques."""
import csv
import json
import platform
import sys
from pathlib import Path
import numpy as np
import scipy
import matplotlib
import numba
from .lattice import build_lattice
from .model import energy_int, SCALE


def validate_results(result_dir):
    root = Path(result_dir)
    summary = json.loads((root/'summary.json').read_text(encoding='utf8'))
    lat = build_lattice(summary['config']['L'])
    all_samples = []
    checks = []
    for file in sorted(root.glob('chain_*.npz')):
        with np.load(file) as a:
            s = a['samples']
            assert np.isfinite(s).all(), 'Valeur non finie dans les échantillons.'
            p = s[:, :, 5:9]
            np.testing.assert_allclose(p.sum(axis=2), 1, atol=2e-13)
            np.testing.assert_allclose(s[:, :, 1], (p[:, :, 0]-.25)/.75, atol=2e-13)
            np.testing.assert_allclose(s[:, :, 2], np.sqrt((4/3)*((p-.25)**2).sum(axis=2)), atol=2e-13)
            b1 = (1-s[:, :, 3])*(lat.n/4)*12*.75
            a1 = (12*(lat.n/4)-b1)/2
            b2 = (1-s[:, :, 4])*(lat.n/4)*6*.75
            expected = -3.813*lat.n+.161*a1+.0063*b2
            np.testing.assert_allclose(s[:, :, 0], expected, atol=2e-9, rtol=0)
            for st in a['final_states']:
                assert st.sum() == lat.n//4
                assert energy_int(st, lat.neighbors) >= round(-3.813*lat.n*SCALE)
            all_samples.append(s)
            checks.append(f'{file.name}: énergie reconstruite par les paires/Cowley, ordre et composition vérifiés pour toutes les mesures')
    all_samples = np.array(all_samples)
    assert len(all_samples) == len(summary['config']['seeds'])
    with (root/'observables.csv').open(encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    np.testing.assert_allclose([float(r['u_eV']) for r in rows],
                               all_samples[:, :, :, 0].mean(axis=(0, 1))/lat.n, atol=1e-10)
    assert all(float(r['cv_kB']) >= 0 for r in rows)
    assert np.all(all_samples[:, :, :, 0] >= -3.813*lat.n-1e-9)
    assert float(rows[0]['q']) > .9, 'Le froid ne rejoint pas un état fortement ordonné.'
    assert float(rows[-1]['q']) < .3, 'La borne chaude est encore fortement ordonnée.'
    assert 0 < np.argmax([float(r['cv_kB']) for r in rows]) < len(rows)-1, 'Maximum thermique sur une borne.'
    checks += ['CSV cohérent avec les archives brutes', 'Bornes analytiques et régimes froid/chaud vérifiés']
    result = {
        'technical_validation': 'passed',
        'critical_diagnostics_ok': summary['critical_diagnostics_ok'],
        'checks': checks,
        'limitations': 'Validation technique et références exactes ; T* reste un maximum de taille finie.',
        'runtime': {'python': platform.python_version(), 'executable': sys.executable,
                    'numpy': np.__version__, 'scipy': scipy.__version__,
                    'matplotlib': matplotlib.__version__, 'numba': numba.__version__}}
    (root/'validation.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf8')
    return result
