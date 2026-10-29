"""Pedagogical figures and report tables, computed from saved trajectories."""
import csv
from datetime import datetime
from itertools import product
import json
import platform
from pathlib import Path
import re
import sys
import numpy as np
from .statistics import bootstrap_peak


def _data(result_dir):
    root = Path(result_dir)
    summary = json.loads((root/'summary.json').read_text(encoding='utf8'))
    with (root/'observables.csv').open(encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    archives = []
    for path in sorted(root.glob('chain_*.npz')):
        with np.load(path) as a:
            archives.append({key: a[key] for key in a.files})
    return root, summary, rows, archives


def record_current_validation(result_dir, test_log):
    """Record a successful current test run and independent reference checks."""
    import scipy
    import numba
    from .lattice import build_lattice, ordered, random_configuration
    from .model import energy_int, delta_reference, observables, SCALE
    from .kernel import delta_fast, sample_observables

    root, summary, _, _ = _data(result_dir)
    match = re.search(r'Ran (\d+) tests?', test_log)
    if not match or not re.search(r'^OK\s*$', test_log, re.MULTILINE):
        raise RuntimeError('A successful test log is required.')
    lat = build_lattice(summary['config']['L'])
    state = ordered(lat)
    ref = observables(state, lat)
    initial_energy = energy_int(state, lat.neighbors)/SCALE
    i = int(np.flatnonzero(state)[0])
    j = int(lat.neighbors[i, 0])
    defect = delta_reference(state, i, j, lat.neighbors)/SCALE
    assert abs(initial_energy+3.813*lat.n) < 1e-9
    np.testing.assert_allclose(ref[1:5], [1, 1, -1/3, 1], atol=1e-13)
    assert abs(defect-.5586) < 1e-12
    rng = np.random.default_rng(984)
    max_gap, max_observable_gap = 0, 0.
    for _ in range(30):
        s = random_configuration(lat, rng)
        a, b = rng.choice(lat.n, 2, replace=False)
        before = energy_int(s, lat.neighbors)
        fast = delta_fast(s, a, b, lat.neighbors)
        max_observable_gap = max(max_observable_gap, float(np.max(abs(
            sample_observables(s, before, lat.neighbors, lat.sublattice)
            -observables(s, lat)))))
        s[a], s[b] = s[b], s[a]
        max_gap = max(max_gap, abs(energy_int(s, lat.neighbors)-before-fast))
    assert max_gap == 0
    assert max_observable_gap < 1e-12
    now = datetime.now().astimezone()
    result = {
        'date': now.date().isoformat(), 'timestamp': now.isoformat(),
        'passed': True, 'test_count': int(match.group(1)),
        'environment': {'python': platform.python_version(), 'os': platform.system(),
                        'numpy': np.__version__, 'scipy': scipy.__version__,
                        'numba': numba.__version__, 'executable': sys.executable},
        'checks': [
            ['Sites et voisinages', rf'${lat.n}$ ; $12+6$', rf'${lat.n}$ ; $12+6$ : OK'],
            ['Énergie ordonnée (eV)', f'${-3.813*lat.n:.3f}$', f'${initial_energy:.3f}$'],
            ['Cowley ordonné', r'$(-1/3,1)$', r'$(-1/3,1)$ : OK'],
            ['Ordre de référence', r'$\eta=q=1$', r'$\eta=q=1$ : OK'],
            ['Défaut Au--Cu (eV)', '$0.5586$', f'${defect:.4f}$'],
            ['Écart énergie locale/globale', '$0$ micro-eV', f'${max_gap}$ : OK'],
            ['Observables compilées/référence', r'Écart $<10^{-12}$', 'OK (30 états)'],
            ['Suite automatique', f'{match.group(1)} tests', 'Tous réussis'],
        ],
    }
    (root/'tests_current.log').write_text(test_log, encoding='utf8')
    (root/'verification_current.json').write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding='utf8')
    return result


def analyze_peak(result_dir, bootstrap=1000):
    root, summary, rows, archives = _data(result_dir)
    samples = np.array([a['samples'] for a in archives])
    t = np.array([float(row['T_K']) for row in rows])
    draws = samples.shape[1]
    requested = max(2, int(np.ceil(2*max(float(r['tau_max_samples']) for r in rows))))
    block = min(requested, max(2, draws//4))
    blocks = sorted(set([max(2, block//2), block, min(draws, 2*block)]))
    cases = []
    for length in blocks:
        peaks = bootstrap_peak(samples[..., 0], t, summary['n_atoms'], length,
                               bootstrap=bootstrap, seed=9421+length)
        values, counts = np.unique(peaks, return_counts=True)
        low, high = np.quantile(peaks, [.025, .975], method='nearest')
        cases.append({'block_samples': length, 'ci95_K': [float(low), float(high)],
                      'se_K': float(np.std(peaks, ddof=1)),
                      'frequencies': {str(float(v)): int(c) for v, c in zip(values, counts)}})
    q = np.array([float(row['q']) for row in rows])
    u = np.array([float(row['u_eV']) for row in rows])
    cv = np.array([float(row['cv_kB']) for row in rows])
    xq = np.array([float(row['Xq_per_eV']) for row in rows])
    result = {'bootstrap': bootstrap, 'seed': 9421, 'draws_per_chain': draws,
              'base_block_samples': block, 'requested_block_samples': requested,
              'base_block_sweeps': block*summary['config']['sample_every'],
              'observed_peak_K': float(t[np.argmax(cv)]), 'cases': cases,
              'indicators_K': {'Maximum de capacité calorifique': float(t[np.argmax(cv)]),
                               'Pente la plus forte de q': float(t[np.argmax(-np.gradient(q, t))]),
                               "Pente la plus forte de l'énergie": float(t[np.argmax(np.gradient(u, t))]),
                               "Maximum des fluctuations d'ordre": float(t[np.argmax(xq)])},
              'scope': 'Grid-valued peak; fixed size and observed chains; synchronized circular blocks.'}
    (root/'peak_uncertainty.json').write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding='utf8')
    return result


def make_pedagogical_figures(result_dir):
    # Import after figures.py has initialized its plotting helpers.
    from .figures import plt, style, save, INK, BLUE, GOLD, RED, GREEN
    from matplotlib.patches import FancyBboxPatch, Polygon
    root, summary, rows, archives = _data(result_dir)
    folder = root/'figures'
    folder.mkdir(exist_ok=True)
    style()
    t = np.array([float(r['T_K']) for r in rows])
    peak = int(np.argmax([float(r['cv_kB']) for r in rows]))
    n = summary['n_atoms']

    # 13: shells centered on an Au in a perfect ordered reference.
    first = np.array([d for d in product((-1, 0, 1), repeat=3)
                      if sum(v*v for v in d) == 2])/2
    second = np.array([tuple(2*s if k == axis else 0 for k in range(3))
                       for axis in range(3) for s in (-1, 1)])/2
    fig = plt.figure(figsize=(11, 4.7))
    for panel, (points, color, title) in enumerate([
        (first, BLUE, 'Première coquille : 12 Cu à a/√2'),
        (second, GOLD, 'Deuxième coquille : 6 Au à a')], 1):
        ax = fig.add_subplot(1, 2, panel, projection='3d')
        for point in points:
            ax.plot(*np.array([[0, 0, 0], point]).T, color=color, alpha=.45)
        ax.scatter(*points.T, c=color, s=85, edgecolors='white', depthshade=False)
        ax.scatter([0], [0], [0], c=GOLD, marker='*', s=240, edgecolors=INK,
                   label='Au de référence', depthshade=False)
        ax.set(xlabel='x/a', ylabel='y/a', zlabel='z/a', title=title,
               xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), zlim=(-1.1, 1.1))
        ax.set_box_aspect((1, 1, 1))
        ax.legend(fontsize=8, loc='upper left')
    fig.suptitle('Deux distances, deux voisinages chimiques dans L1₂', y=.99)
    fig.subplots_adjust(bottom=.15, top=.85, wspace=.12)
    save(fig, folder, '13_coquilles',
         'Géométrie fixe : la coquille dépend de la distance, son occupation dépend de la configuration.')

    # 14: compact algorithm flowchart, without external diagram packages.
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set(xlim=(0, 10), ylim=(0, 11)); ax.axis('off')
    def box(x, y, text, width=5.6, height=.62):
        ax.add_patch(FancyBboxPatch((x-width/2, y-height/2), width, height,
                     boxstyle='round,pad=0.08', facecolor='#f2f6f9', edgecolor=BLUE))
        ax.text(x, y, text, ha='center', va='center', fontsize=10)
    def arrow(a, b, text=None):
        ax.annotate('', xy=b, xytext=a, arrowprops={'arrowstyle': '->', 'color': INK})
        if text:
            ax.text((a[0]+b[0])/2+.12, (a[1]+b[1])/2, text, fontsize=9)
    box(5, 10.2, 'Configuration courante et énergie à température fixée')
    box(5, 9.1, 'Tirer un Au et un Cu ; calculer ΔE sans modifier l’état')
    arrow((5, 9.8), (5, 9.5))
    ax.add_patch(Polygon([(5, 8.65), (8, 7.9), (5, 7.15), (2, 7.9)],
                        facecolor='#fff6df', edgecolor=GOLD))
    ax.text(5, 7.9, 'ΔE ≤ 0 ou u < exp[−ΔE/(kBT)] ?', ha='center', va='center', fontsize=10)
    arrow((5, 8.7), (5, 8.6))
    box(2.6, 6.35, 'Permuter Au et Cu\nE ← E + ΔE', 3.4, .8)
    box(7.4, 6.35, 'Conserver l’état\net l’énergie', 3.4, .8)
    arrow((3.8, 7.6), (2.6, 6.82), 'oui')
    arrow((6.2, 7.6), (7.4, 6.82), 'non')
    box(5, 5.1, 'Compter la tentative, acceptée ou refusée')
    arrow((2.6, 5.9), (4.1, 5.47)); arrow((7.4, 5.9), (5.9, 5.47))
    box(5, 3.95, 'Répéter N tentatives par température : un balayage')
    arrow((5, 4.7), (5, 4.35))
    box(5, 2.8, 'Après le balayage : échanges de répliques si échéance')
    arrow((5, 3.55), (5, 3.2))
    box(5, 1.65, 'Puis mesurer si échéance : E, η, q, α₁, α₂ et pμ')
    arrow((5, 2.4), (5, 2.05))
    ax.text(5, .65, 'Poursuivre les balayages ; contrôler énergie et composition après le segment.\n'
            'Équilibration : conserver les traces ; production : calculer les moyennes.',
            ha='center', fontsize=10)
    fig.suptitle('De la proposition à la mesure', y=.98)
    fig.subplots_adjust(top=.94, bottom=.06)
    save(fig, folder, '14_organigramme')

    # 15: genuine equivalent reference configurations, not MC samples.
    from .lattice import build_lattice, ordered
    small = build_lattice(3)
    mask = np.all(small.xyz < 4, axis=1)
    xyz = small.xyz[mask]/2
    fig = plt.figure(figsize=(12, 4.3))
    for variant in range(4):
        ax = fig.add_subplot(1, 4, variant+1, projection='3d')
        au = ordered(small, variant)[mask].astype(bool)
        ax.scatter(*xyz[~au].T, c=BLUE, s=20, alpha=.28, depthshade=False)
        ax.scatter(*xyz[au].T, c=GOLD, s=65, edgecolors=INK, linewidths=.3, depthshade=False)
        ax.set_title(f'Au sur μ = {variant}\nη = '+('1' if variant == 0 else '−1/3')+' ; q = 1', fontsize=11)
        ax.set(xticks=[], yticks=[], zticks=[])
        ax.set_box_aspect((1, 1, 1)); ax.view_init(23, 32)
    fig.suptitle('Quatre variantes ordonnées équivalentes : Au en or, Cu en bleu', y=.97, fontsize=16)
    fig.subplots_adjust(bottom=.16, top=.77, wspace=.02)
    save(fig, folder, '15_variantes', 'États de référence construits : chaque panneau représente une configuration parfaitement ordonnée.')

    # 16: two planes are necessary to include all four fcc sublattices.
    a = archives[0]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), sharex=True, sharey=True)
    for col, k in enumerate([0, peak, len(t)-1]):
        for row, z in enumerate([0, 1]):
            plane = a['xyz'][:, 2] == z
            points = a['xyz'][plane]/2
            au = a['final_states'][k, plane].astype(bool)
            axes[row, col].scatter(points[:, 0], points[:, 1],
                c=np.where(au, GOLD, BLUE), s=65, edgecolors=INK, linewidths=.3)
            axes[row, col].set(title=f'{t[k]:.0f} K ; z/a = {z/2:g}',
                xlabel='x/a', ylabel='y/a', xlim=(-.3, summary['config']['L']-.2),
                ylim=(-.3, summary['config']['L']-.2))
            axes[row, col].set_aspect('equal')
    fig.suptitle('Coupes des configurations finales : Au en or, Cu en bleu', y=.99, fontsize=16)
    fig.tight_layout(rect=(0, .055, 1, .94))
    save(fig, folder, '16_coupes_2D', 'Deux plans adjacents montrent les quatre sous-réseaux ; un instantané ne remplace pas une moyenne thermique.')

    # 17: cumulative means use production only, at the observed cV maximum.
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.8), sharex=True)
    for c, a in enumerate(archives):
        data = a['samples'][:, peak]
        count = np.arange(1, len(data)+1)
        x = count*summary['config']['sample_every']
        for ax, values in zip(axes, [data[:, 0]/n, data[:, 2]]):
            ax.plot(x, np.cumsum(values)/count, lw=1.4, label=f'Chaîne {c+1}')
            ax.legend(fontsize=9)
    axes[0].set(ylabel='Moyenne cumulée de u (eV/atome)')
    axes[1].set(ylabel='Moyenne cumulée de q', xlabel='Balayages de production depuis la fin de l’équilibration')
    fig.suptitle(f'Accord progressif des moyennes à {t[peak]:.0f} K', y=.99)
    fig.tight_layout(rect=(0, .06, 1, .94))
    save(fig, folder, '17_moyennes_cumulees', 'Les points successifs partagent leurs mesures : ces courbes sont corrélées et ne prouvent pas à elles seules l’équilibre.')

    # 18: synchronized bootstrap and block-length sensitivity.
    analysis = analyze_peak(root)
    fig, ax = plt.subplots(figsize=(10, 4.9))
    for c, case in enumerate(analysis['cases']):
        values = np.array([float(v) for v in case['frequencies']])
        frequencies = np.array(list(case['frequencies'].values()))/analysis['bootstrap']
        ax.bar(values+(c-1)*1.2, frequencies, width=1.1,
               label=f"Blocs de {case['block_samples']} mesures", color=[BLUE, RED, GREEN][c])
    all_values = [float(v) for case in analysis['cases'] for v in case['frequencies']]
    ax.set(xlim=(min(all_values)-10, max(all_values)+10), ylim=(0, 1.05),
           xlabel='Température du maximum sur la grille (K)', ylabel='Fréquence bootstrap')
    ax.set_xticks(sorted(set(all_values)))
    ax.legend(fontsize=10)
    fig.suptitle('Incertitude statistique sur la position du maximum', y=.98)
    fig.tight_layout(rect=(0, .07, 1, .94))
    save(fig, folder, '18_incertitude_pic', '1 000 rééchantillonnages par longueur de bloc ; températures synchronisées dans chaque chaîne ; taille et grille fixées.')


def write_pedagogical_tables(result_dir):
    root, summary, rows, _ = _data(result_dir)
    out = root/'report'
    out.mkdir(exist_ok=True)
    def write(name, lines):
        (out/name).write_text('\n'.join(lines)+'\n', encoding='utf8')
    # Actual current-run metadata, or a clearly dated historical log fallback.
    path = root/'verification_current.json'
    if path.exists():
        v = json.loads(path.read_text(encoding='utf8'))
        env = v['environment']
        label = f"{env['os']}"+r'\newline '+f"Python {env['python']}"
        columns = ''.join(r'>{\raggedright\arraybackslash}p{'+str(width)+r'\linewidth}'
                          for width in (.26, .17, .20, .15, .14))
        lines = [r'\begin{center}\footnotesize\setlength{\tabcolsep}{3pt}',
                 r'\begin{tabular}{'+columns+r'}\toprule',
                 r'Test & Attendu & Résultat & Environnement & Date\\\midrule']
        for test, expected, observed in v['checks']:
            lines.append(f"{test} & {expected} & {observed} & {label} & {v['date']}"+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}',
                  f"Versions : NumPy {env['numpy']}, SciPy {env['scipy']}, Numba {env['numba']}. "
                  r'La validation actuelle est distincte des journaux historiques de la campagne.']
        write('current_validation.tex', lines)
    else:
        write('current_validation.tex', [
            r'La suite est décrite dans \texttt{tests.log}. Aucune nouvelle exécution datée n’est enregistrée ici.'])

    picks = [rows[0]]+[min(rows, key=lambda r: abs(float(r['T_K'])-x)) for x in
                      [1000, summary['T_cv_peak_grid_K'], 1100]]+[rows[-1]]
    lines = [r'\begin{center}\small', r'\begin{tabular}{rrrrrr}\toprule',
             r'$T$ (K) & $u$ (eV/atome) & $q$ & $\alpha_1$ & $\alpha_2$ & $c_V/k_B$\\\midrule']
    for row in picks:
        lines.append(f"{float(row['T_K']):.0f} & {float(row['u_eV']):.6f} & "
                     f"{float(row['q']):.3f} & {float(row['alpha1']):.3f} & "
                     f"{float(row['alpha2']):.3f} & {float(row['cv_kB']):.3f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}']
    write('representative_results.tex', lines)
    if not (root/'peak_uncertainty.json').exists():
        analyze_peak(root)
    analysis = json.loads((root/'peak_uncertainty.json').read_text(encoding='utf8'))
    base = next(c for c in analysis['cases'] if c['block_samples'] == analysis['base_block_samples'])
    low, high = base['ci95_K']
    lines = [f"Sur les archives, {analysis['bootstrap']} tirages par blocs circulaires "
             f"de {analysis['base_block_samples']} mesures "
             f"({analysis['base_block_sweeps']} balayages) donnent un intervalle percentile "
             f"à 95\\,\\% de $[{low:.0f},{high:.0f}]$ K et un écart-type "
             f"de ${base['se_K']:.2f}$ K pour la position du maximum sur la grille.",
             r'\begin{center}\small\begin{tabular}{rrr}\toprule',
             r'Bloc (mesures) & Intervalle percentile 95\,\% (K) & Écart-type (K)\\\midrule']
    for case in analysis['cases']:
        lo, hi = case['ci95_K']
        lines.append(f"{case['block_samples']} & $[{lo:.0f},{hi:.0f}]$ & {case['se_K']:.2f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
              r'\begin{center}\small\begin{tabular}{lr}\toprule',
              r'Indicateur sur la grille & Température (K)\\\midrule']
    for label, temperature in analysis['indicators_K'].items():
        lines.append(f'{label} & {temperature:.0f}'+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}']
    write('peak_uncertainty.tex', lines)
