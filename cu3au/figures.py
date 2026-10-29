"""Figures scientifiques vectorielles (PDF) et PNG 180 dpi, sans réseau."""
import csv
import json
import io
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from .lattice import build_lattice
from .meanfield import g, A, coexistence
from .model import KB

INK = '#183047'
BLUE = '#167d9a'
GOLD = '#dba52d'
RED = '#ce5964'
GREEN = '#319b79'


def style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.titlesize': 13,
        'axes.labelsize': 11, 'figure.titlesize': 19, 'text.color': INK,
        'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.edgecolor': '#b1bec8', 'axes.grid': True, 'grid.alpha': .18,
        'grid.linewidth': .6, 'legend.frameon': False,
        'savefig.facecolor': 'white', 'figure.facecolor': 'white',
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })


def save(fig, folder, name, subtitle=None):
    if subtitle:
        fig.text(.5, .015, subtitle, ha='center', fontsize=9, color='#536779')
    # Finalisation en mémoire puis remplacement atomique : aucun PDF partiel
    # ne peut être consommé par la compilation du rapport.
    for fmt in ('pdf', 'png'):
        buffer = io.BytesIO()
        fig.savefig(buffer, format=fmt, dpi=180, bbox_inches='tight')
        data = buffer.getvalue()
        if fmt == 'pdf' and not data.rstrip().endswith(b'%%EOF'):
            raise RuntimeError(f'PDF non finalisé : {name}')
        temporary = folder/f'{name}.{fmt}.tmp'
        temporary.write_bytes(data)
        if temporary.stat().st_size != len(data):
            raise RuntimeError(f'Écriture incomplète : {temporary}')
        temporary.replace(folder/f'{name}.{fmt}')
    plt.close(fig)


def curves(ax, T, ys, err, color, label):
    ax.plot(T, ys, 'o-', ms=3.5, lw=1.6, color=color, label=label)
    ax.fill_between(T, ys-err, ys+err, color=color, alpha=.16, lw=0)


def make_cowley_figure(result_dir):
    """Figure dédiée à alpha1 et alpha2, depuis le CSV existant."""
    style()
    root = Path(result_dir)
    folder = root/'figures'
    folder.mkdir(exist_ok=True)
    summary = json.loads((root/'summary.json').read_text(encoding='utf8'))
    with (root/'observables.csv').open(encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    col = lambda key: np.array([float(row[key]) for row in rows])
    temperatures = col('T_K')
    peak_temperature = summary['T_cv_peak_grid_K']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharex=True)
    parameters = [
        ('alpha1', RED, -1/3, r'$\alpha_1$', 'Première coquille : 12 voisins', 'lower right'),
        ('alpha2', GREEN, 1, r'$\alpha_2$', 'Deuxième coquille : 6 voisins', 'upper right'),
    ]
    for ax, (key, color, ordered_limit, symbol, title, legend_loc) in zip(axes, parameters):
        curves(ax, temperatures, col(key), col(key+'_se'), color, 'Monte-Carlo')
        ax.axhline(ordered_limit, color=GOLD, ls='--', lw=1.3,
                   label=r'Ordre parfait $L1_2$')
        # ax.axvline(peak_temperature, color=INK, ls=':', lw=1.2,
        #            label=rf'Maximum de $c_V$ : {peak_temperature:.0f} K')
        ax.set(xlabel='Température (K)', ylabel=symbol, title=title)
        ax.legend(loc=legend_loc, fontsize=8.5, frameon=True,
                  facecolor='white', edgecolor='none', framealpha=.95)
        ax.margins(y=.12)
    fig.suptitle('Ordre chimique local : paramètres de Cowley', y=.98)
    fig.tight_layout(rect=(0, .065, 1, .94))
    save(fig, folder, '12_cowley',
         'Bandes : ±1 erreur-type. Les références horizontales correspondent à l’ordre parfait L1₂.')
    return folder/'12_cowley.pdf'


def make_figures(result_dir):
    style()
    root = Path(result_dir)
    folder = root/'figures'
    folder.mkdir(exist_ok=True)
    summary = json.loads((root/'summary.json').read_text(encoding='utf8'))
    with (root/'observables.csv').open(encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    col = lambda k: np.array([float(r[k]) for r in rows])
    T = col('T_K')
    peak = int(np.argmax(col('cv_kB')))
    tc = T[peak]
    archives = [np.load(p) for p in sorted(root.glob('chain_*.npz'))]
    samples = np.array([a['samples'] for a in archives])
    n = summary['n_atoms']
    cfg = summary['config']
    status = ('Diagnostics locaux satisfaits' if summary['critical_diagnostics_ok']
              else 'Au moins un diagnostic local à renforcer')
    pick = sorted(set([max(0, peak-2), peak, min(len(T)-1, peak+2)]))

    # 01 — géométrie d'une maille L12.
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121, projection='3d')
    vertices = np.array(list(__import__('itertools').product([0, 1], repeat=3)))
    faces = np.array([[0, .5, .5], [1, .5, .5], [.5, 0, .5],
                      [.5, 1, .5], [.5, .5, 0], [.5, .5, 1]])
    for i, v in enumerate(vertices):
        for w in vertices[i+1:]:
            if np.sum((v-w)**2) == 1:
                ax.plot(*np.array([v, w]).T, color='#8ba1af', lw=1)
    ax.scatter(*vertices.T, s=130, c=GOLD, edgecolors='white', depthshade=False, label='Au : sommets')
    ax.scatter(*faces.T, s=115, c=BLUE, edgecolors='white', depthshade=False, label='Cu : faces')
    ax.set(xlabel='x/a', ylabel='y/a', zlabel='z/a', title='Maille conventionnelle ordonnée')
    ax.set_box_aspect([1, 1, 1]); ax.view_init(22, 32)
    ax.legend(loc='lower left', bbox_to_anchor=(-.1, -.04), fontsize=10)
    right = fig.add_subplot(122); right.axis('off')
    text = [('$\\mathrm{Cu_3Au}$ · structure $L1_2$', 21, .87),
            ('1 Au + 3 Cu par maille', 16, .71),
            ('8 sommets × 1/8 = 1 Au', 13, .60),
            ('6 faces × 1/2 = 3 Cu', 13, .51),
            (f'Boîte : {cfg["L"]}³ mailles = {n} sites', 16, .35),
            ('Positions figées · occupations variables', 12, .21)]
    for label, size, y in text:
        right.text(.03, y, label, fontsize=size, color=INK, transform=right.transAxes)
    fig.suptitle('Le désordre chimique se produit sur un réseau fixe', y=1.02)
    fig.tight_layout()
    save(fig, folder, '01_reseau', 'Les sphères de bord sont partagées avec les mailles voisines.')

    # 02 — ordre longue et courte portée.
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)
    for ax, key, color, title in zip(axes, ['q', 'alpha1', 'alpha2'], [BLUE, RED, GREEN],
                                   ['Ordre à longue portée', 'Premiers voisins', 'Seconds voisins']):
        curves(ax, T, col(key), col(key+'_se'), color, key if key != 'q' else 'Amplitude q')
        ax.axvline(tc, color=INK, ls=':', alpha=.7)
        ax.set(title=title, xlabel='Température (K)', ylabel={'q': 'q', 'alpha1': r'$\alpha_1$', 'alpha2': r'$\alpha_2$'}[key])
    axes[0].axhline(1, color=GOLD, ls='--', lw=1)
    axes[1].axhline(-1/3, color=GOLD, ls='--', lw=1)
    axes[2].axhline(1, color=GOLD, ls='--', lw=1)
    for ax in axes[1:]:
        ax.axhline(-1/(n-1), color=INK, lw=.8, alpha=.5)
    fig.suptitle('Ordre global et préférences chimiques locales', y=1.02)
    fig.tight_layout()
    save(fig, folder, '02_ordres', 'Bandes : ±1 erreur-type (blocs et dispersion entre chaînes). Pointillés : maximum de Cᵥ sur la grille.')

    # 03 — énergie et réponse thermique.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    curves(axes[0], T, col('u_eV'), col('u_eV_se'), BLUE, 'Monte-Carlo')
    axes[0].axhline(-3.813, color=GOLD, ls='--', label='Fondamental exact')
    axes[0].axhline(summary['random_limit_eV_atom_finite_N'], color=RED, ls=':', label='Mélange aléatoire, T → ∞')
    axes[0].set(xlabel='Température (K)', ylabel='Énergie (eV/atome)', title='Énergie interne')
    axes[0].legend(fontsize=9)
    curves(axes[1], T, col('cv_kB'), col('cv_kB_se'), RED, 'Fluctuations d’énergie')
    axes[1].axvline(tc, color=INK, ls=':', label=f'Maximum sur grille : {tc:.0f} K')
    axes[1].set(xlabel='Température (K)', ylabel=r'$c_V/k_B$ par atome', title='Capacité calorifique configuratoire')
    axes[1].legend(fontsize=9)
    fig.suptitle('L’énergie augmente lorsque l’ordre est détruit', y=1.02)
    for ax in axes:
        ax.set_xlim(650, 1300)
    fig.tight_layout()
    save(fig, folder, '03_thermodynamique', 'Seule la contribution des occupations chimiques est incluse ; le réseau ne vibre pas.')

    # 04 — distributions près du maximum.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colors = [BLUE, RED, GOLD]
    for k, color in zip(pick, colors):
        a = samples[:, :, k, :].reshape(-1, 9)
        for ax, values, xlabel in [(axes[0], a[:, 0]/n, 'Énergie (eV/atome)'),
                                   (axes[1], a[:, 2], 'Amplitude d’ordre q')]:
            ax.hist(values, bins=45, density=True, histtype='step', lw=1.8,
                    color=color, label=f'{T[k]:.0f} K')
            ax.set(xlabel=xlabel, ylabel='Densité de probabilité')
    for ax in axes:
        ax.legend()
    fig.suptitle('Les distributions révèlent les fluctuations et la coexistence éventuelle', y=1.02)
    fig.tight_layout()
    save(fig, folder, '04_distributions', 'Chaînes regroupées pour la visualisation ; une double bosse seule ne constitue pas une preuve thermodynamique.')

    # 05 — diagnostics.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for k, c, name in [('rhat_E', BLUE, 'Énergie'), ('rhat_q', RED, 'Ordre q')]:
        axes[0, 0].plot(T, col(k), 'o-', ms=3, color=c, label=name)
    axes[0, 0].axhline(1.05, color=INK, ls='--', label='Seuil de contrôle 1,05')
    axes[0, 0].set(ylabel='R-hat séparé par rangs', title='Accord entre chaînes')
    axes[0, 0].legend(fontsize=9)
    axes[0, 1].semilogy(T, np.maximum(col('ess_min'), 1), 'o-', color=GREEN, ms=3)
    axes[0, 1].axhline(200, color=INK, ls='--')
    axes[0, 1].set(ylabel='Effectif indépendant estimé', title='Autocorrélation prise en compte')
    middle = (T[:-1]+T[1:])/2
    for i, chain in enumerate(summary['chains']):
        axes[1, 0].plot(middle, chain['swap_acceptance'], 'o-', ms=3, label=f'Chaîne {i+1}')
        axes[1, 1].plot(T, chain['local_acceptance'], 'o-', ms=3, label=f'Chaîne {i+1}')
    axes[1, 0].set(ylabel='Fraction acceptée', title='Échanges entre températures', ylim=(-.02, 1.02))
    axes[1, 1].set(ylabel='Fraction acceptée', title='Échanges Au-Cu', ylim=(-.02, 1.02))
    for ax in axes.flat:
        ax.set_xlabel('Température (K)')
        ax.axvline(tc, color=INK, ls=':', alpha=.4)
    axes[1, 0].legend(fontsize=9)
    fig.suptitle('Contrôler la simulation avant d’interpréter les courbes', y=1.01)
    fig.tight_layout()
    save(fig, folder, '05_diagnostics', status+' · R-hat et ESS sont des indicateurs, pas des garanties absolues d’équilibre.')

    # 06 — traces thermiques.
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for c, a in enumerate(archives):
        burn = a['burn_samples'][:, peak]
        prod = a['samples'][:, peak]
        whole = np.concatenate([burn, prod])
        sweeps = np.arange(1, len(whole)+1)*cfg['sample_every']
        for ax, values in [(axes[0], whole[:, 0]/n), (axes[1], whole[:, 2])]:
            ax.plot(sweeps, values, lw=.65, alpha=.8, label=f'Chaîne {c+1}')
    for ax in axes:
        ax.axvspan(0, cfg['burn_sweeps'], color='#d4dde4', alpha=.55, label='Équilibration écartée')
        ax.axvline(cfg['burn_sweeps'], color=INK, ls='--')
    axes[0].set(ylabel='Énergie (eV/atome)')
    axes[1].set(ylabel='Amplitude q', xlabel='Balayages Monte-Carlo (N échanges proposés chacun)')
    axes[0].legend(ncols=3, fontsize=9)
    fig.suptitle(f'Historique au voisinage du maximum thermique · {tc:.0f} K', y=1.01)
    fig.tight_layout()
    save(fig, folder, '06_traces', 'Une trace appartient à une température fixe : les configurations peuvent y arriver par échange de répliques.')

    # 07 — configurations effectivement obtenues.
    fig = plt.figure(figsize=(15, 5))
    xyz = archives[0]['xyz']/2
    for panel, k in enumerate([0, peak, len(T)-1], 1):
        state = archives[0]['final_states'][k]
        ax = fig.add_subplot(1, 3, panel, projection='3d')
        mask = state == 1
        ax.scatter(*xyz[~mask].T, s=5, c=BLUE, alpha=.10, depthshade=False)
        ax.scatter(*xyz[mask].T, s=15, c=GOLD, alpha=.95, depthshade=False)
        ax.set(title=f'{T[k]:.0f} K', xlabel='x/a', ylabel='y/a', zlabel='z/a')
        ax.set_box_aspect([1, 1, 1]); ax.view_init(21, 33)
        ax.set_xticks([0, cfg['L']/2, cfg['L']]); ax.set_yticks([0, cfg['L']/2, cfg['L']])
    fig.suptitle('Instantanés issus du calcul · Au en or, Cu en bleu', y=.99)
    fig.tight_layout(rect=(0, .09, 1, .93))
    save(fig, folder, '07_configurations', 'Les occupations changent ; les positions atomiques sont identiques sur les trois panneaux.')

    # 08 — champ moyen : deux minima et saut de l'ordre.
    mf = coexistence()
    eta = np.linspace(0, 1, 1001)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for temp, color in [(mf['T_MF_K']-120, BLUE), (mf['T_MF_K'], RED), (mf['T_MF_K']+120, GOLD)]:
        f = -A*eta**2+KB*temp*(g(eta)-g(0))
        axes[0].plot(eta, 1000*f, color=color, label=f'{temp:.0f} K')
    axes[0].axhline(0, color=INK, lw=.8)
    axes[0].set(xlabel=r'$\eta$', ylabel='f(η) − f(0) (meV/atome)', title='Énergie libre de Bragg-Williams')
    axes[0].legend(fontsize=9)
    temps = np.linspace(100, 2400, 400)
    eq = [eta[np.argmin(-A*eta**2+KB*temp*g(eta))] for temp in temps]
    axes[1].plot(temps, eq, color=BLUE, lw=2, label='Champ moyen')
    curves(axes[1], T, col('q'), col('q_se'), RED, 'Monte-Carlo : q')
    axes[1].axvline(mf['T_MF_K'], color=BLUE, ls=':', label=f'T coexistence MF ≈ {mf["T_MF_K"]:.0f} K')
    axes[1].set(xlabel='Température (K)', ylabel='Amplitude d’ordre', title='Une approximation à comparer au calcul')
    axes[1].legend(fontsize=9)
    fig.suptitle('Le champ moyen néglige les corrélations spatiales', y=1.02)
    fig.tight_layout()
    save(fig, folder, '08_champ_moyen', 'Dans la famille à un sous-réseau privilégié, q = η ; les deux courbes comparent des amplitudes.')

    # 09 — quatre sous-réseaux, et eta signé.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    means = samples.mean(axis=(0, 1))
    for sub, color in enumerate([GOLD, BLUE, RED, GREEN]):
        axes[0].plot(T, means[:, 5+sub], 'o-', ms=3, color=color, label=f'Sous-réseau {sub}')
    axes[0].axhline(.25, color=INK, ls=':')
    axes[0].set(xlabel='Température (K)', ylabel='Occupation moyenne par Au', title='Moyennes sur les mesures et les deux chaînes')
    axes[0].legend(fontsize=9)
    for c in range(len(archives)):
        axes[1].plot(T, samples[c, :, :, 1].mean(axis=0), 'o-', ms=3, label=f'η signé · chaîne {c+1}')
    axes[1].plot(T, col('q'), color=INK, lw=2, label='Amplitude q · moyenne')
    axes[1].set(xlabel='Température (K)', ylabel='Paramètre d’ordre', title='Une variante différente peut changer le signe de η')
    axes[1].legend(fontsize=9)
    fig.suptitle('La symétrie entre quatre variantes doit être prise en compte', y=1.02)
    fig.tight_layout()
    save(fig, folder, '09_sous_reseaux', 'Une moyenne sur plusieurs variantes peut égaliser les occupations malgré un ordre instantané fort.')

    # 10 — circulation des répliques ; une sélection évite la surcharge.
    fig, ax = plt.subplots(figsize=(12, 5))
    trace = archives[0]['labels']
    x = np.arange(len(trace))*cfg['sample_every']
    for label in np.linspace(0, len(T)-1, 4, dtype=int):
        ix = np.argmax(trace == label, axis=1)
        ax.plot(x, T[ix], lw=.8, alpha=.85, label=f'Identité {label}')
    ax.set(xlabel='Balayages de production', ylabel='Température de la réplique (K)')
    ax.legend(ncols=4, fontsize=9)
    fig.suptitle('Les configurations voyagent entre les températures', y=1.01)
    fig.tight_layout()
    save(fig, folder, '10_repliques', 'Les répliques chaudes facilitent la sortie des états métastables ; les trajectoires restent corrélées.')

    # 11 — agrandissement de la zone critique, sans interpolation lissée.
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    for ax, key, color, title, ylabel in [
        (axes[0, 0], 'q', BLUE, 'Perte de l’ordre global', 'Amplitude q'),
        (axes[0, 1], 'u_eV', GREEN, 'Énergie interne', 'Énergie (eV/atome)'),
        (axes[1, 0], 'cv_kB', RED, 'Réponse thermique', r'$c_V/k_B$ par atome'),
        (axes[1, 1], 'Xq_per_eV', GOLD, 'Fluctuations de l’amplitude', r'$X_q$ (eV$^{-1}$)')]:
        curves(ax, T, col(key), col(key+'_se'), color, title)
        ax.axvline(tc, color=INK, ls=':')
        ax.axvspan(*summary['T_cv_neighbor_bracket_K'], color=INK, alpha=.06)
        ax.set(xlim=(tc-100, tc+100), title=title, ylabel=ylabel, xlabel='Température (K)')
    fig.suptitle(f'La grille fine résout le voisinage du maximum · {tc:.0f} K', y=1.01)
    fig.tight_layout(rect=(0, .07, 1, .95))
    save(fig, folder, '11_zoom_transition', 'Zone grisée : températures voisines du maximum, pas un intervalle de confiance. Bandes colorées : ±1 erreur-type.')
    make_cowley_figure(root)
    from .pedagogy import make_pedagogical_figures
    make_pedagogical_figures(root)
    for a in archives:
        a.close()
    return sorted(str(p) for p in folder.glob('*.pdf'))
