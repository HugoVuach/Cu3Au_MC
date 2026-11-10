"""Comparaison de tailles sans extrapolation ni interpolation des températures."""
import csv
import json
from pathlib import Path
import numpy as np
from .simulation import save_csv


def compare_sizes(reference, larger):
    from .figures import plt, style, save, measurement_points, BLUE, RED, INK
    roots = [Path(reference), Path(larger)]
    data = []
    summaries = []
    report_rows = []
    for root in roots:
        summary = json.loads((root/'summary.json').read_text(encoding='utf8'))
        summaries.append(summary)
        with (root/'observables.csv').open(encoding='utf8') as stream:
            rows = list(csv.DictReader(stream))
        data.append(rows)
        critical = [rows[i] for i in summary['critical_window_rows']]
        accepted = [np.mean(c['swap_acceptance']) for c in summary['chains']]
        report_rows.append(dict(L=summary['config']['L'], N=summary['n_atoms'],
            T_peak_grid_K=summary['T_cv_peak_grid_K'],
            critical_diagnostics_ok=summary['critical_diagnostics_ok'],
            min_critical_ess=min(float(r['ess_min']) for r in critical),
            max_critical_rhat=max(max(float(r['rhat_E']), float(r['rhat_q'])) for r in critical),
            mean_swap_acceptance=float(np.mean(accepted)),
            total_roundtrips=sum(sum(c['roundtrips']) for c in summary['chains']),
            spatial_capture_every=summary['config'].get('store_states_every', 0)))
    same_grid = summaries[0]['config']['temperatures'] == summaries[1]['config']['temperatures']
    common = sorted(set(summaries[0]['config']['temperatures']) & set(summaries[1]['config']['temperatures']))
    if not common:
        raise ValueError('No common sampled temperatures for size comparison.')
    small, large = (s['config']['L'] for s in summaries)
    folder = roots[1]/f'comparison_L{small}_L{large}'
    folder.mkdir(exist_ok=True)
    style()
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    variables = [('u_eV', 'Énergie (eV/atome)'), ('q', 'Amplitude q'),
                 ('cv_kB', 'cV/kB par atome'), ('alpha1', 'Cowley α₁'), ('alpha2', 'Cowley α₂')]
    for c, (rows, summary, color) in enumerate(zip(data, summaries, [BLUE, RED])):
        t = np.array([float(r['T_K']) for r in rows])
        label = f"L = {summary['config']['L']}"
        for ax, (key, ylabel) in zip(axes.flat, variables):
            y = np.array([float(r[key]) for r in rows])
            err = np.array([float(r[key+'_se']) for r in rows])
            measurement_points(ax, t, y, color, label, err, series=c)
            ax.set(xlim=(1000, 1100), xlabel='Température (K)', ylabel=ylabel)
        rh = np.array([max(float(r['rhat_E']), float(r['rhat_q'])) for r in rows])
        measurement_points(axes[1, 2], t, rh, color, label, series=c)
    axes[1, 2].axhline(1.05, color=INK, ls=':')
    axes[1, 2].set(xlim=(1000, 1100), xlabel='Température (K)', ylabel='Maximum R-hat(E,q)')
    for ax in axes.flat:
        ax.legend(fontsize=9)
    fig.suptitle(f'Comparaison L = {small} et L = {large} : résultats et contrôle d’échantillonnage', fontsize=17)
    fig.tight_layout(rect=(0, .07, 1, .95))
    save(fig, folder, 'comparaison_thermodynamique',
         ('Même grille ; ' if same_grid else 'Grilles distinctes, points effectivement simulés ; ')+
         'examiner ESS, circulation et budgets, sans interpolation ni extrapolation à deux tailles.')
    save_csv(folder/'diagnostics_comparison.csv', report_rows)
    comparison = dict(sizes=report_rows,
        same_temperature_grid=same_grid, common_temperatures_K=common,
        sweep_budgets=[{k:s['config'][k] for k in ['burn_sweeps','sample_sweeps','sample_every']}
                      for s in summaries],
        peak_shift_grid_K=report_rows[1]['T_peak_grid_K']-report_rows[0]['T_peak_grid_K'],
        caveats=['Same sweep count does not imply the same effective sample count.',
                 'Two sizes do not establish a thermodynamic extrapolation.',
                 'Global exchanges, rigid lattice and fixed pair interactions are unchanged.',
                 'Spatial states are saved for larger sizes; dedicated L6 correlation data use a different sampling protocol.',
                 'Temperature grids and sweep budgets may differ; only actually measured points are drawn.'])
    (folder/'summary.json').write_text(json.dumps(comparison, indent=2), encoding='utf8')
    return folder


def write_report_include(reference, larger):
    """Intègre uniquement les résultats terminés et validés, sans nouveau calcul MC."""
    from .reporting import _copy_atomic
    reference, larger = Path(reference), Path(larger)
    target = reference/'report'
    comparison = json.loads((larger/'comparison_L6_L9/summary.json').read_text(encoding='utf8'))
    validation = json.loads((larger/'validation.json').read_text(encoding='utf8'))
    if validation['technical_validation'] != 'passed':
        raise ValueError('L9 technical validation must pass before report integration.')
    summaries = [json.loads((r/'summary.json').read_text(encoding='utf8')) for r in [reference, larger]]
    if summaries[0]['config']['temperatures'] != summaries[1]['config']['temperatures']:
        raise ValueError('Mismatched temperature grids in size report.')
    data = []
    for root in [reference, larger]:
        with (root/'observables.csv').open(encoding='utf8') as stream:
            data.append(list(csv.DictReader(stream)))
    sizes = comparison['sizes']
    cfg = summaries[1]['config']
    for s, summary in zip(sizes, summaries):
        if s['T_peak_grid_K'] != summary['T_cv_peak_grid_K']:
            raise ValueError('Size comparison is stale relative to simulation summary.')
    peaks = [json.loads((root/'peak_uncertainty.json').read_text(encoding='utf8'))
             for root in [reference, larger]]
    base_cases = [next(c for c in p['cases'] if c['block_samples'] == p['base_block_samples']) for p in peaks]
    peak_rows = [next(r for r in rows if float(r['T_K']) == s['T_peak_grid_K'])
                 for rows, s in zip(data, sizes)]
    hot_rows = [rows[-1] for rows in data]
    critical9 = [data[1][i] for i in summaries[1]['critical_window_rows']]
    capture_counts = []
    for path in sorted(larger.glob('chain_*.npz')):
        with np.load(path) as archive:
            capture_counts.append(len(archive['spatial_sample_sweeps']))
    def interval(case):
        return '$['+','.join(f'{x:.0f}' for x in case['ci95_K'])+']$'
    lines = [r'\subsubsection{Périmètre et protocole de la comparaison}',
        r"La campagne L9 est terminée ; les deux archives ont passé les contrôles "
        r"d'énergie, de composition et de cohérence des observables. Cette "
        r"validation technique n'est pas une validation de l'équilibre statistique.",
        f"La boîte contient ${summaries[1]['n_atoms']}$ sites, dont "
        f"${cfg['L']**3}$ Au et ${3*cfg['L']**3}$ Cu. Les "
        f"{len(cfg['temperatures'])} températures, les {cfg['burn_sweeps']} balayages "
        f"d'équilibration et les {cfg['sample_sweeps']} de production sont "
        f"identiques au protocole L6 ; les observables sont mesurées tous les "
        f"{cfg['sample_every']} balayages et les sessions d'échange ont lieu "
        f"tous les {cfg['exchange_every']} balayages. Les deux graines sont "
        f"{cfg['seeds'][0]} et {cfg['seeds'][1]}.",
        r"Les deux calculs utilisent le même Hamiltonien et le même ensemble "
        r"canonique à composition fixée. Les trajectoires L6/L9 ne sont pas "
        r"appariées configuration par configuration. Augmenter la boîte "
        r"augmente le nombre de propositions par balayage ; un budget "
        r"identique de balayages n'implique pas un effectif indépendant identique.",
        r'\subsubsection{Résultats thermodynamiques et bilan des diagnostics}',
        r'\begin{center}\small\begin{tabular}{lrr}\toprule',
        r'Indicateur & $L=6$ & $L=9$\\\midrule',
        f"Sites & {sizes[0]['N']} & {sizes[1]['N']}"+r'\\',
        f"Maximum de $c_V$ sur la grille (K) & {sizes[0]['T_peak_grid_K']:.0f} & {sizes[1]['T_peak_grid_K']:.0f}"+r'\\',
        'Intervalle bootstrap 95\,\% conditionnel (K) & '+interval(base_cases[0])+' & '+interval(base_cases[1])+r'\\',
        f"Hauteur du pic $c_V/k_B$ par atome & {float(peak_rows[0]['cv_kB']):.3f} & {float(peak_rows[1]['cv_kB']):.3f}"+r'\\',
        f"ESS minimal près du pic & {sizes[0]['min_critical_ess']:.1f} & {sizes[1]['min_critical_ess']:.1f}"+r'\\',
        f"Maximum de R-hat$(E,q)$ près du pic & {sizes[0]['max_critical_rhat']:.4f} & {sizes[1]['max_critical_rhat']:.4f}"+r'\\',
        f"Acceptation moyenne des échanges de répliques & {100*sizes[0]['mean_swap_acceptance']:.1f}\,\% & {100*sizes[1]['mean_swap_acceptance']:.1f}\,\%"+r'\\',
        f"Allers-retours froid--chaud--froid (total) & {sizes[0]['total_roundtrips']} & {sizes[1]['total_roundtrips']}"+r'\\',
        r'Diagnostics autour du pic & satisfaits & insuffisants\\',
        r'\bottomrule\end{tabular}\end{center}',
        f"Le maximum se déplace de {comparison['peak_shift_grid_K']:.0f} K, "
        r"soit une seule maille de température. Les deux intervalles "
        r"bootstrap de référence couvrent $[1050,1055]$ K. Ils sont "
        r"conditionnels aux trajectoires observées et ne constituent "
        r"pas un test du déplacement à taille infinie. Pour L9, "
        f"la longueur de bloc de référence vaut {peaks[1]['base_block_samples']} "
        f"mesures, soit {peaks[1]['base_block_sweeps']} balayages, avec "
        f"{peaks[1]['bootstrap']} tirages. Les intervalles avec des blocs "
        r"plus courts peuvent même se réduire à une température : cet "
        r"effet de grille n'ajoute aucune précision physique.",
        r"Le pic plus haut et la variation plus abrupte de $q$ et de l'énergie "
        r"sont compatibles avec une transition moins arrondie dans une boîte "
        r"plus grande. Cependant, leur amplitude reste sensible au temps passé "
        r"dans les différents bassins : ces deux tailles ne mesurent pas "
        r"un exposant, une chaleur latente ou une température critique extrapolée.",
        r'\fig{tailles_01_comparaison}{1}{Résultats L6/L9 sur une grille identique : énergie et chaleur par atome, ordre global, Cowley et R-hat. Symboles distincts, sans interpolation des températures.}',
        r'\begin{center}\small\begin{tabular}{rrrrrl}\toprule',
        r'$T$ (K), L9 & ESS & R-hat max & Blocs/chaîne & Écart demi-séries & Contrôle\\\midrule']
    for row in critical9:
        passed = str(row['diagnostic_ok']).lower() == 'true'
        lines.append(f"{float(row['T_K']):.0f} & {float(row['ess_min']):.1f} & "
                     f"{max(float(row['rhat_E']),float(row['rhat_q'])):.4f} & "
                     f"{float(row['blocks_per_chain']):.0f} & {float(row['stationarity_z']):.2f} & "
                     +('satisfait' if passed else 'insuffisant')+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r"Le critère R-hat$<1.05$ passe, mais il ne suffit pas. À 1050 et "
        r"1055 K, l'ESS est inférieur à 200 et seulement 16 blocs par chaîne "
        r"sont disponibles, contre 20 demandés. À 1060 K, l'ESS échoue "
        r"encore malgré un R-hat proche de 1. L'absence d'aller-retour "
        r"complet ne prouve pas que tout état est gelé, mais signale une "
        r"exploration globale en température moins convaincante. Une "
        r"acceptation moyenne de 61\,\% ne garantit donc pas une circulation "
        r"efficace : les transitions entre bassins peuvent rester rares.",
        r'\subsubsection{Historiques et distributions : indices de coexistence, pas preuve finale}',
        r"À 1055 K, les histogrammes regroupés montrent deux populations "
        r"d'énergie et d'ordre ; les traces révèlent des séjours corrélés "
        r"dans des régimes plus ou moins ordonnés. Cela est compatible "
        r"avec une coexistence à taille finie et explique une partie du "
        r"ralentissement. Le regroupement des chaînes et le nombre limité "
        r"de passages empêchent de garantir les poids d'équilibre des "
        r"deux populations. Les valeurs de variance et les barres d'erreur "
        r"doivent être interprétées avec cette réserve.",
        r'\fig{tailles_02_historiques}{.95}{Traces L9 à 1055 K : équilibration écartée en gris et production ensuite. À température fixe, une configuration peut également changer par échange de répliques ; ce tracé n’est pas une cinétique atomique physique.}',
        r'\fig{tailles_03_distributions}{.95}{Distributions L9 à 1045, 1055 et 1065 K. Deux populations à 1055 K sont visibles ; leur poids reste conditionnel aux états visités.}',
        r'\subsubsection{Ce que L9 améliore et ce qu’il ne corrige pas}',
        f"À 4000 K, $q$ diminue de {float(hot_rows[0]['q']):.4f} "
        f"à {float(hot_rows[1]['q']):.4f}. Le rapport "
        f"${float(hot_rows[1]['q'])/float(hot_rows[0]['q']):.3f}$ est proche "
        f"de $\sqrt{{864/2916}}={np.sqrt(864/2916):.3f}$ : cette tendance "
        r"est compatible avec la diminution du biais positif d'une norme "
        r"dans un système fini, plutôt qu'avec l'apparition d'un ordre "
        r"supplémentaire à haute température. Les paramètres de Cowley "
        r"restent proches entre tailles ; les préférences locales ne "
        r"disparaissent pas avec le seul biais de taille de $q$.",
        r"Le volume disponible et la résolution spatiale augmentent, "
        r"mais le mélange près de la transition se dégrade avec le budget "
        r"utilisé. Le passage à L9 n'efface donc pas les insuffisances "
        r"d'échantillonnage ; il les rend plus visibles. La température "
        r"reste éloignée de celle du matériau réel et les hypothèses "
        r"de réseau rigide, de paires constantes et d'échanges globaux "
        r"restent inchangées.",
        f"Les archives L9 conservent {capture_counts[0]} états par température "
        f"et par chaîne, espacés de {cfg['store_states_every']} balayages, "
        f"soit {sum(capture_counts)*len(cfg['temperatures'])} configurations "
        r"spatiales. Leur présence rend une future analyse de corrélations "
        r"possible, mais aucune courbe $\xi(T)$ ni animation L9 n'est "
        r"calculée dans cette intégration : les résultats spatiaux et les "
        r"GIF de la partie précédente restent explicitement ceux de L6.",
        r"Avant une interprétation quantitative du déplacement de transition, "
        r"il faut allonger la production autour du pic, vérifier la stabilité "
        r"par blocs plus longs et améliorer la circulation des répliques, "
        r"par exemple en densifiant la grille là où se situent les "
        r"goulots d'étranglement. Il faudra ensuite comparer plusieurs "
        r"tailles à précision statistique contrôlée.",
        r"Les données viennent de \texttt{results/production\_L9}, "
        r"et la synthèse comparative de \texttt{comparison\_L6\_L9}. "
        r"Les figures propres à L9 sont intégrées ici sous des noms "
        r"préfixés \texttt{tailles\_}, sans remplacer les figures L6."]
    (target/'bonus_sizes.tex').write_text('\n'.join(lines)+'\n', encoding='utf8')
    files = [(larger/'comparison_L6_L9/comparaison_thermodynamique.pdf', 'tailles_01_comparaison.pdf'),
             (larger/'figures/06_traces.pdf', 'tailles_02_historiques.pdf'),
             (larger/'figures/04_distributions.pdf', 'tailles_03_distributions.pdf')]
    for source, name in files:
        _copy_atomic(source, target/'figures'/name)
    return target/'bonus_sizes.tex'
