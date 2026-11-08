"""Comparative heating/cooling plots and the section 10 report include."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
from .thermal_cycle import load_cycle_archives, schedule


BRANCHES = [('heating', 'Chauffage'), ('cooling', 'Refroidissement')]


def make_cycle_figures(output):
    from .figures import (plt, style, save, RED, BLUE, GOLD, GREEN, INK,
                          measurement_points, transition_window)
    out = Path(output)
    cfg, data = load_cycle_archives(out)
    folder = out/'figures'
    folder.mkdir(exist_ok=True)
    colors = {'heating': RED, 'cooling': BLUE}
    rows = {}
    for branch, _ in BRANCHES:
        with (out/f'observables_{branch}.csv').open(encoding='utf8') as f:
            rows[branch] = list(csv.DictReader(f))
    t = np.asarray(cfg['temperatures'], dtype=float)
    n = 4*cfg['L']**3
    col = lambda branch, key: np.array([float(r[key]) for r in rows[branch]])
    peaks = [t[np.argmax(col(branch, 'cv_kB'))] for branch, _ in BRANCHES]
    # Même fenêtre pour les deux branches, centrée sur le maximum au chauffage.
    zoom_limits = transition_window(peaks[0])
    style()

    def compare(ax, key, ylabel, zoom=False):
        for index, (branch, label) in enumerate(BRANCHES):
            y = col(branch, key)
            error = col(branch, key+'_se') if key+'_se' in rows[branch][0] else None
            measurement_points(ax, t, y, colors[branch], label, error, series=index)
        ax.set(xlabel='Température (K)', ylabel=ylabel)
        if zoom:
            ax.set_xlim(*zoom_limits)
        ax.legend(fontsize=9)

    # B01: global comparison and the temperature window requested for figure 03.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for row in range(2):
        compare(axes[row, 0], 'u_eV', 'Énergie (eV/atome)', zoom=row == 1)
        compare(axes[row, 1], 'cv_kB', r'$c_V^{fluct}/k_B$ par atome', zoom=row == 1)
    axes[0, 0].set_title('Énergie : cycle complet')
    axes[0, 1].set_title('Réponse par fluctuations : cycle complet')
    axes[1, 0].set_title('Énergie : zone de transition')
    axes[1, 1].set_title('Fluctuations : zone de transition')
    fig.suptitle('Cycle thermique : énergie et fluctuations', y=.99)
    fig.tight_layout(rect=(0, .055, 1, .95))
    save(fig, folder, 'bonus_01_thermodynamique',
         'Barres : ±1 erreur-type par branche ; cercles : chauffage, carrés : refroidissement ; stationnarité requise pour les fluctuations.')

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, (key, label) in zip(axes.flat, [
        ('q', 'Amplitude q'), ('eta', 'Ordre signé η'), ('alpha1', 'Cowley α₁'),
        ('alpha2', 'Cowley α₂'), ('Xq_per_eV', 'Fluctuations Xq (eV⁻¹)'), ('binder_q', 'Cumulant Uq')]):
        compare(ax, key, label, zoom=True)
    fig.suptitle('Ordre global et local au chauffage et au refroidissement', y=.99, fontsize=17)
    fig.tight_layout(rect=(0, .06, 1, .95))
    save(fig, folder, 'bonus_02_ordres', 'η dépend de la variante choisie ; q est l’amplitude à privilégier pour comparer les branches.')

    # Ne pas fusionner les variantes ordonnées choisies par les cycles indépendants.
    # Lignes : sens du parcours ; colonnes : cycles, avec la même clé μ partout.
    cycles = len(cfg['seeds'])
    fig, axes = plt.subplots(2, cycles, figsize=(6*cycles, 8),
                             sharex=True, sharey=True, squeeze=False)
    for row, (branch, label) in enumerate(BRANCHES):
        for c in range(cycles):
            ax = axes[row, c]
            occupations = np.array([a['samples'][:, 5:9].mean(axis=0)
                                    for a in data[branch][c]])
            for mu, color in enumerate([GOLD, BLUE, RED, GREEN]):
                measurement_points(ax, t, occupations[:, mu], color,
                                   f'μ = {mu}', series=mu)
            ax.axhline(.25, color=INK, ls=':', lw=1)
            ax.set(xlabel='Température (K)', ylabel='Occupation temporelle moyenne par Au',
                   title=f'{label} · cycle {c+1}', ylim=(-.04, 1.06))
            ax.legend(fontsize=9, ncols=2)
    fig.suptitle('Occupations des sous-réseaux : chaque cycle séparément', y=.99, fontsize=17)
    fig.tight_layout(rect=(0, .085, 1, .95))
    save(fig, folder, 'bonus_03_sous_reseaux',
         'Chaque point : moyenne temporelle d’un seul cycle, sans moyenne entre cycles.\n'
         'Les couleurs et symboles identifient μ ; des changements de variante peuvent subsister dans une moyenne temporelle.')

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, (index, label, divisor) in zip(axes.flat,
            [(0, 'Δu (eV/atome)', n), (2, 'Δq', 1), (3, 'Δα₁', 1), (4, 'Δα₂', 1)]):
        gaps = []
        for c in range(len(cfg['seeds'])):
            gap = np.array([(data['cooling'][c][k]['samples'][:, index].mean()-
                             data['heating'][c][k]['samples'][:, index].mean())/divisor
                            for k in range(len(t))])
            gaps.append(gap)
            measurement_points(ax, t, gap, [BLUE, RED][c % 2], f'Cycle {c+1}', series=c)
        measurement_points(ax, t, np.mean(gaps, axis=0), INK, 'Moyenne appariée', series=len(gaps))
        ax.axhline(0, color=INK, ls=':', lw=1)
        ax.set(xlabel='Température (K)', ylabel=label)
        ax.set_xlim(*zoom_limits)
        ax.legend(fontsize=8)
    fig.suptitle('Mémoire du trajet : refroidissement moins chauffage', y=.99, fontsize=17)
    fig.tight_layout(rect=(0, .06, 1, .95))
    save(fig, folder, 'bonus_04_ecarts', 'Les branches sont appariées ; les symboles distinguent les cycles et leur moyenne, sans test de significativité.')

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    for c in range(len(cfg['seeds'])):
        elapsed = 0
        for branch, label in BRANCHES:
            parts = [data[branch][c][index] for b, index, _ in schedule(cfg) if b == branch]
            series = np.concatenate([np.concatenate([a['burn_samples'], a['samples']]) for a in parts])
            temp = np.concatenate([np.full(len(a['burn_samples'])+len(a['samples']), float(a['temperature_K'])) for a in parts])
            x = elapsed+np.arange(1, len(series)+1)*cfg['sample_every']
            stride = max(1, len(series)//3000)
            for ax, values in zip(axes, [temp, series[:, 0]/n, series[:, 2]]):
                ax.plot(x[::stride], values[::stride], color=colors[branch],
                        ls='-' if c == 0 else '--', lw=.65, alpha=.8,
                        label=f'{label}, cycle {c+1}')
            elapsed = int(x[-1])
    for ax, ylabel in zip(axes, ['Température (K)', 'Énergie (eV/atome)', 'Amplitude q']):
        ax.set_ylabel(ylabel)
    axes[0].legend(ncols=2, fontsize=8)
    axes[-1].set_xlabel('Balayages cumulés du cycle, équilibration comprise')
    fig.suptitle('Historique complet du chauffage puis du refroidissement', y=.99, fontsize=17)
    fig.tight_layout(rect=(0, .06, 1, .95))
    save(fig, folder, 'bonus_05_traces', 'Affichage décimé pour la lisibilité ; toutes les mesures enregistrées, y compris les transitoires, restent dans les archives NPZ.')

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for index, (branch, label) in enumerate(BRANCHES):
        color = colors[branch]
        for ax, values in zip(axes.flat, [
                np.maximum(col(branch, 'rhat_E'), col(branch, 'rhat_q')),
                np.maximum(1, col(branch, 'ess_min')),
                col(branch, 'tau_max_samples'), col(branch, 'local_acceptance')]):
            measurement_points(ax, t, values, color, label, series=index)
    axes[0, 1].set_yscale('log')
    axes[0, 0].axhline(1.05, color=INK, ls='--'); axes[0, 0].legend(fontsize=9)
    axes[0, 1].axhline(200, color=INK, ls='--')
    for ax, label in zip(axes.flat, ['Maximum de R-hat(E,q)', 'ESS minimal estimé', 'τ maximal (mesures)', 'Fraction de mouvements acceptés']):
        ax.set(xlabel='Température (K)', ylabel=label)
    fig.suptitle('Diagnostics séparés des deux branches', y=.99)
    fig.tight_layout(rect=(0, .06, 1, .95))
    save(fig, folder, 'bonus_06_diagnostics', 'Une hystérésis peut refléter une relaxation insuffisante ; des séries constantes ne prouvent pas un mélange efficace.')

    k = int(np.argmin(abs(t-cfg.get('comparison_temperature_K', 1050))))
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), sharex=True, sharey=True)
    for row, (branch, label) in enumerate(BRANCHES):
        a = data[branch][0][k]
        for col_index, z in enumerate([0, 1]):
            plane = a['xyz'][:, 2] == z
            xy = a['xyz'][plane, :2]/2
            au = a['final_state'][plane].astype(bool)
            axes[row, col_index].scatter(xy[:, 0], xy[:, 1], c=np.where(au, GOLD, BLUE),
                                        s=65, edgecolors=INK, linewidths=.3)
            axes[row, col_index].set(title=f'{label} ; z/a = {z/2:g}', xlabel='x/a', ylabel='y/a')
            axes[row, col_index].set_aspect('equal')
    fig.suptitle(f'Instantanés à {t[k]:.0f} K : Au en or, Cu en bleu', y=.99, fontsize=17)
    fig.tight_layout(rect=(0, .06, 1, .94))
    save(fig, folder, 'bonus_07_configurations', 'Premier cycle, même température, deux plans adjacents ; les instantanés ne représentent pas les occupations moyennes.')

    cv_up = col('heating', 'cv_kB'); cv_down = col('cooling', 'cv_kB')
    middle = int(np.argmin(abs(t-(t[np.argmax(cv_up)]+t[np.argmax(cv_down)])/2)))
    picks = sorted(set([max(0, middle-1), middle, min(len(t)-1, middle+1)]))
    fig, axes = plt.subplots(2, len(picks), figsize=(12, 7), squeeze=False)
    for column, index in enumerate(picks):
        for row, (obs, divisor, label) in enumerate([(0, n, 'Énergie (eV/atome)'), (2, 1, 'Amplitude q')]):
            values = {branch: np.concatenate([data[branch][c][index]['samples'][:, obs]/divisor
                       for c in range(len(cfg['seeds']))]) for branch, _ in BRANCHES}
            bins = np.histogram_bin_edges(np.concatenate(list(values.values())), bins=35)
            for branch, legend in BRANCHES:
                axes[row, column].hist(values[branch], bins=bins, density=True,
                    histtype='step', lw=1.6, linestyle='-' if branch == 'heating' else '--',
                    color=colors[branch], label=legend)
            axes[row, column].set(xlabel=label, ylabel='Densité', title=f'{t[index]:.0f} K')
            if column == 0:
                axes[row, column].legend(fontsize=8)
    fig.suptitle('Distributions au voisinage des maxima de fluctuations', y=.99, fontsize=17)
    fig.tight_layout(rect=(0, .06, 1, .95))
    save(fig, folder, 'bonus_08_distributions', 'Cycles regroupés pour la visualisation ; à chaque température les deux branches utilisent les mêmes classes d’histogramme.')
    return sorted(folder.glob('bonus_*.pdf'))


def write_bonus_include(output, reference_dir, complete=False):
    out = Path(output)
    cfg = json.loads((out/'config.json').read_text(encoding='utf8'))
    target = Path(reference_dir)/'report'
    target.mkdir(parents=True, exist_ok=True)
    n = 4*cfg['L']**3
    lines = [r'\subsection{Protocole séquentiel et conservation de la mémoire}',
        f"Le cycle porte sur $L={cfg['L']}$, soit ${n}$ sites, avec "
        f"{len(cfg['seeds'])} trajectoires indépendantes et {len(cfg['temperatures'])} "
        f"températures entre {min(cfg['temperatures']):g} et {max(cfg['temperatures']):g} K.",
        f"À chaque température et dans chaque sens, on écarte {cfg['burn_sweeps']} "
        f"balayages initiaux, puis on mesure pendant {cfg['sample_sweeps']} balayages, "
        f"à raison d'une mesure tous les {cfg['sample_every']} balayages.",
        r"Chaque trajectoire part d'une variante ordonnée. La configuration finale d'un palier "
        r"est transmise au suivant ; après le palier le plus chaud, le programme parcourt les "
        r"températures en sens inverse sans réinitialiser les occupations. Le sommet chaud "
        r"fait l'objet d'un second palier mesuré au début de la branche retour : aucune mesure n'est dupliquée.",
        r"Il n'y a pas d'échange de répliques dans ce bonus : l'objectif est de comparer "
        r"deux histoires thermiques avec le même budget de balayages. Les déplacements "
        r"restent des échanges Au--Cu globaux ; ce cycle ne représente donc pas une vitesse "
        r"de chauffage expérimentale ni une cinétique de diffusion réelle.",
        r'\subsection{Mesures sauvegardées et reprise}',
        r"Chaque archive de palier conserve les séries d'équilibration et de production "
        r"$(E,\eta,q,\alpha_1,\alpha_2,p_0,p_1,p_2,p_3)$, les états initial et final, "
        r"les compteurs d'acceptation et l'état du générateur aléatoire. L'énergie et la "
        r"composition sont vérifiées après chaque segment. Les archives sont écrites "
        r"atomiquement après un palier complet. L'option \texttt{--resume} reprend au "
        r"dernier palier achevé ; un palier interrompu est recalculé depuis son début.",
        r'\par', r'\begin{algocode}{Pseudocode du cycle chauffage--refroidissement}',
        r'\small\begin{verbatim}',
        'POUR chaque graine independante :',
        '    initialiser une configuration ordonnee',
        '    POUR branche dans [chauffage, refroidissement] :',
        '        parcourir les temperatures dans le sens de la branche',
        '        POUR chaque temperature :',
        '            conserver la configuration du palier precedent',
        '            equilibrer par Metropolis a temperature fixee',
        '            mesurer les neuf observables en production',
        '            sauvegarder mesures, etat final et generateur',
        '        FIN POUR',
        '    FIN POUR',
        'FIN POUR',
        'comparer les branches a temperature identique, cycle par cycle',
        r'\end{verbatim}', r'\end{algocode}',
        r'\subsection{Interprétation des différences entre branches}',
        r"On compare principalement $u$, $q$, les deux paramètres de Cowley et les "
        r"fluctuations. L'ordre signé $\eta$ dépend aussi de la variante et ne suffit "
        r"pas à identifier une hystérésis. Les branches sont appariées au sein d'un "
        r"même cycle : $\Delta A_c(T)=A_{c,\downarrow}(T)-A_{c,\uparrow}(T)$.",
        r"Une différence persistante à budget fini peut révéler une métastabilité ou "
        r"une relaxation insuffisante. À équilibre complet, les observables canoniques "
        r"ne dépendent pas du sens de parcours. Une boucle ne prouve donc pas, à elle "
        r"seule, une transition du premier ordre. La capacité calculée par variance "
        r"est notée $c_V^{\mathrm{fluct}}$ : son interprétation comme réponse thermique "
        r"d'équilibre exige des diagnostics satisfaisants."]
    if not complete:
        lines += [r'\begin{attention}',
            r"Le calcul du bonus est en cours ou en attente de sa fin. Les courbes "
            r"et les valeurs ne sont ajoutées qu'après la sauvegarde et la validation "
            r"des deux branches de tous les cycles. Les résultats des questions 1 à 9 "
            r"proviennent de la campagne canonique antérieure, distincte de ce cycle.",
            r'\end{attention}']
    else:
        result = json.loads((out/'summary.json').read_text(encoding='utf8'))
        lines += [r'\subsection{Résultats du cycle terminé}',
            r'\begin{center}\small\begin{tabular}{lrr}\toprule',
            r'Branche & Maximum de $c_V^{\mathrm{fluct}}$ (K) & Paliers avec diagnostic insuffisant\\\midrule']
        for branch, label in BRANCHES:
            lines.append(f"{label} & {result['peaks_cv_grid_K'][branch]:.0f} & "
                         f"{len(result['diagnostic_failed_temperatures_K'][branch])}"+r'\\')
        lines += [r'\bottomrule\end{tabular}\end{center}',
            f"L'écart refroidissement moins chauffage des maxima vaut "
            f"${result['peak_difference_cooling_minus_heating_K']:.0f}$ K sur la grille. "
            r"Il décrit ce protocole, sans fournir une température critique extrapolée.",
            r"Les barres indiquent une erreur-type calculée séparément par branche. "
            r"Dans les comparaisons des moyennes par branche, les cercles creux représentent "
            r"le chauffage et les carrés creux le refroidissement ; "
            r"les vues rapprochées couvrent 100 K autour du maximum au chauffage. "
            r"Pour les différences, les courbes individuelles montrent la dispersion "
            r"entre cycles ; avec seulement deux graines, elle ne constitue pas un "
            r"test robuste de significativité. Les paliers dont les diagnostics "
            r"échouent doivent être prolongés avant une interprétation d'équilibre."]
        _, archives = load_cycle_archives(out)
        lines += [r'\subsection{Occupations des Au : distinguer les cycles et les variantes}',
            r"Le label $\mu\in\{0,1,2,3\}$ repère une position géométrique dans la "
            r"base CFC ; il ne change pas pendant le calcul. Chaque sous-réseau "
            r"contient $L^3$ sites et son occupation vaut "
            r"$p_\mu=N_{\mathrm{Au},\mu}/L^3$. Ainsi $\sum_\mu p_\mu=1$ "
            r"à la composition étudiée. Un état parfaitement ordonné présente "
            r"une occupation égale à 1 et trois égales à 0.",
            r"La figure des occupations sépare les cycles indépendants : "
            r"une colonne par cycle et une ligne par sens du parcours. Chaque point "
            r"est la moyenne temporelle des seules mesures de production de ce cycle "
            r"à cette température. Les mêmes couleurs et symboles désignent les "
            r"mêmes sous-réseaux dans tous les panneaux.",
            f"Au chauffage, à ${min(cfg['temperatures']):g}$ K, les archives donnent :",
            r'\begin{center}\small\begin{tabular}{lrrrrr}\toprule',
            r'Cycle & $\langle p_0\rangle$ & $\langle p_1\rangle$ & $\langle p_2\rangle$ & $\langle p_3\rangle$ & $\langle q\rangle$\\\midrule']
        cold_values = []
        for c in range(len(cfg['seeds'])):
            cold = archives['heating'][c][0]['samples']
            values = np.r_[cold[:, 5:9].mean(axis=0), cold[:, 2].mean()]
            cold_values.append(values)
            lines.append(f'{c+1} & '+' & '.join(f'{v:.3f}' for v in values)+r'\\')
        pooled = np.mean(cold_values, axis=0)
        lines += [r'\midrule Moyenne des cycles & '+' & '.join(f'{v:.3f}' for v in pooled)+r'\\',
            r'\bottomrule\end{tabular}\end{center}',
            r"Moyenner des cycles partis dans des variantes différentes peut "
            r"donc répartir l'occupation moyenne entre plusieurs sous-réseaux, "
            r"même si chacun des états est ordonné. Par exemple, la moyenne de "
            r"$(1,0,0,0)$ et $(0,0,1,0)$ vaut $(1/2,0,1/2,0)$ : "
            r"ce vecteur moyen ne décrit aucune des deux configurations. "
            r"Il ne signale pas une erreur de définition de $\mu$.",
            r"Près de la transition, un même cycle peut aussi changer de variante "
            r"au cours du temps : séparer les cycles ne supprime pas cet effet de "
            r"moyenne temporelle. Pour comparer le degré d'ordre, on calcule $q$ "
            r"sur chaque configuration puis on moyenne les valeurs obtenues. "
            r"En général, $\langle q(\mathbf p)\rangle\ne "
            r"q(\langle\mathbf p\rangle)$ ; une occupation moyenne uniforme "
            r"ne suffit donc pas à conclure au désordre instantané."]
        captions = [
            ('bonus_01_thermodynamique', 'Énergie et réponse par fluctuations, sur le cycle complet et dans la zone de transition.'),
            ('bonus_02_ordres', 'Comparaison de l’ordre global, de l’ordre local et des fluctuations.'),
            ('bonus_03_sous_reseaux', 'Occupations des Au séparées par cycle et par sens du parcours. Chaque point est une moyenne temporelle au sein d’un seul cycle ; couleurs et symboles identifient les sous-réseaux, sans moyenne entre cycles.'),
            ('bonus_04_ecarts', 'Différences appariées refroidissement moins chauffage ; les symboles distinguent les cycles indépendants et leur moyenne.'),
            ('bonus_05_traces', 'Historique incluant les transitoires ; seule la représentation est décimée, pas les archives.'),
            ('bonus_06_diagnostics', 'Autocorrélation, accord des chaînes, effectif utile et acceptation des mouvements.'),
            ('bonus_07_configurations', 'Deux coupes de configurations finales à température identique, dans le premier cycle.'),
            ('bonus_08_distributions', 'Distributions d’énergie et d’ordre au voisinage des maxima, comparées à température égale.')]
        for name, caption in captions:
            # Les historiques denses et les coupes à une seule température
            # restent exportés, sans répétition dans le corps du rapport.
            if name in ('bonus_05_traces', 'bonus_07_configurations'):
                continue
            lines.append(r'\fig{'+name+r'}{.95}{'+caption+'}')
        lines.append(r"Les exports complémentaires \texttt{bonus\_05\_traces.pdf} "
                     r"(historiques denses, transitoires compris) et "
                     r"\texttt{bonus\_07\_configurations.pdf} (coupes à 1050 K) "
                     r"restent dans \texttt{results/thermal\_cycle\_L6/figures}. "
                     r"Les neuf coupes fixes et l'évolution sont examinées dans les GIF, "
                     r"sans répéter ces vues dans le corps du rapport.")
        (target/'figures').mkdir(exist_ok=True)
        for path in (out/'figures').glob('bonus_*.*'):
            if path.suffix in ('.png', '.pdf'):
                destination = target/'figures'/path.name
                temporary = destination.with_suffix(destination.suffix+'.tmp')
                shutil.copy2(path, temporary)
                temporary.replace(destination)
    (out/'report').mkdir(exist_ok=True)
    text = ('\n'.join(lines)+'\n').replace(r'\subsection{', r'\subsubsection{')
    (out/'report'/'bonus_cycle.tex').write_text(text, encoding='utf8')
    (target/'bonus_cycle.tex').write_text(text, encoding='utf8')


def publish_bonus(output, reference_dir, complete=False, compile_pdf=True):
    from .reporting import make_report
    reference = Path(reference_dir).resolve()
    write_bonus_include(output, reference, complete=complete)
    make_report(reference, compile_pdf=False)
    if not compile_pdf:
        return
    root = Path(__file__).resolve().parents[1]
    executable = shutil.which('pdflatex')
    if executable is None:
        raise RuntimeError('pdflatex is unavailable; report sources have been saved.')
    build = Path(output).resolve()/'report'/'build'
    build.mkdir(parents=True, exist_ok=True)
    for iteration in range(3):
        command = [executable, '-interaction=nonstopmode', '-halt-on-error',
                   '-output-directory='+str(build),
                   str(root/'report'/'rapport_template.tex')]
        run = subprocess.run(command, cwd=reference/'report', capture_output=True,
                             text=True, encoding='utf8', errors='replace')
        (Path(output)/'report'/f'compile_pass_{iteration+1}.log').write_text(
            run.stdout+run.stderr, encoding='utf8')
        if run.returncode:
            raise RuntimeError('LaTeX compilation failed: '+run.stdout[-2000:])
        log = (build/'rapport_template.log').read_text(encoding='utf8', errors='replace')
        if iteration >= 1 and 'Rerun' not in log and 'undefined references' not in log:
            break
    # Preserve the last good PDFs if compilation fails or the job is stopped.
    for destination in (root/'report'/'rapport_template.pdf', reference/'report'/'rapport.pdf'):
        temporary = destination.with_suffix('.pdf.tmp')
        shutil.copy2(build/'rapport_template.pdf', temporary)
        temporary.replace(destination)
    shutil.copy2(build/'rapport_template.log', root/'report'/'rapport_template.log')
