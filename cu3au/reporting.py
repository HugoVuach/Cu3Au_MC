"""Rapport LaTeX régénéré exclusivement depuis les données sauvegardées."""
import csv
import json
import re
import shutil
import subprocess
from pathlib import Path
import numpy as np
from .model import KB


def _copy_atomic(source, destination):
    source = Path(source)
    destination = Path(destination)
    # Éviter de remplacer un export inchangé, notamment lorsqu'un lecteur
    # Windows garde son fichier ouvert pendant la compilation du rapport.
    if destination.exists() and source.read_bytes() == destination.read_bytes():
        return str(destination)
    temporary = destination.with_suffix(destination.suffix+'.tmp')
    shutil.copy2(source, temporary)
    temporary.replace(destination)
    return str(destination)


def make_report(result_dir, compile_pdf=True):
    root = Path(result_dir).resolve()
    out = root/'report'
    out.mkdir(exist_ok=True)
    template = Path(__file__).resolve().parents[1]/'report'/'rapport_template.tex'
    shutil.copy2(template, out/'rapport.tex')
    shutil.copytree(root/'figures', out/'figures', dirs_exist_ok=True,
                    copy_function=_copy_atomic)
    summary = json.loads((root/'summary.json').read_text(encoding='utf8'))
    from .pedagogy import write_pedagogical_tables
    write_pedagogical_tables(root)
    cfg = summary['config']
    with (root/'observables.csv').open(encoding='utf8') as f:
        rows = list(csv.DictReader(f))
    test_path = root/'tests_current.log' if (root/'tests_current.log').exists() else root/'tests.log'
    tests = test_path.read_text(encoding='utf8') if test_path.exists() else ''
    found = re.search(r'Ran (\d+) tests?', tests)
    peaks = []
    for file in sorted(root.glob('chain_*.npz')):
        with np.load(file) as a:
            e = a['samples'][:, :, 0]
            t = a['temperatures']
            cv = np.var(e, axis=0)/(summary['n_atoms']*(KB*t)**2)
            peaks.append(float(t[np.argmax(cv)]))
    status = ('critères locaux satisfaits' if summary['critical_diagnostics_ok']
              else 'au moins un critère local non satisfait ; statistiques à renforcer')
    macros = {
        'Lbox': cfg['L'], 'Natoms': summary['n_atoms'],
        'Tpeak': f'{summary["T_cv_peak_grid_K"]:.0f}',
        'Tleft': f'{summary["T_cv_neighbor_bracket_K"][0]:.0f}',
        'Tright': f'{summary["T_cv_neighbor_bracket_K"][1]:.0f}',
        'Ntemps': len(cfg['temperatures']), 'Nchains': len(cfg['seeds']),
        'BurnSweeps': cfg['burn_sweeps'], 'ProductionSweeps': cfg['sample_sweeps'],
        'SampleEvery': cfg['sample_every'], 'ExchangeEvery': cfg['exchange_every'],
        'Ndraws': cfg['sample_sweeps']//cfg['sample_every'],
        'TestCount': found.group(1) if found else 'non exécutés dans ce dossier',
        'DiagnosticStatus': status,
        'RandomLimit': f'{summary["random_limit_eV_atom_finite_N"]:.9f}',
        'ChainPeaks': ', '.join(f'chaîne {i+1} : ${t:.0f}$ K' for i, t in enumerate(peaks))}
    (out/'auto_results.tex').write_text('\n'.join(
        '\\newcommand{\\'+k+'}{'+str(v)+'}' for k, v in macros.items())+'\n', encoding='utf8')
    table = [r'\begin{center}\small',
             r'\begin{tabular}{rrrrrr}\toprule',
             r'$T$ (K) & $q$ & $c_V/k_B$ & ESS & $\widehat R_E$ & $\widehat R_q$\\\midrule']
    for i in summary['critical_window_rows']:
        row = rows[i]
        table.append(f"{float(row['T_K']):.0f} & {float(row['q']):.3f} & "
                     f"{float(row['cv_kB']):.2f} & {float(row['ess_min']):.0f} & "
                     f"{float(row['rhat_E']):.3f} & {float(row['rhat_q']):.3f}"+r'\\')
    table += [r'\bottomrule\end{tabular}\end{center}']
    (out/'critical_table.tex').write_text('\n'.join(table), encoding='utf8')
    low, high = rows[0], rows[-1]
    peak = max(rows, key=lambda r: float(r['cv_kB']))
    paragraph = (
        f"À ${float(low['T_K']):.0f}$ K, le calcul donne "
        f"$u={float(low['u_eV']):.6f}$ eV/atome et $q={float(low['q']):.4f}$. "
        "La comparaison à la borne exacte vérifie l'état ordonné de basse température.\n\n"
        f"Au maximum thermique, à ${float(peak['T_K']):.0f}$ K, "
        f"$q={float(peak['q']):.3f}$ et $c_V/k_B={float(peak['cv_kB']):.2f}$ par atome. "
        "Les historiques et les distributions permettent d'évaluer l'exploration "
        "des configurations dans cette zone.\n\n"
        f"À ${float(high['T_K']):.0f}$ K, on obtient "
        f"$u={float(high['u_eV']):.6f}$ eV/atome, $q={float(high['q']):.3f}$, "
        f"$\\alpha_1={float(high['alpha1']):.3f}$ et $\\alpha_2={float(high['alpha2']):.3f}$. "
        "Ces valeurs sont celles d'une température finie ; elles ne doivent pas "
        "être remplacées par les limites de mélange indépendant.\n\n"
    )
    if summary['critical_diagnostics_ok']:
        paragraph += ("Les critères locaux fixés avant l'analyse sont satisfaits dans la "
                      "fenêtre du maximum. Cela appuie la cohérence statistique du résultat "
                      "à cette taille, tout en laissant ouvertes les erreurs systématiques "
                      "liées à des états métastables non visités.\n")
    else:
        paragraph += ("Au moins un critère de contrôle échoue dans la fenêtre du maximum. "
                      "Le maximum doit donc être présenté comme provisoire : augmenter les "
                      "durées, densifier la grille ou améliorer la circulation des répliques "
                      "avant d'en tirer une conclusion quantitative ferme.\n")
    (out/'discussion_results.tex').write_text(paragraph, encoding='utf8')
    # La synthèse principale reste L6 ; les résultats L9 sont un complément
    # nommé, pour ne pas remplacer silencieusement les macros et figures L6.
    larger = root.parent/'production_L9'
    if root.name == 'production' and (larger/'comparison_L6_L9/summary.json').exists():
        from .size_comparison import write_report_include
        write_report_include(root, larger)
    if compile_pdf:
        executable = shutil.which('pdflatex')
        if not executable:
            raise RuntimeError('pdflatex absent. Installer TeX Live/MiKTeX ou utiliser --no-pdf ; les sources LaTeX sont déjà produites.')
        for iteration in range(2):
            result = subprocess.run([executable, '-interaction=nonstopmode', '-halt-on-error',
                                     'rapport.tex'], cwd=out, capture_output=True, text=True)
            (out/f'compile_pass_{iteration+1}.log').write_text(result.stdout+result.stderr, encoding='utf8')
            if result.returncode:
                raise RuntimeError('Échec de compilation LaTeX :\n'+result.stdout[-4000:])
    return out
