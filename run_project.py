#!/usr/bin/env python3
"""Lanceur complet, à exécuter depuis n'importe quel dossier."""
import argparse
import json
import subprocess
import sys
from pathlib import Path


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description='Cu3Au : tests, Monte-Carlo, figures et rapport LaTeX.')
    parser.add_argument('--preset', choices=['quick', 'production', 'pilot'], default='production')
    parser.add_argument('--config', type=Path, help='Configuration JSON personnalisée.')
    parser.add_argument('--output', type=Path, help='Dossier de résultats, par défaut results/<preset>.')
    parser.add_argument('--reuse', action='store_true', help='Réutiliser les archives si la configuration est identique.')
    parser.add_argument('--no-pdf', action='store_true', help='Produire les sources LaTeX sans compiler le PDF.')
    parser.add_argument('--tests-only', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve() if args.output else base/'results'/args.preset
    out.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, '-m', 'unittest', 'discover', '-s', str(base/'tests'), '-v']
    print('1/6 · Tests physiques et statistiques', flush=True)
    tests = subprocess.run(cmd, cwd=base, capture_output=True, text=True)
    (out/'tests.log').write_text(tests.stdout+tests.stderr, encoding='utf8')
    print(tests.stdout+tests.stderr, flush=True)
    if tests.returncode:
        raise SystemExit('Tests en échec : le calcul est arrêté.')
    if args.tests_only:
        return
    from cu3au.simulation import load_config, campaign
    from cu3au.validation import validate_results
    from cu3au.figures import make_figures
    from cu3au.reporting import make_report
    cfg = load_config(args.config or base/'configs'/f'{args.preset}.json')
    if args.reuse and (out/'summary.json').exists():
        prior = json.loads((out/'summary.json').read_text(encoding='utf8'))
        if prior['config'] != cfg:
            raise SystemExit('Configuration différente des données existantes : choisir un autre dossier ou retirer --reuse.')
        print('2/6 · Réutilisation explicite des archives existantes', flush=True)
    else:
        if (out/'summary.json').exists():
            raise SystemExit('Ce dossier contient déjà une campagne. Utiliser --reuse ou un nouveau --output.')
        if any(out.glob('chain_*.npz')):
            raise SystemExit('Archives partielles présentes sans résumé final. Choisir un nouveau --output pour éviter de mélanger deux campagnes.')
        print('2/6 · Campagne Monte-Carlo', flush=True)
        campaign(cfg, out)
    print('3/6 · Validation indépendante des archives', flush=True)
    validation = validate_results(out)
    print('4/6 · Figures PDF vectorielles et PNG', flush=True)
    figures = make_figures(out)
    print('5/6 · Rapport question par question', flush=True)
    report = make_report(out, compile_pdf=not args.no_pdf)
    print('6/6 · Résumé', flush=True)
    summary = json.loads((out/'summary.json').read_text(encoding='utf8'))
    print(f"Maximum de cV sur grille : {summary['T_cv_peak_grid_K']:.0f} K")
    print(f"Diagnostics locaux : {summary['critical_diagnostics_ok']}")
    print(f'{len(figures)} figures PDF + PNG. Rapport : {report}')
    if not summary['critical_diagnostics_ok']:
        print('ATTENTION : estimation provisoire ; consulter les diagnostics avant une conclusion quantitative.')
    print('Terminé. Validation technique : '+validation['technical_validation'])


if __name__ == '__main__':
    main()
