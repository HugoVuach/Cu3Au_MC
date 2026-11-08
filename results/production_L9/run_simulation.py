"""Run the saved L=9 campaign, validate archives, and generate figures."""
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent
CONFIG = ROOT / 'configs' / 'production_L9.json'
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ['PYTHONUTF8'] = '1'
os.environ.setdefault('MPLCONFIGDIR', str(OUTPUT/'plot_cache'))
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')


def status(state, **details):
    from cu3au.thermal_cycle import write_status
    write_status(OUTPUT, state, config=str(CONFIG), **details)


def main():
    if (OUTPUT / 'summary.json').exists() or any(OUTPUT.glob('chain_*.npz')):
        raise RuntimeError('Existing campaign archives: refusing to overwrite.')
    status('testing', L=9, n_atoms=2916)
    tests = subprocess.run(
        [sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT/'tests'), '-v'],
        capture_output=True, text=True, encoding='utf-8', check=False,
    )
    (OUTPUT/'tests.log').write_text(tests.stdout+tests.stderr, encoding='utf-8')
    print(tests.stdout+tests.stderr, flush=True)
    if tests.returncode:
        raise RuntimeError('Tests failed; simulation was not started.')

    from cu3au.simulation import load_config, campaign
    from cu3au.validation import validate_results
    from cu3au.figures import make_figures

    cfg = load_config(CONFIG)
    status('simulating', L=cfg['L'], n_atoms=4*cfg['L']**3)
    print('Starting L=9 campaign: 2916 sites, 46 temperatures, 2 chains.', flush=True)
    def progress(chain, phase, done, total, captures):
        status('simulating', L=cfg['L'], n_atoms=4*cfg['L']**3,
               chain=chain, chains=len(cfg['seeds']), phase=phase,
               sweeps_done=done, sweeps_target=total, spatial_captures=captures)
    summary = campaign(cfg, OUTPUT, progress=progress)
    status('validating', L=cfg['L'], n_atoms=summary['n_atoms'])
    validate_results(OUTPUT)
    from cu3au.pedagogy import record_current_validation
    record_current_validation(OUTPUT, tests.stdout+tests.stderr)
    status('plotting', L=cfg['L'], n_atoms=summary['n_atoms'])
    figures = make_figures(OUTPUT)
    status('comparing_sizes', L=cfg['L'], n_atoms=summary['n_atoms'])
    from cu3au.size_comparison import compare_sizes
    comparison = compare_sizes(ROOT/'results/production', OUTPUT)
    status('complete', L=cfg['L'], n_atoms=summary['n_atoms'],
           T_peak_K=summary['T_cv_peak_grid_K'],
           critical_diagnostics_ok=summary['critical_diagnostics_ok'],
           figures=len(figures), size_comparison=str(comparison))
    print('Campaign, validation, figures, and L6/L9 comparison complete.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        status('failed', error=traceback.format_exc())
        traceback.print_exc()
        raise SystemExit(1)
