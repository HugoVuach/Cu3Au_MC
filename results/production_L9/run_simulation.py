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
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')


def status(state, **details):
    data = {'state': state, 'pid': os.getpid(), 'config': str(CONFIG), **details}
    temporary = OUTPUT / 'execution.tmp'
    temporary.write_text(json.dumps(data, indent=2), encoding='utf-8')
    temporary.replace(OUTPUT / 'execution.json')


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
    summary = campaign(cfg, OUTPUT)
    status('validating', L=cfg['L'], n_atoms=summary['n_atoms'])
    validate_results(OUTPUT)
    status('plotting', L=cfg['L'], n_atoms=summary['n_atoms'])
    figures = make_figures(OUTPUT)
    status('complete', L=cfg['L'], n_atoms=summary['n_atoms'],
           T_peak_K=summary['T_cv_peak_grid_K'],
           critical_diagnostics_ok=summary['critical_diagnostics_ok'],
           figures=len(figures))
    print('Campaign, validation, and figures complete.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        status('failed', error=traceback.format_exc())
        traceback.print_exc()
        raise SystemExit(1)
