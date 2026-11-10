"""L12 : tests, campagne dédiée, validation, figures et comparaisons L6/L9."""
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent
CONFIG = ROOT/'configs/production_L12.json'
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ['PYTHONUTF8'] = '1'
os.environ.setdefault('MPLCONFIGDIR', str(OUTPUT/'plot_cache'))
sys.stdout.reconfigure(encoding='utf8')
sys.stderr.reconfigure(encoding='utf8')


def status(state, **details):
    from cu3au.thermal_cycle import write_status
    write_status(OUTPUT, state, config=str(CONFIG), **details)


def main():
    if (OUTPUT/'summary.json').exists():
        raise RuntimeError('Campaign already complete; refusing to restart.')
    cfg = json.loads(CONFIG.read_text(encoding='utf8'))
    if (OUTPUT/'config.json').exists():
        if json.loads((OUTPUT/'config.json').read_text(encoding='utf8')) != cfg:
            raise RuntimeError('Frozen run configuration differs; refusing to mix runs.')
    else:
        (OUTPUT/'config.json').write_text(json.dumps(cfg, indent=2)+'\n', encoding='utf8')
    details = dict(L=cfg['L'], n_atoms=4*cfg['L']**3, temperatures=len(cfg['temperatures']))
    status('testing', **details)
    tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT/'tests'), '-v'],
                           capture_output=True, text=True, encoding='utf8')
    log = tests.stdout+tests.stderr
    (OUTPUT/'tests.log').write_text(log, encoding='utf8')
    print(log, flush=True)
    if tests.returncode:
        raise RuntimeError('Tests failed; simulation not started.')
    from cu3au.simulation import load_config, campaign
    from cu3au.validation import validate_results
    from cu3au.figures import make_figures
    from cu3au.pedagogy import record_current_validation
    from cu3au.size_comparison import compare_sizes
    # La copie locale fige la configuration utilisée par la campagne.
    cfg = load_config(OUTPUT/'config.json')
    status('simulating', **details)
    print(f"Starting L={cfg['L']}: {details['n_atoms']} sites, "
          f"{details['temperatures']} temperatures, {len(cfg['seeds'])} chains; "
          f"burn={cfg['burn_sweeps']}, production={cfg['sample_sweeps']}.", flush=True)

    def progress(chain, phase, done, total, captures):
        checkpoints = sorted((OUTPUT/'checkpoints'/f'chain_{chain:02d}').glob('segment_*.npz'))
        status('simulating', **details, chain=chain, chains=len(cfg['seeds']),
               phase=phase, sweeps_done=done, sweeps_target=total, spatial_captures=captures,
               last_checkpoint=str(checkpoints[-1]) if checkpoints else None,
               saved_checkpoints=len(checkpoints))

    summary = campaign(cfg, OUTPUT, progress=progress, checkpoint=True, resume=True)
    status('validating', **details)
    validate_results(OUTPUT)
    record_current_validation(OUTPUT, log)
    status('plotting', **details)
    figures = make_figures(OUTPUT)
    status('comparing_sizes', **details)
    comparisons = [str(compare_sizes(ROOT/'results'/name, OUTPUT))
                   for name in ['production', 'production_L9']]
    status('complete', **details, T_peak_K=summary['T_cv_peak_grid_K'],
           critical_diagnostics_ok=summary['critical_diagnostics_ok'],
           figures=len(figures), size_comparisons=comparisons)
    print('L12 complete: archives validated, figures and both size comparisons generated.', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        status('failed', error=traceback.format_exc())
        traceback.print_exc()
        raise SystemExit(1)
