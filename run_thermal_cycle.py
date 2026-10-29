#!/usr/bin/env python3
"""Run or resume a sequential heating/cooling cycle from a saved config."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=root/'configs/thermal_cycle_L6.json')
    parser.add_argument('--output', type=Path, default=root/'results/thermal_cycle_L6')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--reference-report', type=Path,
                        help='Attach the bonus to this existing equilibrium result directory.')
    parser.add_argument('--no-pdf', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    os.environ['PYTHONUTF8'] = '1'
    os.environ.setdefault('MPLCONFIGDIR', str(out/'plot_cache'))
    sys.stdout.reconfigure(encoding='utf8')
    sys.stderr.reconfigure(encoding='utf8')
    from cu3au.thermal_cycle import (load_cycle_config, config_fingerprint,
                                     run_cycle, write_status)
    try:
        cfg = load_cycle_config(args.config)
        if (out/'config.json').exists():
            prior = json.loads((out/'config.json').read_text(encoding='utf8'))
            if config_fingerprint(prior) != config_fingerprint(cfg):
                raise ValueError('Different configuration in output directory.')
        if any(out.glob('chain_*/step_*.npz')) and not args.resume:
            raise FileExistsError('Existing cycle: use --resume or a new output directory.')
        (out/'config.json').write_text(json.dumps(cfg, indent=2), encoding='utf8')
        write_status(out, 'testing', L=cfg['L'])
        tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover',
                                '-s', str(root/'tests'), '-v'], cwd=root,
                                capture_output=True, text=True, encoding='utf8')
        log = tests.stdout+tests.stderr
        (out/'tests.log').write_text(log, encoding='utf8')
        print(log, flush=True)
        if tests.returncode:
            raise RuntimeError('Tests failed; no production cycle was started.')
        from cu3au.thermal_cycle_reporting import make_cycle_figures, publish_bonus
        if args.reference_report:
            from cu3au.pedagogy import record_current_validation
            record_current_validation(args.reference_report, log)
            write_status(out, 'preparing_report', L=cfg['L'])
            publish_bonus(out, args.reference_report, complete=False,
                          compile_pdf=not args.no_pdf)
        print(f'Starting sequential cycle L={cfg["L"]}: '
              f'{len(cfg["temperatures"])} plateaus per direction, '
              f'{len(cfg["seeds"])} independent cycles.', flush=True)
        result = run_cycle(cfg, out, resume=args.resume)
        write_status(out, 'plotting', L=cfg['L'], n_atoms=result['n_atoms'])
        figures = make_cycle_figures(out)
        if args.reference_report:
            write_status(out, 'publishing_report', L=cfg['L'])
            publish_bonus(out, args.reference_report, complete=True,
                          compile_pdf=not args.no_pdf)
        write_status(out, 'complete', L=cfg['L'], n_atoms=result['n_atoms'],
                     figures=len(figures), peaks_cv_grid_K=result['peaks_cv_grid_K'])
        print('Cycle, validation, figures, and requested report update complete.', flush=True)
    except Exception:
        write_status(out, 'failed', error=traceback.format_exc())
        traceback.print_exc()
        raise SystemExit(1)


if __name__ == '__main__':
    main()
