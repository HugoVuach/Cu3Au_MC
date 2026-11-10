"""Calcule L9 par défaut ; publie le PDF après chaque taille finie.

Relancer ce même fichier reprend après le dernier palier sauvegardé.
Les autres tailles demandent une option explicite, par exemple --sizes 12.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
os.chdir(ROOT)
os.environ['PYTHONUTF8']='1'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'results/production_L12/plot_cache'))


def status(state,**details):
    from datetime import datetime
    out=ROOT/'results/comparison_L6_L9_L12/cycles_execution.json'
    data=dict(state=state,pid=os.getpid(),updated_at=datetime.now().astimezone().isoformat(),**details)
    temporary=out.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(data,indent=2),encoding='utf8');temporary.replace(out)


def main(sizes):
    from cu3au.thermal_cycle import load_cycle_config,run_cycle,write_status
    from cu3au.thermal_cycle_reporting import make_cycle_figures,publish_bonus
    from cu3au.cycle_size_report import completed_cycles
    completed=[r['L'] for r in completed_cycles()]
    plan=ROOT/'results/comparison_L6_L9_L12/cycles_plan.json'
    plan.write_text(json.dumps(dict(requested_sizes=sizes,
        deferred_sizes=[L for L in [9,12] if L not in sizes and L not in completed]),indent=2)+'\n',encoding='utf8')
    for L in sizes:
        cfg=load_cycle_config(ROOT/f'configs/thermal_cycle_L{L}.json')
        out=ROOT/f'results/thermal_cycle_L{L}'
        status('simulating',L=L,planned_sizes=sizes,completed_sizes=completed)
        print(f'Starting/resuming thermal cycle L={L}',flush=True)
        # run_cycle contrôle tous les paliers, y compris ceux déjà archivés.
        def progress(chain,step,path):
            status('simulating',L=L,planned_sizes=sizes,completed_sizes=completed,
                chain=chain+1,step=step+1,completed_plateaus=(chain*2*len(cfg['temperatures'])+step+1),
                total_plateaus=4*len(cfg['temperatures']))
        result=run_cycle(cfg,out,resume=True,progress=progress)
        status('plotting',L=L)
        write_status(out,'plotting',L=L)
        make_cycle_figures(out)
        write_status(out,'complete',L=L,n_atoms=result['n_atoms'])
        status('publishing',L=L)
        publish_bonus(ROOT/'results/thermal_cycle_L6',ROOT/'results/production',complete=True)
        print(f'Cycle L={L} completed and primary PDF updated.',flush=True)
        if L not in completed:
            completed.append(L)
    status('complete',planned_sizes=sizes,completed_sizes=completed)


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf8');sys.stderr.reconfigure(encoding='utf8')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes',nargs='+',type=int,choices=[9,12],default=[9])
    args=parser.parse_args()
    if len(set(args.sizes))!=len(args.sizes):
        parser.error('Each size must appear only once.')
    try:
        main(args.sizes)
    except Exception:
        status('failed',error=traceback.format_exc())
        traceback.print_exc();raise SystemExit(1)
