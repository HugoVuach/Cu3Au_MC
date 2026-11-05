"""Deux plaquettes GIF : cycle thermique archivé et relaxation isotherme dédiée.

Lancement : python -m cu3au.slice_animations
"""
import json
from pathlib import Path

import numpy as np
from PIL import Image

from .crystal_slices import ORIENTATIONS, plane_polygon
from .figures import plt, style, BLUE, GOLD, INK
from .lattice import build_lattice
from .model import energy_int, ordre, SCALE
from .kernel import run_segment


def thermal_history(root):
    """Une trajectoire continue, sans mélanger les deux cycles indépendants."""
    cfg = json.loads((root/'config.json').read_text(encoding='utf8'))
    states, temperatures, branches, sweeps = [], [], [], []
    previous = None
    sources = sorted((root/'chain_01').glob('step_*.npz'))
    for index, path in enumerate(sources):
        with np.load(path) as a:
            if index == 0:
                states.append(a['initial_state'].copy())
                temperatures.append(float(a['temperature_K']))
                branches.append('Départ ordonné')
                sweeps.append(0)
            else:
                np.testing.assert_array_equal(previous, a['initial_state'])
            previous = a['final_state'].copy()
            states.append(previous)
            temperatures.append(float(a['temperature_K']))
            branches.append('Chauffage' if str(a['branch']) == 'heating' else 'Refroidissement')
            sweeps.append((index+1)*(cfg['burn_sweeps']+cfg['sample_sweeps']))
    if len(sources) != 2*len(cfg['temperatures']):
        raise ValueError('Incomplete thermal cycle.')
    return dict(states=np.array(states), temperatures=np.array(temperatures),
                branches=np.array(branches), sweeps=np.array(sweeps), L=cfg['L'])


def isothermal_history(root, output, temperature=1000., total=6000, every=20, seed=20261104):
    """Trempe depuis un état chaud archivé ; aucune équilibration n'est retranchée."""
    target = output/'relaxation_isotherme.npz'
    cfg = json.loads((root/'config.json').read_text(encoding='utf8'))
    source = root/'chain_01'/f"step_{len(cfg['temperatures'])-1:03d}.npz"
    lat = build_lattice(cfg['L'])
    with np.load(source) as a:
        initial = a['final_state'].copy()
        initial_temperature = float(a['temperature_K'])
    if target.exists():
        with np.load(target) as a:
            if (float(a['temperature_K']) != temperature or int(a['total_sweeps']) != total
                    or int(a['every']) != every or int(a['seed']) != seed):
                raise ValueError('Existing isothermal data use another configuration.')
            np.testing.assert_array_equal(a['states'][0], initial)
            return {key: a[key].copy() for key in ['states', 'temperatures', 'branches', 'sweeps', 'L']}
    if total % every:
        raise ValueError('Capture interval must divide total sweeps.')
    states = initial[None, :].copy()
    energies = np.array([energy_int(initial, lat.neighbors)], dtype=np.int64)
    saved, observations = [initial.copy()], []
    rng = np.random.default_rng(seed)
    for done in range(every, total+1, every):
        records, _, _, _, _, _ = run_segment(
            states, energies, np.array([temperature]), lat.neighbors, lat.sublattice,
            every, 0, every, int(rng.integers(1, 2**31-1)), np.zeros(1, dtype=np.int32))
        if (int(states[0].sum()) != lat.n//4 or
                energy_int(states[0], lat.neighbors) != energies[0] or
                energies[0] < round(-3.813*lat.n*SCALE)):
            raise RuntimeError('Composition/energy validation failed.')
        saved.append(states[0].copy())
        observations.append(records[0, 0])
        if done % 1000 == 0:
            print(f'Isothermal: {done}/{total} sweeps; q={records[0, 0, 2]:.3f}', flush=True)
    sweeps = np.arange(0, total+1, every)
    history = dict(states=np.array(saved), temperatures=np.full(len(saved), temperature),
                   branches=np.array(['Trempe puis relaxation']*len(saved)), sweeps=sweeps, L=lat.L)
    np.savez_compressed(target, **history, observations=np.array(observations),
                        xyz=lat.xyz, sublattice=lat.sublattice, temperature_K=temperature,
                        initial_temperature_K=initial_temperature, total_sweeps=total,
                        every=every, seed=seed, source=str(source),
                        rng_state_json=json.dumps(rng.bit_generator.state))
    return history


def render(history, folder, name, title, footnote, duration):
    """Géométrie/cadrage constants ; seules espèces et annotations sont actualisées."""
    from matplotlib.patches import Polygon
    from matplotlib.lines import Line2D
    lat = build_lattice(int(history['L']))
    style()
    fig, axes = plt.subplots(3, 3, figsize=(11, 11), dpi=110)
    artists = []
    metadata = []
    for row, (orientation, normal, horizontal, vertical, equation, xlab, ylab) in enumerate(ORIENTATIONS):
        levels = np.unique(lat.xyz@normal)
        middle = np.argmin(abs(levels-(lat.L-.5)*normal.sum()))
        polygons = []
        for col, level in enumerate(levels[middle-1:middle+2]):
            ax = axes[row, col]
            mask = lat.xyz@normal == level
            points = lat.xyz[mask]/2
            uv = np.column_stack((points@horizontal, points@vertical))
            polygon = plane_polygon(normal, level/2, lat.L, horizontal, vertical)
            polygons.append(polygon)
            ax.add_patch(Polygon(polygon, facecolor='#f3f6f8', edgecolor='#b1bec8', lw=.7, zorder=0))
            cu = ax.scatter([], [], s=15, c=BLUE, edgecolors='white', linewidths=.3)
            au = ax.scatter([], [], s=45, c=GOLD, edgecolors=INK, linewidths=.5)
            ax.set(title=f'({orientation}) · {equation} = {level/2:g} a', xlabel=xlab, ylabel=ylab)
            ax.set_aspect('equal')
            ax.tick_params(labelsize=8)
            ax.xaxis.label.set_size(8)
            ax.yaxis.label.set_size(8)
            ax.title.set_size(10)
            count = ax.text(.02, .98, '', transform=ax.transAxes, va='top', fontsize=8,
                            bbox=dict(facecolor='white', alpha=.85, edgecolor='none'))
            artists.append((mask, uv, cu, au, count))
            metadata.append(dict(orientation=orientation, integer_level=int(level), sites=int(mask.sum())))
        all_uv = np.concatenate(polygons)
        for ax in axes[row]:
            ax.set_xlim(all_uv[:, 0].min()-.3, all_uv[:, 0].max()+.3)
            ax.set_ylim(all_uv[:, 1].min()-.3, all_uv[:, 1].max()+.3)
    fig.suptitle(title, fontsize=17, y=.985)
    status = fig.text(.5, .94, '', ha='center', fontsize=11)
    fig.legend(handles=[Line2D([], [], marker='o', ls='none', color=GOLD, markeredgecolor=INK, label='Au'),
                        Line2D([], [], marker='o', ls='none', color=BLUE, label='Cu')],
               loc='upper center', bbox_to_anchor=(.5, .925), ncols=2, fontsize=10)
    fig.text(.5, .025, footnote, ha='center', fontsize=8, color=INK)
    fig.tight_layout(rect=(.01, .07, .99, .90), h_pad=1.8, w_pad=1.4)
    frames, observables = [], []
    palette = None
    for frame, state in enumerate(history['states']):
        if state.sum() != lat.n//4:
            raise RuntimeError('Incorrect frame composition.')
        _, q, occupations = ordre(state, lat.sublattice)
        observables.append(dict(sweeps=int(history['sweeps'][frame]),
                                temperature_K=float(history['temperatures'][frame]), q=q,
                                occupations=occupations.tolist()))
        for mask, uv, cu, au, count in artists:
            gold = state[mask].astype(bool)
            cu.set_offsets(uv[~gold])
            au.set_offsets(uv[gold])
            count.set_text(f'{gold.sum()} Au / {len(gold)} sites')
        status.set_text(f"{history['branches'][frame]} · T = {history['temperatures'][frame]:g} K · "
                        f"{history['sweeps'][frame]:,} balayages · q = {q:.3f}".replace(',', ' '))
        fig.canvas.draw()
        rgb = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
        if palette is None:
            palette = rgb.quantize(colors=256)
        frames.append(rgb.quantize(palette=palette, dither=Image.Dither.NONE))
        if frame in (0, len(history['states'])-1):
            suffix = 'debut' if frame == 0 else 'fin'
            rgb.save(folder/f'{name}_{suffix}.png')
            fig.savefig(folder/f'{name}_{suffix}.pdf')
        if frame % 50 == 0:
            print(f'{name}: image {frame+1}/{len(history["states"])}', flush=True)
    durations = [duration]*len(frames)
    if name == 'plaquette_relaxation_isotherme':
        # Ralentir les 200 premiers balayages, où se produit l'essentiel de la
        # transformation ; les valeurs de balayage restent affichées exactement.
        durations = [800 if sweep <= 200 else duration for sweep in history['sweeps']]
    durations[0] = 1000
    durations[-1] = 1500
    temporary = folder/f'{name}.gif.tmp'
    frames[0].save(temporary, format='GIF', save_all=True, append_images=frames[1:],
                   loop=0, duration=durations, disposal=2, optimize=False)
    temporary.replace(folder/f'{name}.gif')
    plt.close(fig)
    with Image.open(folder/f'{name}.gif') as gif:
        if gif.n_frames != len(frames):
            raise RuntimeError('GIF frame count mismatch.')
        dimensions = gif.size
    (folder/f'{name}.json').write_text(json.dumps(dict(
        frames=len(frames), dimensions=dimensions, L=lat.L, slices=metadata,
        title=title, note=footnote, observables=observables), indent=2, ensure_ascii=False), encoding='utf8')


def main():
    root = Path('results/thermal_cycle_L6')
    output = root/'visualisation'/'animations'
    output.mkdir(parents=True, exist_ok=True)
    fixed = isothermal_history(root, output)
    render(thermal_history(root), output, 'plaquette_cycle_thermique',
           'Cycle thermique · L = 6 · cycle 1 · neuf coupes fixes',
           'Une image par fin de palier ; les événements entre paliers ne sont pas enregistrés.\n'
           'Balayages cumulés, équilibration comprise ; lecture accélérée, sans temps physique en secondes.', 350)
    render(fixed, output, 'plaquette_relaxation_isotherme',
           'Relaxation après trempe 4000 → 1000 K · L = 6 · neuf coupes fixes',
           'Configurations espacées de 20 balayages à 1000 K ; première image : état chaud avant relaxation.\n'
           'Début ralenti ; échanges Au–Cu globaux ; l’apparition d’ordre ne démontre pas à elle seule une nucléation.', 100)
    print('Both GIFs complete.', flush=True)


if __name__ == '__main__':
    main()
