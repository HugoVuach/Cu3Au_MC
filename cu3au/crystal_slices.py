"""Plaquette de neuf couches atomiques exactes dans un instantané sauvegardé.

Exemple : python -m cu3au.crystal_slices
Les plans sont sélectionnés dans le volume central sans projection de couches.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

from .figures import plt, style, save, BLUE, GOLD, INK
from .model import ordre


ORIENTATIONS = [
    ('001', np.array([0, 0, 1]), np.array([1., 0., 0.]), np.array([0., 1., 0.]),
     'z', 'x/a', 'y/a'),
    ('110', np.array([1, 1, 0]), np.array([1., -1., 0.])/np.sqrt(2),
     np.array([0., 0., 1.]), 'x+y', '(x−y)/(√2 a)', 'z/a'),
    ('111', np.array([1, 1, 1]), np.array([1., -1., 0.])/np.sqrt(2),
     np.array([1., 1., -2.])/np.sqrt(6), 'x+y+z', '(x−y)/(√2 a)', '(x+y−2z)/(√6 a)'),
]


def plane_polygon(normal, offset, length, horizontal, vertical):
    """Intersection du plan avec le cube [0,L]^3 ; coordonnées physiques en a."""
    from itertools import product
    vertices = np.array(list(product((0., float(length)), repeat=3)))
    points = []
    for i, left in enumerate(vertices):
        for right in vertices[i+1:]:
            if np.count_nonzero(left != right) != 1:
                continue
            dl, dr = left@normal-offset, right@normal-offset
            if abs(dl) < 1e-10:
                points.append(left)
            if abs(dr) < 1e-10:
                points.append(right)
            if dl*dr < 0:
                points.append(left+dl/(dl-dr)*(right-left))
    points = np.unique(np.round(points, 12), axis=0)
    projection = np.column_stack((points@horizontal, points@vertical))
    return projection[ConvexHull(projection).vertices]


def make_slice_plate(archive, output):
    from matplotlib.patches import Polygon
    from matplotlib.lines import Line2D
    source = Path(archive)
    folder = Path(output)
    folder.mkdir(parents=True, exist_ok=True)
    with np.load(source) as a:
        xyz = a['xyz'].copy()
        state = a['final_state'].copy()
        sub = a['sublattice'].copy()
        temperature = float(a['temperature_K'])
        branch = str(a['branch'])
        chain = int(a['chain'])
    length = int(xyz.max()//2+1)
    eta, q, p = ordre(state, sub)
    style()
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    metadata = []
    for row, (name, normal, horizontal, vertical, equation, xlab, ylab) in enumerate(ORIENTATIONS):
        # Un plan exact vérifie hX+kY+lZ=S dans les coordonnées entières a/2.
        # Pour (111), seules les valeurs paires de S contiennent des sites CFC.
        levels = np.unique(xyz@normal)
        midpoint = (length-.5)*normal.sum()
        middle = np.argmin(abs(levels-midpoint))
        selected = levels[middle-1:middle+2]
        if len(selected) != 3:
            raise ValueError('Three central occupied layers are required.')
        polygons = []
        for col, level in enumerate(selected):
            ax = axes[row, col]
            mask = xyz@normal == level
            points = xyz[mask]/2
            au = state[mask] == 1
            uv = np.column_stack((points@horizontal, points@vertical))
            polygon = plane_polygon(normal, level/2, length, horizontal, vertical)
            polygons.append(polygon)
            ax.add_patch(Polygon(polygon, facecolor='#f3f6f8', edgecolor='#b1bec8', lw=1, zorder=0))
            ax.scatter(*uv[~au].T, s=24, c=BLUE, edgecolors='white', linewidths=.35, zorder=2)
            ax.scatter(*uv[au].T, s=74, c=GOLD, edgecolors=INK, linewidths=.65, zorder=3)
            ax.set(title=f'({name}) · {equation} = {level/2:g} a', xlabel=xlab, ylabel=ylab)
            ax.set_aspect('equal')
            ax.text(.03, .97, f'{au.sum()} Au / {len(au)} sites', transform=ax.transAxes,
                    va='top', fontsize=10, bbox=dict(facecolor='white', alpha=.9, edgecolor='none'))
            metadata.append(dict(orientation=name, integer_level=int(level),
                                 offset_a=float(level/2), sites=int(mask.sum()), au=int(au.sum()),
                                 sublattices=np.unique(sub[mask]).tolist()))
            assert np.allclose(points@normal, level/2)
        # Même échelle et même cadrage pour les trois couches d'une orientation.
        all_uv = np.concatenate(polygons)
        for ax in axes[row]:
            ax.set_xlim(all_uv[:, 0].min()-.3, all_uv[:, 0].max()+.3)
            ax.set_ylim(all_uv[:, 1].min()-.3, all_uv[:, 1].max()+.3)
    label = 'Refroidissement' if branch == 'cooling' else 'Chauffage'
    # Les archives numérotent les chaînes à partir de zéro.
    cycle = chain+1
    fig.suptitle(f'Coupes cristallographiques d’une même configuration\n'
                 f'{label} · cycle {cycle} · {temperature:g} K · L = {length}',
                 fontsize=22, y=.985)
    handles = [Line2D([], [], marker='o', ls='none', color=GOLD, markeredgecolor=INK, ms=9, label='Au'),
               Line2D([], [], marker='o', ls='none', color=BLUE, ms=6, label='Cu')]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .93), ncols=2)
    fig.tight_layout(rect=(.01, .075, .99, .91), h_pad=2.2, w_pad=2)
    name = 'plaquette_coupes_001_110_111'
    save(fig, folder, name,
         'Neuf couches atomiques exactes : trois plans centraux parallèles par orientation, sans épaisseur ni superposition.\n'
         'Instantané final d’un palier : ces coupes spatiales ne constituent pas une séquence temporelle de nucléation.')
    (folder/f'{name}.json').write_text(json.dumps(dict(
        archive=str(source.resolve()), temperature_K=temperature, branch=branch, cycle=cycle,
        L=length, q=q, occupations=p.tolist(), slices=metadata), indent=2), encoding='utf8')
    return folder/f'{name}.pdf'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path,
                        default=Path('results/thermal_cycle_L6/chain_01/step_072.npz'))
    parser.add_argument('--output', type=Path, default=Path('results/thermal_cycle_L6/visualisation'))
    args = parser.parse_args()
    print(make_slice_plate(args.archive, args.output))
