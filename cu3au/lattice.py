"""Géométrie exacte en unités a/2 ; aucune comparaison de flottants."""
from dataclasses import dataclass
from itertools import product
import numpy as np


@dataclass(frozen=True)
class Lattice:
    L: int
    xyz: np.ndarray
    sublattice: np.ndarray
    neighbors: np.ndarray

    @property
    def n(self):
        return len(self.xyz)


def build_lattice(L=6):
    """Les petites boîtes L<3 sont exclues pour éviter les images multiples."""
    if isinstance(L, bool) or int(L) != L or L < 3:
        raise ValueError("L doit être un entier >= 3.")
    L = int(L)
    basis = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]])
    xyz = np.array([2*np.array(c)+b for c in product(range(L), repeat=3) for b in basis], dtype=np.int32) # On répète une base dans toutes les mailles ( maille CFC conventionnelle en unité de a) : X = 2(i,j,k) + b_(mu)
    sub = np.tile(np.arange(4, dtype=np.int8), L**3)  # Attribuer un sous réseau à chaque site : sub[i] = mu , les sous réseaux coorespondent aux parité des coordonnées
    index = {tuple(x): i for i, x in enumerate(xyz)}  # Construction d'une correspondance entre coordonnées et indices : (X,Y,Z) -> i
    first = [d for d in product((-1, 0, 1), repeat=3) if sum(x*x for x in d) == 2]   # Construire les douze déplacements aux premiers voisins, on examine les 3^3 triplets et on conservet que ceux dont la somme de carré des coordonnées = 2. Ca reprensente avec les signes à 12.
    second = [tuple(2*s if k == axis else 0 for k in range(3)) for axis in range(3) for s in (-1, 1)]  # Construire les 6 déplacements aux second voisins
    shifts = np.array(first+second, dtype=np.int32)  # On regroupe les déplacements, on concatenne la liste des premiers et second voisins
    neighbors = np.array([[index[tuple((x+d) % (2*L))] for d in shifts] for x in xyz], dtype=np.int32) # On construit les indices des voisins par périodicité. 1; on calcul les coordonnée des voisins : X= x+d  2. Ramener le voisin dans sa boite périodique 3. Convertir les coordonnée en indices
    return Lattice(L, xyz, sub, neighbors)   # on renvoie la géométrie finale


def ordered(lat, variant=0):   # on ordonne 0 si  Cu et 1 si Au
    if variant not in range(4):
        raise ValueError("La variante doit appartenir à {0,1,2,3}.")
    return (lat.sublattice == variant).astype(np.int8)


def random_configuration(lat, rng):
    state = np.zeros(lat.n, dtype=np.int8)
    state[rng.choice(lat.n, lat.n//4, replace=False)] = 1
    return state
