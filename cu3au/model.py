"""Deux évaluations indépendantes de l'énergie ; 1=Au, 0=Cu."""
from polars._utils.construction import other
import numpy as np

KB = 8.617333262e-5  # eV/K
SCALE = 1_000_000  # entiers en micro-eV : pas de dérive numérique de l'énergie
J_INT = np.array([[-582000, -689000], [-689000, -635000]], dtype=np.int64)
I_INT = np.array([[0, 6300], [6300, 0]], dtype=np.int64)
E0_PER_ATOM = -3.813

def energy_int(state, neighbors):
    e1 = J_INT[state[:, None], state[neighbors[:, :12]]].sum()
    e2 = I_INT[state[:, None], state[neighbors[:, 12:]]].sum()
    return int(e1+e2)//2


def energy(state, neighbors):
    return energy_int(state, neighbors)/SCALE



def delta_reference(state, i, j, neighbors):
    if i == j or state[i] == state[j]:
        return 0
    a, b = int(state[i]), int(state[j])
    delta = 0
    for site, old, new, other in ((i, a, b, j), (j, b, a, i)):
        for k, v in enumerate(neighbors[site]):
            if v != other:
                table = J_INT if k < 12 else I_INT
                delta += int(table[new, state[v]]-table[old, state[v]])
    return delta


def cowley(state, neighbors):
    mask = state == 1
    out = []
    for sl in (slice(0, 12), slice(12, 18)):
        nn = neighbors[mask, sl]
        p_cu = np.mean(state[nn] == 0)
        out.append(1-p_cu/(1-state.mean()))
    return np.array(out)


def ordre(state, sublattice):
    """eta du sujet"""
    p = np.array([state[sublattice == s].mean() for s in range(4)])
    eta = (p[0]-.25)/.75
    amplitude = np.sqrt((4/3)*np.sum((p-.25)**2))
    return float(eta), float(amplitude), p


def observables(state, lat):
    eta, amplitude, p = ordre(state, lat.sublattice)
    a1, a2 = cowley(state, lat.neighbors)
    return np.r_[energy(state, lat.neighbors), eta, amplitude, a1, a2, p]


def random_energy_limit(n):
    """Espérance EXACTE à T=infini, pour N fini et composition fixée."""
    na, nc = n//4, 3*n//4
    p = np.array([nc*(nc-1), 2*na*nc, na*(na-1)])/float(n*(n-1))
    return (6*np.dot(p, [-.582, -.689, -.635])
            + 3*p[1]*.0063)
