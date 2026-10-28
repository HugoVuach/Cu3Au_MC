"""Noyau Numba. Un balayage = N propositions d'échanges Au-Cu.

La réduction E-E0=161000*N_AuAu(1)+6300*N_AuCu(2), en micro-eV,
est exacte à composition fixée. Le calcul de référence utilise les matrices.
"""
import numpy as np
from numba import njit
from .model import KB, SCALE


@njit(cache=True)
def delta_fast(state, i, j, neighbors):
    if state[i] == state[j]:
        return 0
    if state[i] == 0:
        i, j = j, i
    ai, aj, bi, bj = 0, 0, 0, 0
    for k in range(12):
        v, w = neighbors[i, k], neighbors[j, k]
        if v != j:
            ai += state[v]
        if w != i:
            aj += state[w]
    for k in range(12, 18):
        v, w = neighbors[i, k], neighbors[j, k]
        if v != j:
            bi += state[v]
        if w != i:
            bj += state[w]
    return 161000*(aj-ai)+12600*(bi-bj)


@njit(cache=True)
def sample_observables(state, e_int, neighbors, sub):
    n = len(state)
    p = np.zeros(4)
    b1, b2 = 0, 0
    for i in range(n):
        if state[i] == 1:
            p[sub[i]] += 4/n
            for k in range(12):
                b1 += 1-state[neighbors[i, k]]
            for k in range(12, 18):
                b2 += 1-state[neighbors[i, k]]
    eta = (p[0]-.25)/.75
    q = np.sqrt((4/3)*np.sum((p-.25)**2))
    row = np.empty(9)
    row[0] = e_int/SCALE
    row[1] = eta
    row[2] = q
    row[3] = 1-b1/((n/4)*12*.75)
    row[4] = 1-b2/((n/4)*6*.75)
    row[5:] = p
    return row


@njit(cache=True)
def run_segment(states, energies, temps, neighbors, sub, sweeps,
                exchange_every, sample_every, seed, labels, phase=0):

    np.random.seed(seed)
    r, n = states.shape
    records = np.empty((sweeps//sample_every, r, 9))
    label_trace = np.empty((sweeps//sample_every, r), dtype=np.int32)
    accepted = np.zeros(r, dtype=np.int64)
    swap_ok = np.zeros(r-1, dtype=np.int64)
    swap_n = np.zeros(r-1, dtype=np.int64)
    beta = 1/(KB*temps*SCALE)
    m = 0
    for sweep in range(1, sweeps+1):
        for t in range(r):

            for attempt in range(n):
                i = np.random.randint(n)
                j = np.random.randint(n)
                while states[t, i] == states[t, j]:
                    i = np.random.randint(n)
                    j = np.random.randint(n)
                de = delta_fast(states[t], i, j, neighbors)
                
                if de <= 0 or np.random.random() < np.exp(-beta[t]*de):
                    states[t, i], states[t, j] = states[t, j], states[t, i]
                    energies[t] += de
                    accepted[t] += 1

        if exchange_every > 0 and sweep % exchange_every == 0:
            for t in range(phase % 2, r-1, 2):
                swap_n[t] += 1
                loga = (beta[t]-beta[t+1])*(energies[t]-energies[t+1])
                if loga >= 0 or np.log(max(np.random.random(), 1e-300)) < loga:
                    for i in range(n):
                        states[t, i], states[t+1, i] = states[t+1, i], states[t, i]
                    energies[t], energies[t+1] = energies[t+1], energies[t]
                    labels[t], labels[t+1] = labels[t+1], labels[t]
                    swap_ok[t] += 1
            phase += 1

        if sweep % sample_every == 0:
            for t in range(r):
                records[m, t] = sample_observables(states[t], energies[t], neighbors, sub)
                label_trace[m, t] = labels[t]
            m += 1

    return records, label_trace, accepted, swap_ok, swap_n, phase
