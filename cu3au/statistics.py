"""Erreurs tenant compte de l'autocorrélation et diagnostics entre chaînes."""
import numpy as np
from scipy.stats import rankdata, norm
from .model import KB


def tau_int(x):
    """1/2 + somme des autocorrélations, séquence de paires positives."""
    x = np.asarray(x, dtype=float)
    x = x-x.mean()
    n = len(x)
    if n < 8 or np.dot(x, x) < 1e-24:
        return .5
    f = np.fft.rfft(x, n=2**int(np.ceil(np.log2(2*n))))
    ac = np.fft.irfft(f*f.conjugate())[:n]
    ac /= np.arange(n, 0, -1)
    ac /= ac[0]
    total = 0.
    for k in range(1, n//2-1, 2):
        pair = ac[k]+ac[k+1]
        if pair <= 0:
            break
        total += pair
        if k > 10*(.5+total):
            break
    return max(.5, .5+total)


def split_rhat(chains):
    """R-hat séparé, normalisé par rangs, maximum avec la version repliée."""
    a = np.asarray(chains, dtype=float)
    m, n = a.shape
    half = n//2
    if half < 4:
        return float('nan')
    a = np.concatenate([a[:, :half], a[:, -half:]], axis=0)

    def basic(x):
        ranks = rankdata(x.ravel()).reshape(x.shape)
        z = norm.ppf((ranks-.375)/(x.size+.25))
        w = np.mean(np.var(z, axis=1, ddof=1))
        if w < 1e-24:
            return 1. if np.ptp(z) < 1e-12 else float('inf')
        b = half*np.var(z.mean(axis=1), ddof=1)
        return max(1., np.sqrt(((half-1)*w/half+b/half)/w))
    return max(basic(a), basic(np.abs(a-np.median(a))))


def moments(a, temperature, n):
    """a=(mesures,9). e est totale ; cv/kB est sans dimension."""
    e, eta, q, a1, a2 = a[:, :5].T
    q2 = np.mean(q*q)
    binder = 1-np.mean(q**4)/(3*q2*q2) if q2 > 0 else np.nan
    return np.array([e.mean()/n, eta.mean(), q.mean(), a1.mean(), a2.mean(),
                     np.var(e)/(n*(KB*temperature)**2),
                     n*np.var(q)/(KB*temperature), binder])


NAMES = ['u_eV', 'eta', 'q', 'alpha1', 'alpha2', 'cv_kB', 'Xq_per_eV', 'binder_q']


def bootstrap_peak(energies, temperatures, n, block, bootstrap=1000, seed=9421):
    """Circular blocks, synchronized across temperatures within each chain.

    Return grid-valued cV peak locations. Chains keep their separate thermal
    variances, and each chain receives an independent set of block starts.
    """
    e = np.asarray(energies, dtype=float)
    t = np.asarray(temperatures, dtype=float)
    if (e.ndim != 3 or e.shape[2] != len(t) or n <= 0 or
            not np.all(np.isfinite(e)) or np.any(t <= 0) or
            not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0)):
        raise ValueError('Invalid energies, temperature grid, or system size.')
    chains, draws, nt = e.shape
    if not 1 <= block <= draws or bootstrap < 2:
        raise ValueError('Require 1 <= block <= draws and bootstrap >= 2.')
    rng = np.random.default_rng(seed)
    peaks = np.empty(bootstrap)
    # Center before resampling: absolute binding energy is irrelevant to cV.
    e = e-e.mean(axis=1, keepdims=True)
    for b in range(bootstrap):
        cv = np.zeros(nt)
        for c in range(chains):
            starts = rng.integers(draws, size=int(np.ceil(draws/block)))
            ix = ((starts[:, None]+np.arange(block)) % draws).ravel()[:draws]
            cv += np.var(e[c, ix, :], axis=0)/(n*(KB*t)**2)
        peaks[b] = t[np.argmax(cv/chains)]
    return peaks


def summarize(samples, temperatures, n, seed=318, bootstrap=200):
    """samples=(chaînes, mesures, températures,9).

Les variances thermiques sont calculées séparément dans chaque chaîne : on
ne confond pas un désaccord entre chaînes avec une fluctuation d'équilibre.
Bootstrap par blocs circulaires de longueur >= 2*tau, par chaîne.
"""
    rng = np.random.default_rng(seed)
    chains, draws, nt, _ = samples.shape
    rows = []
    for t, temperature in enumerate(temperatures):
        a = samples[:, :, t, :]
        taus = np.array([[tau_int(c[:, k]) for k in (0, 2)] for c in a])
        tau = float(taus.max())
        block = min(max(2, int(np.ceil(2*tau))), max(2, draws//4))
        est = np.mean([moments(c, temperature, n) for c in a], axis=0)
        boot = np.empty((bootstrap, len(NAMES)))
        for b in range(bootstrap):
            per_chain = []
            for c in a:
                starts = rng.integers(0, draws, size=int(np.ceil(draws/block)))
                ix = ((starts[:, None]+np.arange(block)) % draws).ravel()[:draws]
                per_chain.append(moments(c[ix], temperature, n))
            boot[b] = np.mean(per_chain, axis=0)
        errors = np.std(boot, axis=0, ddof=1)
        # L'écart entre chaînes indépendantes doit aussi apparaître dans les erreurs.
        if chains > 1:
            chain_est = np.array([moments(c, temperature, n) for c in a])
            errors = np.maximum(errors, np.std(chain_est, axis=0, ddof=1)/np.sqrt(chains))
        row = {'T_K': float(temperature)}
        row.update({key: float(v) for key, v in zip(NAMES, est)})
        row.update({key+'_se': float(v) for key, v in zip(NAMES, errors)})
        row['tau_max_samples'] = tau
        row['ess_min'] = float(np.min(np.sum(draws/(2*taus), axis=0)))
        row['rhat_E'] = split_rhat(a[:, :, 0])
        row['rhat_q'] = split_rhat(a[:, :, 2])
        row['block_samples'] = block
        row['blocks_per_chain'] = draws//block
        row['stationarity_z'] = float(max(
            abs(c[:draws//2, 0].mean()-c[draws//2:, 0].mean()) /
            max(np.std(c[:, 0], ddof=1)*np.sqrt(8*tau/draws), 1e-12)
            for c in a))
        row['diagnostic_ok'] = bool(row['ess_min'] >= 200 and
                                     row['rhat_E'] < 1.05 and row['rhat_q'] < 1.05 and
                                     row['blocks_per_chain'] >= 20 and
                                     row['stationarity_z'] < 3)
        rows.append(row)
    return rows
