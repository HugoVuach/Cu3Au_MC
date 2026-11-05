"""Ordre L1₂ local par maille, spectres et longueur de second moment.

Les composantes sont les enveloppes des trois ondes Qx,Qy,Qz.
La FFT est faite sur les centres de mailles de pas a, pas sur la concentration k=0.
Les spectres connexes soustraient la moyenne d'ensemble de CHAQUE chaîne,
jamais la moyenne spatiale de chaque configuration (qui annulerait le mode zéro).
"""
import numpy as np

SIGNS = np.array([[1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])


def local_order(states, L):
    states = np.asarray(states)
    if states.shape[-1] != 4*L**3:
        raise ValueError('Expected 4 L^3 sites, ordered by cell then basis.')
    cells = states.reshape(states.shape[:-1]+(L, L, L, 4))
    return np.einsum('...m,am->...a', cells-.25, SIGNS)


def order_fft(states, L):
    return np.fft.fftn(local_order(states, L), axes=(-4, -3, -2))


def spectra(states, L):
    """states=(samples,sites). Spectres / (3 V), avec V=L^3."""
    f = order_fft(states, L)
    raw = np.mean(abs(f)**2, axis=0).sum(axis=-1)/(3*L**3)
    coherent = (abs(f.mean(axis=0))**2).sum(axis=-1)/(3*L**3)
    connected = raw-coherent
    # Un petit négatif d'arrondi n'est pas une puissance physique négative.
    if connected.min() < -1e-10:
        raise RuntimeError('Negative connected power beyond rounding tolerance.')
    return raw, np.maximum(connected, 0), coherent


def within_variant_spectrum(states, L, minimum=8):
    """Variance conditionnelle à la variante dominante, avec poids de fréquence.

    Diagnostic de fluctuations intra-variante, distinct de la covariance canonique.
    Utile seulement en phase nettement ordonnée ; le conditionnement biaise la
    physique près du désordre et ne doit pas y être appelé longueur d'équilibre.
    """
    f = order_fft(states, L)
    occupations = states.reshape(len(states), L**3, 4).mean(axis=1)
    variant = occupations.argmax(axis=1)
    result = np.zeros((L, L, L))
    used = 0
    for v in range(4):
        subset = f[variant == v]
        if len(subset) < minimum:
            continue
        power = np.maximum(np.mean(abs(subset)**2, axis=0)-abs(subset.mean(axis=0))**2, 0)
        result += len(subset)*power.sum(axis=-1)/(3*L**3)
        used += len(subset)
    return result/used if used else result*np.nan, used/len(states)


def neighboring_modes(spectrum):
    return np.mean([spectrum[tuple(sign if j == axis else 0 for j in range(3))]
                    for axis in range(3) for sign in [-1, 1]])


def second_moment_length(s0, s1, L):
    """xi/a de l'enveloppe sur le réseau de mailles ; NaN si indéterminable.

    Aucun écrêtage à L et aucun remplacement d'une racine négative par zéro.
    """
    if not np.isfinite(s0+s1) or s1 <= 1e-14 or s0 < s1:
        return float('nan')
    return float(np.sqrt(s0/s1-1)/(2*np.sin(np.pi/L)))


def radial_correlation(spectrum):
    """Wiener-Khinchin : C(r)=IFFT(S), puis moyenne des déplacements de même r².

    Avec S=|FFT(psi)|²/(3V), IFFT(S) = <psi(R).psi(R+r)>_R / 3.
    """
    L = spectrum.shape[0]
    correlation = np.fft.ifftn(spectrum).real
    d = np.minimum(np.arange(L), L-np.arange(L))
    r2 = sum(v*v for v in np.meshgrid(d, d, d, indexing='ij'))
    shells = np.unique(r2)
    return np.sqrt(shells), np.array([correlation[r2 == s].mean() for s in shells])


def mode_amplitudes(states, L):
    f = order_fft(states, L)
    indices = [(0, 0, 0)]+[tuple(sign if j == axis else 0 for j in range(3))
                           for axis in range(3) for sign in [-1, 1]]
    return np.stack([f[:, i, j, k, :] for i, j, k in indices], axis=1)


def length_from_amplitudes(amplitudes, L, connected=True):
    power = np.mean(abs(amplitudes)**2, axis=0)
    if connected:
        power -= abs(amplitudes.mean(axis=0))**2
    power = power.sum(axis=-1)/(3*L**3)
    return second_moment_length(float(power[0]), float(power[1:].mean()), L)


def time_lengths(states, L):
    """Hors équilibre : longueurs spectrales indicatives, jamais xi d'équilibre."""
    f = order_fft(states, L)
    power = (abs(f)**2).sum(axis=-1)/(3*L**3)
    qgrid = 2*np.pi*np.fft.fftfreq(L)
    grids = np.meshgrid(qgrid, qgrid, qgrid, indexing='ij')
    khat2 = sum(4*np.sin(k/2)**2 for k in grids)
    results = []
    for s in power:
        s0, s1 = float(s[0, 0, 0]), float(neighboring_modes(s))
        diffuse = s.copy()
        diffuse[0, 0, 0] = 0
        denominator = float(np.sum(khat2*diffuse))
        width = np.sqrt(diffuse.sum()/denominator) if denominator > 1e-14 else np.nan
        results.append([second_moment_length(s0, s1, L), width,
                        s0/s.sum() if s.sum() else np.nan, s0, s1])
    return np.asarray(results)
