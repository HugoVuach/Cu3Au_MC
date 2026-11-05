"""Analyse des configurations spatiales, erreurs par blocs et figures du bonus 2."""
import json
from pathlib import Path
import numpy as np
from .correlation import (spectra, mode_amplitudes, neighboring_modes,
                          second_moment_length, radial_correlation, time_lengths,
                          within_variant_spectrum)
from .statistics import tau_int, split_rhat
from .simulation import save_csv


def mode_powers(amplitudes, L):
    raw = np.mean(abs(amplitudes)**2, axis=0)
    connected = np.maximum(raw-abs(amplitudes.mean(axis=0))**2, 0)
    return raw.sum(axis=-1)/(3*L**3), connected.sum(axis=-1)/(3*L**3)


def estimate_length(powers, L):
    return second_moment_length(float(powers[0]), float(powers[1:].mean()), L)


def bootstrap_lengths(amplitudes, L, block, repetitions, seed):
    rng = np.random.default_rng(seed)
    n = amplitudes[0].shape[0]
    estimates = []
    for _ in range(repetitions):
        powers = []
        for a in amplitudes:
            starts = rng.integers(n, size=int(np.ceil(n/block)))
            indices = ((starts[:, None]+np.arange(block)) % n).ravel()[:n]
            powers.append(mode_powers(a[indices], L)[1])
        estimates.append(estimate_length(np.mean(powers, axis=0), L))
    finite = np.asarray(estimates)[np.isfinite(estimates)]
    if len(finite) < 2:
        return [float('nan'), float('nan')], float('nan'), len(finite)/repetitions
    return np.quantile(finite, [.025, .975]).tolist(), float(np.std(finite, ddof=1)), len(finite)/repetitions


def analyze(output='results/correlation_L6'):
    out = Path(output)
    cfg = json.loads((out/'config.json').read_text(encoding='utf8'))
    L = cfg['L']
    archives = [np.load(p) for p in sorted(out.glob('chain_*.npz'))]
    if len(archives) != len(cfg['seeds']):
        raise ValueError('Missing correlation configuration archives.')
    for a in archives:
        if str(a['config_json']) != json.dumps(cfg, sort_keys=True):
            raise ValueError('Archive configuration mismatch.')
    rows, raw_spectra, connected_spectra = [], [], []
    amplitude_series = []
    for index, temperature in enumerate(cfg['temperatures']):
        raw, connected, amplitudes, conditional, conditional_fractions = [], [], [], [], []
        for a in archives:
            states = a['states'][:, index]
            if not np.all(states.sum(axis=1) == L**3):
                raise ValueError('Incorrect composition in correlation samples.')
            r, c, _ = spectra(states, L)
            raw.append(r)
            connected.append(c)
            amplitudes.append(mode_amplitudes(states, L))
            conditional_spectrum, fraction = within_variant_spectrum(states, L)
            conditional.append(conditional_spectrum)
            conditional_fractions.append(fraction)
        amplitude_series.append(amplitudes)
        raw = np.mean(raw, axis=0)
        connected = np.mean(connected, axis=0)
        raw_spectra.append(raw)
        connected_spectra.append(connected)
        powers = [np.mean(abs(a)**2, axis=-1)/(L**3) for a in amplitudes]
        taus = [tau_int(v) for a, p in zip(archives, powers) for v in [
            a['samples'][:, index, 0], a['samples'][:, index, 2], p[:, 0], p[:, 1:].mean(axis=1)]]
        taus += [tau_int(a[:, 0, component].real) for a in amplitudes for component in range(3)]
        n = len(amplitudes[0])
        block = min(max(2, int(np.ceil(2*max(taus)))), n//4)
        ci, se, fraction = bootstrap_lengths(amplitudes, L, block, cfg['bootstrap'], 731+index)
        q = np.array([a['samples'][:, index, 2] for a in archives])
        conditional_mean = np.mean(conditional, axis=0)
        conditional_xi = second_moment_length(conditional_mean[0, 0, 0], neighboring_modes(conditional_mean), L)
        # Aucun conditionnement publié comme longueur d'équilibre au voisinage
        # de la transition. On ne montre ce diagnostic que pour q>0.8 et T<=1000 K.
        conditional_applicable = temperature <= 1000 and q.mean() > .8 and min(conditional_fractions) >= .9
        energies = np.array([a['samples'][:, index, 0] for a in archives])
        rhat_q, rhat_e = split_rhat(q), split_rhat(energies)
        rhat_s0 = split_rhat(np.array([p[:, 0] for p in powers]))
        # Le mode le plus sensible à l'absence de visites entre variantes.
        rhat_components = max(split_rhat(np.array([a[:, 0, component].real
                                                  for a in amplitudes])) for component in range(3))
        xi = second_moment_length(connected[0, 0, 0], neighboring_modes(connected), L)
        xi_raw = second_moment_length(raw[0, 0, 0], neighboring_modes(raw), L)
        ess = len(archives)*n/(2*max(taus))
        flags = []
        if max(rhat_q, rhat_e, rhat_s0) > 1.05:
            flags.append('chain_disagreement')
        if rhat_components > 1.05:
            flags.append('variant_mixing_unresolved')
        if n/block < 8 or ess < 200:
            flags.append('insufficient_effective_samples')
        if fraction < .95 or not np.isfinite(xi):
            flags.append('unstable_ratio')
        if xi >= L/2:
            flags.append('finite_size_dominated')
        coherent_fraction = float((raw[0, 0, 0]-connected[0, 0, 0])/raw[0, 0, 0])
        row = dict(T_K=temperature, xi_connected_a=xi, xi_raw_a=xi_raw,
            xi_within_variant_a=conditional_xi if conditional_applicable else float('nan'),
            within_variant_coverage=min(conditional_fractions),
            ci95_low_a=ci[0], ci95_high_a=ci[1], bootstrap_se_a=se,
            valid_bootstrap_fraction=fraction, block_samples=block, tau_max_samples=max(taus),
            effective_samples=ess, rhat_q=rhat_q, rhat_E=rhat_e, rhat_S0=rhat_s0,
            rhat_order_components=rhat_components, q_mean=float(q.mean()),
            coherent_fraction_zero_mode=coherent_fraction,
            S0_connected=float(connected[0, 0, 0]), S1_connected=float(neighboring_modes(connected)),
            diagnostics_ok=not flags, flags=';'.join(flags))
        rows.append(row)
        print('T=', temperature, 'xi/a=', round(xi, 3), 'flags=', row['flags'], flush=True)
    save_csv(out/'correlation_vs_temperature.csv', rows)
    np.savez_compressed(out/'spectra.npz', temperatures=cfg['temperatures'],
                        raw=np.array(raw_spectra), connected=np.array(connected_spectra))
    # Sensibilité à la longueur des blocs autour de la transition.
    sensitivity = []
    for i in [int(np.argmin(abs(np.array(cfg['temperatures'])-t))) for t in [1040, 1050, 1060]]:
        base = rows[i]['block_samples']
        for block in sorted(set([max(2, base//2), base, min(n//2, 2*base)])):
            ci, se, fraction = bootstrap_lengths(amplitude_series[i], L, block, cfg['bootstrap'], 1731+i+block)
            sensitivity.append(dict(T_K=cfg['temperatures'][i], block_samples=block,
                                    ci95_low_a=ci[0], ci95_high_a=ci[1], se_a=se,
                                    valid_fraction=fraction))
    save_csv(out/'block_sensitivity.csv', sensitivity)
    with np.load('results/thermal_cycle_L6/visualisation/animations/relaxation_isotherme.npz') as a:
        temporal = time_lengths(a['states'], L)
        time_rows = []
        for i, (sweep, values) in enumerate(zip(a['sweeps'], temporal)):
            ell, width, fraction, s0, s1 = values
            time_rows.append(dict(sweeps=int(sweep), temperature_K=float(a['temperature_K']),
                ell_apparent_a=float(ell), ell_diffuse_width_a=float(width),
                zero_mode_fraction=float(fraction), S0=float(s0), S1=float(s1),
                box_dominated=bool(ell >= L/2), ratio_defined=bool(np.isfinite(ell))))
        save_csv(out/'lengths_vs_time.csv', time_rows)
    summary = dict(config=cfg, snapshots_per_chain_temperature=n,
                   failed_temperatures_K=[r['T_K'] for r in rows if not r['diagnostics_ok']],
                   connected_definition='mean of per-chain <|F|^2>-|<F>|^2; no per-frame spatial centering',
                   spatial_field='three L12 cell envelopes, conventional cell step a',
                   initial_time_q=float(np.sqrt(temporal[0, 3]/L**3)),
                   final_time_q=float(np.sqrt(temporal[-1, 3]/L**3)),
                   time_box_dominated_count=int(np.count_nonzero(temporal[:, 0] >= L/2)),
                   temporal_frames=len(temporal),
                   bootstrap_repetitions=cfg['bootstrap'])
    def clean(value):
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items()}
        if isinstance(value, list):
            return [clean(v) for v in value]
        if isinstance(value, float) and not np.isfinite(value):
            return None
        return value
    (out/'summary.json').write_text(json.dumps(clean(summary), indent=2), encoding='utf8')
    make_figures(out, rows, raw_spectra, connected_spectra, time_rows)
    for a in archives:
        a.close()
    return summary


def make_figures(out, rows, raw, connected, time_rows):
    from .figures import plt, style, save, BLUE, RED, GOLD, GREEN, INK, measurement_points
    style()
    folder = out/'figures'
    folder.mkdir(exist_ok=True)
    T = np.array([r['T_K'] for r in rows])
    L = json.loads((out/'config.json').read_text())['L']
    col = lambda k: np.array([r[k] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax in axes:
        measurement_points(ax, T, col('xi_raw_a'), GOLD, 'Spectre brut', series=1)
        measurement_points(ax, T, col('xi_connected_a'), BLUE, 'Partie connexe')
        measurement_points(ax, T, col('xi_within_variant_a'), GREEN, 'Intra-variante (phase froide)', series=2)
        for t, r in zip(T, rows):
            if np.isfinite(r['ci95_low_a']+r['ci95_high_a']):
                ax.vlines(t, r['ci95_low_a'], r['ci95_high_a'], color=BLUE, alpha=.55, lw=.9)
            if not r['diagnostics_ok']:
                ax.plot(t, r['xi_connected_a'], 'x', color=RED, ms=8)
        ax.axhline(L/2, color=INK, ls=':', label='L/2 : repère de taille finie')
        ax.set(xlabel='Température (K)', ylabel='Longueur / a')
        ax.legend(fontsize=8)
    axes[0].set_title('Estimation canonique, L = 6')
    axes[1].set(xlim=(1000, 1100), title='Voisinage de la transition')
    visible = (T >= 1000) & (T <= 1100)
    ymax = np.nanmax(np.r_[col('xi_raw_a')[visible], col('ci95_high_a')[visible]])
    axes[1].set_ylim(-.3, 1.08*ymax)
    fig.tight_layout(rect=(0, .08, 1, .98))
    save(fig, folder, 'correlation_01_xi',
         'Traits verticaux : bootstrap 95 % conditionnel aux chaînes ; croix : réserve statistique ou taille finie ; pas de barre pour le diagnostic intra-variante.')

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for t, color in zip([900, 1040, 1050, 1060, 1300], [BLUE, GREEN, RED, GOLD, INK]):
        i = int(np.argmin(abs(T-t)))
        radius, cr = radial_correlation(connected[i])
        scale = cr[0]
        measurement_points(axes[0], radius, cr/scale if scale > 0 else cr, color, f'{T[i]:g} K', series=list(T).index(T[i]) % 6)
        axis_spectrum = np.array([np.mean([connected[i][tuple(k if j == axis else 0 for j in range(3))]
                                          for axis in range(3)]) for k in range(L//2+1)])
        measurement_points(axes[1], np.arange(L//2+1)*2*np.pi/L, axis_spectrum, color, f'{T[i]:g} K', series=list(T).index(T[i]) % 6)
    axes[0].axhline(0, color=INK, lw=.6)
    axes[0].set(xlabel='Distance minimale r/a', ylabel='C connecté(r) / C connecté(0)', title='Espace réel : moyennes par coquille')
    axes[1].set(xlabel='Décalage |k| a autour des ondes d’ordre', ylabel='S connecté(k)', title='Fourier : moyenne des trois axes')
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout(rect=(0, .08, 1, .98))
    save(fig, folder, 'correlation_02_real_fourier', 'FFT du champ d’ordre par maille ; aucune hypothèse de décroissance exponentielle ajustée à ces quelques distances.')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    values = np.array([[r[k] for r in time_rows] for k in ['sweeps', 'ell_apparent_a', 'ell_diffuse_width_a', 'zero_mode_fraction']])
    sweep, ell, width, fraction = values
    for ax in axes[0]:
        measurement_points(ax, sweep, ell, BLUE, 'Longueur apparente brute')
        measurement_points(ax, sweep, width, GREEN, 'Largeur diffuse, mode zéro exclu', series=1)
        ax.axhline(L/2, color=INK, ls=':', label='L/2')
        ax.set(xlabel='Balayages depuis la trempe', ylabel='Longueur / a')
        ax.legend(fontsize=8)
    axes[0, 0].set_title('Relaxation à 1000 K : 6000 balayages')
    axes[0, 1].set(xlim=(0, 200), title='Début de la transformation')
    measurement_points(axes[1, 0], sweep, fraction, GOLD, 'Fraction du mode uniforme')
    axes[1, 0].set(xlabel='Balayages', ylabel='S(0) / somme des S(k)', ylim=(0, 1.02))
    measurement_points(axes[1, 1], T, col('coherent_fraction_zero_mode'), RED, 'Part cohérente, campagne canonique')
    axes[1, 1].set(xlabel='Température (K)', ylabel='Part cohérente du mode zéro', ylim=(0, 1.02))
    fig.tight_layout(rect=(0, .075, 1, .97))
    save(fig, folder, 'correlation_03_relaxation', 'Hors équilibre : ni longueur de fluctuation canonique, ni rayon de germe ; le mode uniforme peut dominer la longueur apparente.')

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    for key, label, color in [('rhat_q', 'q', BLUE), ('rhat_S0', 'S(0) brut', RED), ('rhat_order_components', 'Composantes signées', GREEN)]:
        measurement_points(axes[0], T, col(key), color, label)
    axes[0].set_yscale('log')
    axes[0].axhline(1.05, color=INK, ls=':')
    axes[0].set(ylabel='R-hat (échelle logarithmique)', title='Accord des deux chaînes')
    axes[0].legend(fontsize=8)
    measurement_points(axes[1], T, col('effective_samples'), BLUE)
    axes[1].axhline(200, color=INK, ls=':')
    axes[1].set(ylabel='ESS conservateur estimé', title='Autocorrélation E, q et puissances')
    measurement_points(axes[2], T, col('valid_bootstrap_fraction'), GREEN)
    axes[2].set(ylabel='Fraction des ratios définis', ylim=(0, 1.05), title='Stabilité du bootstrap')
    for ax in axes:
        ax.set_xlabel('Température (K)')
    fig.tight_layout(rect=(0, .075, 1, .98))
    save(fig, folder, 'correlation_04_diagnostics', 'Un désaccord des composantes signées peut signaler des variantes différentes ; un bon accord de q ne prouve pas leur exploration.')


if __name__ == '__main__':
    analyze()
