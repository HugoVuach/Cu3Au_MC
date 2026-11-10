"""Analyse croisée L6/L9/L12 depuis les archives, sans nouvelle dynamique.

Les grilles restent distinctes et les statistiques non linéaires sont d'abord
calculées dans chaque chaîne. Les états spatiaux ne sont jamais assimilés à
autant d'observations indépendantes. CLI : python -m cu3au.cross_analysis.
"""
import csv
import json
from pathlib import Path
import numpy as np

from .model import KB
from .simulation import save_csv
from .statistics import tau_int, split_rhat


ROOT = Path(__file__).resolve().parents[1]
SOURCES = ['production', 'production_L9', 'production_L12']


def read_rows(path):
    with Path(path).open(encoding='utf8') as stream:
        return list(csv.DictReader(stream))


def column(rows, key):
    return np.array([float(r[key]) for r in rows])


def energy_channels(samples, n):
    """Excès exacts E-E0, en eV totaux : AuAu(1) et AuCu(2)."""
    return (.161*n*(.375+1.125*samples[..., 3]),
            .0063*n*1.125*(1-samples[..., 4]))


def covariance_budget(samples, temperature, n):
    """Variance par chaîne, puis moyenne ; covariance comprise, diviseur M."""
    e1, e2 = energy_channels(samples, n)
    d1 = e1-e1.mean(axis=1, keepdims=True)
    d2 = e2-e2.mean(axis=1, keepdims=True)
    scale = n*(KB*np.asarray(temperature))**2
    return np.array([(d1*d1).mean(axis=(0, 1))/scale,
                     (d2*d2).mean(axis=(0, 1))/scale,
                     2*(d1*d2).mean(axis=(0, 1))/scale])


def correlation(x, y):
    if min(np.std(x), np.std(y)) < 1e-12:
        return float('nan')
    return float(np.corrcoef(x, y)[0, 1])


def phase_statistics(samples, threshold=.5):
    """Classement descriptif q>seuil, sans identification de phase rigoureuse.

    Les changements sont comptés à température fixe ; avec échanges de
    répliques ils ne correspondent pas à des événements cinétiques physiques.
    """
    high = samples[..., 2] > threshold
    return high.mean(axis=1), np.sum(high[:, 1:] != high[:, :-1], axis=1)


def load_campaigns():
    campaigns = []
    for name in SOURCES:
        path = ROOT/'results'/name
        summary = json.loads((path/'summary.json').read_text(encoding='utf8'))
        validation = json.loads((path/'validation.json').read_text(encoding='utf8'))
        if validation['technical_validation'] != 'passed':
            raise ValueError(f'Unvalidated campaign: {path}')
        cfg = summary['config']
        series = []
        for file in sorted(path.glob('chain_*.npz')):
            with np.load(file) as a:
                if not np.array_equal(a['temperatures'], cfg['temperatures']):
                    raise ValueError('Archive temperature grid mismatch.')
                series.append(a['samples'])
        if len(series) != len(cfg['seeds']):
            raise ValueError('Missing independent chain.')
        samples = np.stack(series)
        rows = read_rows(path/'observables.csv')
        if not np.array_equal(column(rows, 'T_K'), cfg['temperatures']):
            raise ValueError('CSV temperature grid mismatch.')
        campaigns.append(dict(path=path, summary=summary, cfg=cfg, rows=rows,
                              L=cfg['L'], N=summary['n_atoms'], samples=samples))
    return campaigns


def analyze_archives(campaigns, out):
    all_rows, chain_rows, overview, sensitivity = [], [], [], []
    for run in campaigns:
        a, n, L = run['samples'], run['N'], run['L']
        t = np.asarray(run['cfg']['temperatures'])
        cvparts = covariance_budget(a, t, n)
        e1, e2 = energy_channels(a, n)
        reconstruction = np.max(abs(a[..., 0]-(-3.813*n+e1+e2)))
        if reconstruction > 1e-7:
            raise ValueError('Energy/Cowley identity failed.')
        if not np.allclose(cvparts.sum(axis=0), column(run['rows'], 'cv_kB'), atol=1e-8):
            raise ValueError('Covariance budget does not reproduce cV.')
        fractions, changes = phase_statistics(a)
        derivative = np.gradient(column(run['rows'], 'u_eV'), t)/KB
        run_rows = []
        for j, T in enumerate(t):
            b = a[:, :, j]
            q_of_mean = np.sqrt(4/3*np.sum((b[..., 5:9].mean(axis=1)-.25)**2, axis=1))
            tau = np.array([[tau_int(c[:, k]) for k in range(9)] for c in b])
            row = dict(L=L, T_K=float(T), cv_1=float(cvparts[0, j]),
                       cv_2=float(cvparts[1, j]), cv_cross=float(cvparts[2, j]),
                       cv_derivative=float(derivative[j]),
                       q_mean=float(b[..., 2].mean()), q_of_temporal_mean=float(q_of_mean.mean()),
                       q_sqrtN=float(b[..., 2].mean()*np.sqrt(n)),
                       q2=float(np.mean(b[..., 2]**2)), q2_random_exact=3/(n-1),
                       rho_E_q=float(np.nanmean([correlation(c[:, 0], c[:, 2]) for c in b]))
                       if np.var(b[..., 0]) > 1e-20 else float('nan'),
                       tau_E=float(tau[:, 0].max()), tau_q=float(tau[:, 2].max()),
                       tau_alpha1=float(tau[:, 3].max()), tau_alpha2=float(tau[:, 4].max()),
                       tau_occupations=float(tau[:, 5:].max()),
                       rhat_occupations=max(split_rhat(b[..., k]) for k in range(5, 9)),
                       phase_weight_difference=float(abs(fractions[0, j]-fractions[1, j])),
                       q_half_changes=int(changes[:, j].sum()),
                       ordered_fraction=float(fractions[:, j].mean()),
                       cv_between_chain_means=float(np.var(b[..., 0].mean(axis=1))/(n*(KB*T)**2)),
                       local_acceptance=float(np.mean([c['local_acceptance'][j]
                           for c in run['summary']['chains']])),
                       energy_identity_max_error_eV=float(reconstruction))
            run_rows.append(row)
            for c, chain in enumerate(b):
                occ = chain[:, 5:9]
                # La variante dominante est ambiguë lorsque q est petit.
                ordered = chain[:, 2] > .8
                dominant = occ.argmax(axis=1)
                cr = dict(L=L, chain=c+1, T_K=float(T), q=float(chain[:, 2].mean()),
                          q_of_mean=float(q_of_mean[c]),
                          u_eV=float(chain[:, 0].mean()/n),
                          cv_kB=float(np.var(chain[:, 0])/(n*(KB*T)**2)),
                          ordered_fraction=float(fractions[c, j]),
                          q_half_changes=int(changes[c, j]), ordered_samples=int(ordered.sum()))
                cr.update({f'p{k}': float(occ[:, k].mean()) for k in range(4)})
                cr.update({f'variant{k}_given_q_gt_08': float(np.mean(dominant[ordered] == k))
                           if ordered.any() else float('nan') for k in range(4)})
                chain_rows.append(cr)
            if 1040 <= T <= 1065:
                base = int(run['rows'][j]['block_samples'])
                for block in sorted({max(2, base//2), base, min(a.shape[1]//2, 2*base)}):
                    # Moyennes de blocs disjoints : diagnostic de stabilité de
                    # l'erreur, pas un nouveau bootstrap ni un IC de coexistence.
                    blocks = b.shape[1]//block
                    means = b[:, :blocks*block, 2].reshape(len(b), blocks, block).mean(axis=2)
                    sensitivity.append(dict(L=L, T_K=float(T), block_samples=block,
                        blocks_per_chain=blocks,
                        q_batch_se=float(np.sqrt(np.sum(np.var(means, axis=1, ddof=1)/blocks))/len(b))))
        run['cross'] = run_rows
        all_rows.extend(run_rows)
        peak = int(np.argmax(column(run['rows'], 'cv_kB')))
        critical = [run['rows'][i] for i in run['summary']['critical_window_rows']]
        rates = np.mean([c['swap_acceptance'] for c in run['summary']['chains']], axis=0)
        run['swap_rates'] = rates
        cfg = run['cfg']
        overview.append(dict(L=L, N=n, temperatures=len(t), burn=cfg['burn_sweeps'],
            production=cfg['sample_sweeps'], peak_T=float(t[peak]),
            peak_cv=float(column(run['rows'], 'cv_kB')[peak]),
            peak_q=float(column(run['rows'], 'q')[peak]),
            min_ess=min(float(r['ess_min']) for r in critical),
            max_rhat=max(max(float(r[k]) for k in ['rhat_E','rhat_q']) for r in critical),
            critical_ok=run['summary']['critical_diagnostics_ok'],
            min_swap=float(rates.min()), mean_swap=float(rates.mean()),
            weakest_link_low=float(t[np.argmin(rates)]),
            weakest_link_high=float(t[np.argmin(rates)+1]),
            roundtrips=int(sum(sum(c['roundtrips']) for c in run['summary']['chains'])),
            hot_q=float(a[:, :, -1, 2].mean()),
            peak_weight_difference=run_rows[peak]['phase_weight_difference']))
        print(f'Analyzed L={L}: peak {t[peak]}, max energy residual {reconstruction:.3g}', flush=True)
    save_csv(out/'all_observables_cross.csv', all_rows)
    save_csv(out/'chain_variant_statistics.csv', chain_rows)
    save_csv(out/'block_sensitivity.csv', sensitivity)
    save_csv(out/'overview.csv', overview)
    return overview


def analyze_spatial(campaigns, out):
    """Nouvelle exploitation des états L9/L12 ; L6 a sa campagne dédiée.

    On soustrait la moyenne complexe par chaîne, puis on moyenne les
    puissances, avant de prendre le rapport. Le bootstrap vise les séries
    capturées, avec bloc dicté aussi par les composantes signées lentes.
    """
    from .correlation import mode_amplitudes
    from .correlation_analysis import mode_powers, estimate_length, bootstrap_lengths
    targets = [900, 980, 1000, 1020, 1030, 1040, 1045, 1050, 1052.5, 1055, 1060, 1080, 1100, 1300, 4000]
    rows = []
    for run in campaigns[1:]:
        L, t = run['L'], run['cfg']['temperatures']
        # Une seule décompression par archive, pas une par température.
        states, times = [], []
        for file in sorted(run['path'].glob('chain_*.npz')):
            with np.load(file) as a:
                states.append(a['spatial_states'])
                times.append(a['spatial_sample_sweeps'])
        for T in targets:
            if T not in t:
                continue
            j = t.index(T)
            amps = [mode_amplitudes(s[:, j], L) for s in states]
            powers = [mode_powers(f, L) for f in amps]
            raw = np.mean([p[0] for p in powers], axis=0)
            conn = np.mean([p[1] for p in powers], axis=0)
            tau = max(tau_int(x) for f in amps for x in
                      [np.sum(abs(f[:, 0])**2, axis=-1),
                       np.mean(np.sum(abs(f[:, 1:])**2, axis=-1), axis=-1),
                       *[f[:, 0, k].real for k in range(3)]])
            draws = len(amps[0])
            requested = int(np.ceil(2*tau))
            block = min(max(2, requested), draws//4)
            ci, se, valid = bootstrap_lengths(amps, L, block, 300, 14231+L+j)
            xi = estimate_length(conn, L)
            rhat = max(split_rhat(np.stack([f[:, 0, k].real for f in amps])) for k in range(3))
            q2 = np.mean([np.sum(abs(f[:, 0])**2, axis=-1)/(3*L**6) for f in amps])
            # S(0)=V <q²> : vérification directe des captures vs les séries.
            indices = times[0]//run['cfg']['sample_every']-1
            saved_q2 = np.mean(run['samples'][:, indices, j, 2]**2)
            if not np.isclose(q2, saved_q2, atol=1e-10):
                raise ValueError('Spatial states are not synchronized with observations.')
            flags = []
            if xi >= L/2:
                flags.append('finite_size')
            if draws/tau < 200 or draws//block < 20:
                flags.append('low_ESS_or_few_blocks')
            if rhat > 1.05:
                flags.append('signed_component_disagreement')
            if requested > block:
                flags.append('block_capped')
            if valid < .95 or not np.isfinite(xi):
                flags.append('unstable_ratio')
            rows.append(dict(L=L, T_K=T, captures_per_chain=draws,
                capture_every=run['cfg']['store_states_every'],
                xi_connected_a=xi, xi_raw_a=estimate_length(raw, L),
                ci95_low_a=ci[0], ci95_high_a=ci[1], bootstrap_se_a=se,
                valid_bootstrap_fraction=valid, block_samples=block,
                requested_block_samples=requested, effective_samples=draws/tau,
                rhat_order_components=rhat, q2=float(q2),
                S0_connected=float(conn[0]), S1_connected=float(conn[1:].mean()),
                flags=';'.join(flags)))
            print(f'Spatial L={L}, T={T}: xi/a={xi:.3g} {flags}', flush=True)
        del states
    save_csv(out/'spatial_comparison.csv', rows)
    return rows


def make_plots(campaigns, spatial, out):
    from .figures import plt, style, save, measurement_points, BLUE, RED, GREEN, GOLD, INK
    style()
    colors = [BLUE, RED, GREEN]
    figures = []

    def finish(fig, name, title, note):
        fig.suptitle(title, fontsize=16)
        fig.tight_layout(rect=(0, .07, 1, .94))
        save(fig, out, name, note)
        figures.append(name)

    def draw(ax, key, label, limits=(1035, 1070), cross=False, errors=True):
        visible = []
        for i, (run, color) in enumerate(zip(campaigns, colors)):
            rows = run['cross'] if cross else run['rows']
            t = column(rows, 'T_K')
            y = column(rows, key)
            err = column(rows, key+'_se') if errors and key+'_se' in rows[0] else None
            measurement_points(ax, t, y, color,
                               f"L = {run['L']}", err, series=i)
            mask = np.isfinite(y)
            if limits:
                mask &= (t >= limits[0]) & (t <= limits[1])
            visible.extend(y[mask] if err is None else np.r_[y[mask]-err[mask],y[mask]+err[mask]])
        ax.set(xlabel='Température (K)', ylabel=label)
        if limits:
            ax.set_xlim(*limits)
        finite = np.asarray(visible)[np.isfinite(visible)]
        if len(finite):
            low, high = finite.min(), finite.max()
            margin = max((high-low)*.08, abs(high)*.005, 1e-6)
            ax.set_ylim(low-margin,high+margin)
        ax.legend(fontsize=8)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    for ax, (key, label) in zip(axes.flat, [('q','Ordre global q'),('u_eV','u (eV/atome)'),
            ('cv_kB','cV/kB par atome'),('alpha1','Cowley α₁'),('alpha2','Cowley α₂'),('Xq_per_eV','Xq (eV⁻¹)')]):
        draw(ax, key, label)
    finish(fig, 'croise_01_transition', 'L6, L9 et L12 : même figure, trois tailles',
           'Symboles creux, barres ±1 erreur-type ; grilles réelles, sans interpolation. Près du pic, L9/L12 restent sous-échantillonnés.')

    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    for ax, key, label, limits, cross in [
        (axes[0,0], 'q', 'Ordre global q', (650,1300), False),
        (axes[0,1], 'q_sqrtN', 'q √N', (1100,4000), True),
        (axes[1,0], 'binder_q', 'Cumulant Uq (convention scalaire)', (1035,1070), False),
        (axes[1,1], 'rho_E_q', 'Corrélation temporelle E–q', (1000,1100), True)]:
        draw(ax,key,label,limits,cross)
    finish(fig, 'croise_02_ordres', 'Ordre global, biais de taille et fluctuations complémentaires',
           'q est une norme positive ; un cumulant négatif signale des distributions larges, sans prouver à lui seul une transition de premier ordre.')

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    for ax,key,label,cross in [(axes[0,0],'ess_min','ESS estimé (E,q)',False),
        (axes[0,1],'rhat_q','R-hat(q)',False),(axes[0,2],'blocks_per_chain','Blocs par chaîne',False),
        (axes[1,0],'tau_alpha1','τ(α₁), en mesures',True),
        (axes[1,1],'tau_occupations','Maximum τ(pμ), en mesures',True)]:
        draw(ax,key,label,cross=cross,errors=False)
        ax.set_yscale('log')
    axes[0,0].axhline(200,color=INK,ls=':',lw=1)
    axes[0,1].axhline(1.05,color=INK,ls=':',lw=1)
    axes[0,2].axhline(20,color=INK,ls=':',lw=1)
    for i,(run,color) in enumerate(zip(campaigns,colors)):
        t=np.asarray(run['cfg']['temperatures'])
        rates=np.mean([c['swap_acceptance'] for c in run['summary']['chains']],axis=0)
        measurement_points(axes[1,2], (t[:-1]+t[1:])/2,rates,color,f"L = {run['L']}",series=i)
    axes[1,2].set(xlabel='Milieu du lien en température (K)',ylabel='Acceptation des échanges',xlim=(950,1120))
    axes[1,2].legend(fontsize=8)
    finish(fig,'croise_03_diagnostics','Pourquoi une boîte plus grande ne garantit pas un meilleur équilibre',
           'τ(pμ) suit aussi les variantes ; la moyenne des taux d’échange peut masquer un lien presque infranchissable.')

    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    for i,(run,color) in enumerate(zip(campaigns,colors)):
        j=int(np.argmax(column(run['rows'],'cv_kB')))
        T=run['cfg']['temperatures'][j]
        b=run['samples'][:,:,j]
        steps=np.arange(1,b.shape[1]+1)*run['cfg']['sample_every']
        bins_energy=np.histogram_bin_edges(b[...,0].ravel()/run['N'],bins=40)
        for c in range(len(b)):
            axes[i,0].plot(steps,b[c,:,2],lw=.5,color=[BLUE,RED][c],label=f'Chaîne {c+1}')
            axes[i,1].hist(b[c,:,2],bins=np.linspace(0,1,41),density=True,histtype='step',color=[BLUE,RED][c])
            axes[i,2].hist(b[c,:,0]/run['N'],bins=bins_energy,density=True,histtype='step',color=[BLUE,RED][c])
        axes[i,0].set(title=f"L={run['L']} · {T:g} K",ylabel='q',xlabel='Balayages de production')
        axes[i,0].legend(fontsize=8)
        axes[i,1].set(xlabel='q',ylabel='Densité')
        axes[i,2].set(xlabel='u (eV/atome)',ylabel='Densité')
    finish(fig,'croise_04_populations','Au maximum de chaque taille : traces et distributions par chaîne',
           'Températures des lignes différentes ; les histogrammes ne sont pas fusionnés. Les passages à T fixe incluent les échanges de répliques.')

    fig,axes=plt.subplots(2,3,figsize=(13,7.5))
    for i,(run,color) in enumerate(zip(campaigns,colors)):
        t=column(run['cross'],'T_K')
        for k,(key,label) in enumerate([('cv_1','Var(E₁)'),('cv_2','Var(E₂)'),('cv_cross','2 Cov(E₁,E₂)')]):
            measurement_points(axes[0,i],t,column(run['cross'],key),[BLUE,GOLD,RED][k],label,series=k)
        axes[0,i].set(title=f"L={run['L']} : budget exact de cV",xlim=(1035,1070),xlabel='T (K)',ylabel='Contribution à cV/kB')
        axes[0,i].legend(fontsize=8)
        measurement_points(axes[1,i],t,column(run['rows'],'cv_kB'),color,'Variance de E',series=0)
        measurement_points(axes[1,i],t,column(run['cross'],'cv_derivative'),INK,'Différence finie de u(T)/kB',series=1)
        axes[1,i].set(xlim=(1035,1070),xlabel='T (K)',ylabel='cV/kB par atome')
        axes[1,i].legend(fontsize=8)
    finish(fig,'croise_05_energie','Relier les paires, Cowley et les fluctuations d’énergie',
           'Les covariances sont indispensables. La dérivée est une différence finie sur la grille, sans barre d’erreur ni valeur de test statistique.')

    cycle=ROOT/'results/thermal_cycle_L6'
    fig,axes=plt.subplots(2,3,figsize=(13,7.5))
    variables=[('q','Ordre q'),('u_eV','u (eV/atome)'),('cv_kB','Fluctuations cV/kB'),
               ('alpha1','Cowley α₁'),('alpha2','Cowley α₂'),('Xq_per_eV','Xq (eV⁻¹)')]
    for ax,(key,label) in zip(axes.flat,variables):
        draw(ax,key,label,limits=(1000,1100))
        for k,(branch,label_b,color) in enumerate([('heating','L6 chauffage',GOLD),('cooling','L6 refroidissement',INK)]):
            rows=read_rows(cycle/f'observables_{branch}.csv')
            measurement_points(ax,column(rows,'T_K'),column(rows,key),color,label_b,column(rows,key+'_se'),series=k+3)
        ax.legend(fontsize=7)
    finish(fig,'croise_06_protocoles','Trois productions canoniques et les deux branches disponibles du cycle L6',
           'Cette plaquette conserve le cycle L6 comme référence ; les cycles des autres tailles sont comparés séparément après leur achèvement.')

    fig,axes=plt.subplots(2,3,figsize=(13,7.5))
    for i,run in enumerate(campaigns):
        t=run['cfg']['temperatures'];a=run['samples']
        for c,color in enumerate([BLUE,RED]):
            meanq=a[c,:,:,2].mean(axis=0)
            qm=np.sqrt(4/3*np.sum((a[c,:,:,5:9].mean(axis=0)-.25)**2,axis=-1))
            measurement_points(axes[0,i],t,meanq,color,f'⟨q⟩, chaîne {c+1}',series=c)
            measurement_points(axes[0,i],t,qm,color,f'q(⟨p⟩), chaîne {c+1}',series=c+2)
            measurement_points(axes[1,i],t,(a[c,:,:,2]>.5).mean(axis=0),color,f'Chaîne {c+1}',series=c)
        axes[0,i].set(title=f"L={run['L']} : moyenne et norme",xlim=(900,1100),xlabel='T (K)',ylabel='Amplitude')
        axes[1,i].set(xlim=(1035,1070),ylim=(-.03,1.03),xlabel='T (K)',ylabel='Fraction des mesures avec q > 0,5')
        for ax in axes[:,i]:ax.legend(fontsize=7)
    finish(fig,'croise_07_variantes','Variantes, amplitude d’ordre et poids des populations',
           'q > 0,5 est un seuil descriptif commun, pas une frontière de phase ajustée. Une moyenne temporelle peut effacer les variantes.')

    fig,axes=plt.subplots(1,3,figsize=(13,4.6))
    spatial6=read_rows(ROOT/'results/correlation_L6/correlation_vs_temperature.csv')
    for i,(L,color) in enumerate(zip([6,9,12],colors)):
        rows=spatial6 if L==6 else [r for r in spatial if r['L']==L]
        t=column(rows,'T_K');xi=column(rows,'xi_connected_a')
        label=f'L={L}'+(' (campagne dédiée)' if L==6 else ' (captures production)')
        measurement_points(axes[0],t,xi,color,label,series=i)
        measurement_points(axes[1],t,xi/L,color,label,series=i)
        err=column(rows,'bootstrap_se_a')
        measurement_points(axes[2],t,xi,color,label,err,series=i)
    for ax in axes:
        ax.set(xlabel='T (K)',ylabel='ξ connexe / a',xlim=(980,1100));ax.legend(fontsize=7)
    axes[0].set_yscale('log');axes[0].set_title('Grandes valeurs : ordre et mélange de variantes')
    axes[1].set(ylabel='ξ connexe / (La)',title='Rapport à la taille de la boîte')
    axes[1].axhline(.5,color=INK,ls=':')
    axes[2].set(xlim=(1060,1320),ylim=(0,5),title='Phase chaude : barres conditionnelles')
    finish(fig,'croise_08_spatial','Première comparaison spatiale L6/L9/L12',
           'Même définition du champ, protocoles et captures différents ; ξ n’est pas un rayon de germe. Les grands rapports et erreurs sont à contrôler.')
    return figures


def main():
    out=ROOT/'results/comparison_L6_L9_L12'
    out.mkdir(exist_ok=True)
    campaigns=load_campaigns()
    overview=analyze_archives(campaigns,out)
    spatial=analyze_spatial(campaigns,out)
    figures=make_plots(campaigns,spatial,out)
    from .cycle_size_report import completed_cycles
    payload=dict(sizes=overview, figures=figures, sources=SOURCES,
                 heating_cooling_completed_sizes=[r['L'] for r in completed_cycles()],
                 caveats=['L12 has a denser grid and twice the sweep budget.',
                          'Technical validation is not equilibrium validation.',
                          'L6 spatial correlation campaign has a separate protocol.'])
    (out/'summary.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    from .cross_report import write_include
    write_include(campaigns,overview,spatial,figures,out)
    print('Cross analysis, figures and LaTeX include complete.',flush=True)


if __name__=='__main__':
    main()
