"""Comparaisons des seuls cycles entièrement terminés et validés."""
import json
from pathlib import Path
import numpy as np
from .cross_analysis import ROOT, read_rows, column


def completed_cycles():
    runs=[]
    for L in [6,9,12]:
        out=ROOT/'results'/f'thermal_cycle_L{L}'
        if not (out/'summary.json').exists():
            continue
        summary=json.loads((out/'summary.json').read_text(encoding='utf8'))
        if 'passed' not in summary.get('validation',''):
            raise ValueError(f'Unvalidated cycle L{L}')
        cfg=summary['config']
        if cfg['L']!=L or cfg['exchange_every']!=0:
            raise ValueError('Incompatible cycle metadata.')
        runs.append(dict(L=L,path=out,cfg=cfg,summary=summary,
            branches={branch:read_rows(out/f'observables_{branch}.csv')
                      for branch in ['heating','cooling']}))
    return runs


def write_cycle_size_include():
    from .reporting import _copy_atomic
    from .figures import plt, style, save, BLUE, RED, GREEN
    runs=completed_cycles()
    out=ROOT/'results/comparison_L6_L9_L12'
    target=ROOT/'results/production/report'
    lines=[r'\subsection{Partie 4 -- Chauffage et refroidissement : comparaison des tailles}',
        r'\label{sec:cycles-tailles}',
        r'Les cycles séquentiels additionnels conservent le protocole L6 : '
        r'46 températures dans chaque sens, deux trajectoires indépendantes, '
        r'6000 balayages écartés puis 30000 mesurés par palier, pas 10 et '
        r'aucun échange de répliques. Les états sont transmis sans '
        r'réinitialisation du chauffage au refroidissement. Chaque palier '
        r'achevé est sauvegardé pour permettre une reprise exacte. '
        r'Les tailles sont comparées à même budget de balayages, '
        r'ce qui ne garantit pas un même effectif indépendant.']
    plan_path=out/'cycles_plan.json'
    planned=json.loads(plan_path.read_text(encoding='utf8'))['requested_sizes'] if plan_path.exists() else [9]
    missing=[L for L in planned if L not in [r['L'] for r in runs]]
    if missing:
        lines += [r'\begin{attention}\raggedright',
            'Les cycles '+', '.join(f'L{L}' for L in missing)+
            r' sont programmés dans le lancement additionnel ; leurs branches '
            r'ne sont pas encore disponibles dans cette version. Les '
            r'courbes sont ajoutées uniquement après achèvement des '
            r'deux branches et contrôle des archives. Le suivi détaillé '
            r'est dans \texttt{results/comparison\_L6\_L9\_L12/cycles\_execution.json}.',
            r'\par\end{attention}']
    deferred=[L for L in [9,12] if L not in planned and L not in [r['L'] for r in runs]]
    if deferred:
        lines.append('Le cycle '+', '.join(f'L{L}' for L in deferred)+
            r' est différé à la demande de l’utilisateur et ne sera pas lancé '
            r'automatiquement. Les productions canoniques L6/L9/L12 '
            r'déjà terminées restent comparées dans la partie précédente.')
    if len(runs)>=2:
        ref=runs[0]['cfg']
        for run in runs[1:]:
            for key in ['temperatures','seeds','burn_sweeps','sample_sweeps','sample_every','exchange_every']:
                if run['cfg'][key]!=ref[key]:
                    raise ValueError(f'Different cycle protocol for {key}.')
        style()
        variables=[('q','Ordre q'),('u_eV','u (eV/atome)'),('cv_kB','cV fluct / kB'),
                   ('alpha1','Cowley α₁'),('alpha2','Cowley α₂'),('Xq_per_eV','Xq (eV⁻¹)')]
        for name,limits,title in [('croise_09_cycles',(1000,1100),'Chauffage et refroidissement : tailles superposées'),
                                  ('croise_10_cycles_global',(1000,1100),'Cycles thermiques : zoom 1000–1100 K')]:
            fig,axes=plt.subplots(2,3,figsize=(13,7.5))
            for ax,(key,label) in zip(axes.flat,variables):
                bounds=[]
                for run,color in zip(runs,[BLUE,RED,GREEN]):
                    for branch,marker,lab in [('heating','^','chauffage'),('cooling','v','refroidissement')]:
                        rows=run['branches'][branch];t=column(rows,'T_K');y=column(rows,key);err=column(rows,key+'_se')
                        ax.errorbar(t,y,yerr=err,fmt=marker,ls='none',mfc='none',ms=5,
                            color=color,capsize=1.5,elinewidth=.7,label=f"L{run['L']} {lab}")
                        mask=(t>=limits[0])&(t<=limits[1]);bounds.extend(np.r_[y[mask]-err[mask],y[mask]+err[mask]])
                lo,hi=np.min(bounds),np.max(bounds);pad=max((hi-lo)*.07,1e-7)
                ax.set(xlim=limits,ylim=(lo-pad,hi+pad),xlabel='T (K)',ylabel=label)
                ax.legend(fontsize=7,ncols=2)
            fig.suptitle(title,fontsize=16);fig.tight_layout(rect=(0,.065,1,.94))
            save(fig,out,name,'Couleur = taille ; triangle vers le haut = chauffage, vers le bas = refroidissement. Barres conditionnelles ±1 erreur-type.')
            _copy_atomic(out/f'{name}.pdf',target/'figures'/f'{name}.pdf')
        fig,axes=plt.subplots(2,3,figsize=(13,7.5))
        for run,color in zip(runs,[BLUE,RED,GREEN]):
            paired=read_rows(run['path']/'paired_differences.csv')
            for ax,key,label in zip(axes.flat[:4],['q','u_eV','alpha1','alpha2'],
                    ['Δq','Δu (eV/atome)','Δα₁','Δα₂']):
                ax.errorbar(column(paired,'T_K'),column(paired,'delta_'+key),
                    yerr=column(paired,'delta_'+key+'_between_cycles_se'),
                    fmt='o' if run['L']==6 else 's',ls='none',mfc='none',ms=5,
                    color=color,capsize=2,label=f"L{run['L']}")
                ax.axhline(0,color='#183047',ls=':',lw=.8)
                ax.set(xlabel='T (K)',ylabel=label,xlim=(1035,1070));ax.legend(fontsize=8)
            for branch,marker,label in [('heating','^','chauffage'),('cooling','v','retour')]:
                rows=run['branches'][branch];t=column(rows,'T_K')
                for ax,y,ylabel in [(axes[1,1],column(rows,'ess_min'),'ESS (E,q)'),
                    (axes[1,2],np.maximum(column(rows,'rhat_E'),column(rows,'rhat_q')),'Maximum R-hat(E,q)')]:
                    ax.plot(t,y,marker=marker,ls='none',mfc='none',ms=5,color=color,
                        label=f"L{run['L']} {label}")
                    ax.set(xlabel='T (K)',ylabel=ylabel,xlim=(1035,1070))
                axes[1,1].set_yscale('log')
        axes[1,1].axhline(200,color='#183047',ls=':',lw=1)
        axes[1,2].axhline(1.05,color='#183047',ls=':',lw=1)
        # Les bornes verticales doivent suivre le zoom, sans les R-hat
        # des séries froides presque constantes hors de la fenêtre.
        critical_rhat=[max(float(r['rhat_E']),float(r['rhat_q'])) for run in runs
            for rows in run['branches'].values() for r in rows if 1035<=float(r['T_K'])<=1070]
        axes[1,2].set_ylim(.99,max(critical_rhat)+.03)
        for ax in axes[1,1:]:ax.legend(fontsize=7,ncols=2)
        fig.suptitle('Écarts appariés des branches et fiabilité statistique',fontsize=16)
        fig.tight_layout(rect=(0,.065,1,.94))
        save(fig,out,'croise_11_cycles_diagnostics',
            'Δ = refroidissement − chauffage ; barres des écarts : dispersion entre deux cycles seulement. Les seuils de diagnostic ne prouvent pas l’équilibre.')
        _copy_atomic(out/'croise_11_cycles_diagnostics.pdf',target/'figures/croise_11_cycles_diagnostics.pdf')
        # Les deux exports ont désormais la même fenêtre demandée : une seule
        # plaquette est reproduite dans le rapport pour éviter un doublon.
        size_labels=', '.join(f'L{r["L"]}' for r in runs)
        lines += [r'\fig{croise_10_cycles_global}{1}{Superposition des cycles chauffage/refroidissement '+size_labels+r' entre 1000 et 1100 K : couleurs pour les tailles et triangles opposés pour les deux sens du parcours.}',
                  r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
                  r'Taille & Pic chauffage (K) & Pic retour (K) & Paliers insuffisants aller & retour\\\midrule']
        gaps=[]
        for run in runs:
            s=run['summary'];heat=run['branches']['heating'];cool=run['branches']['cooling']
            if not np.array_equal(column(heat,'T_K'),column(cool,'T_K')):
                raise ValueError('Cannot compare different branch grids.')
            peaks=s['peaks_cv_grid_K'];bad=s['diagnostic_failed_temperatures_K']
            lines.append(f"L{run['L']} & {peaks['heating']:g} & {peaks['cooling']:g} & "
                         f"{len(bad['heating'])} & {len(bad['cooling'])}"+r'\\')
            difference=column(cool,'q')-column(heat,'q');j=int(np.argmax(abs(difference)))
            gaps.append(f"L{run['L']} : $\\Delta q={difference[j]:.3f}$ à {float(heat[j]['T_K']):g} K")
        lines += [r'\bottomrule\end{tabular}\end{center}',
            r'Le plus grand écart de moyenne $\Delta q=q_\downarrow-q_\uparrow$ '
            r'sur les températures mesurées vaut '+ '; '.join(gaps)+'.',
            r'Ces extrêmes sont sélectionnés après inspection de la grille et '
            r'n’ont pas un statut de test statistique. Les différences appariées '
            r'par trajectoire restent disponibles dans chaque '
            r'\texttt{paired\_differences.csv}. Les deux trajectoires d’un '
            r'même volume sont indépendantes ; leurs deux branches ne le sont '
            r'pas, et les volumes ne sont pas appariés configuration par configuration.',
            r'Une évolution de l’écart des branches avec L peut associer une '
            r'barrière plus coûteuse et une relaxation moins complète. Les '
            r'paliers signalés insuffisants empêchent de transformer ces boucles '
            r'en coexistence d’équilibre, en spinodales expérimentales ou en '
            r'taux de nucléation. L’acceptation globale Au--Cu ne constitue '
            r'pas une dynamique de diffusion locale. La comparaison aux '
            r'productions canoniques se fait avec ces réserves de protocole.']
        lines += [r'\fig{croise_11_cycles_diagnostics}{1}{Écarts appariés entre branches pour chaque taille et contrôles ESS/R-hat. La dispersion entre deux trajectoires ne constitue pas un test robuste de significativité.}']
        lines += detailed_cycle_analysis(runs,out,target)
        lines += detailed_L12_analysis(runs,out,target)
    text='\n\n'.join(lines)+'\n'
    (out/'bonus_size_cycles.tex').write_text(text,encoding='utf8')
    _copy_atomic(out/'bonus_size_cycles.tex',target/'bonus_size_cycles.tex')
    summary_path=out/'summary.json'
    if summary_path.exists():
        summary=json.loads(summary_path.read_text(encoding='utf8'))
        summary['heating_cooling_completed_sizes']=[r['L'] for r in runs]
        summary['deferred_thermal_cycle_sizes']=deferred
        summary['cycle_figures']=[name for name in ['croise_09_cycles','croise_10_cycles_global',
            'croise_11_cycles_diagnostics','croise_12_cycle_L9_reference',
            'croise_13_cycle_L12_reference','croise_14_cycle_L12_populations'] if (out/f'{name}.pdf').exists()]
        summary_path.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    return runs


def detailed_cycle_analysis(runs,out,target):
    """Lecture des mesures L9 finies, sans rééchantillonner la dynamique."""
    from .figures import plt,style,save,measurement_points,BLUE,RED,GREEN
    from .reporting import _copy_atomic
    from .thermal_cycle import load_cycle_archives
    from .simulation import save_csv
    selected=next((r for r in runs if r['L']==9),None)
    if selected is None:
        return []
    equilibrium=read_rows(ROOT/'results/production_L9/observables.csv')
    heat=selected['branches']['heating'];cool=selected['branches']['cooling']
    style();fig,axes=plt.subplots(2,3,figsize=(13,7.5))
    for ax,key,label in zip(axes.flat,
            ['q','u_eV','cv_kB','alpha1','alpha2','binder_q'],
            ['Ordre q','u (eV/atome)','cV fluct / kB','Cowley α₁','Cowley α₂','Cumulant Uq']):
        bounds=[]
        for i,(rows,color,name) in enumerate([(equilibrium,BLUE,'Canonique avec répliques'),
                (heat,RED,'Chauffage séquentiel'),(cool,GREEN,'Refroidissement séquentiel')]):
            t=column(rows,'T_K');y=column(rows,key);err=column(rows,key+'_se')
            measurement_points(ax,t,y,color,name,err,series=i)
            mask=(t>=1035)&(t<=1070);bounds.extend(np.r_[y[mask]-err[mask],y[mask]+err[mask]])
        lo,hi=np.min(bounds),np.max(bounds);pad=max((hi-lo)*.08,1e-6)
        ax.set(xlabel='T (K)',ylabel=label,xlim=(1035,1070),ylim=(lo-pad,hi+pad))
        ax.legend(fontsize=7)
    fig.suptitle('L9 : trois protocoles comparés à taille fixée',fontsize=16)
    fig.tight_layout(rect=(0,.065,1,.94))
    save(fig,out,'croise_12_cycle_L9_reference',
         'Même Hamiltonien, même taille et même grille ; protocoles distincts. La production canonique elle-même reste insuffisante près du pic.')
    _copy_atomic(out/'croise_12_cycle_L9_reference.pdf',target/'figures/croise_12_cycle_L9_reference.pdf')
    _copy_atomic(selected['path']/'figures/bonus_03_sous_reseaux.pdf',target/'figures/cycle_L9_sous_reseaux.pdf')
    lines=[r'\subsubsection{L9 terminé : distinguer accord des maxima et équilibre des branches}',
        r'Les 184 paliers L9 sont achevés. Les contrôles des énergies, '
        r'de composition, des neuf observables et de continuité entre '
        r'paliers passent. Le maximum des fluctuations est à 1055 K '
        r'dans les deux sens, comme dans la production canonique L9. '
        r'Cet accord porte sur une case de la grille de 5 K ; il '
        r'ne démontre ni disparition de toute mémoire ni équilibre complet.',
        r'\fig{croise_12_cycle_L9_reference}{1}{À L9 fixé, superposition de la production canonique et des branches séquentielles. Les différences de moyennes et de réponses se lisent avec leurs erreurs et les diagnostics.}',
        r'\begin{center}\small\begin{tabular}{llrrrr}\toprule',
        r'Protocole L9 & $T$ (K) & $\langle q\rangle$ & SE$(q)$ & $c_V^{\rm fluct}/k_B$ & ESS\\\midrule']
    comparison=[]
    for T in [1050,1055]:
        for name,rows in [('Canonique',equilibrium),('Chauffage',heat),('Retour',cool)]:
            r=next(x for x in rows if float(x['T_K'])==T)
            lines.append(f"{name} & {T} & {float(r['q']):.3f} & {float(r['q_se']):.3f} & "
                         f"{float(r['cv_kB']):.2f} & {float(r['ess_min']):.1f}"+r'\\')
            comparison.append(dict(L=9,T_K=T,protocol=name,**{k:float(r[k]) for k in
                ['q','q_se','u_eV','cv_kB','alpha1','alpha2','ess_min','rhat_E','rhat_q','blocks_per_chain']}))
    save_csv(out/'cycle_L9_protocol_comparison.csv',comparison)
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'À 1050 K, le refroidissement L9 reste en moyenne plus ordonné '
        r'que le chauffage : $\Delta q=+0.0796$. Mais l’erreur de dispersion '
        r'entre les deux cycles appariés vaut 0.1538, supérieure à cet '
        r'écart. L6 donne au même point un écart négatif, '
        r'$-0.0684$, avec dispersion 0.0557. Le changement de signe '
        r'interdit d’affirmer que l’hystérésis augmente régulièrement '
        r'avec le volume. Il peut refléter les temps de résidence, '
        r'les passages aléatoires entre populations et le faible nombre '
        r'de répétitions. On ne déduit aucune largeur d’hystérésis '
        r'significative de ces seuls extrêmes.',
        r'Le lien avec l’énergie est cohérent : à 1050 K en L9, '
        r'le retour a un $u$ inférieur de 0.001924 eV/atome, '
        r'$\Delta\alpha_1=-0.00861$ et $\Delta\alpha_2=+0.05137$. '
        r'La population plus ordonnée favorise donc les voisinages '
        r'de $L1_2$ et abaisse l’énergie. L’identité énergie--Cowley '
        r'expliquée dans la partie précédente relie ces écarts : '
        r'ils ne constituent pas des confirmations indépendantes.',
        r'\subsubsection{Pourquoi les contrôles critiques ne permettent pas une conclusion plus forte}',
        r'\begin{center}\small\begin{tabular}{llrrrr}\toprule',
        r'Branche L9 & $T$ (K) & ESS & R-hat max & Blocs/chaîne & Écart demi-séries\\\midrule']
    for name,rows in [('Chauffage',heat),('Retour',cool)]:
        for T in [1045,1050,1055,1060]:
            r=next(x for x in rows if float(x['T_K'])==T)
            lines.append(f"{name} & {T} & {float(r['ess_min']):.1f} & "
                f"{max(float(r['rhat_E']),float(r['rhat_q'])):.3f} & {int(r['blocks_per_chain'])} & "
                f"{float(r['stationarity_z']):.2f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Au chauffage à 1050 K, R-hat atteint environ 1.29 ; '
        r'au maximum de 1055 K, l’ESS est inférieur à 15 dans les '
        r'deux branches. Inversement, un ESS additionné supérieur à 200 '
        r'au retour à 1060 K coexiste avec R-hat supérieur à 1.08 : '
        r'un critère isolé ne suffit pas. La présence de seulement '
        r'cinq températures signalées par branche L9, contre neuf '
        r'et dix pour L6, ne mesure pas une meilleure précision '
        r'au maximum ; la fenêtre affectée est plus étroite mais '
        r'le ralentissement y est sévère. La production canonique '
        r'L9 reste elle aussi insuffisante près de son maximum.',
        r'\subsubsection{Trajectoires individuelles, variantes et régimes froid/chaud}']
    cfg,data=load_cycle_archives(selected['path'])
    j=cfg['temperatures'].index(1050)
    lines += [r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
        r'À 1050 K, L9 & $q_\uparrow$ & $q_\downarrow$ & $\Delta q$ & Fraction $q>0.5$ aller/retour\\\midrule']
    for c in range(2):
        h=data['heating'][c][j]['samples'][:,2];k=data['cooling'][c][j]['samples'][:,2]
        lines.append(f"Cycle {c+1} & {h.mean():.3f} & {k.mean():.3f} & {k.mean()-h.mean():+.3f} & "
                     f"{(h>.5).mean():.3f}/{(k>.5).mean():.3f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Le premier cycle donne ainsi un écart de signe opposé au second. '
        r'Le seuil $q>0.5$ sert uniquement à décrire les populations '
        r'visitées ; il n’est pas une frontière de phase ajustée. '
        r'Leur pondération explique une grande partie du désaccord '
        r'de l’ordre moyen, sans mettre en cause le label géométrique $\mu$.',
        r'\fig{cycle_L9_sous_reseaux}{1}{Occupations L9 séparées par trajectoire et par sens. Les moyennes ne doivent pas être fusionnées en un instantané ; au retour, la seconde trajectoire sélectionne une autre variante.}',
        r'À 300 K, la première trajectoire retrouve la variante 0 ; '
        r'la seconde part de la variante 2 et finit sur la variante 3. '
        r'Dans les deux cas $q=1$ à la précision numérique : '
        r'changer de sous-réseau privilégié n’est pas perdre l’ordre. '
        r'Le diagnostic insuffisant à 500 K appartient à un autre '
        r'régime : la moyenne de $q$ reste supérieure à 0.9995 dans chaque trajectoire, '
        r'avec seulement six ou sept valeurs distinctes d’énergie. '
        r'La statistique des excitations rares est limitée ; ce '
        r'signalement ne localise pas une seconde transition.',
        r'À 1300 K, les branches L9 coïncident pratiquement '
        r'($q\simeq0.0624$, $\alpha_1\simeq-0.152$, '
        r'$\alpha_2\simeq0.153$), avec des diagnostics favorables. '
        r'À 4000 K, $q\simeq0.0354$, contre environ 0.0657 pour L6 : '
        r'la diminution suit le biais de norme finie décrit dans '
        r'la partie 3. Les ordres locaux restent non nuls. '
        r'Les zones éloignées du pic sont donc plus reproductibles '
        r'entre protocoles que la pondération des populations critiques.',
        r'\begin{oral}{Cycle L9 -- Quelle conclusion retenir après le calcul ?}',
        r'\begin{itemize}',
        r'\item \textbf{Un résultat terminé n’est pas un équilibre garanti.} '
        r'Les archives sont complètes et cohérentes ; les deux branches '
        r'placent leur maximum à 1055 K. Mais les longues corrélations '
        r'et les différences entre trajectoires limitent les valeurs '
        r'quantitatives au pic. Une barre d’erreur conditionnelle '
        r'ne corrige pas un temps de résidence mal exploré.',
        r'\item \textbf{La mémoire apparente doit être répétée.} À 1050 K, '
        r'les deux trajectoires donnent des écarts de signes opposés. '
        r'Il faut davantage de trajectoires et des paliers critiques '
        r'plus longs avant d’attribuer un sens et une largeur robustes '
        r'à la boucle. Ce protocole à échanges globaux n’établit '
        r'pas une vitesse physique de nucléation.',
        r'\end{itemize}\end{oral}']
    return lines


def detailed_L12_analysis(runs,out,target):
    """Explique les résultats L12 finaux, en distinguant mémoire et équilibre."""
    from .figures import plt,style,save,measurement_points,BLUE,RED,GREEN
    from .reporting import _copy_atomic
    from .thermal_cycle import load_cycle_archives
    from .simulation import save_csv
    run=next((r for r in runs if r['L']==12),None)
    if run is None:
        return []
    equilibrium=read_rows(ROOT/'results/production_L12/observables.csv')
    cfg,data=load_cycle_archives(run['path'])
    style();fig,axes=plt.subplots(2,3,figsize=(13,7.5))
    for ax,key,label in zip(axes.flat,['q','u_eV','cv_kB','alpha1','alpha2','binder_q'],
            ['Ordre q','u (eV/atome)','cV fluct / kB','Cowley α₁','Cowley α₂','Cumulant Uq']):
        bounds=[]
        for i,(rows,color,name) in enumerate([(equilibrium,BLUE,'Production avec répliques'),
                (run['branches']['heating'],RED,'Chauffage'),(run['branches']['cooling'],GREEN,'Refroidissement')]):
            t=column(rows,'T_K');y=column(rows,key);err=column(rows,key+'_se')
            measurement_points(ax,t,y,color,name,err,series=i)
            mask=(t>=1035)&(t<=1070);bounds.extend(np.r_[y[mask]-err[mask],y[mask]+err[mask]])
        lo,hi=np.min(bounds),np.max(bounds);pad=max((hi-lo)*.08,1e-6)
        ax.set(xlim=(1035,1070),ylim=(lo-pad,hi+pad),xlabel='T (K)',ylabel=label)
        ax.legend(fontsize=7)
    fig.suptitle('L12 : mémoire thermique et production canonique',fontsize=16)
    fig.tight_layout(rect=(0,.065,1,.94))
    save(fig,out,'croise_13_cycle_L12_reference',
        'La production à répliques a une grille plus fine et deux fois plus de balayages ; un bon diagnostic dans une seule population ne prouve pas son poids d’équilibre.')
    _copy_atomic(out/'croise_13_cycle_L12_reference.pdf',target/'figures/croise_13_cycle_L12_reference.pdf')
    fig,axes=plt.subplots(2,3,figsize=(13,7.5))
    statistics=[]
    for row,T in enumerate([1050,1055]):
        j=cfg['temperatures'].index(T)
        for col,branch in enumerate(['heating','cooling']):
            for c,color in enumerate([BLUE,RED]):
                a=data[branch][c][j]['samples'];q=a[:,2]
                steps=np.arange(1,len(q)+1)*cfg['sample_every']
                axes[row,col].plot(steps,q,color=color,lw=.6,label=f'Trajectoire {c+1}')
                axes[row,2].hist(q,bins=np.linspace(0,1,41),density=True,histtype='step',
                    color=color,ls='-' if branch=='heating' else '--',
                    label=f'{"Chauffage" if branch=="heating" else "Retour"} {c+1}')
                statistics.append(dict(L=12,T_K=T,branch=branch,trajectory=c+1,
                    q=float(q.mean()),fraction_q_gt_05=float(np.mean(q>.5)),
                    u_eV=float(a[:,0].mean()/(4*12**3))))
            axes[row,col].set(xlabel='Balayages de production',ylabel='q',ylim=(-.02,1.02),
                title=f'{T} K · '+('chauffage' if branch=='heating' else 'refroidissement'))
            axes[row,col].legend(fontsize=7)
        axes[row,2].set(xlabel='q',ylabel='Densité',title=f'{T} K · distributions séparées')
        axes[row,2].legend(fontsize=6.5,ncols=2)
    fig.suptitle('L12 : des poids de populations très différents selon le trajet',fontsize=16)
    fig.tight_layout(rect=(0,.065,1,.94))
    save(fig,out,'croise_14_cycle_L12_populations',
        'Séries temporelles sans échange de répliques ; les balayages ne sont pas un temps physique. Les histogrammes gardent les trajectoires et branches séparées.')
    _copy_atomic(out/'croise_14_cycle_L12_populations.pdf',target/'figures/croise_14_cycle_L12_populations.pdf')
    save_csv(out/'cycle_L12_population_statistics.csv',statistics)
    lines=[r'\subsubsection{L12 achevé : une mémoire plus marquée, des poids d’équilibre encore incertains}',
        r'Les 184 paliers L12 sont achevés et validés techniquement. '
        r'Le maximum de variance est à 1055 K au chauffage et 1050 K '
        r'au refroidissement. La production avec échanges de répliques '
        r'présente son maximum à 1052.5 K, avec une grille deux fois '
        r'plus fine et deux fois plus de balayages. Ces trois nombres '
        r'désignent des critères observés dans des protocoles différents. '
        r'Le milieu des deux maxima de branche n’est pas un nouvel '
        r'estimateur validé de la transition, même s’il coïncide ici '
        r'avec le maximum de la production à répliques.',
        r'\fig{croise_13_cycle_L12_reference}{1}{L12 à taille fixée : comparaison de la production à répliques et des deux branches. Les points à 1052.5 K appartiennent uniquement à la grille canonique.}',
        r'\begin{center}\small\begin{tabular}{llrrrr}\toprule',
        r'Protocole L12 & $T$ (K) & $\langle q\rangle$ & SE$(q)$ & ESS & R-hat max\\\midrule']
    for T in [1050,1055]:
        for name,rows in [('Répliques',equilibrium),('Chauffage',run['branches']['heating']),
                          ('Retour',run['branches']['cooling'])]:
            r=next(x for x in rows if float(x['T_K'])==T)
            lines.append(f"{name} & {T} & {float(r['q']):.3f} & {float(r['q_se']):.3f} & "
                f"{float(r['ess_min']):.1f} & {max(float(r['rhat_E']),float(r['rhat_q'])):.3f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'À 1050 K, le chauffage reste presque entièrement ordonné '
        r'($\langle q\rangle=0.772$), tandis que le retour mélange '
        r'des populations et donne 0.302. Les deux trajectoires '
        r'apparient chacune leurs branches :',
        r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
        r'L12, 1050 K & $q_\uparrow$ & $q_\downarrow$ & $\Delta q$ & Fraction $q>0.5$ aller/retour\\\midrule']
    for c in [1,2]:
        h=next(r for r in statistics if r['T_K']==1050 and r['branch']=='heating' and r['trajectory']==c)
        k=next(r for r in statistics if r['T_K']==1050 and r['branch']=='cooling' and r['trajectory']==c)
        lines.append(f"Trajectoire {c} & {h['q']:.3f} & {k['q']:.3f} & {k['q']-h['q']:+.3f} & "
            f"{h['fraction_q_gt_05']:.3f}/{k['fraction_q_gt_05']:.3f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Les deux écarts L12 ont le même signe, contrairement au cas '
        r'L9 à cette température. Leur moyenne est $\Delta q=-0.4695$, '
        r'avec une erreur de dispersion entre seulement deux trajectoires '
        r'de 0.1781. C’est une mémoire de parcours nette dans les '
        r'trajectoires observées ; cela ne donne ni un intervalle de '
        r'confiance fiable à 95\,\%, ni une loi de taille. Avec deux '
        r'répétitions, il ne faut pas transformer automatiquement '
        r'« moyenne $\pm2$ SE » en intervalle de confiance.',
        r'L’écart énergétique est cohérent : au retour à 1050 K, '
        r'$\Delta u=+0.01048$ eV/atome, '
        r'$\Delta\alpha_1=+0.04695$ et $\Delta\alpha_2=-0.27893$. '
        r'L’état moyen de retour est donc moins ordonné et plus '
        r'coûteux énergétiquement. Cette différence de moyennes '
        r'hors équilibre n’est pas une chaleur latente mesurée.',
        r'\fig{croise_14_cycle_L12_populations}{1}{Traces et distributions L12 à 1050 et 1055 K. Un écart de moyennes peut provenir de temps de résidence différents dans les populations ; les données ne sont pas fusionnées.}',
        r'Au chauffage à 1055 K, une trajectoire donne '
        r'$\langle q\rangle\simeq0.536$, l’autre 0.084. R-hat '
        r'atteint environ 1.44. Au retour à 1050 K, l’ESS est '
        r'proche de 20 et R-hat de 1.45. Dans ces deux paliers, '
        r'le choix des blocs est plafonné à 750 mesures, avec '
        r'seulement quatre blocs par chaîne : l’incertitude '
        r'conditionnelle elle-même reste fragile.',
        r'Inversement, le chauffage à 1050 K et le retour à 1055 K '
        r'passent les diagnostics locaux, tout en restant presque '
        r'entièrement dans des populations opposées. Une petite '
        r'erreur et R-hat proche de 1 signifient alors que les '
        r'trajectoires se ressemblent \emph{dans la population visitée}. '
        r'Ils ne prouvent pas que son poids canonique soit correct. '
        r'Le passage de neuf/dix températures signalées en L6 '
        r'à seulement deux par branche en L12 ne prouve donc '
        r'pas une meilleure exploration de la transition.',
        r'À 300 K, les deux trajectoires retrouvent $q=1$ : '
        r'la première passe de la variante 0 à la variante 3, '
        r'la seconde conserve la variante 2. Au chaud, les branches '
        r'redeviennent proches et $q\simeq0.023$ à 4000 K, '
        r'ce qui complète la tendance de taille finie observée en '
        r'L6/L9. Les GIF restent ceux de L6, et les longueurs '
        r'L12 de la partie 3 viennent des captures canoniques, '
        r'pas d’un nouveau GIF ou d’une analyse spatiale du cycle.',
        r'\begin{oral}{Bilan final -- Que permettent réellement les trois tailles ?}',
        r'\begin{itemize}',
        r'\item \textbf{Résultats solides dans le modèle.} Le fondamental '
        r'et sa normalisation sont établis ; l’ordre global diminue '
        r'au chaud, les préférences locales persistent et le résidu '
        r'de norme décroît avec le volume. Les protocoles situent '
        r'la transformation dans la même zone de 1050--1055 K.',
        r'\item \textbf{Limite principale restante.} Les grandes boîtes '
        r'ressèrent la transition et accentuent la dépendance au '
        r'parcours, mais les poids des populations ne sont pas '
        r'assez bien contrôlés au pic. Il faut davantage de '
        r'trajectoires, des paliers critiques plus longs et '
        r'des vérifications du mélange avant une extrapolation '
        r'de température ou une largeur d’hystérésis quantitative.',
        r'\item \textbf{Portée physique.} Augmenter L ne rétablit '
        r'ni vibrations, ni dilatation, ni diffusion par lacunes. '
        r'Une longueur spectrale n’est pas un rayon de germe ; '
        r'un balayage n’est pas une seconde ; les maxima du '
        r'modèle ne sont pas ceux du matériau réel.',
        r'\end{itemize}\end{oral}']
    return lines
