"""Sources du bonus 2, analyse des GIF et publication du rapport complet."""
import csv
import json
from pathlib import Path
import shutil
import numpy as np


def read_rows(path):
    with Path(path).open(encoding='utf8') as stream:
        return list(csv.DictReader(stream))


def number(value, digits=3):
    x = float(value)
    return f'{x:.{digits}f}' if np.isfinite(x) else r'\textemdash'


def atomic_copy(source, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix+'.tmp')
    shutil.copy2(source, temporary)
    temporary.replace(destination)


def animation_figure(animation_dir, folder):
    from .figures import plt, style, save, BLUE, RED, GOLD, measurement_points
    style()
    thermal = json.loads((animation_dir/'plaquette_cycle_thermique.json').read_text(encoding='utf8'))
    fixed = json.loads((animation_dir/'plaquette_relaxation_isotherme.json').read_text(encoding='utf8'))
    t = np.array([r['sweeps'] for r in thermal['observables']])
    q = np.array([r['q'] for r in thermal['observables']])
    x = np.array([r['sweeps'] for r in fixed['observables']])
    y = np.array([r['q'] for r in fixed['observables']])
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    for sl, color, label, series in [(slice(0, 47), RED, 'Chauffage', 0), (slice(47, None), BLUE, 'Refroidissement', 1)]:
        measurement_points(axes[0], t[sl]/1e6, q[sl], color, label, series=series)
    axes[0].set(xlabel='Millions de balayages cumulés', ylabel='Ordre global q', title='Cycle thermique : fins de palier')
    axes[0].legend(fontsize=8)
    for ax in axes[1:]:
        measurement_points(ax, x, y, GOLD)
        ax.set(xlabel='Balayages depuis la trempe', ylabel='Ordre global q', ylim=(0, 1.02))
    axes[1].set_title('Relaxation à 1000 K')
    axes[2].set(xlim=(0, 200), title='Les 200 premiers balayages')
    fig.tight_layout(rect=(0, .075, 1, .98))
    save(fig, folder, 'animations_01_evolution',
         'Les instants sont ceux des configurations des GIF ; aucune interpolation temporelle des positions ou des espèces.')
    return thermal, fixed


def write_includes(output='results/correlation_L6', reference='results/production'):
    root = Path(__file__).resolve().parents[1]
    out = Path(output)
    target = Path(reference)/'report'
    target.mkdir(parents=True, exist_ok=True)
    folder = out/'figures'
    folder.mkdir(exist_ok=True)
    animations = root/'results/thermal_cycle_L6/visualisation/animations'
    thermal, fixed = animation_figure(animations, folder)
    summary = json.loads((out/'summary.json').read_text(encoding='utf8'))
    cfg = summary['config']
    rows = read_rows(out/'correlation_vs_temperature.csv')
    time_rows = read_rows(out/'lengths_vs_time.csv')
    thresholds = {threshold: next(r['sweeps'] for r in fixed['observables'] if r['q'] >= threshold)
                  for threshold in [.8, .9]}
    gif_lines = [r'\subsubsection{Deux plaquettes animées : protocole et accès aux fichiers}',
        r"Les fichiers GIF accompagnent le PDF dans le sous-dossier \texttt{animations}. "
        r"Le PDF présente des images fixes et des liens vers les animations externes ; "
        r"il ne joue pas lui-même les GIF. Les liens supposent de conserver ce sous-dossier "
        r"à côté du rapport.",
        r'\begin{itemize}',
        r'\item \href{animations/plaquette_cycle_thermique.gif}{Ouvrir la plaquette du cycle chauffage--refroidissement}.',
        r'\item \href{animations/plaquette_relaxation_isotherme.gif}{Ouvrir la plaquette de relaxation à température fixe}.',
        r'\end{itemize}',
        r"Chaque plaquette contient neuf couches atomiques exactes : trois couches "
        r"parallèles pour chacune des orientations $(001)$, $(110)$ et $(111)$. "
        r"Les positions des coupes et le cadrage restent fixes ; Au est en or et Cu "
        r"en bleu. Pour $(001)$, on suit $z/a=2,2.5,3$ ; pour $(110)$, "
        r"$(x+y)/a=5,5.5,6$ ; pour $(111)$, $(x+y+z)/a=7,8,9$. "
        r"Ces trois orientations sont des vues complémentaires du même volume, "
        r"pas trois simulations distinctes.",
        f"Le GIF thermique comporte {thermal['frames']} images : l'état initial "
        r"puis les 92 fins de palier du premier cycle, de 300 à 4000 K puis en sens "
        r"inverse. Les événements à l'intérieur d'un palier ne sont pas sauvegardés "
        r"dans cette campagne historique. Le compteur cumule équilibration et production.",
        f"Le GIF isotherme comporte {fixed['frames']} états. Il part de l'état final "
        r"à 4000 K du chauffage du cycle 1, puis impose immédiatement 1000 K. "
        r"Un calcul dédié de 6000 balayages, graine 20261104, sauvegarde les "
        r"occupations tous les 20 balayages. Aucun transitoire n'est éliminé. "
        r"La première image est l'état chaud juste après le changement de température, "
        r"avant toute relaxation ; l'étiquette 1000 K indique la consigne. "
        r"Les 200 premiers balayages sont ralentis à la lecture, sans modifier "
        r"les valeurs du compteur.",
        r'\subsubsection{Ce que montrent les animations et ce qu’elles ne démontrent pas}',
        r"Le cycle thermique rend visible le passage entre un motif ordonné et "
        r"une occupation désordonnée. La variante retrouvée au refroidissement "
        r"peut différer de celle du départ : un changement de variante n'est "
        r"pas, à lui seul, une perte d'ordre. Une seule couche peut paraître pauvre "
        r"en Au parce que le sous-réseau privilégié se trouve dans une couche voisine.",
        f"Dans la relaxation, $q$ passe de ${fixed['observables'][0]['q']:.3f}$ "
        f"à ${fixed['observables'][-1]['q']:.3f}$ entre 0 et 6000 balayages. "
        f"Les premiers états enregistrés avec $q\geq0.8$ et $q\geq0.9$ "
        f"apparaissent respectivement à {thresholds[.8]} et {thresholds[.9]} balayages "
        r"(résolution de capture : 20 balayages). L'évolution ultérieure fluctue "
        r"et n'est pas strictement monotone.",
        r"Cette trajectoire montre une mise en ordre rapide de la petite boîte. "
        r"Elle ne suffit pas à identifier un rayon critique, un temps d'attente "
        r"de nucléation ou une loi de croissance de domaines. Il faudrait définir "
        r"un domaine local, suivre ses interfaces et répéter les trempes avec "
        r"plusieurs graines et tailles. Les échanges Au--Cu sont globaux : "
        r"les balayages ne constituent pas une durée physique ni une cinétique "
        r"de diffusion atomique calibrée.",
        r'\fig{animations_01_evolution}{1}{Ordre global des configurations présentes dans les deux GIF et agrandissement du début de la relaxation.}',
        r'\fig{animation_isotherme_debut}{.93}{Neuf coupes de l’état initial chaud, avant relaxation à la nouvelle consigne de 1000 K.}',
        r'\fig{animation_isotherme_fin}{.93}{Les mêmes neuf coupes après 6000 balayages à 1000 K. Leur comparaison indique une mise en ordre ; les images seules ne démontrent pas une nucléation.}']

    correlation_lines = [r'\subsubsection{Pourquoi mesurer une longueur spatiale ?}',
        r"L'amplitude $q$ décrit l'ordre global, tandis que $\alpha_1$ et $\alpha_2$ "
        r"décrivent les deux premières coquilles. La longueur de corrélation $\xi$ "
        r"caractérise la portée spatiale des corrélations de l'ordre chimique. "
        r"Elle n'est pas automatiquement le rayon d'un germe ni la taille "
        r"géométrique d'un domaine.",
        r'\subsubsection{Champ local : conserver les quatre variantes \texorpdfstring{L1$_2$}{L12}}',
        r"On numérote les mailles conventionnelles par $\mathbf R$, de pas $a$, "
        r"et leurs quatre sites par $\mu$. Avec $n_{\mathbf R\mu}=1$ pour Au "
        r"et 0 pour Cu, on définit trois composantes signées :",
        r'\[\psi_\alpha(\mathbf R)=\sum_{\mu=0}^{3}\epsilon_{\alpha\mu}'
        r'(n_{\mathbf R\mu}-1/4),\qquad'
        r'\epsilon=\begin{pmatrix}1&1&-1&-1\\1&-1&1&-1\\1&-1&-1&1\end{pmatrix}.\]',
        r"Dans les quatre états parfaitement ordonnés, le champ vaut "
        r"$(1,1,1)$, $(1,-1,-1)$, $(-1,1,-1)$ ou $(-1,-1,1)$. "
        r"Le signe distingue donc la variante, ce que ne ferait pas une "
        r"amplitude locale toujours positive. Avec $V=L^3$ et "
        r"$\mathbf m=V^{-1}\sum_{\mathbf R}\boldsymbol\psi(\mathbf R)$, "
        r"on retrouve exactement $q^2=|\mathbf m|^2/3$.",
        r"Ce champ représente les enveloppes, par maille, des trois ondes "
        r"d'ordre $\mathbf Q_x=(2\pi/a)(1,0,0)$, $\mathbf Q_y=(2\pi/a)(0,1,0)$ "
        r"et $\mathbf Q_z=(2\pi/a)(0,0,1)$. À décalage nul, ses composantes "
        r"sont exactement les amplitudes de ces ondes. À décalage non nul, "
        r"on étudie la FFT du champ par maille : cette convention d'enveloppe "
        r"est explicitement distincte d'un calcul de diffraction atomique "
        r"incluant les phases de position internes à la base.",
        r"La concentration totale étant fixée, son mode uniforme ne fluctue "
        r"pas. Le mode $\mathbf k=0$ utilisé ci-dessous est celui du champ "
        r"d'ordre, pas celui de la concentration brute.",
        r'\subsubsection{Corrélation réelle, FFT et second moment}',
        r'\[F_\alpha(\mathbf k)=\sum_{\mathbf R}\psi_\alpha(\mathbf R)e^{-i\mathbf k\cdot\mathbf R},'
        r'\qquad S_{\rm brut}(\mathbf k)=\frac1{3V}\sum_\alpha\langle|F_\alpha(\mathbf k)|^2\rangle.\]',
        r'\[S_{\rm c}(\mathbf k)=\frac1{3V}\sum_\alpha\left['
        r'\langle|F_\alpha(\mathbf k)|^2\rangle-|\langle F_\alpha(\mathbf k)\rangle|^2\right].\]',
        r"La soustraction se fait sur la moyenne temporelle complexe de chaque "
        r"chaîne, puis les spectres des chaînes sont moyennés. On ne soustrait "
        r"jamais la moyenne spatiale de chaque configuration : cela annulerait "
        r"artificiellement le mode dont on veut mesurer les fluctuations. "
        r"La FFT inverse donne exactement, sur la boîte périodique,",
        r'\[C_{\rm c}(\mathbf r)=\frac1{3V}\sum_{\mathbf R}\left['
        r'\langle\boldsymbol\psi(\mathbf R)\cdot\boldsymbol\psi(\mathbf R+\mathbf r)\rangle'
        r'-\langle\boldsymbol\psi(\mathbf R)\rangle\cdot\langle\boldsymbol\psi(\mathbf R+\mathbf r)\rangle\right].\]',
        r"Dans la limite de volume infini et pour une corrélation isotrope "
        r"sommable, la longueur de second moment vérifie",
        r'\[\xi_{\rm 2nd}^{\,2}=\frac{\sum_{\mathbf r}|\mathbf r|^2C_{\rm c}(\mathbf r)}'
        r'{6\sum_{\mathbf r}C_{\rm c}(\mathbf r)}.\]',
        r"Le mot moment désigne ici une pondération par la distance au carré, "
        r"pas la variance d'énergie de $c_V$. Sur une boîte périodique, une "
        r"somme naïve en $r^2$ dépend du choix des images. Le code emploie "
        r"l'estimateur spectral fini, avec $\overline S_1$ la moyenne sur les "
        r"six vecteurs $\pm(2\pi/(La))\mathbf e_j$ :",
        r'\[\boxed{\frac{\xi_{2,L}}a=\frac1{2\sin(\pi/L)}'
        r'\sqrt{\frac{S_{\rm c}(0)}{\overline S_1}-1}}.\]',
        r"Le facteur sinus correspond au réseau de mailles de pas $a$. "
        r"Cet estimateur compare le pic d'ordre et ses premiers voisins "
        r"réciproques ; une FFT complète n'est pas indispensable pour ces "
        r"sept modes, mais elle permet aussi de reconstruire $C(\mathbf r)$. "
        r"Il ne prouve pas une forme lorentzienne à tous les vecteurs d'onde. "
        r"La longueur de second moment n'est pas strictement identique à "
        r"la longueur de décroissance exponentielle. On n'ajuste donc pas "
        r"une exponentielle sur les quelques distances accessibles ici.",
        r'\subsubsection{Retirer une moyenne ne supprime pas nécessairement le mélange de variantes}',
        r"La part cohérente $|\langle F\rangle|^2$ est sauvegardée séparément. "
        r"Cependant, si la chaîne visite plusieurs variantes ordonnées, "
        r"$\langle\boldsymbol\psi\rangle$ peut être proche de zéro alors "
        r"que $q$ reste élevé. Le spectre connexe canonique contient alors "
        r"encore la fluctuation collective entre variantes : une grande "
        r"valeur de $\xi_{2,L}$ ne décrit pas les petits défauts d'une "
        r"phase pure.",
        r"Un diagnostic supplémentaire calcule la variance de Fourier "
        r"conditionnellement au sous-réseau majoritaire, puis moyenne ces "
        r"variances avec les fréquences des variantes. Ce diagnostic "
        r"intra-variante enlève la variation entre leurs moyennes. Il est "
        r"affiché uniquement pour $T\leq1000$ K, $\langle q\rangle>0.8$ "
        r"et au moins 90\,\% de couverture par des groupes de huit états "
        r"ou plus. Il s'agit d'un conditionnement explicite, pas d'un "
        r"remplacement de la corrélation canonique ; il serait biaisant "
        r"dans le régime désordonné. Aucune incertitude statistique "
        r"n'est attribuée ici à ce seul diagnostic.",
        r'\subsubsection{Campagne dédiée, incertitudes et contrôles}',
        f"La configuration \texttt{{configs/correlation\_L6.json}} utilise "
        f"$L={cfg['L']}$, {len(cfg['temperatures'])} températures entre "
        f"{min(cfg['temperatures'])} et {max(cfg['temperatures'])} K, "
        f"deux chaînes de graines {cfg['seeds'][0]} et {cfg['seeds'][1]}. "
        f"Après {cfg['burn_sweeps']} balayages d'équilibration, "
        f"{cfg['sample_sweeps']} balayages produisent "
        f"{summary['snapshots_per_chain_temperature']} configurations par chaîne "
        f"et température, espacées de {cfg['sample_every']} balayages. "
        f"Les échanges de répliques ont lieu tous les {cfg['exchange_every']} balayages.",
        r"Les états de départ sont les configurations finales de la campagne "
        r"canonique antérieure à la température disponible la plus proche ; les nouvelles "
        r"mesures proviennent uniquement de la nouvelle production. Contrairement "
        r"aux anciennes archives contenant surtout des observables, les "
        r"nouvelles archives conservent chaque état atomique échantillonné.",
        f"Les intervalles sont les percentiles 2.5 et 97.5 d'un bootstrap "
        f"circulaire de {cfg['bootstrap']} répétitions, par blocs contigus dans "
        r"chaque chaîne, avec recalcul de la moyenne complexe et du rapport "
        r"à chaque répétition. La longueur des blocs vaut environ deux fois "
        r"le plus grand temps d'autocorrélation estimé parmi $E$, $q$, les "
        r"puissances des modes et les composantes signées. Une sensibilité "
        r"à demi/double longueur est conservée autour de la transition. "
        r"Ces intervalles sont conditionnels aux trajectoires observées "
        r"et ne corrigent ni une phase non visitée ni la taille finie.",
        r"Les diagnostics comprennent R-hat de $E$, $q$, $S(0)$ et des "
        r"composantes signées, l'effectif utile, les échanges acceptés "
        r"et la fraction de rapports bootstrap définis. Une racine "
        r"négative ou un dénominateur nul produit une valeur manquante, "
        r"jamais un zéro fabriqué. Une longueur supérieure à $L/2$ est "
        r"signalée sans être écrêtée : ce seuil est un repère prudent, "
        r"pas une borne mathématique de l'estimateur.",
        r'\begin{algocode}{Pseudocode : mesures spatiales et analyse}',
        r'\small\begin{verbatim}',
        'POUR chaque chaine et temperature :',
        '    equilibrer, puis capturer les etats pendant la production',
        '    verifier composition et energie exacte a chaque capture',
        '    construire les trois champs signes psi sur les mailles',
        '    calculer F = FFT(psi)',
        '    accumuler <abs(F)^2> et <F> complexe separement',
        '    S_connexe = (<abs(F)^2> - abs(<F>)^2) / (3*V)',
        '    conserver aussi S_brut et la part coherente',
        'moyenner les spectres des chaines, puis calculer le rapport',
        'reconstruire C = IFFT(S_connexe)',
        'reechantillonner des blocs, recalculer les moyennes et xi',
        'signaler ratios indefinis, diagnostics et taille finie',
        r'\end{verbatim}', r'\end{algocode}',
        r'\subsubsection{Résultats à taille finie et comparaison des définitions}',
        r'\begin{center}\scriptsize\begin{tabular}{rrrrrl}\toprule',
        r'$T$ (K) & $\xi_{\rm brut}/a$ & $\xi_{\rm c}/a$ & IC 95\,\% de $\xi_{\rm c}/a$ & Intra-variante$/a$ & Réserve\\\midrule']
    for r in rows:
        flag = r['flags']
        status = 'taille finie' if flag == 'finite_size_dominated' else ('voir CSV' if flag else 'aucune détectée')
        correlation_lines.append(f"{r['T_K']} & {number(r['xi_raw_a'])} & "
            f"{number(r['xi_connected_a'])} & $[{number(r['ci95_low_a'])},{number(r['ci95_high_a'])}]$ & "
            f"{number(r['xi_within_variant_a'])} & {status}"+r'\\')
    correlation_lines += [r'\bottomrule\end{tabular}\end{center}',
        r"Dans la phase froide, la grande longueur canonique et la longueur "
        r"intra-variante beaucoup plus courte répondent à des questions "
        r"différentes : cohérence de l'ordre à l'échelle de la boîte, ou "
        r"fluctuations conditionnelles autour d'une variante. L'augmentation "
        r"du rapport spectral ne constitue pas la mesure d'un domaine "
        r"géométrique plus grand que le volume.",
        r"Au-dessus de la transition, la longueur devient plus courte. "
        r"Avec seulement six mailles par direction, le premier décalage "
        r"accessible est $k_{\min}a=2\pi/6\simeq1.047$ : le régime "
        r"des très petits vecteurs d'onde est peu résolu. Cette campagne "
        r"ne détermine ni un exposant critique, ni une divergence, ni "
        r"une température critique extrapolée. Il faudrait plusieurs tailles "
        r"et étudier notamment $\xi/L$ en conservant la même définition.",
        r'\fig{correlation_01_xi}{1}{Longueurs spectrales brute et connexe, et diagnostic conditionnel intra-variante à basse température. Les grandes valeurs doivent être interprétées avec les réserves de taille finie.}',
        r'\fig{correlation_02_real_fourier}{1}{Corrélations réelles connexes obtenues par FFT inverse et spectres des enveloppes d’ordre. Un plateau à basse température peut subsister du fait des fluctuations collectives entre variantes.}',
        r'\fig{correlation_04_diagnostics}{1}{Diagnostics de la campagne spatiale ; l’accord statistique ne supprime pas les effets de taille finie.}',
        r'\subsubsection{Après la trempe : distinguer croissance de l’ordre et longueur de fluctuation}',
        r"Les 301 configurations du GIF isotherme permettent une analyse "
        r"spatiale instantanée à 1000 K. On calcule le même rapport de "
        r"puissances sans moyenne d'ensemble, noté $\ell_{\rm app}(t)$. "
        r"Ce n'est pas $\xi_{\rm c}(T)$ : le système relaxe et aucune "
        r"moyenne d'ensemble à temps fixé sur plusieurs trempes n'est disponible.",
        r"On enregistre aussi la fraction $f_0=S(0)/\sum_{\mathbf k}S(\mathbf k)$ "
        r"et une longueur de largeur diffuse, en excluant explicitement "
        r"le mode uniforme :",
        r'\[\frac{\ell_{\rm diff}}a=\left['
        r'\frac{\sum_{\mathbf k\ne0}S(\mathbf k)}'
        r'{\sum_{\mathbf k\ne0}\widehat k_a^{\,2}S(\mathbf k)}'
        r'\right]^{1/2},\qquad \widehat k_a^{\,2}=4\sum_j\sin^2(k_ja/2).\]',
        r"Cette largeur est un second moment dans l'espace réciproque, "
        r"distinct du second moment en espace réel qui définit $\xi$. "
        r"Elle sonde le spectre diffus restant, incluant les défauts et les "
        r"courtes distances, et n'est pas un rayon de domaine.",
        f"Dans cette trajectoire, $\ell_{{\rm app}}>L/2$ pour "
        f"{summary['time_box_dominated_count']} des {summary['temporal_frames']} "
        r"images. Dès 40 balayages, sa valeur est dominée par l'ordre "
        r"qui se forme à l'échelle de la boîte ; à l'instant initial, "
        r"le rapport ne définit pas de longueur réelle. Ces observations "
        r"empêchent d'interpréter une courbe croissante comme une loi "
        r"de croissance de germes. Elles expliquent pourquoi les GIF "
        r"doivent être accompagnés d'une analyse quantitative et de "
        r"plusieurs répétitions.",
        r'\fig{correlation_03_relaxation}{1}{Longueur apparente, largeur diffuse et poids de l’ordre uniforme pendant la relaxation, avec zoom initial ; comparaison à la part cohérente de la campagne canonique.}',
        r'\subsubsection{Fichiers, vérifications et références de méthode}',
        r"Le code est dans \texttt{correlation.py}, \texttt{correlation\_campaign.py}, "
        r"\texttt{correlation\_analysis.py} et \texttt{correlation\_report.py}. "
        r"Les configurations, CSV, spectres complets et diagnostics sont dans "
        r"\texttt{results/correlation\_L6}. Les tests vérifient les quatre "
        r"variantes, $S(0)=V\langle q^2\rangle$, FFT contre somme directe, "
        r"FFT inverse contre corrélation de paires, valeurs indéfinies et "
        r"séparation d'un mélange de variantes parfaitement ordonnées.",
        r"Pour la définition de second moment sur un réseau périodique : "
        r"\href{https://arxiv.org/abs/hep-lat/9112002}{Edwards et al., Multi-Grid Monte Carlo III, équations (3.14)--(3.18)}. "
        r"Pour les trois composantes d'ordre Cu$_3$Au et les domaines : "
        r"\href{https://www.homepages.ucl.ac.uk/~ucapikr/harald3.pdf}{Reichert et al., Europhysics Letters 53, 570--576 (2001)}. "
        r"Ce dernier travail porte sur une surface expérimentale : ses "
        r"mécanismes cinétiques ne sont pas assimilés à notre boîte périodique rigide."]

    for name, lines in [('bonus_animations.tex', gif_lines), ('bonus_correlation.tex', correlation_lines)]:
        (out/name).write_text('\n'.join(lines)+'\n', encoding='utf8')
        atomic_copy(out/name, target/name)
    for path in folder.glob('*.pdf'):
        atomic_copy(path, target/'figures'/path.name)
    for suffix in ['debut', 'fin']:
        atomic_copy(animations/f'plaquette_relaxation_isotherme_{suffix}.pdf',
                    target/'figures'/f'animation_isotherme_{suffix}.pdf')
    for name in ['plaquette_cycle_thermique.gif', 'plaquette_relaxation_isotherme.gif']:
        for parent in [target, root/'report']:
            atomic_copy(animations/name, parent/'animations'/name)


if __name__ == '__main__':
    write_includes()
    from .thermal_cycle_reporting import publish_bonus
    publish_bonus('results/thermal_cycle_L6', 'results/production', complete=True)
    print('Correlation bonus and GIF analysis published; full report compiled.', flush=True)
