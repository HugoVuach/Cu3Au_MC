"""Synthèse pédagogique et critique des trois productions sauvegardées."""
import json
from pathlib import Path
import numpy as np

from .cross_analysis import column, read_rows, ROOT


def fmt(x, digits=3):
    if not np.isfinite(float(x)):
        return r'non défini'
    return f'{float(x):.{digits}f}'


def write_include(campaigns, overview, spatial, figures, out):
    from .reporting import _copy_atomic
    target=ROOT/'results/production/report'
    lines=[
        r'\subsubsection{Question commune : que change réellement l’augmentation de taille ?}',
        r'Les trois productions sont terminées et techniquement validées. '
        r'Les contrôles de composition, d’énergie et d’identité des observables '
        r'passent pour $L=6,9,12$ ; les diagnostics statistiques près du maximum '
        r'passent pour L6 mais restent insuffisants pour L9 et L12. '
        r'La comparaison ci-dessous exploite les archives brutes, y compris des '
        r'observables qui n’étaient pas représentées dans les anciennes figures. '
        r'Elle remplace la précédente synthèse à deux tailles ; les exports '
        r'L6/L9 et L6/L12 restent disponibles sans être reproduits une seconde fois.',
        r'\begin{center}\small\begin{tabular}{lrrr}\toprule',
        r'Protocole & $L=6$ & $L=9$ & $L=12$\\\midrule']
    for label,key,digits in [('Nombre de sites','N',0),('Températures simulées','temperatures',0),
        ('Balayages écartés par chaîne','burn',0),('Balayages de production','production',0),
        ('Maximum sur la grille (K)','peak_T',1),('Hauteur du maximum $c_V/k_B$','peak_cv',2),
        ('ESS minimal dans la fenêtre du pic','min_ess',1),('R-hat maximal $(E,q)$','max_rhat',3),
        ('Allers-retours froid--chaud--froid','roundtrips',0)]:
        lines.append(label+' & '+' & '.join(fmt(s[key],digits) for s in overview)+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Les deux premières campagnes utilisent 46 températures et le même '
        r'budget. L12 utilise 90 températures, un pas de 2.5 K près du pic '
        r'et deux fois plus de balayages. Les observables sont mesurées tous '
        r'les 10 balayages, les sessions d’échange tous les 5, avec deux chaînes '
        r'par taille. L12 ne constitue donc pas un contrôle où seule la taille '
        r'change. Les points sont les températures effectivement simulées ; '
        r'aucune interpolation ne crée une mesure entre deux points.',
        r'\begin{oral}{1. Comment comparer trois tailles sans confondre volume et précision ?}',
        r'\begin{itemize}',
        r'\item \textbf{Commencer par les quantités intensives.} Doubler le '
        r'côté de L6 à L12 multiplie le nombre d’atomes par huit, mais '
        r'ne doit pas multiplier l’énergie \emph{par atome} par huit. '
        r'On compare donc $u=E/N$, $c_V=C_V/N$ et les paramètres normalisés. '
        r'Un balayage effectue en revanche $N$ propositions : un même nombre '
        r'de balayages représente davantage de calcul, sans assurer davantage '
        r'd’états indépendants.',
        r'\item \textbf{Lire les contrôles avant de conclure.} Une barre d’erreur '
        r'petite est conditionnelle aux bassins visités. L’augmentation du volume '
        r'peut diminuer l’arrondi de la transition tout en augmentant les séjours '
        r'métastables. Un désaccord des chaînes interdit de présenter le '
        r'déplacement du maximum comme une mesure précise de la limite infinie.',
        r'\end{itemize}\end{oral}',
        r'\subsubsection{Ordre global et local : augmentation sous la transition, diminution au chaud}',
        r'\fig{croise_01_transition}{1}{Superposition des trois tailles : amplitude d’ordre, énergie, capacité calorifique, Cowley et fluctuations de q. La fenêtre 1035--1070 K rend visible la transition sans traits reliant les mesures.}',
        r'À 1045 K, avant le pic, les amplitudes moyennes sont respectivement '
        +', '.join('$'+fmt(next(r for r in run['rows'] if float(r['T_K'])==1045)['q'])+'$' for run in campaigns)
        +r'. L’augmentation de $q$ avec la taille est compatible avec une phase '
        r'ordonnée moins contaminée par l’arrondi de taille finie. À 1000 K, '
        r'les trois valeurs restent proches de 0.88--0.89 : augmenter la boîte '
        r'ne crée donc pas un ordre supplémentaire à toutes les températures. '
        r'À 1055 K, L6 est encore dans une zone large de mélange, alors que L12 '
        r'est surtout dans la population peu ordonnée. Comparer les courbes à '
        r'température fixe et comparer leurs valeurs au maximum propre de chaque '
        r'taille répondent ainsi à deux questions différentes.',
        r'Les deux paramètres de Cowley se rapprochent dans les phases '
        r'éloignées de la transition, alors que $q$ présente une forte dépendance '
        r'en taille au chaud. À 1300 K, $\alpha_1\simeq-0.152$ et '
        r'$\alpha_2\simeq0.153$ pour les trois tailles : les préférences '
        r'chimiques locales persistent après la perte de l’ordre global. '
        r'Ce constat associe deux échelles spatiales distinctes ; il ne faut '
        r'pas appeler les paramètres de Cowley des fractions de phase.',
        r'\fig{croise_02_ordres}{1}{Ordre sur une gamme élargie, amplitude multipliée par la racine de N au chaud, cumulant et corrélation énergie--ordre. Les deux derniers diagnostics complètent les courbes moyennes.}',
        r'Dans un mélange strictement aléatoire à composition fixée, le '
        r'nombre d’Au d’un sous-réseau suit une loi hypergéométrique. '
        r'Elle donne $\operatorname{Var}(p_\mu)=9/[16(N-1)]$, donc '
        r'$\langle q^2\rangle=3/(N-1)$. On attend une amplitude résiduelle '
        r'd’ordre $N^{-1/2}$, même sans ordre macroscopique. '
        r'À 4000 K, les valeurs mesurées de $\langle q\rangle$ sont '
        +', '.join('$'+fmt(s['hot_q'],5)+'$' for s in overview)+r' ; '
        r'les produits $\sqrt N\langle q\rangle$ valent '
        +', '.join('$'+fmt(s['hot_q']*np.sqrt(s['N']),3)+'$' for s in overview)
        +r'. Cette quasi-superposition soutient l’interprétation en fluctuations '
        r'finies. Ce n’est pas un ajustement d’exposant critique : 4000 K '
        r'est encore une température finie, avec des corrélations locales.',
        r'\subsubsection{Réponses, coexistence et rôle des poids des populations}',
        r'Le maximum de $c_V$ est successivement à 1050, 1055 et 1052.5 K : '
        r'il ne se déplace pas de façon monotone avec L. Sa hauteur augmente '
        r'fortement, ce qui est compatible avec un arrondi moindre et avec '
        r'une variance due à des populations d’énergie différentes. '
        r'Il serait prématuré d’en déduire une loi $c_V^{\max}\propto L^3$, '
        r'une chaleur latente ou un ordre thermodynamique de transition : '
        r'la résolution de grille et surtout les poids d’équilibre ne sont '
        r'pas contrôlés avec une précision comparable.',
        r'\begin{center}\small\begin{tabular}{lrrr}\toprule',
        r'Critère sur la grille (K) & L6 & L9 & L12\\\midrule']
    for key,label,kind in [('cv_kB',r'Maximum $c_V$', 'max'),
        ('Xq_per_eV',r'Maximum $X_q$','max'),('binder_q',r'Minimum $U_q$ près du pic','min')]:
        values=[]
        for run in campaigns:
            rows=[r for r in run['rows'] if 1030<=float(r['T_K'])<=1100]
            chosen=(max if kind=='max' else min)(rows,key=lambda r:float(r[key]))
            values.append(fmt(chosen['T_K'],1))
        lines.append(label+' & '+' & '.join(values)+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Les maxima de réponse ne sont pas plusieurs mesures indépendantes : '
        r'ils utilisent les mêmes configurations. Le minimum négatif de $U_q$, '
        r'particulièrement marqué pour L12, renseigne sur une distribution '
        r'large ou très asymétrique de $q$. La convention '
        r'$1-\langle q^4\rangle/(3\langle q^2\rangle^2)$ est scalaire '
        r'alors que $q$ est une norme à trois composantes ; sa limite chaude '
        r'n’est pas zéro. Un creux négatif ne suffit donc pas, seul, à établir '
        r'une coexistence thermodynamique.',
        r'\fig{croise_04_populations}{1}{Traces et distributions séparées par chaîne, au maximum propre à chaque taille. Les deux populations deviennent plus distinctes, mais les temps de séjour et leurs poids diffèrent entre chaînes.}',
        r'Pour chiffrer cette difficulté, on classe à titre descriptif les '
        r'mesures selon $q>0.5$. Ce seuil commun n’est pas une séparatrice '
        r'de phase déterminée par une énergie libre. Au maximum propre à '
        r'chaque taille, on obtient :',
        r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
        r'Taille & $T$ (K) & Fraction chaîne 1 & Fraction chaîne 2 & Changements\\\midrule']
    for run,s in zip(campaigns,overview):
        j=int(np.argmax(column(run['rows'],'cv_kB')))
        a=run['samples'][:,:,j,2]>.5
        lines.append(f"L{run['L']} & {s['peak_T']:g} & {a[0].mean():.3f} & {a[1].mean():.3f} & "
                     f"{np.sum(a[:,1:]!=a[:,:-1])}"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'À L12, les fractions des deux chaînes diffèrent de '
        f"{100*overview[2]['peak_weight_difference']:.1f} points de pourcentage. "
        r'C’est un désaccord majeur sur les populations visitées, malgré '
        r'la présence de passages dans les deux séries. Un histogramme '
        r'fusionné aurait masqué cette asymétrie.',
        r'Le comptage de changements dépend du seuil et du pas de mesure, '
        r'et inclut des changements par échange de répliques. Il ne mesure '
        r'ni un taux de nucléation ni un temps de basculement atomique. '
        r'Un désaccord de fractions explique simultanément le désaccord de '
        r'$\langle q\rangle$, celui de l’énergie et l’instabilité des variances.',
        r'\subsubsection{Audit statistique : ce que L12 n’a pas réparé}',
        r'\fig{croise_03_diagnostics}{1}{Diagnostics croisés : effectif utile, R-hat, nombre de blocs, lenteur de Cowley et des sous-réseaux, acceptation par lien en température. La moyenne des acceptations ne résume pas la circulation.}',
        r'\begin{center}\small\begin{tabular}{lrrrrl}\toprule',
        r'Taille, $T$ (K) & ESS & R-hat max & Blocs & Écart demi-séries & Statut\\\midrule']
    for run in campaigns[1:]:
        for j in run['summary']['critical_window_rows']:
            r=run['rows'][j]
            lines.append(f"L{run['L']}, {float(r['T_K']):g} & {float(r['ess_min']):.1f} & "
                f"{max(float(r['rhat_E']),float(r['rhat_q'])):.3f} & {int(r['blocks_per_chain'])} & "
                f"{float(r['stationarity_z']):.2f} & "
                +('satisfait' if r['diagnostic_ok']=='True' else 'insuffisant')+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'L12 résout mieux les températures voisines de la transition et '
        r'augmente les séries, mais à 1052.5 K le maximum de R-hat dépasse '
        r'1.26 et l’ESS minimal reste inférieur à 45. À 1050 K, un ESS '
        r'additionné apparemment élevé coexiste avec seulement sept blocs '
        r'par chaîne : une chaîne rapide peut masquer la lenteur de l’autre '
        r'dans une somme d’ESS. Le contrôle de la chaîne la plus lente et '
        r'le désaccord entre chaînes restent indispensables.',
        r'La moyenne des acceptations est proche de 61\,\% pour L9 et L12, '
        r'mais aucun aller-retour complet froid--chaud--froid n’est enregistré. '
        r'Les liens les moins efficaces sont :',
        r'\begin{center}\small\begin{tabular}{lrrr}\toprule',
        r'Taille & Lien (K) & Acceptation minimale & Acceptation moyenne\\\midrule']
    for s in overview:
        lines.append(f"L{s['L']} & {s['weakest_link_low']:g}--{s['weakest_link_high']:g} & "
            f"{100*s['min_swap']:.3f}\,\% & {100*s['mean_swap']:.1f}\,\%"+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'La densification de L12 a donc réellement amélioré le plus '
        r'mauvais lien : environ 29\,\% d’acceptation minimale, contre '
        r'0.37\,\% pour L9. Le maintien de séjours lents et l’absence '
        r'd’aller-retour montrent cependant qu’un chevauchement suffisant '
        r'entre températures voisines ne supprime pas les barrières entre '
        r'populations. Avec davantage de températures, il faut aussi '
        r'davantage de déplacements le long de l’échelle pour la parcourir.',
        r'Le manque d’allers-retours est un indicateur global, sans être '
        r'une preuve que toutes les moyennes à chaque température soient '
        r'fausses. Il faut allonger les fenêtres critiques, vérifier les '
        r'temps de résidence et adapter la grille aux liens réellement '
        r'limitants. Les trois tailles actuelles ne justifient aucune '
        r'extrapolation quantitative vers $L\to\infty$.',
        r'La nouvelle analyse estime aussi $\tau_{\alpha_1}$, '
        r'$\tau_{\alpha_2}$ et $\max_\mu\tau_{p_\mu}$ : ces colonnes '
        r'n’étaient pas toutes tracées auparavant. Une occupation signée '
        r'peut relaxer plus lentement que la norme $q$ lorsque la boîte '
        r'change de variante. Un R-hat des occupations élevé au froid '
        r'peut simplement signaler des variantes différentes dans les '
        r'chaînes ; il limite l’estimation de l’occupation canonique '
        r'moyenne, sans invalider automatiquement une énergie invariante. '
        r'Les erreurs des moyennes de $q$ ont été recalculées pour des blocs '
        r'disjoints de longueurs moitié, référence et double dans '
        r'\texttt{block\_sensitivity.csv}. Ce contrôle reste conditionnel '
        r'aux trajectoires et devient fragile lorsqu’il ne reste que '
        r'quelques blocs.']
    lines += [r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
              r'Au maximum de chaque taille & $b$ & SE$_{b/2}(q)$ & SE$_b(q)$ & SE$_{2b}(q)$\\\midrule']
    sensitivity=read_rows(out/'block_sensitivity.csv')
    for run,s in zip(campaigns,overview):
        cases=[r for r in sensitivity if int(r['L'])==run['L'] and float(r['T_K'])==s['peak_T']]
        cases.sort(key=lambda r:int(r['block_samples']))
        lines.append(f"L{run['L']} & {cases[1]['block_samples']} & "+
                     ' & '.join(fmt(r['q_batch_se'],4) for r in cases)+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Ces erreurs par moyennes de blocs n’incluent pas une correction '
        r'de dispersion entre chaînes. Elles augmentent encore avec la '
        r'longueur des blocs, surtout pour L12 ; on ne peut donc pas '
        r'affirmer une stabilité complète. Les barres publiées conservent '
        r'par ailleurs le maximum entre bootstrap et dispersion des chaînes.',
        r'\begin{center}\small\begin{tabular}{lrrr}\toprule',
        r'Position du pic : bootstrap conditionnel & Bloc (mesures) & IC 95\,\% (K) & SE (K)\\\midrule']
    for run in campaigns:
        p=json.loads((run['path']/'peak_uncertainty.json').read_text(encoding='utf8'))
        base=next(c for c in p['cases'] if c['block_samples']==p['base_block_samples'])
        low,high=base['ci95_K']
        lines.append(f"L{run['L']} & {p['base_block_samples']} & $[{low:g},{high:g}]$ & "
                     +fmt(base['se_K'],2)+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Le bootstrap du maximum rééchantillonne les températures '
        r'conjointement dans chaque chaîne. Ses intervalles se chevauchent, '
        r'mais ils ignorent les bassins manqués et le désaccord de poids '
        r'd’équilibre. Le résultat L12 n’acquiert donc pas une précision '
        r'physique de 2.5 K simplement parce que sa grille est plus fine.',
        r'\subsubsection{Analyse croisée exacte : Cowley, défauts et capacité calorifique}',
        r'Le comptage des paires relie exactement l’énergie et les deux '
        r'paramètres de Cowley. En notant $A_1$ le nombre de paires Au--Au '
        r'de première coquille et $B_2$ celui des paires Au--Cu de seconde coquille,',
        r'\[\frac{A_1}{N}=\frac38+\frac98\alpha_1,\qquad '
        r'\frac{B_2}{N}=\frac98(1-\alpha_2),\]',
        r'\[u=-3.813+0.161\left(\frac38+\frac98\alpha_1\right)'
        r'+0.0063\frac98(1-\alpha_2)\quad\mathrm{eV/atome}.\]',
        r'Cette identité est vérifiée sur toutes les mesures des trois '
        r'productions : les accords entre énergie et Cowley ne sont donc '
        r'pas trois preuves indépendantes de transition. Ils constituent '
        r'un contrôle de cohérence du calcul et relient les défauts '
        r'locaux au coût énergétique.',
        r'Pour $E_1=0.161A_1$ et $E_2=0.0063B_2$, les fluctuations satisfont',
        r'\[\frac{c_V}{k_B}=\frac{\operatorname{Var}(E_1)'
        r'+\operatorname{Var}(E_2)+2\operatorname{Cov}(E_1,E_2)}'
        r'{N(k_BT)^2}.\]',
        r'La covariance est indispensable, même si la seconde interaction '
        r'est petite : les deux défauts se développent sur les mêmes '
        r'configurations. Le budget ci-dessous est calculé dans chaque '
        r'chaîne avant moyenne et reproduit le $c_V$ archivé.',
        r'\fig{croise_05_energie}{1}{Décomposition exacte des fluctuations entre première coquille, seconde coquille et covariance. En bas, confrontation descriptive à la pente discrète de l’énergie moyenne.}',
        r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
        r'Au pic & Var$(E_1)$ & Var$(E_2)$ & $2$ Cov$(E_1,E_2)$ & Somme\\\midrule']
    for run in campaigns:
        j=int(np.argmax(column(run['rows'],'cv_kB')));r=run['cross'][j]
        parts=[r[k] for k in ['cv_1','cv_2','cv_cross']]
        lines.append(f"L{run['L']} & "+' & '.join(fmt(v,3) for v in parts+[sum(parts)])+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'Les valeurs sont sans dimension, en unités de $c_V/k_B$ par atome. '
        r'La comparaison à $\Delta u/(k_B\Delta T)$ sert de contrôle '
        r'descriptif : une dérivée discrète arrondit un pic étroit et utilise '
        r'des moyennes corrélées entre températures. Un désaccord au pic '
        r'peut venir de cette discrétisation, des poids des populations ou '
        r'de l’équilibration ; il ne démontre pas, seul, une erreur de '
        r'formule de capacité calorifique.',
        r'Un contrôle supplémentaire mesure le terme artificiel qui '
        r'serait ajouté si l’on fusionnait les chaînes avant de calculer '
        r'la variance : $\operatorname{Var}_c(\overline E_c)/[N(k_BT)^2]$. '
        r'Ce terme n’est pas inclus dans les résultats thermiques. '
        r'Au pic propre à chaque taille, il vaudrait respectivement '
        +', '.join('$'+fmt(run['cross'][int(np.argmax(column(run['rows'],'cv_kB')))]['cv_between_chain_means'],3)+'$'
                   for run in campaigns)+r'. Une divergence des chaînes '
        r'aurait ainsi gonflé artificiellement le pic L12 si cette '
        r'précaution méthodologique n’avait pas été prise.',
        r'L’acceptation des échanges Au--Cu locaux, distincte de celle des '
        r'échanges de répliques, vaut au maximum respectivement '
        +', '.join('$'+fmt(100*run['cross'][int(np.argmax(column(run['rows'],'cv_kB')))]['local_acceptance'],2)+r'\,\%$'
                   for run in campaigns)+r'. Un mouvement accepté '
        r'n’est pas une information indépendante : ces taux ne remplacent '
        r'ni le temps d’autocorrélation ni le nombre de passages entre bassins.',
        r'\subsubsection{Variantes et occupations : relire les bosses sans inventer d’ordre}',
        r'\fig{croise_07_variantes}{1}{Moyenne de la norme q et norme des occupations moyennes, séparées par chaîne et par taille ; fractions descriptives des deux populations autour de la transition.}',
        r'Une augmentation d’une occupation particulière peut signaler '
        r'un changement de variante et non une augmentation de l’ordre '
        r'global. Les tableaux sauvegardent maintenant $p_0,\ldots,p_3$ '
        r'par chaîne et les fréquences des variantes dominantes conditionnées '
        r'à $q>0.8$. À faible ordre, l’argument maximal de $p_\mu$ '
        r'est ambigu ; on ne l’interprète donc pas comme une variante '
        r'physiquement bien formée. L’écart entre $\langle q\rangle$ et '
        r'$q(\langle\mathbf p\rangle)$ est un diagnostic de moyenne, '
        r'pas une mesure directe du nombre de domaines.',
        r'\subsubsection{Superposition avec le chauffage et le refroidissement réellement disponibles}',
        r'\fig{croise_06_protocoles}{1}{Superposition des trois productions canoniques et des branches de chauffage/refroidissement L6. Chaque symbole indique aussi le protocole ; aucune branche L9 ou L12 n’est fabriquée à partir d’une production canonique.}',
        r'Une production par échanges de répliques mesure des températures '
        r'fixes simultanément ; ce n’est pas une rampe de chauffage, même '
        r'si son CSV est trié par température croissante. Le cycle séquentiel '
        r'propage au contraire le dernier état d’un palier vers le suivant, '
        r'puis revient depuis le véritable état chaud. Les deux branches '
        r'sont appariées à l’intérieur d’un cycle et portent sa mémoire.',
        r'Les archives achevées de chauffage/refroidissement existent '
        r'uniquement pour L6 à la date de cette analyse. La figure permet '
        r'donc de comparer cinq séries effectivement disponibles : trois '
        r'productions canoniques et les deux branches L6. Une superposition '
        r'des branches pour chacune des trois tailles nécessitera deux '
        r'campagnes de cycle dédiées. L’écart des branches L6 reste un effet '
        r'de protocole à durée finie : l’énergie et les ordres peuvent être '
        r'comparés, mais leur variance ne devient un $c_V$ d’équilibre '
        r'que si le palier est stationnaire. Une hystérésis observée ne '
        r'fournit donc pas directement deux températures de transition.',
        r'\subsubsection{Données spatiales jusque-là inexploitées : première comparaison des trois tailles}',
        r'Les productions L9 et L12 conservent chacune 300 configurations '
        r'par température et par chaîne, séparées respectivement de 100 '
        r'et 200 balayages. Les états sélectionnés à 15 températures '
        r'cibles sont analysés ici avec le même champ local à trois '
        r'composantes que dans la partie 2. L6 utilise la campagne spatiale '
        r'dédiée déjà décrite : les budgets ne sont pas identiques. '
        r'La comparaison illustre une tendance, sans constituer une étude '
        r'de convergence à précision statistique égale.',
        r'On calcule les amplitudes aux sept modes $0$ et '
        r'$\pm k_{\min}\mathbf e_j$, on soustrait les moyennes complexes '
        r'dans chaque chaîne, puis on moyenne les puissances connexes '
        r'avant de former le rapport de second moment. Les captures sont '
        r'vérifiées contre les mesures de $q$ au même balayage via '
        r'$S(0)=L^3\langle q^2\rangle$. Les incertitudes utilisent '
        r'300 tirages de blocs circulaires avec longueur choisie aussi '
        r'd’après les trois composantes signées lentes. Le nombre de blocs, '
        r'les racines indéfinies et les désaccords de chaînes sont signalés '
        r'dans \texttt{spatial\_comparison.csv}.',
        r'\fig{croise_08_spatial}{1}{Comparaison des longueurs connexes : valeurs absolues, rapport à la boîte et régime chaud avec erreurs-types conditionnelles. L6 vient d’une campagne dédiée ; L9 et L12 de captures des productions principales.}',
        r'\begin{center}\small\begin{tabular}{lrrrr}\toprule',
        r'Taille et température & $\xi/a$ & ESS spatial & R-hat signé max & Limites détectées\\\midrule']
    for L in [9,12]:
        for T in [1000,1050,1060,1300]:
            r=next(x for x in spatial if x['L']==L and x['T_K']==T)
            # Statut court ; toutes les causes restent dans le CSV.
            flags=[]
            if 'finite_size' in r['flags']:flags.append('boîte')
            if 'low_ESS' in r['flags']:flags.append('statistique')
            if 'disagreement' in r['flags']:flags.append('variantes')
            if 'capped' in r['flags']:flags.append('blocs bornés')
            if 'unstable' in r['flags']:flags.append('rapport')
            lines.append(f"L{L}, {T} K & {fmt(r['xi_connected_a'],2)} & "
                f"{fmt(r['effective_samples'],1)} & {fmt(r['rhat_order_components'],2)} & "
                +(', '.join(flags) if flags else 'aucune détectée')+r'\\')
    lines += [r'\bottomrule\end{tabular}\end{center}',
        r'La longueur ne croît pas nécessairement avec L : à 1000 K, '
        r'la valeur connexe L12 est inférieure à celle de L9, malgré '
        r'la boîte plus grande. Les réserves sur les composantes signées '
        r'montrent que cette comparaison mélange des fluctuations '
        r'intra-variante et des visites collectives de variantes. '
        r'Ce n’est pas une mesure d’un domaine qui rétrécirait en '
        r'augmentant le volume. À 1060 K, L12 a simultanément un $q$ '
        r'plus faible et une longueur plus courte que L6/L9 : sa '
        r'population est déjà surtout peu ordonnée, tandis que '
        r'les autres boîtes restent davantage dans la queue de '
        r'transition. Il faut donc croiser $\xi$, les distributions '
        r'et les poids, plutôt que comparer seulement trois longueurs.',
        r'Une grande valeur froide reflète aussi l’ordre uniforme et la '
        r'visite incomplète des variantes. Un $\xi$ supérieur au côté '
        r'de la boîte n’est pas le rayon d’un domaine plus grand que '
        r'le système. À haute température, les longueurs beaucoup plus '
        r'courtes et les paramètres de Cowley proches entre tailles '
        r'sont compatibles avec un ordre local de portée finie. '
        r'On ne lit pas un croisement universel de $\xi/(La)$ dans ces '
        r'quelques points, et l’on n’ajuste ni exposant critique ni '
        r'loi de croissance des germes.',
        r'\subsubsection{Bilan croisé et données complémentaires accessibles}',
        r'\begin{itemize}',
        r'\item \textbf{Concordance physique.} Le fondamental reste le même '
        r'par atome pour toutes les tailles. La chute de $q$, la hausse de '
        r'l’énergie, les changements de Cowley, les distributions et les '
        r'réponses situent tous une zone proche de 1050--1055 K. Cette '
        r'concordance est physiquement cohérente, mais les identités '
        r'énergie--Cowley et les données partagées empêchent de la compter '
        r'comme autant de preuves statistiques indépendantes.',
        r'\item \textbf{Ce que le volume améliore.} L’ordre résiduel au chaud '
        r'diminue comme une amplitude de taille finie ; la variation près '
        r'de la transition se resserre et la résolution réciproque '
        r'$k_{\min}a=2\pi/L$ s’améliore de 1.047 à 0.698 puis 0.524. '
        r'Les préférences locales loin du pic sont beaucoup moins sensibles '
        r'au volume que l’ordre global.',
        r'\item \textbf{Ce que le volume n’améliore pas ici.} Les passages '
        r'entre populations deviennent trop lents pour le budget disponible. '
        r'L12 conserve un désaccord important des chaînes malgré la grille '
        r'densifiée et la production doublée. Le maximum de 1052.5 K '
        r'reste provisoire. Le réseau rigide, les paires fixes, l’absence '
        r'de phonons et la dynamique d’échanges globaux restent identiques : '
        r'aucune de ces omissions n’est corrigée en augmentant L.',
        r'\item \textbf{Ce qu’il faut ensuite mesurer.} Allonger les séries '
        r'critiques et vérifier leurs poids de population ; comparer '
        r'des tailles à diagnostics satisfaits ; répéter les cycles '
        r'L6/L9/L12 achevés avec davantage de graines et des paliers '
        r'critiques plus longs pour tester la mémoire du parcours ; répéter '
        r'plusieurs trempes avec une dynamique locale si l’objectif '
        r'devient la nucléation physique. Une différence entre deux '
        r'courbes seule ne permet pas de choisir entre ces mécanismes.',
        r'\end{itemize}',
        r'{\raggedright Les fichiers de \texttt{results/comparison\_L6\_L9\_L12} '
        r'contiennent : \texttt{overview.csv} (budgets et diagnostics), '
        r'\texttt{all\_observables\_cross.csv} (moments, covariances, '
        r'autocorrélations et contrôle exact énergie--Cowley à toutes '
        r'les températures), \texttt{chain\_variant\_statistics.csv} '
        r'(occupations et populations par chaîne), '
        r'\texttt{block\_sensitivity.csv} (stabilité des erreurs), '
        r'\texttt{spatial\_comparison.csv} (modes de Fourier et limites). '
        r'Les neuf observables instantanées, $X_q$, le cumulant, les '
        r'identités et acceptations de répliques, les erreurs, les '
        r'demi-séries et les captures spatiales ont ainsi été examinés. '
        r'Les historiques bruts d’équilibration gardent leur rôle de '
        r'diagnostic, sans être incorporés aux moyennes de production.\par}',
        r'\begin{oral}{1. Quelle conclusion défendre à partir de toutes les données ?}',
        r'\begin{itemize}',
        r'\item \textbf{Dire ce qui est robuste et ce qui reste conditionnel.} '
        r'La structure froide et son énergie sont vérifiées analytiquement. '
        r'Au chaud, diminuer le biais de $q$ en augmentant N et conserver '
        r'un Cowley non nul sont deux observations compatibles : l’ordre '
        r'global disparaît tandis que des préférences de voisinage survivent. '
        r'Près du pic, la cohérence entre ordres et énergie est forte, '
        r'mais le poids des deux populations n’est pas assez stable dans '
        r'les grandes boîtes. On rapporte donc des maxima observés, '
        r'pas une température expérimentale ni une limite infinie précise.',
        r'\item \textbf{Ne pas transformer une amélioration géométrique en '
        r'validation complète.} L12 laisse davantage de place aux fluctuations '
        r'et résout mieux Fourier, mais il ralentit aussi l’exploration. '
        r'La question suivante n’est pas seulement « encore plus grand ? » : '
        r'c’est « davantage de passages indépendants, avec quel protocole '
        r'et quels diagnostics ? ». La réponse guide les calculs utiles.',
        r'\end{itemize}\end{oral}']
    from .cycle_size_report import completed_cycles
    completed=[r['L'] for r in completed_cycles()]
    if len(completed)>1:
        for i,line in enumerate(lines):
            if line.startswith('Les archives achevées de chauffage/refroidissement'):
                lines[i]=(r'La figure précédente garde les trois productions canoniques '
                    r'et le cycle L6 comme références de protocole. Les cycles '
                    +', '.join(f'L{L}' for L in completed)+
                    r' sont maintenant achevés ; leurs branches de chauffage et '
                    r'de refroidissement sont superposées en '
                    r'section~\ref{sec:cycles-tailles}, avec leurs diagnostics propres. '
                    r'Les fluctuations d’un palier ne peuvent être interprétées '
                    r'comme une capacité calorifique d’équilibre que si '
                    r'la stationnarité et l’exploration sont suffisantes.')
    text='\n\n'.join(lines)+'\n'
    (out/'bonus_sizes.tex').write_text(text,encoding='utf8')
    _copy_atomic(out/'bonus_sizes.tex',target/'bonus_sizes.tex')
    for name in figures:
        _copy_atomic(out/f'{name}.pdf',target/'figures'/f'{name}.pdf')
    return target/'bonus_sizes.tex'


def refresh_include():
    """Régénère le texte depuis l’analyse existante, sans refaire les FFT."""
    from .cross_analysis import load_campaigns
    out=ROOT/'results/comparison_L6_L9_L12'
    summary=json.loads((out/'summary.json').read_text(encoding='utf8'))
    campaigns=load_campaigns()
    cross=read_rows(out/'all_observables_cross.csv')
    for run in campaigns:
        run['cross']=[{k:float(v) for k,v in r.items()} for r in cross if int(r['L'])==run['L']]
    spatial=[{k:(v if k=='flags' else float(v)) for k,v in r.items()}
             for r in read_rows(out/'spatial_comparison.csv')]
    return write_include(campaigns,summary['sizes'],spatial,summary['figures'],out)
