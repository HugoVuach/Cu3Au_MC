from pathlib import Path
import re

folder = Path(__file__).resolve().parents[1]
path = folder / 'presentation_epuree.tex'
source = path.read_text(encoding='utf8')
# Title, subtitle and exactly two spoken-presentation prompts per slide.
slides = [
    ("", "",
     r"Nous simulons les échanges entre Au et Cu à composition fixée pour suivre la mise en ordre chimique.",
     r"Nous comparons les observables, les cycles thermiques et trois tailles de boîte : $L=6,9,12$."),
    ("Question 1 -- Construire le réseau et accepter les échanges atomiques",
     r"Hamiltonien de paires, variation locale d'énergie et règle de Metropolis",
     r"Le réseau CFC et les voisins restent fixes ; seuls les occupants Au et Cu sont échangés.",
     r"Metropolis accepte parfois une hausse d'énergie pour échantillonner les excitations thermiques sans changer la composition."),
    ("Question 2 -- Mesurer l'ordre chimique dans les deux premières coquilles",
     r"De la probabilité conditionnelle aux paramètres de Cowley $\alpha_1$ et $\alpha_2$",
     r"On compte les voisins Cu autour des Au, puis on compare cette probabilité à la concentration globale en Cu.",
     r"Un paramètre négatif indique un excès de voisins Cu ; l'état $L1_2$ parfait donne $\alpha_1=-1/3$ et $\alpha_2=1$."),
    ("Question 3 -- Distinguer l'ordre global et les quatre variantes",
     r"Paramètre signé $\eta$, occupations des sous-réseaux et amplitude invariante $q$",
     r"$\eta$ dépend du sous-réseau de référence, tandis que $q=1$ reconnaît les quatre variantes parfaitement ordonnées.",
     r"Une moyenne uniforme des occupations peut masquer des changements de variante : on calcule $q$ avant de le moyenner."),
    ("Question 4 -- Vérifier le programme avec des références exactes",
     r"Géométrie, énergie locale, observables et reprise des trajectoires : 34 tests réussis",
     r"Les calculs locaux sont confrontés à un recalcul global et à des états dont les valeurs sont connues exactement.",
     r"Ces tests vérifient le fonctionnement du code ; la qualité de l'équilibre se contrôle séparément sur les séries simulées."),
    ("Question 5 -- Calculer l'énergie de l'état ordonné à température nulle",
     r"Comptage des paires et passage de la boîte entière à l'atome ou à la formule Cu$_3$Au",
     r"On divise l'énergie totale par 864 pour obtenir l'énergie par atome, ou par 216 pour l'énergie par formule.",
     r"La positivité de $E-E_0$ montre que $L1_2$ est un minimum global pour le Hamiltonien et la composition considérés."),
    ("Question 6 -- Échantillonner les configurations à chaque température",
     r"Deux chaînes indépendantes et échanges de répliques pour la production canonique L6",
     r"Après l'équilibration, les mesures sont classées par température, après les éventuels échanges de configurations entre répliques.",
     r"L'accord des moyennes cumulées aide au diagnostic, mais ne prouve pas à lui seul que tous les états ont été explorés."),
    ("Question 7 -- Relier les configurations aux grandeurs thermiques",
     r"Énergie, ordre et réponses calculées à partir des fluctuations de production",
     r"La capacité calorifique utilise la variance de l'énergie totale ; Metropolis fournit les configurations et leurs poids thermiques.",
     r"Les variances sont calculées dans chaque chaîne avant moyenne pour éviter d'ajouter artificiellement leur désaccord."),
    ("Question 8 -- Repérer la transition et l'incertitude sur le maximum",
     r"Maximum de capacité calorifique et bootstrap par blocs sur la grille L6",
     r"Le maximum de $c_V$ se trouve à 1050 K en L6 ; les autres indicateurs situent la transformation dans la même zone.",
     r"L'intervalle bootstrap concerne la case du maximum sur cette grille et ne corrige ni la taille finie ni les états non visités."),
    ("Question 9 -- Comparer les résultats au champ moyen et au matériau réel",
     r"Ce que décrit le modèle de paires et ce que ses approximations limitent",
     r"Le modèle reproduit la perte d'ordre global et la persistance de préférences chimiques locales à haute température.",
     r"Le réseau rigide, les interactions simplifiées et l'absence de diffusion par lacunes limitent l'accord quantitatif avec le matériau réel."),
    ("Annexe A1 -- Suivre la mémoire d'un cycle de chauffage et refroidissement",
     r"Bonus 1 : paliers séquentiels sans échange de répliques, pour deux trajectoires L6",
     r"Chaque palier reprend l'état final du précédent, puis le refroidissement repart du véritable état chaud.",
     r"Un écart entre les branches révèle une mémoire à durée finie ; il ne définit pas automatiquement une hystérésis d'équilibre."),
    ("Annexe A2 -- Lire les occupations des Au sans confondre les variantes",
     r"Bonus 1 : séparation des deux cycles et des branches de chauffage/refroidissement",
     r"Chaque panneau suit une trajectoire et un sens du cycle ; couleurs et symboles gardent le même indice de sous-réseau.",
     r"Une bosse d'occupation peut traduire un changement de variante ; elle doit être confrontée à $q$ et aux autres occupations."),
    ("Annexe A3 -- Visualiser la mise en ordre sur neuf coupes atomiques",
     r"Bonus 1 : trois couches par orientation $(001)$, $(110)$ et $(111)$, avec accès aux GIF",
     r"Les neuf vues montrent des couches fixes du même volume, avec Au en or et Cu en bleu.",
     r"Cliquer sur les plaquettes ouvre la relaxation isotherme ; la formule bleue ouvre le cycle à température variable."),
    ("Annexe A4 -- Mesurer la mise en ordre après une trempe à 1000 K",
     r"Bonus 1--2 : ordre global, poids du mode uniforme et longueur diffuse pendant la relaxation",
     r"Cette trajectoire montre une mise en ordre rapide de la petite boîte, avec des captures tous les 20 balayages.",
     r"La longueur apparente est dominée par l'ordre uniforme ; une seule trempe ne donne ni rayon critique ni vitesse physique de croissance."),
    ("Annexe A5 -- Extraire une longueur de corrélation par Fourier",
     r"Bonus 2 : champ d'ordre signé, spectre connexe et estimateur de second moment",
     r"Le champ à trois composantes distingue les quatre variantes ; son mode uniforme mesure l'ordre, pas la concentration totale fixée.",
     r"Le rapport entre le centre du pic et les premiers modes voisins définit une longueur de second moment, distincte d'un rayon de domaine."),
    ("Annexe A6 -- Interpréter la longueur de corrélation en fonction de $T$",
     r"Bonus 2 : comparaison des spectres connexes, des corrélations réelles et du diagnostic intra-variante",
     r"Au froid, les changements collectifs de variante peuvent produire une grande longueur connexe, bien supérieure à la boîte.",
     r"Le diagnostic intra-variante et les contrôles statistiques empêchent d'assimiler cette valeur à une taille géométrique de domaine."),
    ("Annexe A7 -- Comparer l'ordre et la transition pour trois tailles",
     r"Bonus 3 : superposition L6/L9/L12 et amplitude résiduelle de l'ordre au chaud",
     r"À 1045 K, l'ordre est plus marqué dans les grandes boîtes ; au chaud, le résidu de $q$ diminue avec le volume.",
     r"La quasi-superposition de $\sqrt N\langle q\rangle$ soutient un effet de fluctuations finies, sans constituer une mesure d'exposant critique."),
    ("Annexe A8 -- Vérifier si une boîte plus grande améliore l'équilibre",
     r"Bonus 3 : effectif utile, désaccord des chaînes et circulation des répliques",
     r"Les grandes tailles rendent le pic plus marqué, mais les deux chaînes L12 ne visitent pas les populations avec les mêmes poids.",
     r"La grille plus fine et la production doublée ne suppriment pas ce défaut : augmenter $L$ ne garantit pas une meilleure convergence."),
    ("Annexe A9 -- Relier les défauts locaux à l'énergie et à la capacité",
     r"Bonus 3 : identité énergie--Cowley et rôle de la covariance entre coquilles",
     r"L'énergie et les paramètres de Cowley sont reliés par une identité exacte : leur accord n'est pas une preuve indépendante de transition.",
     r"La covariance contribue fortement au pic de capacité et doit être conservée, même si la seconde interaction est faible."),
    ("Annexe A10 -- Comparer la portée spatiale de l'ordre entre tailles",
     r"Bonus 3 : longueurs $\xi/a$, rapports $\xi/(La)$ et diagnostics des composantes signées",
     r"Au chaud, les longueurs courtes sont plus proches entre tailles ; près de la transition, les différences de populations pèsent davantage.",
     r"Les captures et budgets diffèrent entre campagnes, et les composantes signées restent lentes : la comparaison n'est pas une extrapolation contrôlée."),
    ("Annexe A11 -- Superposer les cycles thermiques des trois tailles",
     r"Bonus 4 : chauffage/refroidissement entre 1000 et 1100 K et écarts appariés à 1050 K",
     r"Les branches d'un cycle sont appariées ; les écarts sont calculés trajectoire par trajectoire avant leur moyenne.",
     r"Avec seulement deux cycles et des diagnostics critiques insuffisants, aucune loi monotone d'hystérésis avec la taille n'est établie."),
    ("Annexe A12 -- Expliquer la forte mémoire thermique observée en L12",
     r"Bonus 4 : confronter moyennes de branche, séries temporelles et populations individuelles",
     r"À 1050 K, le retour est moins ordonné que l'aller, mais les deux trajectoires de retour restent très différentes.",
     r"Un diagnostic local favorable peut correspondre à un seul bassin visité ; l'écart d'énergie mesuré n'est donc pas une chaleur latente."),
    ("Annexe A13 -- Organiser les mouvements, les balayages et les mesures",
     r"Du mouvement élémentaire aux deux protocoles de température et à l'analyse spatiale",
     r"Un balayage contient $N$ propositions, refus compris ; les mesures suivent un calendrier fixé après les échanges de répliques éventuels.",
     r"Le cycle transmet l'état au palier suivant, tandis que la production à répliques échantillonne plusieurs températures simultanément."),
    ("Annexe A14 -- Récapituler les paramètres physiques et numériques",
     r"Même Hamiltonien pour les trois tailles, avec des budgets canoniques et cycliques distincts",
     r"Passer de L6 à L12 multiplie le nombre de sites par huit, sans modifier la composition ni les énergies de paire.",
     r"Un même nombre de balayages ne fournit pas nécessairement le même nombre de configurations indépendantes."),
]

pattern = re.compile(r'(\\begin\{frame\}(?:\[[^\]]*\])?(?:\{[^\n]*\})?)(.*?)(\\end\{frame\})', re.S)
matches = list(pattern.finditer(source))
assert len(matches) == len(slides) == 24
assert r'\deuxpoints{' not in source
counter = 0

def enrich(match):
    global counter
    title, subtitle, first, second = slides[counter]
    body = match.group(2)
    if counter == 0:
        start = match.group(1)
        body = body.replace(r'{\fontsize{15}{18}\selectfont Monte-Carlo\par}',
                            r"{\fontsize{13}{16}\selectfont Simulation Monte-Carlo d'un alliage sur réseau rigide\par}")
        body += '\\vspace{7mm}\n'
    else:
        start = r'\begin{frame}{' + title + '}\n' + r'\framesubtitle{' + subtitle + '}'
    body = re.sub(r'\\(fpic|gifpic)\{([^}]+)\}\{(\d+)mm\}', resize, body)
    if counter == 6:
        body = body.replace('{53mm}', '{46mm}').replace('{44mm}', '{36mm}')
    counter += 1
    return start + body + r'\deuxpoints{' + first + '}{' + second + '}\n' + match.group(3)

def resize(match):
    height = int(match.group(3))
    if height >= 90:
        height = 83
    elif height >= 80:
        height = 73
    return '\\' + match.group(1) + '{' + match.group(2) + '}{' + str(height) + 'mm}'

source = pattern.sub(enrich, source)
anchor = r'\setbeamerfont{frametitle}{size=\fontsize{15}{17}\selectfont,series=\bfseries}'
source = source.replace(anchor, anchor + '\n' + r'\setbeamerfont{framesubtitle}{size=\fontsize{8.5}{10}\selectfont}'
                        + '\n' + r'\setbeamerfont{itemize/enumerate body}{size=\fontsize{9}{10.7}\selectfont}'
                        + '\n' + r'\setbeamertemplate{itemize item}{\textcolor{teal}{\textbullet}}')
source = source.replace(r'\begin{document}',
                        r'\newcommand{\deuxpoints}[2]{\par\smallskip\begin{itemize}\setlength{\itemsep}{2pt}\setlength{\parskip}{0pt}\item #1\item #2\end{itemize}}'
                        + '\n\n' + r'\begin{document}')
path.write_text(source, encoding='utf8')
print('Updated 24 slides with titles, subtitles and two bullet points each.')
