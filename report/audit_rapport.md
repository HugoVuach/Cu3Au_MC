# Audit du rapport Cu₃Au par rapport à l’énoncé

Audit du 28 octobre 2026, fondé sur la nouvelle version de `subject.md`, qui contient les neuf questions originales.

## Conclusion

**Les neuf questions sont traitées. Je n’ai pas trouvé de question entièrement omise, ni d’erreur dans les valeurs analytiques principales vérifiées.** Le réseau, les paramètres de paire, la composition, les deux coquilles, Metropolis, les fonctions `cowley` et `ordre`, les résultats thermiques et la discussion correspondent aux demandes du sujet.

Deux réserves empêchent cependant de dire que tout est définitivement validé :

- **Question 4 :** les archives contiennent un compte rendu de 15 tests réussis, mais leur relance dans l’environnement Python actuel échoue avant l’exécution des tests physiques, car Numba manque. Une dépendance Polars superflue, non déclarée, compromet également une installation propre du projet.
- **Question 8 :** le résultat obtenu est un maximum de capacité calorifique à **1050 K pour la boîte de 864 sites**, sur une grille discrète. Il répond à la demande d’une estimation numérique de la transition, mais ne détermine pas une température critique à la limite thermodynamique. Le rapport expose déjà cette limite correctement.

L’amélioration prioritaire est pédagogique : montrer visuellement les coquilles et les variantes, expliciter les comptages et les unités, puis guider la lecture des résultats. Ajouter davantage de théorie avancée n’est pas indispensable pour répondre au sujet.

## Périmètre et vérifications effectuées

Documents examinés :

- `subject.md`, en particulier les questions 1 à 9, lignes 67 à 101 ;
- `report/rapport_template.tex` et son PDF actuel de **33 pages** ;
- les fichiers injectés `auto_results.tex`, `critical_table.tex`, `discussion_results.tex` ;
- les sources de géométrie, énergie, noyau Monte-Carlo, statistiques et figures ;
- les deux archives de production, `observables.csv`, `summary.json`, les journaux historiques de tests et de validation ;
- les figures, notamment les sous-réseaux, les configurations et les distributions.

Les archives contiennent un tableau de forme `(2, 3000, 46, 9)` : deux chaînes, 3000 mesures de production, 46 températures et neuf observables instantanées. J’ai recalculé les moyennes et contrôlé leur correspondance avec le CSV :

| Contrôle sur les données existantes | Résultat |
|---|---|
| Énergie moyenne par atome : archives contre CSV | Écart maximal inférieur à $10^{-13}$ eV/atome |
| Moyennes de η, q, α₁ et α₂ : archives contre CSV | Écarts maximaux inférieurs à $4\times10^{-14}$ |
| Énergie de chaque mesure reconstruite par les comptages de paires déduits des α | Écart maximal inférieur à $10^{-12}$ eV |
| Somme des quatre occupations de sous-réseaux | Accord avec 1 à l’arrondi flottant près |

Ces contrôles confirment la cohérence des résultats archivés. Ils ne remplacent ni l’exécution complète des tests du programme actuel, ni une preuve d’équilibration parfaite. Aucune nouvelle campagne Monte-Carlo n’a été exécutée pour cet audit.

## Vérification question par question

Les lignes indiquées ci-dessous concernent `report/rapport_template.tex` dans son état au moment de l’audit.

| Question du sujet | Réponse présente | Évaluation | Ajout conseillé |
|---|---|---|---|
| **1. Écrire le programme Monte-Carlo : réseau, voisins, énergie, échanges et Metropolis.** | Section 1, ligne 91 ; géométrie entière, PBC, 12 + 6 voisins, énergie de paires, ΔE local, balance détaillée ; annexe Algorithme, ligne 676, et pseudocode, ligne 1258. | **Couverte.** Les coordonnées modulo 12 en unités de a/2 sont équivalentes aux PBC modulo 6 mailles demandées. Les échanges conservent les 216 Au et 648 Cu. | Un dessin des deux coquilles et un organigramme court reliant proposition, décision, balayage et mesure. |
| **2. Écrire la fonction `cowley`.** | Section 2, ligne 197 ; fonction nommée `cowley` dans `model.py`, et version compilée dans `kernel.py`. Formule, comptage des paires, valeurs ordonnées et correction de taille finie. Figure dédiée `12_cowley`. | **Couverte.** Le comptage correspond à la probabilité conditionnelle demandée. | Un exemple chiffré explicitant la différence entre probabilité conditionnelle et paramètre α. |
| **3. Écrire la fonction `ordre` avec le paramètre η défini sur les sommets.** | Section 3, ligne 250 ; fonction `ordre` dans `model.py`. Le η du sujet est conservé ; q est ajouté pour reconnaître les quatre variantes. Figure `09_sous_reseaux`. | **Couverte, mais la figure peut tromper un débutant.** Les occupations moyennes proches de 1/4 à basse température ne signifient pas nécessairement un désordre instantané. | Quatre petites configurations ordonnées et un exemple distinguant moyenne des occupations et moyenne de q. |
| **4. Tester le programme complet.** | Section 4, ligne 284 ; 15 tests historiques, défaut d’antisites, comparaison ΔE local/global, contrôles de composition et d’énergie par segment, validation des archives. | **Documentée et étayée historiquement ; relance locale incomplète.** Les modules de tests ne se chargent pas dans l’environnement actuel faute de Numba. | Un tableau « test / valeur attendue / résultat / environnement / date », et une relance après remise en état de l’installation. |
| **5. Déterminer l’énergie de cohésion à T = 0 depuis L1₂.** | Section 5, ligne 331 ; comptage des paires, E₀, u₀, énergie par formule Cu₃Au et preuve du minimum global. | **Couverte et démonstration solide.** Valeurs vérifiées : −3294,432 eV pour la boîte, −3,813 eV/atome et −15,252 eV par unité Cu₃Au. | Montrer les deux divisions et écrire explicitement « par unité de formule » ; tableau des paires dans les deux coquilles. |
| **6. Générer des configurations d’équilibre à plusieurs températures.** | Section 6, ligne 371 ; équilibration/production, protocole, départs, échanges de répliques, configurations finales et traces de circulation. | **Couverte avec diagnostics d’équilibre.** Les diagnostics favorables rendent les résultats plausibles sans garantir l’exploration de tous les états. | Des coupes 2D lisibles des configurations et une comparaison des moyennes cumulées des deux chaînes. |
| **7. Calculer les grandeurs physiques en fonction de T.** | Section 7, ligne 436 ; énergie, η, q, α₁, α₂, cᵥ, fluctuations d’ordre, erreurs et limites de haute température. | **Couverte, au-delà du minimum explicitement demandé.** Les unités de cᵥ sont cohérentes avec les énergies totales utilisées. | Expliquer d’abord la chaîne « configurations → énergies mesurées → variance → cᵥ » ; ajouter un tableau de quelques températures représentatives. |
| **8. Déterminer la température de transition ordre–désordre.** | Section 8, ligne 526 ; maximum à 1050 K, voisins de grille 1045 et 1055 K, accord des chaînes, distributions, zoom et comparaison au champ moyen. | **Couverte comme estimation à taille finie.** L’intervalle des voisins n’est justement pas présenté comme un intervalle de confiance. | Quantifier l’incertitude statistique sur la position du maximum, et résumer la concordance des indices de transition. Plusieurs tailles seraient un approfondissement. |
| **9. Discuter les résultats.** | Section 9, ligne 614 ; résultats chiffrés injectés, compétition énergie–entropie, ordre local/global, limites du réseau rigide et de taille finie. | **Couverte.** La discussion distingue correctement modèle et matériau réel. | Une conclusion courte répondant directement : état froid, zone de transition, état chaud, portée de la valeur 1050 K. |

Les échanges de répliques, q, le cumulant de Binder et le champ moyen sont des compléments utiles. **Ils ne sont pas des exigences explicites de l’énoncé et ne doivent pas masquer les fonctions `cowley`, `ordre` et le Metropolis de base.** Les deux fonctions demandées existent bien ; le noyau de production emploie des calculs compilés équivalents.

## Points à corriger ou à clarifier en priorité

### 1. Rendre la validation du programme reproductible

`results/production/tests.log` indique historiquement « Ran 15 tests ... OK ». Cette preuve doit être conservée avec son contexte d’exécution, sans la confondre avec un nouveau test sur le poste actuel.

La tentative `python -m unittest discover -s tests -v` effectuée pendant la préparation de l’audit a échoué au chargement des deux modules : `ModuleNotFoundError: No module named 'numba'`. Il s’agit de deux erreurs de chargement, **pas de deux tests physiques ayant donné un résultat faux**. Le Python du dossier `.venv` ne dispose pas non plus des dépendances scientifiques nécessaires.

De plus, `cu3au/model.py:2` importe `other` depuis un module privé de Polars. Cet import est inutilisé et Polars ne figure pas dans les dépendances déclarées. Une installation neuve peut donc échouer même après installation des bibliothèques scientifiques prévues.

**Proposition :** supprimer cet import superflu, préparer un environnement compatible avec les versions déclarées, puis relancer les 15 tests et la validation des archives. Vérifier aussi l’annonce « Python 3.11 ou ultérieur » avec les versions de dépendances effectivement installables. L’audit n’a ni installé de dépendances ni modifié le code de simulation.

### 2. Identifier sans ambiguïté le rapport à remettre

Le rapport récent est `report/rapport_template.pdf`, avec 33 pages. `results/production/report/rapport.pdf` contient encore 17 pages et ne reprend pas les derniers ajouts. Le README dirige toujours le lecteur vers ce second fichier.

**Proposition :** désigner un rapport final unique et synchroniser les instructions de compilation. Les mentions de « onze figures » et de « 17 pages » dans les documents de livraison sont historiques ; la livraison actuelle compte douze figures et un rapport enrichi. Ne pas écraser silencieusement les traces historiques de validation : dater les nouvelles informations.

### 3. Préciser la fréquence des échanges entre répliques

La section 6 parle d’une tentative sur chaque lien concerné tous les 5 balayages, tandis que l’annexe décrit l’alternance des deux parités. La formulation mérite d’être rendue explicite :

> Une session d’échange a lieu tous les 5 balayages. Les liens pairs et impairs sont alternés ; un lien donné est donc testé tous les 10 balayages.

Cela évite de confondre fréquence des sessions et fréquence des tentatives sur un lien particulier.

### 4. Mieux justifier les incertitudes de certaines observables

La longueur des blocs est choisie à partir des temps d’autocorrélation de E et q. Cela ne démontre pas automatiquement que ce choix couvre la corrélation temporelle de η, des occupations et des deux α. Les changements de variante peuvent notamment être lents même lorsque E et q changent peu.

**Proposition :** mesurer aussi les temps de corrélation des α et du sous-réseau dominant, ou montrer la stabilité des erreurs lorsque la longueur des blocs augmente. Ce point demande une analyse supplémentaire des archives ; il ne prouve pas que les erreurs publiées sont fausses.

À 300 K, E et q sont constants dans les mesures : le code renvoie alors τ = 1/2, ESS = 6000 et R-hat = 1. Il faut expliquer que ces nombres reflètent ici le traitement des séries constantes et **ne mesurent pas une exploration complète de toutes les excitations ou variantes**. Une erreur nulle sur cᵥ dans une série sans fluctuations observées n’est pas une preuve que la capacité calorifique physique est exactement nulle à cette température finie.

## Ajouts de figures recommandés

| Priorité | Proposition | Ce que le lecteur comprendrait | Données ou travail nécessaires |
|---|---|---|---|
| **P1** | **Schéma des coquilles autour d’un Au**, avec les 12 premiers voisins et les 6 seconds, distances et couleurs distinctes. | Pourquoi il existe deux α et deux familles d’interactions ; comment compter les voisins. | Géométrie déjà disponible ; aucune nouvelle simulation. |
| **P1** | **Les quatre variantes L1₂**, puis une configuration réellement désordonnée, avec un tableau des pμ, η et q. | Pourquoi η peut être négatif dans un état parfaitement ordonné et pourquoi une moyenne peut annuler η. | Illustration construite à partir des états de référence existants. |
| **P1** | **Compléter `09_sous_reseaux` par une lecture instantanée** : pμ à un instant froid et fréquence de visite des variantes, ou trace de `argmax(pμ)`. | La différence entre « 25 % en moyenne » et « mélange à chaque instant ». | Les pμ enregistrés suffisent ; aucun nouveau calcul Monte-Carlo. |
| **P1** | **Coupes 2D des configurations**, à froid, près de 1050 K et à chaud. | La disparition de la périodicité et l’apparition d’antisites plus facilement que dans un nuage 3D dense. | Configurations finales déjà archivées ; préciser qu’il s’agit d’instantanés. |
| **P2** | **Schéma d’un échange Au–Cu**, avant/après, avec les liaisons modifiées et la liaison mutuelle inchangée. | D’où vient ΔE et pourquoi on exclut le partenaire dans les sommes locales. | Illustration du défaut dont l’énergie est déjà calculée. |
| **P2** | **Moyennes cumulées de u et q**, séparées par chaîne, près du maximum. | Pourquoi on écarte une équilibration et comment vérifier que les moyennes se stabilisent. | Séries existantes ; comparer également plusieurs fenêtres de production. |
| **P2** | **Zoom conjoint sur q, α₁ et α₂** autour de la transition. | L’ordre local subsiste lorsque l’ordre global devient faible. | CSV existant. Garder masqués dans la figure Cowley les repères que l’utilisateur a demandé de masquer. |
| **P2** | **Histogrammes par chaîne** autour du maximum, en complément des histogrammes regroupés. | Vérifier que les deux bosses ne proviennent pas seulement de trajectoires explorant des régions différentes. | Archives existantes ; classes communes et incertitudes par blocs si nécessaire. |
| **P3** | **Comparaison de plusieurs tailles L** : q(T), cᵥ(T), déplacement du maximum et cumulants. | Distinguer le résultat de la boîte étudiée de la limite thermodynamique. | Nouvelles campagnes ; approfondissement, non condition minimale de conformité au sujet. |

Il n’est pas nécessaire d’ajouter tous ces dessins. Les quatre propositions P1 apporteraient le plus de clarté. Le graphique des configurations 3D a aussi un défaut de mise en page visible : son texte inférieur se superpose aux axes. Il faut réserver davantage de marge avant de l’utiliser dans des diapositives.

Les figures `12_cowley`, `03_thermodynamique`, `06_traces`, `11_zoom_transition` et `04_distributions` couvrent déjà les résultats principaux. La priorité est de les expliquer plutôt que de multiplier les courbes similaires.

## Démonstrations et exemples pédagogiques à ajouter

### Question 2 : passer d’un comptage à une probabilité, puis à α

Ajouter un exemple encadré :

\[
N_{\mathrm{Au}}=216,\quad z_1=12,\quad B_1=1944
\quad\Longrightarrow\quad
p_1^{\mathrm{Cu|Au}}=\frac{1944}{216\times12}=0.75
\quad\Longrightarrow\quad\alpha_1=0.
\]

Puis comparer avec l’ordre parfait : $B_1=2592$, donc $p_1=1$ et α₁ = −1/3. Insister sur les deux normalisations successives : nombre de voisins examinés, puis concentration globale. Pour α₂, montrer que les six seconds voisins d’un Au sont tous Au dans L1₂, d’où $B_2=0$ et α₂ = 1.

### Question 3 : montrer explicitement l’effet de la moyenne

Les quatre configurations parfaitement ordonnées ont les occupations `(1,0,0,0)`, `(0,1,0,0)`, `(0,0,1,0)` et `(0,0,0,1)`. Leur moyenne uniforme vaut `(1/4,1/4,1/4,1/4)` alors que chacune a q = 1. Présenter côte à côte :

\[
\left\langle q(p)\right\rangle=1,
\qquad q\!\left(\langle p\rangle\right)=0.
\]

Le rapport contient déjà l’argument. Un dessin et ce calcul très court le rendraient immédiatement compréhensible. Préciser que la figure actuelle moyenne aussi sur les deux chaînes.

### Question 5 : expliciter le changement d’unité

Afficher les étapes qui sont actuellement condensées :

\[
u_0=\frac{E_0}{864}=-3.813\ \mathrm{eV/atome},\qquad
E_0^{\mathrm{Cu_3Au}}=4u_0=\frac{E_0}{216}
=-15.252\ \text{eV/unité de formule}.
\]

Expliquer qu’une unité Cu₃Au contient quatre atomes. Conserver la distinction entre énergie négative du modèle et énergie positive nécessaire à la dissociation, sans assimiler automatiquement ce modèle de paires à une mesure expérimentale de cohésion.

### Question 7 : introduire cᵥ avant la formule de fluctuations

Introduire d’abord sa signification physique et ses unités :

\[
C_V=\left(\frac{\partial\langle E\rangle}{\partial T}\right)_V,
\qquad c_V=C_V/N.
\]

Puis montrer en quelques lignes comment le poids canonique conduit à la variance de E. Le rapport donne déjà la dérivation essentielle ; il lui manque surtout une entrée pédagogique et le lien explicite avec la série mesurée. Expliquer que Metropolis produit les configurations, tandis que le calcul de cᵥ utilise leurs énergies. La proportion de mouvements acceptés n’est pas la capacité calorifique.

### Question 8 : justifier le diagnostic de transition et sa précision

Ajouter un tableau récapitulatif associant la baisse de q, le maximum de cᵥ, l’évolution de u et les distributions. Les valeurs déjà archivées sont :

| T (K) | q | α₁ | α₂ | cᵥ/kB |
|---:|---:|---:|---:|---:|
| 300 | 1,000 | −0,333 | 1,000 | 0,000 observé |
| 1000 | 0,881 | −0,283 | 0,803 | 1,041 |
| 1050 | 0,575 | −0,233 | 0,527 | 5,873 |
| 1100 | 0,192 | −0,180 | 0,242 | 0,694 |
| 4000 | 0,065 | −0,067 | 0,026 | 0,028 |

Au maximum, l’erreur-type publiée sur cᵥ/kB est d’environ 0,102. **Cette erreur sur la hauteur du pic n’est pas une erreur en kelvins sur sa position.** Pour estimer cette dernière, on peut rééchantillonner les archives par blocs, en conservant les dépendances entre températures d’une même chaîne, et rechercher le maximum dans chaque tirage. Il faudra toujours distinguer cette incertitude statistique, la résolution de grille et le biais de taille finie.

Les échanges de répliques ne totalisent que 4 puis 3 allers-retours complets froid–chaud–froid en production. Ce résultat mérite d’être cité près du diagnostic d’équilibre : il ne réfute pas les résultats, mais invite à ne pas transformer des R-hat favorables en garantie absolue. Les histogrammes montrent une structure compatible avec une coexistence ; le rapport a raison de ne pas en faire une preuve définitive de premier ordre à partir d’une seule taille.

## Organisation pédagogique proposée

1. Ajouter un petit glossaire au début : site, maille, sous-réseau μ, coquille, configuration, réplique, chaîne, balayage, équilibration et production.
2. Garder dans chaque question une réponse directe, un exemple et la figure la plus utile. Placer les détails d’implémentation dans l’annexe déjà écrite.
3. Présenter les résultats dans l’ordre : structure froide → occupations et ordre → évolution avec T → capacité calorifique → estimation de transition → limites.
4. Conserver le champ moyen comme comparaison complémentaire clairement identifiée : il n’est pas requis explicitement par le sujet. Son calcul est développé ; si le lecteur le trouve trop dense, en déplacer les détails en annexe plutôt que les supprimer.
5. Terminer par une conclusion de quelques phrases donnant les valeurs clés et ce qu’elles permettent réellement d’affirmer.

La bibliographie du sujet et celle du rapport ne sont pas identiques. Ce n’est pas une question laissée sans réponse : les références de Metropolis et des échanges de répliques sont pertinentes pour les méthodes employées. Si une comparaison à la littérature ou à l’expérience est ajoutée, elle devra citer une source effectivement consultée et distinguer les résultats de ce modèle des données du matériau réel. Aucun nouveau résultat expérimental n’est supposé ou introduit dans cet audit.

## Ordre de travail recommandé

- **D’abord :** clarifier les unités de Q5, le sens des moyennes de Q3, les fréquences d’échange et la provenance historique des tests ; identifier le PDF final.
- **Ensuite :** ajouter les dessins des coquilles et des variantes, améliorer la figure 09 et les coupes des configurations ; insérer le glossaire et les exemples courts.
- **Puis :** remettre l’environnement de tests en état et confirmer la reproductibilité du code actuel ; examiner la stabilité des erreurs par blocs et du pic à partir des archives.
- **En approfondissement :** prolonger les simulations, ajouter des graines et comparer plusieurs tailles si l’objectif est une détermination thermodynamique plus précise.

Cet audit décrit l'état initial et les recommandations au moment de sa rédaction.

## Suivi de mise à jour : ajouts des neuf questions intégrés

Les ajouts de la colonne « Ajout conseillé » ont depuis été intégrés au rapport L=6 :

| Question | Mise à jour réalisée |
|---|---|
| 1 | Schéma des coquilles `13_coquilles` et organigramme `14_organigramme` ; glossaire. |
| 2 | Exemple chiffré distinguant comptage, probabilité conditionnelle et α. |
| 3 | Quatre configurations de référence dans `15_variantes`, tableau des pμ, η et q, et démonstration de l'effet de la moyenne. |
| 4 | Environnement réparé, 17 tests réussis, contrôles indépendants actuels et tableau daté dans le rapport. |
| 5 | Comptage des paires dans les deux coquilles et divisions explicites par 864 et 216. |
| 6 | Coupes dans deux plans `16_coupes_2D` et comparaison des moyennes cumulées `17_moyennes_cumulees`. |
| 7 | Définition physique de cᵥ, chaîne du calcul depuis Metropolis et tableau de températures représentatives. |
| 8 | Bootstrap synchronisé sur les températures : 1000 tirages pour des blocs de 18, 37 et 74 mesures ; intervalle percentile [1050,1055] K ; concordance des indicateurs et figure `18_incertitude_pic`. |
| 9 | Conclusion physique courte sur les états froid/chaud et la portée du maximum de 1050 K. |

La validation actuelle figure dans `results/production/verification_current.json` et `validation_current.json`. Les journaux historiques sont conservés. Aucune campagne L=9 n'a été relancée ; l'étude de plusieurs tailles reste un approfondissement distinct.
