# Analyse croisée L6/L9/L12 — 9 novembre 2026

Le rapport principal reste `report/rapport_template.pdf`. Les questions 1–9
conservent L6 comme référence et la partie 10.3 intègre désormais les trois
productions terminées. La partie 10.4 distingue explicitement les cycles
séquentiels encore en calcul. Aucun résultat manquant n'est interpolé ou inventé.

## Résultats et interprétation

- Maxima observés de la chaleur configuratoire : 1050, 1055 et 1052,5 K.
  L12 possède 90 températures et deux fois le budget de balayages ; la
  comparaison n'isole donc pas parfaitement le seul effet de taille.
- À 1045 K, q augmente de 0,652 à 0,784 puis 0,795, tandis qu'à 1000 K les
  amplitudes restent proches de 0,88–0,89. L'augmentation d'ordre avec la
  taille dépend du régime ; elle n'est pas uniforme à toute température.
- À 4000 K, q√N vaut 1,922, 1,924 et 1,926 : les petits résidus d'ordre sont
  compatibles avec une norme de fluctuations finies. Cowley reste non nul,
  montrant que ce résultat ne signifie pas indépendance chimique complète.
- Le minimum d'acceptation des échanges passe d'environ 0,37 % en L9 à 29 %
  en L12. Pourtant aucun aller-retour complet n'est enregistré dans ces deux
  campagnes ; la densification améliore le chevauchement, sans résoudre
  l'exploration des populations critiques.
- Au maximum L12, les fractions descriptives q > 0,5 diffèrent de 66,6 points
  entre les chaînes. R-hat atteint environ 1,268 et l'ESS minimal vaut 44,7.
  Les distributions sont affichées séparément. L12 reste provisoire au pic.
- L'énergie est reconstruite depuis Cowley sur toutes les mesures, avec
  un résidu maximal inférieur à 8e-12 eV. La décomposition de cV conserve
  les deux variances et leur covariance. Les accords ne sont pas présentés
  comme des preuves statistiques indépendantes.
- Les captures spatiales L9 et L12 sont exploitées à 14 et 15 températures.
  Les rapports de Fourier connexes, leurs erreurs par blocs et leurs réserves
  sont comparés à la campagne spatiale dédiée L6. La grande longueur froide
  n'est pas interprétée comme un rayon de germe ; les mélanges de variantes,
  les composantes lentes et les rapports indéfinis restent signalés.

## Figures et données supplémentaires

Huit nouvelles plaquettes PDF/PNG sont dans
`results/comparison_L6_L9_L12` : transition, ordre et cumulant, diagnostics,
traces/distributions par chaîne, décomposition énergétique, comparaison des
protocoles, variantes/populations et longueurs spatiales. Elles remplacent
les trois reproductions de l'ancienne synthèse L6/L9.

Les CSV sauvegardent aussi les autocorrélations de toutes les observables,
les occupations par chaîne, les variantes conditionnelles, les acceptations
locales, le terme de variance entre moyennes de chaînes (exclu de cV), les
poids des populations et la sensibilité des erreurs à la longueur des blocs.

Les encadrés des neuf questions ont été développés en listes explicatives,
avec des exemples numériques et des distinctions entre code exact,
échantillonnage et approximation physique. Les encadrés des deux premiers
bonus reprennent les mêmes conventions. Le corps reste en 10 points.

Deux défauts de texte généré ont été corrigés : un échappement Python
transformait `\\texttt` en tabulation, et la colonne de réserve spatiale ne
détaillait pas les causes multiples. Les listes utilisent les packages
LaTeX déjà installés, sans nouvelle dépendance.

## Cycles additionnels

Après la demande de restriction à L9, `run_size_cycles.py` lance **L9 seul**
par défaut, avec les mêmes 46 températures,
deux graines, 6000 balayages écartés et 30000 mesurés par palier que L6.
Chaque palier achevé est sauvegardé. Relancer ce script reprend les paliers
valides. Le suivi est `results/comparison_L6_L9_L12/cycles_execution.json`
et les suivis détaillés sont dans `results/thermal_cycle_L9/execution.json`
Le cycle L12 est différé et ne sera pas lancé automatiquement. Sa configuration
reste disponible pour un éventuel lancement ultérieur explicite.

Le premier processus a disparu après 33 paliers sauvegardés. Les journaux
ont été conservés dans `cycle_launches/20261109_191531`, puis L9 a été
repris avec `run_size_cycles.py --sizes 9`, sans refaire les paliers achevés.

Après chaque taille complète : validation intégrée au cycle, figures dédiées,
superpositions des branches disponibles et recompilation automatique du PDF
principal. Les deux branches par taille seront comparées en symboles creux,
avec une vue de transition et une vue élargie. Les données et le rapport L6
restent présents pendant ces nouveaux calculs.

## Vérification

34 tests passent, dont trois tests supplémentaires sur l'analyse croisée.
Les contrôles numériques des budgets de covariance et de synchronisation
des captures avec les observables passent. L'audit final des figures et de
la compilation est sauvegardé dans `report_audit.json`, à côté des CSV.

## Cycle L9 achevé et analyse finale

Les 184 paliers L9 sont terminés et validés. Le lanceur s'est achevé le
9 novembre à 21 h 04, après sa première publication automatique du PDF.
Le cycle L12 reste différé, conformément à la demande de l'utilisateur.

La mise à jour suivante développe l'analyse : superpositions L6/L9 des
deux branches, écarts appariés et diagnostics, comparaison à la production
canonique L9, occupations séparées des sous-réseaux L9. Les huit figures
propres au cycle L9 sont régénérées. Les valeurs de référence de L6 et les
comparaisons canoniques L6/L9/L12 sont conservées.

Les deux maxima de branche L9 sont à 1055 K. À 1050 K, les écarts de q
des deux trajectoires ont des signes opposés : l'écart moyen +0,0796 a
une dispersion entre cycles de 0,1538. Il n'établit pas une hystérésis
reproductible. R-hat atteint 1,29 au chauffage à 1050 K ; au pic, les
effectifs utiles des deux branches sont inférieurs à 15. Le défaut de
diagnostic à 500 K relève des excitations rares, avec ordre presque parfait,
et n'est pas interprété comme une autre transition. La seconde trajectoire
passe de la variante 2 à la variante 3 au retour, tout en retrouvant q=1.

Les données de confrontation des trois protocoles L9 à 1050 et 1055 K
sont exportées dans `cycle_L9_protocol_comparison.csv`.
