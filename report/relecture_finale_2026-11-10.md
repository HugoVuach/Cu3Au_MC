# Relecture finale du rapport — 10 novembre 2026

Cette revue porte sur les questions 1 à 9, les quatre parties du bonus,
les inclusions générées, les algorithmes et leur pseudocode. Les neuf
réponses du sujet sont présentes. Les archives de simulation ne sont pas
modifiées et aucun nouveau calcul Monte-Carlo n'est lancé.

## Corrections et clarifications

- La couverture nomme explicitement la production canonique L6 à laquelle
  se rapportent les valeurs principales. Un guide distingue les quatre
  parties du bonus et leurs protocoles.
- Le vocabulaire distingue chaîne, réplique, domaine et sous-réseau,
  ainsi que la structure L1₂ et la taille L = 12. μ est un indice géométrique,
  pas un potentiel chimique. Les longueurs sont exprimées en unités de a.
- Dans la variation locale d'énergie, V désigne les listes de voisins et
  U les énergies de paire, pour éviter de donner deux sens au même symbole.
  β et le nombre total de paires M₁ sont explicités.
- La nullité de Cowley n'est plus présentée comme une preuve de mélange
  entièrement indépendant ; un faible q ne suffit pas non plus à exclure
  des régions localement ordonnées.
- Les définitions de SE, ESS, R-hat et du bootstrap précèdent les seuils de
  diagnostic. Les unités de Xq et des rapports sans dimension sont précisées.
  La référence gaussienne 4/9 du cumulant d'une norme à trois composantes
  explique sa limite chaude sans l'imposer aux données.
- Les rapports de longueur sont écrits ξ/(La), pour séparer la longueur
  physique du nombre sans dimension de mailles L.
- La légende des cycles nomme désormais les trois tailles effectivement
  représentées. La fenêtre demandée 1000–1100 K est conservée. La seconde
  exportation identique reste absente du corps du rapport.

## Point scientifique principal complété : le cycle L12

La précédente publication automatique superposait bien les trois tailles,
mais développait surtout l'analyse du cycle L9. La version finale ajoute
la confrontation des trois protocoles L12, leurs traces et distributions
séparées et une synthèse critique.

À 1050 K, q moyen vaut environ 0,772 au chauffage et 0,302 au retour ; les
écarts appariés des deux trajectoires sont négatifs. L'écart moyen est
−0,4695, avec une erreur de dispersion entre deux trajectoires de 0,1781.
C'est une mémoire de parcours dans les données disponibles, sans intervalle
de confiance à 95 % établi ni loi universelle de taille.

Le maximum de variance L12 se situe à 1055 K au chauffage et à 1050 K au
retour, contre 1052,5 K dans la production à répliques. Le milieu des
maxima des branches n'est pas présenté comme un nouvel estimateur validé.
Les paliers critiques ont R-hat proche de 1,44–1,45 et seulement quatre
blocs dans la chaîne la plus lente. Inversement, deux branches peuvent
passer des contrôles locaux en restant dans des populations opposées.
Un faible nombre de températures signalées ne prouve donc pas un meilleur
équilibre global.

Deux plaquettes supplémentaires sont intégrées au PDF : comparaison des
protocoles L12 et traces/distributions par trajectoire. Leurs chiffres sont
exportés dans `cycle_L12_population_statistics.csv`. Les GIF et les longueurs
spatiales restent attribués aux campagnes qui les ont réellement produits.

## Références vérifiées

- Paras et Allanore (2021), doi:10.1103/PhysRevResearch.3.023239 : référence
  de transition à 662 K et caractère du premier ordre du matériau réel.
- Benisek et Dachs (2015), doi:10.1016/j.jallcom.2014.12.215 : contributions
  calorimétriques et dépendance au protocole (pic à 675 K pour l'échantillon
  le plus ordonné, réaction vers 660 K pour les échantillons initialement
  désordonnés). Cette nuance complète la comparaison au modèle.
- Kessler et al., arXiv:cond-mat/0207611 : dynamique atome–lacune et
  scénarios dépendant de la profondeur de trempe ; aucune calibration de
  temps physique n'est transposée à nos échanges globaux.
- Edwards et al., arXiv:hep-lat/9112002, équations 3.14–3.18 : définition
  spectrale de second moment et distinction avec la longueur exponentielle.
- Reichert et al., EPL 53, 570–576 (2001) : trois composantes d'ondes
  d'ordre et concentration ; étude de surface distinguée de notre boîte
  périodique en volume.

## Vérification et limites

Les 34 tests sont relancés avec succès le 10 novembre. Le journal est
`results/comparison_L6_L9_L12/review_tests_20261110.log`. Les vérifications
finales de compilation, de figures, de renvois et de présence des GIF sont
enregistrées dans `results/comparison_L6_L9_L12/report_audit.json`.

La relecture clarifie la portée des résultats ; elle ne supprime pas
l'échantillonnage insuffisant au pic, ne donne pas de cinétique physique
aux balayages et ne convertit pas le modèle rigide en prédiction quantitative
de la température expérimentale. Ces limites restent affichées dans le PDF.
