# Revue du rapport — 8 novembre 2026

Le rapport compilé passe de **55 à 44 pages**, avec un corps en 10 points,
des marges de 19 mm, des espaces réduits et des figures limitées en hauteur.
Les démonstrations, les commentaires des algorithmes et le pseudocode sont conservés.

## Audit des figures

Avant révision : 33 appels de figures, 33 noms distincts, aucune image strictement
identique après comparaison des PDF rendus. Les répétitions perçues étaient
principalement des représentations différentes des mêmes données.

Après révision : **26 figures distinctes**, aucune référence répétée. Sept
représentations sont retirées du corps du PDF, sans supprimer leurs fichiers :

| Export non reproduit | Vue ou renvoi conservé |
|---|---|
| `07_configurations` | `16_coupes_2D`, plus lisible pour les mêmes états canoniques |
| `02_ordres` | Cowley dédié, occupations et ordre global, puis zoom de transition |
| `03_thermodynamique` | Énergie et capacité déjà réunies dans `11_zoom_transition` |
| `bonus_05_traces` | Export dense conservé ; historique des états du GIF analysé séparément |
| `bonus_07_configurations` | Export à 1050 K conservé ; coupes multiples dans les GIF |
| `animation_isotherme_debut` | Plaquette fixe conservée à côté du GIF |
| `animation_isotherme_fin` | Plaquette fixe conservée à côté du GIF |

Les graphiques canoniques et les comparaisons chauffage/refroidissement sont
maintenus lorsqu'ils répondent à des protocoles ou questions différents.
Les zooms et vues globales explicitement réunis en panneaux ne sont pas
considérés comme des fichiers dupliqués.

## Consolidation scientifique

- Tableau initial distinguant production canonique, cycle thermique, trempe et
  campagne spatiale : budgets, capture des états, échanges de répliques et usages.
- Séparation entre exactitude du Hamiltonien discret, incertitude Monte-Carlo
  conditionnelle et erreur de modèle ; un test de code ne prouve pas l'équilibre.
- Discussion sourcée du matériau réel : transition autour de 660–663 K, ordre
  local persistant, calorimétrie avec contributions absentes du réseau rigide,
  diffusion locale et scénarios d'ordering dépendant de la trempe.
- Analyse explicite de l'écart du maximum L6 à l'expérience, sans l'attribuer
  artificiellement à une seule cause et sans promettre que L9 le corrigera.
- Distinction des variantes, du biais positif de la norme q dans une petite
  boîte, des fluctuations intra-variante et des longueurs apparentes de Fourier.
- Explication de l'absence de temps physique dans les GIF et de l'impossibilité
  de conclure à un rayon critique ou à une loi de nucléation avec une trajectoire.
- Correction d'une séquence d'échappement Python qui rendait mal le symbole
  `ell_app` dans le texte LaTeX généré.

## Vérifications

**28 tests réussis**, dont deux nouveaux tests de capture spatiale : correspondance
des états sauvegardés avec les neuf observables, respect des intervalles même
lorsqu'ils ne coïncident pas avec les anciens segments, et validation des paramètres.
Compilation LaTeX finale : aucune référence indéfinie, aucun avertissement et
aucune boîte débordante. Les liens relatifs des deux GIF pointent sur des fichiers présents.

## L9 : relance après finalisation du rapport

La campagne est lancée dans `results/production_L9` avec sa configuration dédiée :
2916 sites, 729 Au, 2187 Cu, 46 températures, deux chaînes, 6000 balayages
d'équilibration et 30000 de production par chaîne. Les anciennes traces du calcul
interrompu sont archivées dans `results/production_L9/launches`.

L'option `store_states_every=100` ajoute 300 captures atomiques par température
et par chaîne, sans modifier le Hamiltonien ni l'intervalle des observables.
Ces états permettent de poursuivre l'analyse spatiale sans refaire une production.
Le suivi indique chaîne, phase, balayages et nombre de captures. Les archives de
chaînes sont publiées atomiquement ; il n'existe pas de reprise en milieu de chaîne.

Après achèvement : validation, figures L9 et comparaison automatique aux résultats
L6 dans `comparison_L6_L9` (courbes par atome, maximum, R-hat, ESS, échanges et
allers-retours). Aucun résultat L9 anticipé n'est ajouté au rapport L6 finalisé.
La comparaison de deux tailles reste exploratoire et ne remplace ni une étude de
convergence statistique, ni une extrapolation thermodynamique.

## Intégration après achèvement de L9

Les résultats achevés sont maintenant ajoutés à la partie 3 du bonus dans le
PDF principal, avec trois figures propres à L9 ou comparatives, un tableau
L6/L9 et le détail des diagnostics près du maximum. Le texte ne présente plus
L9 comme un calcul à venir. Le décalage de 5 K est qualifié de provisoire :
ESS trop faible, nombre de blocs insuffisant à deux températures et aucun
aller-retour complet. Les fichiers L6 et leurs macros ne sont pas remplacés.
