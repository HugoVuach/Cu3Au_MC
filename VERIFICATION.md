# Vérification de la livraison

Exécution réalisée le 21 octobre 2026 sous Linux, Python 3.12.14.
Versions utilisées : NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8,
Numba 0.68.0. Le code livré ne nécessite pas de GPU.

## Vérifications réussies

- **15 tests automatiques** physiques, probabilistes et statistiques.
- Géométrie CFC : distances périodiques, 12 et 6 voisins distincts et réciproques.
- Quatre variantes L12 d'énergie identique, ordre invariant égal à un.
- Énergie totale ordonnée : **−3294.432 eV** ; par atome : **−3.813 eV**.
- Échange Au–Cu premiers voisins depuis L12 : **+0.5586 eV**.
- Comparaison entre variations locales et recalculs globaux de l'énergie,
  y compris pour des partenaires voisins et non voisins, sur L=3 et L=6.
- Identité de borne exacte et minimum global vérifiés.
- Balance détaillée de Metropolis dans un ensemble restreint énuméré exactement.
- Signe et balance détaillée de l'échange entre répliques.
- Reproductibilité d'un court calcul à graine fixée.
- Limite aléatoire à haute température, avec correction de composition finie.
- Autocorrélation testée sur une série AR(1) de temps de corrélation connu.
- R-hat testé sur des chaînes concordantes et volontairement séparées.
- Unités de capacité calorifique et rejet d'une grille invalide.

## Campagne de production réellement exécutée

Configuration : `configs/production.json`.

| Paramètre | Valeur |
|---|---:|
| Taille L / sites | 6 / 864 |
| Au / Cu | 216 / 648 |
| Températures | 46 |
| Chaînes indépendantes | 2 |
| Graines | 20261021 et 20261022 |
| Équilibration par chaîne | 6 000 balayages |
| Production par chaîne | 30 000 balayages |
| Intervalle entre mesures | 10 balayages |
| Mesures par chaîne et par température | 3 000 |
| Propositions d'échanges Au–Cu, équilibration comprise | 2 861 568 000 |
| Durée simulation et statistiques dans cet environnement | environ 525 secondes |

À la fin de chaque segment, les énergies suivies et recalculées concordent
**exactement en arithmétique entière**, et la composition reste constante.
L'identité reliant énergie et comptage de paires est vérifiée pour toutes les
mesures sauvegardées. Les valeurs CSV sont contrôlées contre les archives brutes.

## Résultat et diagnostics près du maximum

Le maximum de capacité calorifique apparaît à **1 050 K dans chacune des deux
chaînes**, ainsi que dans leur moyenne. Les voisins sur la grille sont 1 045
et 1 055 K. Ces bornes ne constituent pas un intervalle de confiance.

| T (K) | q moyen | cV/kB par atome | ESS minimal estimé | R-hat énergie | R-hat q |
|---:|---:|---:|---:|---:|---:|
| 1040 | 0.715 | 4.35 | 2114 | 1.0033 | 1.0045 |
| 1045 | 0.652 | 5.31 | 1476 | 1.0033 | 1.0037 |
| 1050 | 0.575 | 5.87 | 1191 | 1.0026 | 1.0030 |
| 1055 | 0.491 | 5.69 | 1253 | 1.0038 | 1.0029 |
| 1060 | 0.417 | 4.93 | 1688 | 1.0015 | 1.0016 |

Les critères locaux fixés dans le code sont satisfaits dans cette fenêtre.
Les identités de répliques effectuent au total 4 puis 3 trajets complets
froid–chaud–froid pendant la production. Leur nombre reste un diagnostic
complémentaire, pas une preuve d'exploration parfaite.

À 300 K : u = −3.813 eV/atome, q = 1.
À 4 000 K : u ≈ −3.757820 eV/atome, q ≈ 0.0654.
Les corrélations locales à 4 000 K ne sont pas forcées à zéro.

## Chaîne de production des livrables

- Régénération depuis les archives : tests, validation, 11 figures PDF/PNG,
  sources LaTeX et compilation PDF, avec succès.
- Nouveau lancement depuis un dossier vide, protocole court : simulation,
  statistiques, tests, validation, 11 figures et compilation PDF, avec succès.
- Le protocole court signale des diagnostics critiques insuffisants : il sert
  uniquement à tester le fonctionnement ; sa transition n'est pas celle retenue.
- Rapport principal : **17 pages**, compilé par deux passages de pdflatex.
- Inspection visuelle des 17 pages ; absence de débordements LaTeX signalés.
- Chaque PDF de figure a été ouvert par `pdfinfo` sans erreur.
- L'archive finale est vérifiée par contrôle ZIP et un manifeste SHA-256
  permet de contrôler les fichiers extraits.

## Portée de la validation

Les valeurs de référence à T=0 et le coût du défaut sont exacts pour ce
modèle. Le maximum à 1 050 K est un résultat Monte-Carlo de **taille finie**.
Il ne constitue ni une extrapolation thermodynamique, ni une température
expérimentale de Cu3Au. Les diagnostics favorables et les erreurs par blocs
ne couvrent pas toutes les erreurs systématiques possibles.

La portabilité Windows/macOS est prévue par l'emploi de Python, NumPy et
Numba, mais la livraison a été exécutée et validée sous Linux. Les variations
de plateforme peuvent modifier les trajectoires aléatoires ; on attend un
accord statistique, pas nécessairement des fichiers identiques bit à bit.
