# Relecture du rapport final — 11 novembre 2026

Rapport de Hugo Vuacheux et Mazarine Lambert. La revue couvre les neuf
questions de `subject.md`, les quatre parties du bonus, les inclusions
numériques, les figures et la méthodologie.

## Bilan scientifique

Les neuf réponses demandées sont présentes. Les conventions d'énergie,
les comptages de paires, le test de Metropolis, les échanges de répliques,
les paramètres de Cowley et d'ordre sont cohérents avec le sujet et le code.
L'identité énergie–Cowley, la composition et les relations entre occupations,
eta et q sont de nouveau vérifiées sur toutes les mesures des productions
L6/L9/L12. Les trois validations techniques passent.

La distinction entre validation du code et exploration statistique reste
essentielle : les critères locaux près du pic passent pour L6, mais pas
pour L9/L12. L'ESS minimal vaut environ 1191, 37 et 45, et le R-hat maximal
en fenêtre critique environ 1,005, 1,026 et 1,268. Les maxima observés sont
1050, 1055 et 1052,5 K. Ils ne constituent ni une extrapolation à volume
infini ni une température du matériau réel.

Les distributions et la mémoire du parcours sont compatibles avec un
premier ordre, mais ne le démontrent pas pour le Hamiltonien spatial.
Le caractère du matériau réel et la prédiction du champ moyen sont des
informations distinctes. Les longueurs spectrales ne sont pas des rayons
de domaines, et les balayages ne représentent pas un temps expérimental.

## Modifications apportées

- Ajout des deux auteurs sur la couverture et dans les métadonnées PDF ;
  introduction et conclusion formulées comme le travail du binôme.
- Ajout d'une conclusion générale après les résultats des quatre bonus,
  avec les résultats établis, les limites et quatre priorités de travail.
- Explicitation des tests permettant de distinguer premier ordre et
  transition continue : populations convergées, histogrammes d'énergie,
  séparation énergétique persistante et lois de taille à précision comparable.
  La capacité par atome et la capacité totale ne sont pas confondues.
- Vérification des entrées de la FFT, des modes autorisés et des indices
  entiers ; distinction entre intensité instantanée, moyenne signée et
  spectre connexe.
- Contrôle des quatre panneaux de relaxation : les trois premiers suivent
  la trempe ; le dernier vient de la campagne canonique en température.
  Les figures précédentes sont évacuées avant la partie sur la trempe.
- Distinction entre cohésion et formation : dans ce modèle, l'énergie de
  formation vaut −0,2415 eV/atome, soit −0,966 eV par formule, par rapport
  aux cristaux CFC purs de Cu et Au décrits par les mêmes paires.
- Correction d'une perspective obsolète demandant de réaliser les cycles
  L9/L12 : ils sont achevés ; il faut désormais les répéter et les prolonger.
- Correction rédactionnelle de « ressèrent » et formulation plus prudente
  sur la dépendance des cycles au parcours.
- Pour la version à rendre, retrait des 18 encadrés d'approfondissement
  pédagogique, y compris le guide de présentation orale ; suppression de
  l'annexe « Algorithme » et de son pseudocode. Les renvois vers l'annexe
  supprimée sont corrigés. La méthode demeure exposée dans les questions.
  Les blocs pédagogiques sont aussi filtrés lors de la publication des
  inclusions générées, afin qu'une régénération ne les réintroduise pas.

Les textes générés ont été corrigés dans leurs générateurs Python et
resynchronisés, pour conserver les modifications lors d'une régénération.
Les simulations de production et leurs données brutes ne sont pas relancées.

## Prolongements, par ordre de priorité

1. Prolonger les séries critiques, multiplier les graines et contrôler les
   proportions de séjour ainsi que les chaînes les plus lentes.
2. Comparer plusieurs tailles à qualité d'échantillonnage comparable avant
   de conclure sur la nature de la transition et sa limite thermodynamique.
3. Répéter les cycles à plusieurs durées de palier ; suivre des domaines
   locaux et leurs interfaces dans plusieurs trempes. Attribuer aussi une
   incertitude au diagnostic spatial intra-variante.
4. Pour une comparaison quantitative au matériau : enrichir et valider les
   interactions, traiter les relaxations/vibrations et le volume ; pour la
   cinétique, employer une dynamique locale et des barrières calibrées.

La référence méthodologique ajoutée est Challa, Landau et Binder,
*Finite-size effects at temperature-driven first-order transitions*,
Phys. Rev. B 34, 1841–1852 (1986),
https://doi.org/10.1103/PhysRevB.34.1841.

## Vérifications reproductibles

- 34 tests réussis : `results/comparison_L6_L9_L12/review_tests_20261111.log`.
- Contrôles numériques : `validation_review_20261111.json` dans chacune
  des trois productions canoniques.
- Compilation, renvois, citations, fichiers inclus, figures et liens GIF :
  `results/comparison_L6_L9_L12/report_audit_20261111.json`.
- Le PDF principal et sa copie dans `results/production/report/rapport.pdf`
  sont synchronisés après validation du document compilé.
