# Présentation du projet Cu₃Au

## Version épurée pour l'oral

`presentation_epuree.pdf` et `presentation_epuree.tex` proposent une version
séparée de 24 diapositives : titres explicites, un sous-titre et deux phrases
en puces par diapositive, avec figures, formules et paramètres. La question 8
conserve une seule puce après le retrait de la figure du bootstrap et de son commentaire.
Les encadrés « À retenir » restent supprimés. Les titres et pieds de figures
contenant des commentaires sont masqués dans les copies de `figures_epurees/`.
Les axes, légendes, données et symboles sont conservés. Les deux GIF restent
accessibles sur la diapositive A3 (clic sur les plaquettes pour la relaxation ;
clic sur la formule du cycle pour le chauffage/refroidissement).

```powershell
python rapport/compiler_epuree.py
```

La version détaillée ci-dessous est réorganisée séparément.

## Version détaillée

- `presentation.pdf` : diaporama à projeter, au format 16:9 (240 × 135 mm).
- `presentation.tex` : source LaTeX Beamer modifiable, police de corps de 9 points, légendes et références plus petites.
- `figures/` : figures PDF vectorielles copiées du rapport et des plaquettes.
- `animations/` : les deux GIF, accessibles depuis les liens de l'annexe A3.
- `verification.json` : contrôles de compilation, pages des questions et provenance des figures.

La version détaillée comporte **25 diapositives** : une couverture, deux
diapositives pour la question 1, puis les questions 2, 3, 5, 6, 7, 8 et 9.
La question 4 (tests) termine les annexes, après les quatorze annexes complémentaires.
Les diapositives des questions 1–9 se réfèrent à la production canonique L6.
Les tailles L9 et L12 sont explicitement comparées dans les annexes.

| Pages | Contenu |
|---|---|
| 1 | Titre, cadre et parcours de lecture |
| 2–3 | Question 1 : réseau/Hamiltonien, puis organigramme de l'algorithme |
| 4–5 | Questions 2 et 3 : Cowley et ordre à longue portée |
| 6–10 | Questions 5 à 9 |
| 11–13 (A1–A3) | Cycle thermique, occupations des sous-réseaux, plaquettes et GIF |
| 14 (A4) | Trempe isotherme et limites de l'interprétation cinétique |
| 15–16 (A5–A6) | Longueur de corrélation : définition, FFT, résultats et réserves |
| 17–20 (A7–A10) | Ordre, diagnostics, énergie et corrélations pour L6/L9/L12 |
| 21–22 (A11–A12) | Chauffage/refroidissement des trois tailles et analyse L12 |
| 23 (A13) | Boucles de l'algorithme complet en pseudocode |
| 24 (A14) | Références, provenance des données et vocabulaire |
| 25 | Annexe : question 4, vérifications du programme |

Les figures ne sont pas répétées à l'identique dans le diaporama. Les formules
et le texte sont éditables dans le fichier `.tex`. Les nombreuses courbes
proviennent des archives utilisées par le rapport principal ; elles ne sont
pas recalculées par une nouvelle simulation.

## Recompiler

Depuis la racine du projet :

```powershell
python rapport/compiler.py
```

Le script nécessite `pdflatex` (MiKTeX) et PyMuPDF dans le Python utilisé.
Il actualise les copies des figures/GIF depuis les résultats existants,
compile et contrôle le PDF avant de publier `presentation.pdf`.
Les fichiers temporaires et journaux restent dans `build/`.

Pour déplacer la présentation avec ses animations, copier le dossier `rapport`
en conservant `animations` à côté du PDF. Le lecteur PDF doit autoriser
l'ouverture des fichiers externes. Le PDF affiche des images fixes ; les
animations s'ouvrent dans un lecteur de GIF.
