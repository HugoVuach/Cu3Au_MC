# Présentation du projet Cu₃Au

- `presentation.pdf` : diaporama à projeter, au format 16:9 (240 × 135 mm).
- `presentation.tex` : source LaTeX Beamer modifiable, police de corps de 9 points, légendes et références plus petites.
- `figures/` : figures PDF vectorielles copiées du rapport et des plaquettes.
- `animations/` : les deux GIF, accessibles depuis les liens de l'annexe A3.
- `verification.json` : contrôles de compilation, pages des questions et provenance des figures.

Le diaporama comporte **24 diapositives** : une couverture, exactement une diapositive
pour chacune des neuf questions du sujet, puis quatorze annexes.
Les diapositives des questions 1–9 se réfèrent à la production canonique L6.
Les tailles L9 et L12 sont explicitement comparées dans les annexes.

| Pages | Contenu |
|---|---|
| 1 | Titre, cadre et parcours de lecture |
| 2–10 | Questions 1 à 9, dans l'ordre du sujet |
| 11–13 (A1–A3) | Cycle thermique, occupations des sous-réseaux, plaquettes et GIF |
| 14 (A4) | Trempe isotherme et limites de l'interprétation cinétique |
| 15–16 (A5–A6) | Longueur de corrélation : définition, FFT, résultats et réserves |
| 17–20 (A7–A10) | Ordre, diagnostics, énergie et corrélations pour L6/L9/L12 |
| 21–22 (A11–A12) | Chauffage/refroidissement des trois tailles et analyse L12 |
| 23 (A13) | Boucles de l'algorithme complet en pseudocode |
| 24 (A14) | Références, provenance des données et vocabulaire |

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
