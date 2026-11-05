# Cu3Au — transition ordre–désordre par Monte-Carlo

Projet complet répondant aux neuf questions de l'énoncé : réseau CFC rigide,
864 sites, composition Cu3Au conservée, interactions de premiers et seconds
voisins, Metropolis, Cowley, ordre à longue portée, énergie de cohésion,
configurations thermiques, courbes, transition et discussion.

Le dossier livré comprend les résultats réellement calculés, dix-huit figures
scientifiques en PDF vectoriel et PNG, les archives des échantillons et un
rapport français LaTeX/PDF avec démonstrations et encadrés pour l'oral.

## Lire les résultats déjà fournis

- `report/rapport_template.pdf` : rapport enrichi compilé, à lire en priorité.
- `results/production/report/rapport.pdf` : copie synchronisée du rapport enrichi.
- `results/production/report/rapport.tex` : source du rapport, modifiable.
- `results/production/report/auto_results.tex` : nombres injectés depuis la simulation.
- `results/production/figures/` : figures séparées, adaptées à des diapositives.
- `results/production/observables.csv` : mesures, erreurs et diagnostics par température.
- `results/production/summary.json` : protocole, maximum thermique, taux d'acceptation.
- `results/production/validation.json` : vérifications des résultats.
- `results/production/tests.log` : résultats des tests automatiques.
- `results/production/tests_current.log` et `verification_current.json` : tests et contrôles datés de la mise à jour.
- `results/production/peak_uncertainty.json` : bootstrap conjoint du maximum de capacité calorifique et sensibilité à la longueur des blocs.
- `results/production/chain_01.npz`, `chain_02.npz` : données brutes et configurations.
- `VERIFICATION.md` : compte rendu de l'exécution livrée.

## Installation

Python **3.11 ou ultérieur**, recommandé : **Python 3.12**. Les versions de
dépendances utilisées pendant la validation sont figées dans `requirements.txt`.

```bash
python -m venv .venv
```

Activation Windows PowerShell :

```powershell
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Activation macOS/Linux :

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
```

La compilation du rapport nécessite `pdflatex`, fourni par TeX Live ou MiKTeX,
avec `amsmath`, `geometry`, `graphicx`, `booktabs`, `lmodern`, `microtype`,
`fancyhdr`, `hyperref` et `tcolorbox`. Le français de Babel est utilisé s'il
est disponible. Les figures sont produites par Matplotlib sans LaTeX externe.
La première exécution compile le noyau Numba, ce qui ajoute un court délai.

## Exécuter de A à Z

Depuis le dossier du projet :

```bash
python run_project.py --preset production --output results/mon_execution
```

Le lanceur enchaîne les tests, les deux campagnes indépendantes, les
statistiques par blocs, les vérifications indépendantes, les dix-huit figures et
la compilation du rapport. Il affiche l'avancement. Compter plusieurs minutes
pour la production, selon le processeur. Les échanges de répliques sont
calculés dans un noyau compilé ; ce projet ne nécessite pas de GPU.

Pour vérifier rapidement l'installation :

```bash
python run_project.py --preset quick --output results/essai --no-pdf
```

Ce protocole court est un test de fonctionnement : il n'est pas recommandé
pour conclure quantitativement sur la transition.

Pour régénérer les figures et le rapport depuis les données livrées :

```bash
python run_project.py --preset production --reuse
```

Pour produire les sources LaTeX sans installation LaTeX : ajouter `--no-pdf`.
Les données et figures restent complètes. Une compilation PDF demandée mais
impossible fait échouer la commande avec un message explicite.

Pour les tests uniquement :

```bash
python run_project.py --tests-only
```

Une configuration personnalisée est un JSON de même structure que
`configs/production.json` :

```bash
python run_project.py --config configs/production.json --output results/autre_run
```

Pour compiler manuellement le rapport déjà livré :

```bash
cd results/production/report
pdflatex -interaction=nonstopmode -halt-on-error rapport.tex
pdflatex -interaction=nonstopmode -halt-on-error rapport.tex
```

Le fichier `report/rapport_template.tex` à la racine est le **gabarit**, utilisé
par le programme. Pour un rapport autonome, employer le dossier généré
`results/production/report/`, qui inclut les données LaTeX et les figures.

## Bonus 10 : cycle chauffage puis refroidissement

Le module `cu3au/thermal_cycle.py` suit une configuration de palier en palier,
depuis un état ordonné jusqu'à 4000 K puis jusqu'à 300 K. Les deux trajectoires
indépendantes utilisent la configuration `configs/thermal_cycle_L6.json` :
L=6, 46 températures par sens, 6000 balayages initiaux puis 30000 balayages
de production à chaque palier, mesures tous les 10 balayages. Les échanges
de répliques sont désactivés pour conserver l'histoire du cycle.

```powershell
.\.venv\Scripts\python.exe run_thermal_cycle.py --config configs/thermal_cycle_L6.json --output results/thermal_cycle_L6 --reference-report results/production
```

Après une interruption, ajouter `--resume` à cette commande. La reprise
utilise le dernier palier entièrement sauvegardé et rejoue un éventuel
palier partiel ; elle ne réinitialise pas toute la trajectoire.

Fichiers produits dans `results/thermal_cycle_L6/` :

- `run.log` et `execution.json` : progression et état du processus ;
- `chain_01/step_000.npz`, etc. : toutes les mesures enregistrées en
  équilibration et production, états initial/final et générateur aléatoire ;
- `observables_heating.csv`, `observables_cooling.csv` : moyennes, erreurs,
  occupations, acceptations et diagnostics, classés par température croissante ;
- `paired_differences.csv` : écarts refroidissement moins chauffage,
  calculés cycle par cycle ;
- `summary.json` : maxima, diagnostics et portée des résultats ;
- `figures/bonus_01_thermodynamique.pdf` à `bonus_08_distributions.pdf`,
  avec leurs versions PNG : comparaison des branches et diagnostics.

Avec `--reference-report results/production`, la section 10 présente d'abord
le protocole et signale le calcul en cours. Après validation des deux branches,
les huit figures et la synthèse chiffrée sont insérées automatiquement et
les deux copies du PDF sont recompilées. Les derniers PDF valides sont
conservés si une compilation échoue. `--no-pdf` produit les sources sans
compilation ; omettre `--reference-report` laisse le rapport principal intact.

Une différence entre les branches est un effet de protocole à durée finie.
Elle ne démontre pas à elle seule une transition du premier ordre. Les
fluctuations ne deviennent une réponse thermique d'équilibre qu'après
vérification des diagnostics. Les deux branches d'un cycle sont appariées ;
la dispersion entre deux cycles ne constitue pas une preuve de significativité.

## Conventions scientifiques

- `state[i] = 1` : Au ; `0` : Cu. Toujours N/4 Au, 3N/4 Cu.
- Coordonnées entières en unités a/2 ; conditions périodiques modulo 2L.
- Un balayage = N propositions **Au–Cu**. Les tirages de même espèce sont retirés.
- Énergie exacte stockée en **micro-eV entiers**, exponentielles en flottants.
- Mouvements globaux de permutation : échantillonnage d'équilibre, pas diffusion physique.
- `eta` est exactement le paramètre du sujet, référencé au sous-réseau 0.
- `q` est une amplitude reconnaissant les quatre variantes de L12.
- `cv_kB` est la capacité calorifique **par atome en unités de kB**.
- Les erreurs sont des **erreurs-types**, estimées par bootstrap de blocs et
  au minimum par la dispersion des chaînes ; elles ne couvrent pas tous les biais.
- Le maximum de Cv est celui de la **grille et de la boîte finie**. Le voisinage
  des deux températures adjacentes n'est pas un intervalle de confiance.
- La comparaison de champ moyen (~1712 K) est une approximation distincte.

Le projet emploie Metropolis canonique avec échanges de répliques pour aider
l'équilibrage. Chaque température garde la distribution canonique cible. La
grille de production est resserrée près de la transition repérée pendant le
pilote, sans imposer de valeur à la simulation.

## Valeurs analytiques de contrôle

| Test | Valeur |
|---|---:|
| Sites / Au / Cu | 864 / 216 / 648 |
| Coordination premier / second voisinage | 12 / 6 |
| Paires premier / second voisinage | 5184 / 2592 |
| Énergie ordonnée totale | −3294.432 eV |
| Énergie ordonnée par atome | −3.813 eV |
| Cowley dans L12 | alpha1 = −1/3 ; alpha2 = 1 |
| Ordre dans la variante de référence | eta = 1 ; q = 1 |
| Échange Au–Cu premiers voisins depuis L12 | +0.5586 eV |
| Cowley aléatoire à N fixé | −1/(N−1) |

L'identité exacte `E − E0 = 0.161 N_AuAu(1) + 0.0063 N_AuCu(2)` démontre
le minimum global et sert de contrôle indépendant de toutes les mesures.

## Lire les archives NPZ

```python
import numpy as np
d = np.load('results/production/chain_01.npz')
print(d.files)
print(d['samples'].shape)  # mesures, températures, 9 observables
```

Colonnes des neuf observables : énergie totale (eV), eta, q, alpha1, alpha2,
p0, p1, p2, p3 (fractions d'Au dans chacun des quatre sous-réseaux).
`burn_samples` contient les mesures d'équilibration écartées ; `labels`
suit les identités de répliques par température ; `final_states` contient
les configurations finales. `xyz` est en unités de a/2.

## Diagnostics et limites

La validation technique ne se confond pas avec la preuve d'équilibre. Le code
publie R-hat par rangs sur énergie et q, temps d'autocorrélation, ESS,
comparaison entre demi-séries, taux d'acceptation des échanges et trajectoires
des répliques. Il signale explicitement les contrôles non satisfaits.
Pour renforcer un résultat, allonger burn/production dans un nouveau JSON,
ajouter des graines, densifier les températures si leurs échanges sont rares,
et comparer plusieurs tailles L>=3. Une étude d'une seule taille ne détermine
pas la limite thermodynamique ni, à elle seule, l'ordre de la transition.

Le réseau est rigide : pas de phonons, de dilatation, de relaxation locale,
de lacunes ni de paramètres de paire dépendant de T. Les températures
obtenues concernent le modèle fourni, sans calibration expérimentale ajoutée.

## Références

- Metropolis et al. (1953), doi:10.1063/1.1699114.
- Hukushima et Nemoto (1996), doi:10.1143/JPSJ.65.1604.
- Binder (1987), doi:10.1088/0034-4885/50/7/001.

Les figures sont générées à partir des données sauvegardées, sans image ni
donnée téléchargée pendant l'exécution. Le programme ne requiert aucun accès
réseau après installation des dépendances.

## Plaquettes animées des coupes cristallographiques

```powershell
.\.venv\Scripts\python.exe -m cu3au.slice_animations
```

Les deux GIF sont écrits dans `results/thermal_cycle_L6/visualisation/animations`.
Chaque image montre trois couches atomiques exactes pour chacune des orientations
(001), (110), (111), avec plans, cadrage et couleurs Au/Cu fixes.

- `plaquette_cycle_thermique.gif` : cycle 1 archivé, état initial puis 92 fins de
  palier, du chauffage au refroidissement. Le compteur inclut équilibration et
  production. Les configurations intermédiaires entre paliers ne sont pas disponibles.
- `plaquette_relaxation_isotherme.gif` : calcul dédié L=6 après trempe de l'état
  archivé à 4000 K vers 1000 K, graine 20261104, 6000 balayages, captures tous les
  20 balayages sans suppression du transitoire. Les 200 premiers balayages sont
  ralentis dans le GIF ; le compteur indique toujours la valeur exacte.

Les états et observables de relaxation sont conservés dans
`relaxation_isotherme.npz`. Un second lancement réutilise cette archive après
vérification des paramètres. Les JSON accompagnant les GIF décrivent les plans,
les températures, les balayages et q à chaque image. Des vues initiales/finales
PNG et PDF sont également exportées. Le temps Monte-Carlo n'est pas un temps
physique ; ces images d'occupation seules ne démontrent pas une nucléation.

## Bonus, partie 2 : corrélations spatiales et second moment

```powershell
.\.venv\Scripts\python.exe -m cu3au.correlation_campaign
.\.venv\Scripts\python.exe -m cu3au.correlation_analysis
.\.venv\Scripts\python.exe -m cu3au.correlation_report
```

La configuration `configs/correlation_L6.json` définit une campagne canonique
dédiée : L=6, 14 températures, deux graines indépendantes, 3000 balayages
d'équilibration et 12000 de production ; états capturés tous les 20 balayages
et échanges de répliques tous les 5 balayages. Chaque chaîne complète est
sauvegardée atomiquement et réutilisée uniquement si sa configuration correspond.
Les états de départ viennent de l'ancienne production à la température disponible
la plus proche, puis sont rééquilibrés.

`cu3au/correlation.py` construit trois champs d'ordre signés par maille, adaptés
aux quatre variantes L1₂. La FFT porte sur ce champ d'ordre de pas a ; son mode
zéro n'est pas celui de la concentration totale, fixée. Le module distingue
spectres brut, connexe et part cohérente, reconstruit les corrélations réelles
et estime la longueur de second moment sur la boîte périodique.

Les résultats sont conservés dans `results/correlation_L6` :

- `chain_01.npz`, `chain_02.npz` : états atomiques, observables, étiquettes et
  compteurs d'échanges pour toutes les températures.
- `correlation_vs_temperature.csv` : longueurs brute/connexe, diagnostic
  intra-variante froid, bootstrap par blocs, R-hat, ESS et réserves explicites.
- `spectra.npz` : spectres complets brut et connexe.
- `block_sensitivity.csv` : sensibilité du bootstrap aux longueurs de blocs.
- `lengths_vs_time.csv` : analyse des états du GIF isotherme, longueur apparente,
  largeur diffuse hors mode uniforme, poids du mode uniforme.
- `figures/` : quatre figures de corrélations et une analyse des historiques des GIF.

Le bootstrap réestime la moyenne complexe dans chaque chaîne et chaque répétition.
Les racines négatives/dénominateurs nuls donnent des valeurs manquantes ; les
grandes longueurs ne sont pas écrêtées. La visite des variantes peut maintenir
une forte puissance canonique au mode zéro même après soustraction de la moyenne.
Le diagnostic intra-variante est conditionnel, limité au régime froid, et n'est
pas substitué à la définition canonique. Les longueurs temporelles ne sont ni
des longueurs de fluctuation d'équilibre ni des rayons géométriques de germes.

La dernière commande ajoute l'analyse des GIF, le protocole et les résultats de
corrélation à la section 10 puis recompile `report/rapport_template.pdf` et sa
copie de production. Elle requiert les deux GIF et `relaxation_isotherme.npz`
générés par `cu3au.slice_animations`. Les GIF sont également copiés dans un
sous-dossier `animations` à côté de chaque PDF afin de conserver des liens relatifs.
