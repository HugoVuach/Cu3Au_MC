# Cu3Au — transition ordre–désordre par Monte-Carlo

Projet complet répondant aux neuf questions de l'énoncé : réseau CFC rigide,
864 sites, composition Cu3Au conservée, interactions de premiers et seconds
voisins, Metropolis, Cowley, ordre à longue portée, énergie de cohésion,
configurations thermiques, courbes, transition et discussion.

Le dossier livré comprend les résultats réellement calculés, onze figures
scientifiques en PDF vectoriel et PNG, les archives des échantillons et un
rapport français LaTeX/PDF avec démonstrations et encadrés pour l'oral.

## Lire les résultats déjà fournis

- `results/production/report/rapport.pdf` : rapport compilé.
- `results/production/report/rapport.tex` : source du rapport, modifiable.
- `results/production/report/auto_results.tex` : nombres injectés depuis la simulation.
- `results/production/figures/` : figures séparées, adaptées à des diapositives.
- `results/production/observables.csv` : mesures, erreurs et diagnostics par température.
- `results/production/summary.json` : protocole, maximum thermique, taux d'acceptation.
- `results/production/validation.json` : vérifications des résultats.
- `results/production/tests.log` : résultats des tests automatiques.
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
statistiques par blocs, les vérifications indépendantes, les onze figures et
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
