# Projet : Temps de traversée par effet tunnel

**Mondir Boughrous & Tom Jaffrain — Groupe MI3_projet_I**

---

## Objectif

Déterminer le temps que met une particule quantique (modélisée par un paquet d'ondes gaussien 1D) pour franchir une barrière rectangulaire de potentiel par **effet tunnel**, numériquement ET analytiquement.

## Fichiers

- utils.py : constantes et fonctions utiles (hbar=1, m=1, barriere, coef T...)
- OndePLane1d3I.py : partie 1, ondes planes 1D et superposition
- PaquetOndeGauss1d3I.py : partie 2, paquet d'ondes gaussien
- schrodinger_solver.py : partie 3, solveur Schrodinger (Euler explicite)
- tunnel_effect.py : partie 4, effet tunnel et temps de traversee
- main.py : lance tout ou une partie

## Unités réduites

Pour simplifier les calculs, on utilise les unités réduites :

- hbar = 1
- m = 1
- k0 = 2.0
- a0 = 3.0
- E0 = hbar^2 k0^2 / (2m) = 2.0
- v_g = hbar k0 / m = 2.0

## Utilisation

```bash
python3 main.py          # menu interactif
python3 OndePLane1d3I.py # ou lancer une partie directement
```

## Dépendances

Python 3, numpy, matplotlib.

```bash
pip install numpy matplotlib
```

## Contenu

Partie 1 : onde plane 1D Psi(x,t) = A exp(i(kx - wt)), superposition de
3 ondes planes -> battements, enveloppe |A|[1 + cos(Dk x/2)].

Partie 2 : paquet d'ondes gaussien, solution analytique particule libre.
A t=0 : Psi(x,0) = (2/(pi a^2))^(1/4) exp(i k0 x) exp(-x^2/a^2).
Dispersion au cours du temps, norme conservee.

Partie 3 : solveur Schrodinger 1D par Euler explicite, derivees par
differences finies. Validation avec la solution analytique pour V=0.

Partie 4 : mesure du temps de traversee tau0 (sans barriere) et
tau_t (avec barriere, effet tunnel). Coefficient de transmission T.
Influence de la largeur a et de la hauteur V0.

## Remarques

- Le schema d'Euler explicite est stable si dt < 2 m dx^2 / hbar.
- Les simulations peuvent prendre quelques minutes.
