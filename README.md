# Projet : Temps de traversée par effet tunnel

**Mondir Boughrous & Tom Jaffrain — Groupe MI3_projet_I**

---

## Objectif

Determiner le temps que met une particule quantique (modelisee par un
paquet d'ondes gaussien 1D) pour franchir une barriere rectangulaire
de potentiel par effet tunnel, numeriquement et analytiquement.

## Fichiers

- utils.py : constantes et fonctions utiles (hbar=1, m=1, barriere, coef T...)
- OndePLane1d3I.py : partie 1, ondes planes 1D et superposition
- PaquetOndeGauss1d3I.py : partie 2, paquet d'ondes gaussien
- schrodinger_solver.py : partie 3, solveur Schrodinger (Euler explicite)
- tunnel_effect.py : partie 4, effet tunnel et temps de traversee
- main.py : lance tout ou une partie

## Unites reduites

Pour simplifier les calculs, on utilise les unites reduites :

| Grandeur | Valeur |
|----------|--------|
| hbar (constante de Planck reduite) | 1 |
| m (masse) | 1 |
| k0 (nombre d'onde central) | 2.0 |
| a0 (largeur initiale du paquet) | 3.0 |
| E0 = hbar^2 k0^2 / (2m) | 2.0 |
| v_g = hbar k0 / m | 2.0 |

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

Partie 1 : onde plane 1D et superposition de 3 ondes planes, on voit
les battements.

Partie 2 : paquet d'ondes gaussien, solution analytique pour une
particule libre. Le paquet se disperse avec le temps mais la norme
reste egale a 1.

Partie 3 : solveur Schrodinger 1D par Euler explicite, derivees par
differences finies. Comparaison avec la solution analytique pour V=0
(particule libre).

Partie 4 : mesure du temps de traversee tau0 sans barriere et tau_t
avec barriere (effet tunnel). Calcul du coefficient de transmission T.
Etude de l'influence de la largeur a et de la hauteur V0.

## Remarques

- Le schema d'Euler est stable si dt est assez petit devant dx^2.
- Les simulations peuvent prendre quelques minutes.
