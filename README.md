# Projet : Temps de traversée par effet tunnel

**Mondir Boughrous & Tom Jaffrain — Groupe MI3_projet_I**

---

## Objectif

Déterminer le temps que met une particule quantique (modélisée par un paquet d'ondes gaussien 1D) pour franchir une barrière rectangulaire de potentiel par **effet tunnel**, numériquement ET analytiquement.

## Structure du projet

```
projet_tunnel/
├── utils.py                  # Constantes partagées, fonctions utilitaires
├── OndePLane1d3I.py          # Partie 1 : Onde plane 1D + superposition
├── PaquetOndeGauss1d3I.py    # Partie 2 : Paquet d'ondes gaussien analytique
├── schrodinger_solver.py     # Partie 3 : Solveur Schrödinger numérique
├── tunnel_effect.py          # Partie 4 : Effet tunnel, temps de traversée
├── main.py                   # Point d'entrée (lance tout ou une partie)
├── prof/                     # Codes originaux du professeur (référence)
│   ├── OndePlane.py
│   └── PaquetOndes.py
├── graphs/                   # Graphiques générés (créé automatiquement)
└── README.md
```

## Unités réduites

Pour simplifier les calculs, on utilise les unités réduites :

| Grandeur | Valeur |
|----------|--------|
| ℏ (constante de Planck réduite) | 1 |
| m (masse) | 1 |
| k₀ (nombre d'onde central) | 2.0 |
| a₀ (largeur initiale du paquet) | 3.0 |
| E₀ = ℏ²k₀²/(2m) | 2.0 |
| v_g = ℏk₀/m | 2.0 |

## Utilisation

### Lancer une partie spécifique

```bash
python OndePLane1d3I.py        # Partie 1 : ondes planes
python PaquetOndeGauss1d3I.py  # Partie 2 : paquet gaussien
python schrodinger_solver.py   # Partie 3 : solveur Schrödinger
python tunnel_effect.py        # Partie 4 : effet tunnel
```

### Lancer tout le projet

```bash
python main.py
```

## Dépendances

- Python ≥ 3.10
- NumPy
- Matplotlib

```bash
pip install numpy matplotlib
```

## Résumé des résultats

### Partie 1 — Ondes planes
- Implémentation d'une onde plane 1D : Ψ(x,t) = A·e^(i(kx−ωt))
- Superposition de 3 ondes planes → phénomène de **battements**
- Enveloppe : |A|[1 + cos(Δk·x/2)]

### Partie 2 — Paquet d'ondes gaussien
- Solution analytique pour une particule libre
- À t=0 : Ψ(x,0) = (2/(πa²))^(1/4)·e^(ik₀x)·e^(−x²/a²)
- Dispersion du paquet au cours du temps
- Conservation de la norme vérifiée

### Partie 3 — Solveur Schrödinger
- Schéma d'Euler explicite pour l'équation de Schrödinger dépendante du temps
- Dérivées par différences finies
- Validation : comparaison numérique/analytique pour V=0

### Partie 4 — Effet tunnel
- Mesure numérique de τ₀ (temps sans barrière) et τ_t (temps tunnel)
- Coefficient de transmission T (formule analytique)
- Étude de l'influence de la largeur a et de la hauteur V₀
- Comparaison avec les prédictions classiques

## Remarques

- Le schéma d'Euler explicite est conditionnellement stable : Δt < 2mΔx²/ℏ
- Les simulations peuvent prendre quelques minutes selon la résolution choisie
