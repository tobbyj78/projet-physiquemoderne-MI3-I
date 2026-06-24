# utils.py — constantes et fonctions utiles
# unites reduites : hbar = 1, m = 1

import numpy as np

# --- constantes ---
HBAR = 1.0
MASS = 1.0

# --- parametres du paquet d'ondes ---
K0 = 2.0
A0 = 3.0
X0 = -10.0

# --- grandeurs derivees ---
E0 = (HBAR**2 * K0**2) / (2.0 * MASS)
VG = (HBAR * K0) / MASS
VP = (HBAR * K0) / (2.0 * MASS)
LAMBDA = (2.0 * np.pi) / K0


def barriere(x, V0, a):
    # V(x)=V0 pour 0<x<a, 0 sinon
    V = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] > 0.0 and x[i] < a:
            V[i] = V0
    return V


def pos_max(densite, x):
    return x[np.argmax(densite)]


def coef_transmission(E, V0, a):
    # coefficient T pour une barriere rectangulaire
    if E >= V0:
        kp = np.sqrt(2.0 * MASS * (E - V0)) / HBAR
        k = np.sqrt(2.0 * MASS * E) / HBAR
        arg = kp * a
        denom = 1.0 + (V0**2 * np.sin(arg)**2) / (4.0 * E * (E - V0))
    else:
        kappa = np.sqrt(2.0 * MASS * (V0 - E)) / HBAR
        arg = kappa * a
        denom = 1.0 + (V0**2 * np.sinh(arg)**2) / (4.0 * E * (V0 - E))
    return 1.0 / denom


def info():
    print("=" * 50)
    print("  Projet : temps de traversee par effet tunnel")
    print("=" * 50)
    print("  hbar = %.1f,  m = %.1f" % (HBAR, MASS))
    print("  k0 = %.1f,  a0 = %.1f,  x0 = %.1f" % (K0, A0, X0))
    print("  E0 = %.1f,  v_g = %.1f" % (E0, VG))
    print("=" * 50)
