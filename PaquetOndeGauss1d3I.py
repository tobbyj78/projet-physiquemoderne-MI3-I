# PaquetOndeGauss1d3I.py — paquet d'ondes gaussien (solution analytique)
# Projet tunnel — MI3

import numpy as np
import matplotlib.pyplot as plt

HBAR = 1.0
MASS = 1.0
K0 = 2.0
A0 = 3.0


def gauss_wp(k0, a, x, t):
    # Psi(x,t) pour un paquet d'ondes gaussien libre 1D
    a2 = a * a
    xc = (HBAR * k0 * t) / MASS
    denom = 1.0 + 1.0j * (2.0 * HBAR * t) / (MASS * a2)
    pref = (2.0 / (np.pi * a2))**0.25

    if np.abs(t) < 1e-15:
        # t = 0 : formule simplifiee pour eviter 0/0
        psi = pref * np.exp(1.0j * k0 * x) * np.exp(-x**2 / a2)
    else:
        norm_t = pref / np.sqrt(denom)
        phase = 1.0j * k0 * x - 1.0j * (HBAR * k0**2 * t) / (2.0 * MASS)
        gauss = np.exp(-((x - xc)**2) / (a2 * denom))
        psi = norm_t * np.exp(phase) * gauss

    return psi


def pos_analytique(t, k0=K0):
    return (HBAR * k0 * t) / MASS


def largeur_analytique(t, a=A0):
    a2 = a * a
    f = (2.0 * HBAR * t) / (MASS * a2)
    return (a / 2.0) * np.sqrt(1.0 + f**2)


def hauteur_analytique(t, a=A0):
    a2 = a * a
    f = (2.0 * HBAR * t) / (MASS * a2)
    return np.sqrt(2.0 / (np.pi * a2)) / np.sqrt(1.0 + f**2)


def plot_paquet_t0():
    a = A0
    x = np.linspace(-3*a, 3*a, 500)
    psi = gauss_wp(K0, a, x, 0.0)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Paquet d'ondes gaussien a t=0", fontsize=14)

    ax1.plot(x, np.real(psi), color="#3266ad", lw=1.5)
    ax1.set_ylabel("Re(Psi)")
    ax1.grid(True, alpha=0.4)
    ax1.axhline(0, color="gray", lw=0.5, ls="--")

    ax2.plot(x, np.imag(psi), color="#c44e52", lw=1.5)
    ax2.set_ylabel("Im(Psi)")
    ax2.grid(True, alpha=0.4)
    ax2.axhline(0, color="gray", lw=0.5, ls="--")

    dens = np.abs(psi)**2
    ax3.plot(x, dens, color="#1d9e75", lw=2.0)
    ax3.fill_between(x, 0, dens, color="#1d9e75", alpha=0.2)
    ax3.set_ylabel("|Psi|^2")
    ax3.set_xlabel("x")
    ax3.grid(True, alpha=0.4)

    for ax in (ax1, ax2, ax3):
        ax.set_xlim(x[0], x[-1])

    plt.tight_layout()
    plt.savefig("graphs/paquet_gaussien_t0.png", dpi=150)
    plt.show()


def plot_dispersion():
    a = A0
    temps = [0.0, 2.0, 5.0, 10.0, 15.0]
    x_max = pos_analytique(temps[-1]) + 5 * largeur_analytique(temps[-1])
    x = np.linspace(-3*a, x_max, 1000)

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle("Dispersion du paquet gaussien", fontsize=14)
    couleurs = plt.cm.viridis(np.linspace(0.15, 0.90, len(temps)))

    for t_i, c in zip(temps, couleurs):
        psi = gauss_wp(K0, a, x, t_i)
        dens = np.abs(psi)**2
        xc = pos_analytique(t_i)
        dx = largeur_analytique(t_i)
        ax.plot(x, dens, color=c, lw=1.5,
                label="t=%.0f  (xc=%.1f, dx=%.2f)" % (t_i, xc, dx))

    ax.set_xlabel("x")
    ax.set_ylabel("|Psi(x,t)|^2")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.4)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(bottom=0)

    plt.savefig("graphs/paquet_gaussien_dispersion.png", dpi=150)
    plt.show()


def verif_norme():
    a = A0
    temps = np.linspace(0, 20, 11)
    x = np.linspace(-50, 80, 2000)

    print("\nVerification norme :")
    print("  %5s  %12s" % ("t", "norme"))
    for t in temps:
        psi = gauss_wp(K0, a, x, t)
        norme = np.trapezoid(np.abs(psi)**2, x)
        print("  %5.1f  %.8f" % (t, norme))
    print("  -> la norme est conservee (~1)")


if __name__ == "__main__":
    import os
    os.makedirs("graphs")
    print("=== Partie 2 : Paquet d'ondes gaussien ===")
    verif_norme()
    plot_paquet_t0()
    plot_dispersion()
    print("Termine.")
