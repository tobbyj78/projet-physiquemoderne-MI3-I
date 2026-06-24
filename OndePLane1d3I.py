# OndePLane1d3I.py — ondes planes 1D et superposition
# Projet tunnel — MI3

import numpy as np
import matplotlib.pyplot as plt

HBAR = 1.0
MASS = 1.0


def onde_plane(amp, k, omega, x, t):
    # Psi(x,t) = A * exp(i(kx - wt))
    return amp * np.exp(1.0j * (k * x - omega * t))


def omega(k, dispersion="schrodinger"):
    # relation de dispersion : schrodinger -> w = hbar k^2 / (2m)
    if dispersion == "schrodinger":
        return (HBAR * k**2) / (2.0 * MASS)
    else:
        return k


def plot_onde_unique():
    k0 = 2.0
    lam = 2.0 * np.pi / k0
    x = np.linspace(0, 5 * lam, 500)
    w0 = omega(k0)
    psi = onde_plane(1.0, k0, w0, x, 0.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle("Onde plane 1D", fontsize=13)

    ax1.plot(x, np.real(psi), color="#3266ad", lw=1.5)
    ax1.set_ylabel("Re(Psi)")
    ax1.grid(True, alpha=0.4)
    ax1.axhline(0, color="gray", lw=0.5, ls="--")

    ax2.plot(x, np.imag(psi), color="#c44e52", lw=1.5)
    ax2.set_ylabel("Im(Psi)")
    ax2.set_xlabel("x")
    ax2.grid(True, alpha=0.4)
    ax2.axhline(0, color="gray", lw=0.5, ls="--")

    for ax in (ax1, ax2):
        ax.set_xlim(x[0], x[-1])

    plt.savefig("graphs/onde_plane_unique.png", dpi=150)
    plt.show()


def plot_superposition():
    k0 = 2.0
    dk = k0 / 4.0
    A = 1.0
    t = 0.0

    ks = [k0 - dk/2, k0, k0 + dk/2]
    amps = [A/2, A, A/2]
    omegas = [omega(k) for k in ks]

    L = np.pi / dk
    x = np.linspace(-L, L, 800)

    psi_tot = np.zeros(len(x), dtype=complex)
    psi_indiv = []
    for k, amp, w in zip(ks, amps, omegas):
        p = onde_plane(amp, k, w, x, t)
        psi_indiv.append(p)
        psi_tot += p

    # enveloppe : |A| * |1 + cos(dk * x / 2)|
    env = np.abs(A) * np.abs(1.0 + np.cos(dk * x / 2.0))
    densite = np.abs(psi_tot)**2

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle("Superposition de 3 ondes planes", fontsize=14)

    couleurs = ["#3266ad", "#1d9e75", "#c44e52"]
    for i, (p, k) in enumerate(zip(psi_indiv, ks)):
        ax1.plot(x, np.real(p), color=couleurs[i], lw=0.8, alpha=0.7,
                 label="k%d = %.3f" % (i+1, k))
    ax1.plot(x, np.real(psi_tot), color="crimson", lw=2.0, label="somme")
    ax1.set_ylabel("Re(Psi)")
    ax1.legend(fontsize=8, ncol=2)
    ax1.grid(True, alpha=0.4)
    ax1.axhline(0, color="gray", lw=0.5, ls="--")

    ax2.plot(x, densite, color="crimson", lw=2.0, label="|Psi|^2")
    ax2.plot(x, env**2, color="gray", lw=1.5, ls="--",
             label="|A|^2 [1+cos(dk x/2)]^2")
    ax2.set_ylabel("|Psi|^2")
    ax2.set_xlabel("x")
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.4)

    for ax in (ax1, ax2):
        ax.set_xlim(x[0], x[-1])

    plt.tight_layout()
    plt.savefig("graphs/superposition_3waves.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    import os
    os.makedirs("graphs")
    print("=== Partie 1 : Ondes planes ===")
    plot_onde_unique()
    plot_superposition()
    print("Termine.")
