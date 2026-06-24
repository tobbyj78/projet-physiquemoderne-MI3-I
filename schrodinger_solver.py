# schrodinger_solver.py — solveur Schrodinger 1D (Euler explicite)
# Projet tunnel — MI3

import numpy as np
import matplotlib.pyplot as plt
import utils
from PaquetOndeGauss1d3I import gauss_wp


def derive1(f, dx):
    # derivee premiere : (f_{i+1} - f_i) / dx
    n = len(f)
    df = np.zeros(n, dtype=complex)
    for i in range(n - 1):
        df[i] = (f[i + 1] - f[i]) / dx
    df[-1] = df[-2]
    return df


def derive2(f, dx):
    # derivee seconde : (f_{i+1} - 2f_i + f_{i-1}) / dx^2
    n = len(f)
    d2f = np.zeros(n, dtype=complex)
    for i in range(1, n - 1):
        d2f[i] = (f[i + 1] - 2.0 * f[i] + f[i - 1]) / (dx * dx)
    d2f[0] = d2f[1]
    d2f[-1] = d2f[-2]
    return d2f


def valider_derivees():
    print("\n-- Validation derivees --")
    x = np.linspace(0, 5, 100)
    dx = x[1] - x[0]
    f = x**2

    df_num = derive1(f, dx)
    d2f_num = derive2(f, dx)

    df_exact = 2.0 * x
    d2f_exact = 2.0 * np.ones(len(x))

    mask = np.abs(df_exact) > 1e-10
    err1 = np.max(np.abs(df_num[mask] - df_exact[mask]) / np.abs(df_exact[mask]))
    err2 = np.max(np.abs(d2f_num[1:-1] - d2f_exact[1:-1]))

    print("  erreur relative max f'  = %.6f" % err1)
    print("  erreur relative max f'' = %.6f" % err2)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle("Validation derivees — f(x)=x^2", fontsize=13)

    ax1.plot(x, df_exact, "k-", lw=2, label="exact")
    ax1.plot(x, df_num, "r--", lw=1.5, label="numerique")
    ax1.set_ylabel("f'(x)")
    ax1.legend()
    ax1.grid(True, alpha=0.4)

    ax2.plot(x, d2f_exact, "k-", lw=2, label="exact")
    ax2.plot(x, d2f_num, "r--", lw=1.5, label="numerique")
    ax2.set_ylabel("f''(x)")
    ax2.set_xlabel("x")
    ax2.legend()
    ax2.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig("graphs/validation_derivees.png", dpi=150)
    plt.show()


def pas_euler(psi, dx, dt, V):
    # Psi^{n+1} = Psi^n + (i*hbar*dt/2m)*d^2Psi/dx^2 - (i*dt/hbar)*V*Psi^n
    d2psi = derive2(psi, dx)
    kin = (1.0j * utils.HBAR * dt) / (2.0 * utils.MASS) * d2psi
    pot = (-1.0j * dt / utils.HBAR) * V * psi
    return psi + kin + pot


def solve(psi0, x, temps, V=None, stock=1):
    nx = len(x)
    nt = len(temps)
    dx = x[1] - x[0]
    dt = temps[1] - temps[0]

    if V is None:
        V = np.zeros(nx)

    n_frames = (nt + stock - 1) // stock
    hist = np.zeros((n_frames, nx), dtype=complex)
    hist[0] = psi0.copy()

    psi = psi0.copy()
    idx = 1

    for n in range(1, nt):
        psi = pas_euler(psi, dx, dt, V)
        if n % stock == 0 and idx < n_frames:
            hist[idx] = psi.copy()
            idx += 1

    return hist


def valider_solveur():
    print("\n-- Validation solveur Schrodinger (V=0) --")

    dx = 0.1
    x = np.arange(-20, 40 + dx, dx)
    dt = 0.001
    tmax = 10.0
    temps = np.arange(0, tmax + dt, dt)
    stock = 100

    print("  nx=%d, dx=%.2f, dt=%.4f, tmax=%.0f" % (len(x), dx, dt, tmax))

    psi0 = gauss_wp(utils.K0, utils.A0, x - utils.X0, 0.0)

    print("  resolution en cours...")
    hist = solve(psi0, x, temps, V=None, stock=stock)

    psi_ana = gauss_wp(utils.K0, utils.A0, x - utils.X0, tmax)
    dens_num = np.abs(hist[-1])**2
    dens_ana = np.abs(psi_ana)**2
    err = np.max(np.abs(dens_num - dens_ana)) / np.max(dens_ana)
    print("  erreur relative max sur |Psi|^2 a t=%.0f : %.6f" % (tmax, err))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle("Validation : numerique vs analytique (V=0)", fontsize=13)

    for idx, t_val in [(0, 0.0), (len(hist)//2, tmax/2), (len(hist)-1, tmax)]:
        psi_a = gauss_wp(utils.K0, utils.A0, x - utils.X0, t_val)
        psi_n = hist[idx]
        ax1.plot(x, np.abs(psi_a)**2, lw=1.5, alpha=0.7,
                 label="analytique t=%.1f" % t_val)
        ax1.plot(x, np.abs(psi_n)**2, "--", lw=1.2, alpha=0.7,
                 label="numerique t=%.1f" % t_val)

    ax1.set_ylabel("|Psi(x,t)|^2")
    ax1.legend(fontsize=8, ncol=2)
    ax1.grid(True, alpha=0.4)

    ax2.plot(x, np.abs(dens_num - dens_ana), color="crimson", lw=1.5)
    ax2.set_ylabel("| |Psi_num|^2 - |Psi_ana|^2 |")
    ax2.set_xlabel("x")
    ax2.grid(True, alpha=0.4)

    for ax in (ax1, ax2):
        ax.set_xlim(x[0], x[-1])

    plt.tight_layout()
    plt.savefig("graphs/validation_schrodinger.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    import os
    os.makedirs("graphs")
    print("=== Partie 3 : Solveur Schrodinger ===")
    valider_derivees()
    valider_solveur()
    print("Termine.")
