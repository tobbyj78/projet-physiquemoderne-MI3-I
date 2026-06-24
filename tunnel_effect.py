# tunnel_effect.py — temps de traversee par effet tunnel
# Projet tunnel — MI3

import numpy as np
import matplotlib.pyplot as plt
import utils
from PaquetOndeGauss1d3I import gauss_wp
from schrodinger_solver import solve


def mesurer_tau(hist_dens, x, t_stock, x_in, x_out):
    # mesure le temps pour aller de x_in a x_out en suivant le max
    nf = len(t_stock)
    pos = np.zeros(nf)
    for i in range(nf):
        pos[i] = utils.pos_max(hist_dens[i], x)

    i_in = np.argmax(pos >= x_in)
    i_out = np.argmax(pos >= x_out)
    t_in = t_stock[i_in]
    t_out = t_stock[i_out]
    tau = t_out - t_in

    return tau, t_in, t_out, pos


def simuler(V0, a, t_max=30.0, dx=0.05, dt=0.0005):
    x_min = -15.0
    x_max = a + 25.0
    x = np.arange(x_min, x_max + dx, dx)

    temps = np.arange(0.0, t_max + dt, dt)
    stock = max(1, int(0.02 / dt))
    t_stock = temps[::stock]

    V = utils.barriere(x, V0, a)
    psi0 = gauss_wp(utils.K0, utils.A0, x - utils.X0, 0.0)

    print("  simulation : V0=%.1f, a=%.1f, t_max=%.1f" % (V0, a, t_max))
    print("  nx=%d, dx=%.3f, dt=%.5f, frames=%d" % (len(x), dx, dt, len(t_stock)))
    print("  resolution...")

    hist_psi = solve(psi0, x, temps, V=V, stock=stock)
    hist_dens = np.abs(hist_psi)**2

    tau, t_in, t_out, pos = mesurer_tau(hist_dens, x, t_stock, 0.0, a)

    print("  t_entree(x=0)=%.4f  t_sortie(x=%.1f)=%.4f" % (t_in, a, t_out))
    print("  tau = %.4f" % tau)

    return {
        "x": x, "t_stock": t_stock, "hist_dens": hist_dens,
        "V": V, "V0": V0, "a": a,
        "tau": tau, "t_in": t_in, "t_out": t_out, "pos": pos,
    }


def calculer_tau0(a, dx=0.05, dt=0.0005):
    print("\n-- 4.1.b : tau0 sans barriere --")
    res = simuler(V0=0.0, a=a, t_max=20.0, dx=dx, dt=dt)
    tau0_num = res["tau"]
    tau0_th = a / utils.VG
    ecart = abs(tau0_num - tau0_th) / tau0_th * 100.0
    print("\n  tau0 numerique = %.4f" % tau0_num)
    print("  tau0 theorique = %.4f" % tau0_th)
    print("  ecart relatif  = %.2f %%" % ecart)
    return tau0_num, tau0_th, res


def calculer_taut(V0, a, dx=0.05, dt=0.0005):
    print("\n-- 4.1.c : taut avec barriere (V0=%.1f, a=%.1f) --" % (V0, a))
    res = simuler(V0=V0, a=a, t_max=30.0, dx=dx, dt=dt)
    print("\n  taut numerique = %.4f" % res["tau"])
    return res["tau"], res


def etude_a(V0, liste_a, dx=0.05, dt=0.0005):
    print("\n-- 4.1.d : influence de a (V0=%.1f) --" % V0)
    tau0_num, tau0_th, taut_num = [], [], []

    for a in liste_a:
        print("\n  -> a = %.1f" % a)
        r0 = simuler(V0=0.0, a=a, t_max=20.0 + a/utils.VG, dx=dx, dt=dt)
        tau0_num.append(r0["tau"])
        tau0_th.append(a / utils.VG)

        rt = simuler(V0=V0, a=a, t_max=30.0 + a/utils.VG, dx=dx, dt=dt)
        taut_num.append(rt["tau"])

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle("Influence de la largeur a (V0=%.1f)" % V0, fontsize=13)
    ax.plot(liste_a, tau0_th, "k--", lw=1.5, label="tau0 theorique (a/v_g)")
    ax.plot(liste_a, tau0_num, "ko-", lw=1.5, ms=6, label="tau0 numerique")
    ax.plot(liste_a, taut_num, "rs-", lw=1.5, ms=6, label="taut numerique")
    ax.set_xlabel("Largeur a")
    ax.set_ylabel("Temps de traversee")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.4)
    plt.savefig("graphs/influence_a.png", dpi=150)
    plt.show()

    return tau0_num, tau0_th, taut_num


def etude_V0(a, liste_V0, dx=0.05, dt=0.0005):
    print("\n-- 4.1.e : influence de V0 (a=%.1f) --" % a)
    taut_num, T_vals = [], []

    for V0 in liste_V0:
        print("\n  -> V0 = %.1f" % V0)
        rt = simuler(V0=V0, a=a, t_max=30.0, dx=dx, dt=dt)
        taut_num.append(rt["tau"])
        T = utils.coef_transmission(utils.E0, V0, a)
        T_vals.append(T)
        print("    T = %.6e" % T)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Influence de la hauteur V0 (a=%.1f)" % a, fontsize=13)

    ax1.plot(liste_V0, taut_num, "rs-", lw=1.5, ms=6)
    ax1.axhline(a / utils.VG, color="gray", ls="--", lw=1,
                label="tau0 = %.2f" % (a / utils.VG))
    ax1.axvline(utils.E0, color="blue", ls=":", lw=1.2,
                label="E0 = %.1f" % utils.E0)
    ax1.set_ylabel("taut (temps tunnel)")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.4)

    ax2.plot(liste_V0, T_vals, "b^-", lw=1.5, ms=6)
    ax2.axvline(utils.E0, color="blue", ls=":", lw=1.2)
    ax2.set_ylabel("Coefficient T")
    ax2.set_xlabel("Hauteur V0")
    ax2.set_yscale("log")
    ax2.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig("graphs/influence_V0.png", dpi=150)
    plt.show()

    return taut_num, T_vals


def visu_tunnel(V0, a):
    print("\n-- Visualisation tunnel (V0=%.1f, a=%.1f) --" % (V0, a))
    res = simuler(V0=V0, a=a, t_max=25.0, dx=0.04, dt=0.0004)

    x = res["x"]
    V = res["V"]
    t_stock = res["t_stock"]
    hist_dens = res["hist_dens"]
    nt = len(t_stock)

    indices = [0, int(0.25*nt), int(0.40*nt),
               int(0.50*nt), int(0.60*nt), nt-1]
    indices = [min(i, nt-1) for i in indices]

    fig, axes = plt.subplots(len(indices), 1, figsize=(12, 10), sharex=True)
    fig.suptitle("Traversee par effet tunnel — V0=%.1f, a=%.1f (E0=%.1f)" %
                 (V0, a, utils.E0), fontsize=14)
    couleurs = plt.cm.viridis(np.linspace(0.2, 0.9, len(indices)))

    for ax, idx, c in zip(axes, indices, couleurs):
        t_val = t_stock[idx]
        dens = hist_dens[idx]
        Vmax = np.max(V)
        if Vmax > 0:
            Vaff = V / Vmax * np.max(dens) * 1.1
        else:
            Vaff = np.zeros(len(x))

        ax.fill_between(x, 0, Vaff, color="gray", alpha=0.3,
                        label="Barriere" if idx == indices[0] else "")
        ax.plot(x, dens, color=c, lw=1.5)
        ax.set_ylabel("|Psi|^2", fontsize=9)
        ax.text(0.98, 0.85, "t = %.1f" % t_val, transform=ax.transAxes,
                fontsize=9, ha="right", color=c, fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        if idx == indices[0]:
            ax.legend(fontsize=8, loc="upper left")

    axes[-1].set_xlabel("x")
    plt.tight_layout()
    plt.savefig("graphs/instantanés_tunnel.png", dpi=150)
    plt.show()


def analyse_transmission(E, V0, a):
    T = utils.coef_transmission(E, V0, a)
    print("\n-- Analyse transmission --")
    print("  E=%.1f, V0=%.1f, a=%.1f" % (E, V0, a))
    if E >= V0:
        print("  regime E >= V0")
    else:
        kappa = np.sqrt(2.0 * utils.MASS * (V0 - E)) / utils.HBAR
        print("  regime E < V0 (tunnel), kappa=%.4f, kappa*a=%.4f" % (kappa, kappa*a))
    print("  T = %.6e" % T)
    return T


if __name__ == "__main__":
    import os
    os.makedirs("graphs")
    print("=== Partie 4 : Effet tunnel ===")

    V0_def = 5.0
    a_def = 5.0

    analyse_transmission(utils.E0, V0_def, a_def)

    tau0_num, tau0_th, _ = calculer_tau0(a=a_def)
    taut_num, _ = calculer_taut(V0=V0_def, a=a_def)

    print("\n-- Comparaison tau0 / taut --")
    print("  tau0 = %.4f,  taut = %.4f" % (tau0_num, taut_num))
    print("  rapport taut/tau0 = %.4f" % (taut_num / tau0_num))

    visu_tunnel(V0=V0_def, a=a_def)

    etude_a(V0=V0_def, liste_a=[2, 3, 4, 5, 6, 8])
    etude_V0(a=a_def, liste_V0=[1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])

    print("\nTermine.")
