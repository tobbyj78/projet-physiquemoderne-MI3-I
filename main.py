# main.py — lanceur du projet
# Projet tunnel — MI3

import os
import utils
import OndePLane1d3I
import PaquetOndeGauss1d3I as p2
import schrodinger_solver as s3
import tunnel_effect as t4

os.makedirs("graphs")


def partie1():
    print("\n--- Partie 1 : Ondes planes ---")
    OndePLane1d3I.plot_onde_unique()
    OndePLane1d3I.plot_superposition()


def partie2():
    print("\n--- Partie 2 : Paquet gaussien ---")
    p2.verif_norme()
    p2.plot_paquet_t0()
    p2.plot_dispersion()


def partie3():
    print("\n--- Partie 3 : Solveur Schrodinger ---")
    s3.valider_derivees()
    s3.valider_solveur()


def partie4():
    print("\n--- Partie 4 : Effet tunnel ---")
    V0_def = 5.0
    a_def = 5.0

    t4.analyse_transmission(utils.E0, V0_def, a_def)

    tau0, tau0_th, _ = t4.calculer_tau0(a=a_def)
    taut, _ = t4.calculer_taut(V0=V0_def, a=a_def)

    print("\n-- Comparaison --")
    print("  tau0 = %.4f,  taut = %.4f,  taut/tau0 = %.4f" %
          (tau0, taut, taut/tau0))

    t4.visu_tunnel(V0=V0_def, a=a_def)
    t4.etude_a(V0=V0_def, liste_a=[2, 3, 4, 5, 6, 8])
    t4.etude_V0(a=a_def, liste_V0=[1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0])


if __name__ == "__main__":
    print()
    print("  Projet : Temps de traversee par effet tunnel")
    print("  Mondir Boughrous & Tom Jaffrain — MI3_projet_I")
    print()
    print("  1. Partie 1 — Ondes planes")
    print("  2. Partie 2 — Paquet d'ondes gaussien")
    print("  3. Partie 3 — Solveur Schrodinger")
    print("  4. Partie 4 — Effet tunnel")
    print("  5. TOUT lancer")
    print("  0. Quitter")
    print()

    try:
        c = input("Choix : ").strip()
    except (EOFError, KeyboardInterrupt):
        c = "0"

    if c == "1":
        partie1()
    elif c == "2":
        partie2()
    elif c == "3":
        partie3()
    elif c == "4":
        partie4()
    elif c == "5":
        utils.info()
        partie1()
        partie2()
        partie3()
        partie4()
        print("\n=== Termine. Graphiques dans graphs/ ===")
    else:
        print("Au revoir.")
