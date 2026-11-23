"""
Objectives:
Calculate the possible electronic transitions from n (initial state) to n+2
(final state) as a function of the magnetic field B.


Here, gamma3 is taken into account. Otherwise, gamma4 would lead to transitions
between n and n+1 or n-1.


The possible transitions will then be compared with the peaks observed in the
polarization operator at fixed magnetic fields in order to identify the
corresponding transitions.

Remark:
For n -> n+2 transitions, the initial Landau level starts at n = 0.
"""





import numpy as np
import matplotlib.pyplot as plt

from fonctions.hamiltonien import energies_seules
from fonctions.fermi_energy import tableaux_energies
from fonctions.parameters import Unitesphy


def nombre_branches(n):
    if n == 0:
        return 1
    elif n == 1:
        return 3
    else:
        return 4


def calcul_delta_E_n_nplus2(B_values, para_swm, para_num, para_phy):
    """
        Calculate the transition energies from Landau level n to n-2
        as a function of the magnetic field B.
    
        The function computes the energy branches for each Landau level,
        determines the Fermi energy for each magnetic field, and identifies
        the allowed transitions from occupied initial states to unoccupied
        final states.
    
        Returns:
            transitions: Transition energies for each pair of branches,
                        as a function of B and q.
            valeurs_q: The sampled q values.
        """
    valeurs_q = np.linspace(-0.5, 0.5, para_num.Nq)

    # Calculate energies
    E2 = {}

    for n in range(para_num.n_max + 1):
        for lambd1 in range(nombre_branches(n)):
            E2[n, lambd1] = np.zeros(
                (para_num.N_B, para_num.Nq)
            )

    for i, B in enumerate(B_values):
        for j, q in enumerate(valeurs_q):
            for n in range(para_num.n_max + 1):

                energies = energies_seules(n, B, q, para_swm, para_phy)

                for lambd1 in range(len(energies)):
                    E2[n, lambd1][i, j] = energies[lambd1]

    # Calculate Fermi energy
    EF = np.zeros(para_num.N_B)

    for i, B in enumerate(B_values):

        energies_totales = tableaux_energies(
            B, para_swm, para_num, para_phy
        )

        N = len(energies_totales)

        EF[i] = 0.5 * (
            energies_totales[N // 2 - 1]
            + energies_totales[N // 2]
        )

    # Calculate transitions n -> n+2
    transitions2 = {}

    for n in range(0, para_num.n_max - 1):
        for lambd1 in range(nombre_branches(n)):
            for lambd2 in range(nombre_branches(n + 2)):

                E_initial = E2[n, lambd1]
                E_final = E2[n + 2, lambd2]

                difference = E_final - E_initial

                condition = (
                    (E_initial < EF[:, None])
                    & (E_final > EF[:, None])
                )

                transitions2[n, lambd1, lambd2] = np.where(
                    condition, difference, 0
                )

    return transitions2, valeurs_q


def plot_delta_E_n_nplus2(B_values, para_swm, para_num, para_phy):
    """
         Plot the electronic transition energies from n to n-2 as a function
         of the magnetic field B.
    
        The function represents the allowed transitions for each Landau level
        and each pair of energy branches. It also displays the reference
        phonon energy omega_0 to identify possible resonances.
        """
    transitions2, valeurs_q = calcul_delta_E_n_nplus2(
        B_values, para_swm, para_num, para_phy
    )

    unites_phy = Unitesphy()

    noms_branches = {
        0: ["0"],
        1: ["--", "0", "++"],
    }

    for n in range(2, para_num.n_max + 1):
        noms_branches[n] = ["--", "-", "+", "++"]

    # Transitions n -> n+2
    for n in range(0, para_num.n_max - 1):

        plt.figure(figsize=(10, 6))

        nombre_transitions = (
            nombre_branches(n) * nombre_branches(n + 2)
        )

        couleurs = plt.cm.tab10(
            np.linspace(0, 1, nombre_transitions)
        )

        indice_couleur = 0

        for lambd1 in range(nombre_branches(n)):
            for lambd2 in range(nombre_branches(n + 2)):

                matrice = transitions2[n, lambd1, lambd2]
                couleur = couleurs[indice_couleur]
                transition_visible = False

                for j in range(len(valeurs_q)):

                    energie_cm = (
                        matrice[:, j] * unites_phy.eV_to_cm1
                    )

                    energie_cm = np.where(
                        (matrice[:, j] > 0)
                        & (energie_cm > 0)
                        & (energie_cm <= 250),
                        energie_cm,
                        np.nan
                    )

                    if np.any(np.isfinite(energie_cm)):
                        transition_visible = True

                        plt.plot(
                            B_values,
                            energie_cm,
                            color=couleur,
                            alpha=0.4,
                            linewidth=1
                        )

                if transition_visible:
                    etiquette = (
                        f"({n}, {noms_branches[n][lambd1]}) "
                        f"→ ({n+2}, {noms_branches[n+2][lambd2]})"
                    )

                    plt.plot(
                        [],
                        [],
                        color=couleur,
                        linewidth=2,
                        label=etiquette
                    )

                indice_couleur += 1

        plt.axhline(
            y=44,
            color="green",
            linestyle="-",
            linewidth=1.5,
            label=r"$\omega_0 = 44\, \mathrm{cm}^{-1}$"
        )

        plt.xlabel("Magnetic field B (T)")
        plt.ylabel("Transition energy (cm⁻¹)")
        plt.xlim(0, 1.8)
        plt.xticks(np.arange(0, 1.81, 0.1))
        plt.ylim(0, 111)
        plt.yticks(np.arange(0, 112, 5))
        plt.title(f"Electronic transitions {n} → {n+2}")
        plt.grid(True, alpha=0.25)
        plt.legend(fontsize=9, loc="best")
        plt.tight_layout()
        plt.show()
