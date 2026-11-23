"""
Objectives:

Build an array containing all eigenenergies of the system.

The energies are collected for all Landau levels n and all q values,
then sorted in ascending order.

The Fermi energy is determined from the middle of the sorted energy array.
 """


import numpy as np
from fonctions.hamiltonien import energies_seules


def tableaux_energies(B, para_swm, para_num,para_phy): 

    valeurs_q = np.linspace(
        -0.5, 0.5, para_num.Nq
    )

    toutes_energies = []

    for n in range(para_num.n_max + 1):

        for q in valeurs_q:

            energies = energies_seules(
                n, B, q, para_swm,para_phy
            )

            toutes_energies.extend(energies)

    toutes_energies = np.sort(
        np.array(toutes_energies)
    )

    return toutes_energies

# fonction qui donne l'energie de fermi
def energie_fermi(B, para_swm, para_num,para_phy):

    energies_totales = tableaux_energies(
        B, para_swm, para_num,para_phy
    )

    N = len(energies_totales)

    if N < 2:
        raise ValueError(
            "Il faut au moins deux niveaux d'energie."
        )

    if N % 2 == 0:
        
    # Even number of energy levels
        EF = 0.5 * (
        energies_totales[N // 2 - 1]
        + energies_totales[N // 2]
    )
        
    else:
    # Odd number of energy levels
        EF = energies_totales[N // 2]

    return EF
