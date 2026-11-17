
import numpy as np
from fonctions.hamiltonien import energies_seules

#Fonction qui cree un tableau de taille(nombre totale des etats). 
#On met dans ce tableau tout les energies propres par ordre croissant.
#but: calcule enrgie de Fermi (milieu de ce tableau).
def tableaux_energies(B, para_swm, para_num): 

    valeurs_q = np.linspace(
        -0.5, 0.5, para_num.Nq
    )

    toutes_energies = []

    for n in range(para_num.n_max + 1):

        for q in valeurs_q:

            energies = energies_seules(
                n, B, q, para_swm
            )

            toutes_energies.extend(energies)

    toutes_energies = np.sort(
        np.array(toutes_energies)
    )

    return toutes_energies

# fonctionqui donne l'energie de fermi
def energie_fermi(B, para_swm, para_num):

    energies_totales = tableaux_energies(
        B, para_swm, para_num
    )

    N = len(energies_totales)

    if N < 2:
        raise ValueError(
            "Il faut au moins deux niveaux d'énergie."
        )

    EF = 0.5 * (
        energies_totales[N // 2 - 1]
        + energies_totales[N // 2]
    )

    return EF
