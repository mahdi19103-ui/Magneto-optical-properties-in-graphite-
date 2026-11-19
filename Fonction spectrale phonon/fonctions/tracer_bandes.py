

import numpy as np
import matplotlib.pyplot as plt

from fonctions.hamiltonien import energies_seules


def tracer_structure_bandes(
    B, para_num, para_swm, para_phy,
    xlim=(-0.5, 0.5),
    ylim=(-0.5, 0.6),
    taille=(6, 5),
    pas_x=0.1,
    pas_y=0.2,
    Ef=None
):

    plt.figure(figsize=taille)

    valeurs_q = np.linspace(-0.5, 0.5, para_num.Nq)

    for n in range(para_num.n_max + 1):

        energies_n = []

        for q in valeurs_q:
            energies = energies_seules(n, B, q, para_swm, para_phy)
            energies_n.append(energies)

        energies_n = np.array(energies_n)

        # Couleurs des branches
        for i in range(energies_n.shape[1]):

            if n == 0:
                couleur = "magenta"

            elif n == 1:
                if i == 0:
                    couleur = "black"
                elif i == 1:
                    couleur = "magenta"
                else:
                    couleur = "black"

            else:
                if i == 0:
                    couleur = "black"
                elif i == 1:
                    couleur = "blue"
                elif i == 2:
                    couleur = "red"
                else:
                    couleur = "black"

            plt.plot(valeurs_q, energies_n[:, i], color=couleur)

    # Niveau de Fermi
    if Ef is not None:
        plt.axhline(y=Ef, color="green", linestyle="-")

    # Axes
    plt.xlabel(r'$k_z$', fontsize=19)
    plt.ylabel('Energy (eV)', fontsize=17)

    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.xticks(
        np.arange(xlim[0], xlim[1] + pas_x / 2, pas_x),
        fontsize=10
    )

    plt.yticks(
        np.arange(ylim[0], ylim[1] + pas_y / 2, pas_y),
        fontsize=10
    )

    plt.grid(False)
    plt.tight_layout()
    plt.show()




