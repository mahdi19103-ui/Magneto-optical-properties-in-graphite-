
import numpy as np
import matplotlib.pyplot as plt

from fonctions.hamiltonien import energies_seules


def tracer_structure_bandes(
    B, para_num, para_swm,
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
            energies = energies_seules(n, B, q, para_swm)
            energies_n.append(energies)

        energies_n = np.array(energies_n)

        # -----------------------------
        # Couleurs des branches
        # -----------------------------

        for i in range(energies_n.shape[1]):

            # n = 0 : bande centrale
            if n == 0:
                couleur = "magenta"

            # n = 1 : trois branches [--, 0, ++]
            elif n == 1:

                if i == 0:
                    couleur = "black"    # Haute énergie négative

                elif i == 1:
                    couleur = "magenta"   # Bande centrale

                else:
                    couleur = "black"    # Haute énergie positive

            # n >= 2 : quatre branches [--, -, +, ++]
            else:

                if i == 0:
                    couleur = "black"    # Haute énergie négative

                elif i == 1:
                    couleur = "blue"     # Basse énergie négative

                elif i == 2:
                    couleur = "red"      # Basse énergie positive

                else:
                    couleur = "black"    # Haute énergie positive

            plt.plot(
                valeurs_q,
                energies_n[:, i],
                color=couleur
            )

    # Niveau de Fermi
    if Ef is not None:
        plt.axhline(
            y=Ef,
            color="green",
            linestyle="-"
        )

    # -----------------------------
    # Axes et présentation
    # -----------------------------

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



