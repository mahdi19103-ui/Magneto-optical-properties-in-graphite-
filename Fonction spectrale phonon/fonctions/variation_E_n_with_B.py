
import numpy as np
import matplotlib.pyplot as plt

from fonctions.hamiltonien import energies_seules


def tracer_niveaux_Landau(
    q,
    para_swm,para_phy,
    n_max=10,
    B_min=0,
    B_max=10,
    pas_B=2,
    E_min=None,
    E_max=None
):

    # Valeurs du champ magnétique
    B_values = np.arange(B_min, B_max + pas_B, pas_B)

    # Créer la figure
    fig, ax = plt.subplots(figsize=(3, 3))

    # Calculer les niveaux pour chaque n
    for n in range(n_max + 1):

        energies_n = []

        for B in B_values:

            energies = energies_seules(
                n=n,
                B=B,
                q=q,
                para_swm=para_swm,para_phy=para_phy
            )

            energies_n.append(energies)

        energies_n = np.array(energies_n)

        # Couleurs des branches de basse énergie
        if n == 0:
            couleurs = ["magenta"]

        elif n == 1:
            couleurs = ["blue", "magenta", "red"]

        else:
            couleurs = ["blue", "red"]

        # Tracer uniquement les branches de basse énergie
        if n >= 2:
            energies_basse_energie = energies_n[:, 1:3]
        else:
            energies_basse_energie = energies_n

        for j in range(energies_basse_energie.shape[1]):

            ax.plot(
                B_values,
                energies_basse_energie[:, j],
                color=couleurs[j],
                linewidth=1.5
            )

    # Axes
    ax.set_xlabel(r"$B\ (\mathrm{T})$", fontsize=16)
    ax.set_ylabel(r"$E\ (\mathrm{eV})$", fontsize=18)

    """# Placer l'axe E à droite
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")"""

    ax.tick_params(axis='y', labelsize=16)
    # Limites des axes
    ax.set_xlim(B_min, B_max)

    if E_min is not None and E_max is not None:
        ax.set_ylim(E_min, E_max)

# Grille
    ax.grid(True, alpha=0.3)

