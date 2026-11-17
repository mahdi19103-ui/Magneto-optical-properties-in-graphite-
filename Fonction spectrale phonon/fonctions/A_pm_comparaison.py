
import numpy as np
import matplotlib.pyplot as plt

from fonctions.A_pm import A_pm


def trace_A_pm_comparaison(
    B_values,
    eta_cm,
    eta_ph_cm,
    valeurs_q,
    omega_cm_values,
    param_phys,
    para_swm,
    para_num,
    unites_phy,
    x_lim=(0, 80),
    y_lim=None,
    pas_x=20,
    interacting=True
):

    # Création de la figure
    plt.figure(figsize=(8, 5))

    # Couleurs pour les différentes courbes
    couleurs = ["red", "blue", "green"]

    # Calcul et tracé pour chaque champ magnétique
    for B, couleur in zip(B_values, couleurs):

        A_values = A_pm(
            B=B,
            eta_cm=eta_cm,
            eta_ph_cm=eta_ph_cm,
            valeurs_q=valeurs_q,
            omega_cm_values=omega_cm_values,
            param_phys=param_phys,
            para_swm=para_swm,
            para_num=para_num,
            unites_phy=unites_phy,
            interacting=interacting
        )

        plt.plot(
            omega_cm_values,
            A_values,
            color=couleur,
            linewidth=2,
            label=rf"$B={B}$ T"
       
        )


    # Ligne verticale à 44 cm⁻¹
    plt.axvline(
    x=44,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label=r"$\omega_0=44\ \mathrm{cm}^{-1}$"
    )

    # Mise en forme
    plt.xlabel(r"$\omega$ (cm$^{-1}$)", fontsize=20)
    plt.ylabel(r"$A_{+-}(\omega)$", fontsize=18)

    plt.xlim(x_lim)

    plt.xticks(
        np.arange(x_lim[0], x_lim[1] + pas_x, pas_x),fontsize=18
    )

    if y_lim is not None:
        plt.ylim(y_lim)

    plt.legend(fontsize=14)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
