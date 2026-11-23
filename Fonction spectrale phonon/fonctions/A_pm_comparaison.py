
import numpy as np
import matplotlib.pyplot as plt

from fonctions.A_pm import A_pm


def plot_A_pm_comparison(
    B_values,
    eta_cm,
    eta_ph_cm,
    q_values,
    omega_cm_values,
    para_swm,
    para_num,
    para_phy,
    unites_phy,
    x_lim=(0, 80),
    y_lim=None,
    pas_x=20,
    interacting=True
):

    # Create the figure
    plt.figure(figsize=(8, 5))

    # Colors for the different curves
    couleurs = ["red", "blue", "green"]

    # Calculate and plot for each magnetic field
    for B, couleur in zip(B_values, couleurs):

        A_pm_values = A_pm(
            B=B,
            eta_cm=eta_cm,
            eta_ph_cm=eta_ph_cm,
            q_values=q_values,
            omega_cm_values=omega_cm_values,
            para_swm=para_swm,
            para_num=para_num,
            para_phy=para_phy,
            unites_phy=unites_phy,
            interacting=interacting
        )

        plt.plot(
            omega_cm_values,
            A_pm_values,
            color=couleur,
            linewidth=2,
            label=rf"$B={B}$ T"
        )

    # Reference vertical line at 44 cm^-1
    plt.axvline(
        x=44,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label=r"$\omega_0=44\ \mathrm{cm}^{-1}$"
    )

    # Formatting
    plt.xlabel(r"$\omega$ (cm$^{-1}$)", fontsize=20)
    plt.ylabel(r"$A_{+-}(\omega)$", fontsize=18)

    plt.xlim(x_lim)

    plt.xticks(
        np.arange(x_lim[0], x_lim[1] + pas_x, pas_x),
        fontsize=18
    )

    if y_lim is not None:
        plt.ylim(y_lim)

    plt.legend(fontsize=14)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

