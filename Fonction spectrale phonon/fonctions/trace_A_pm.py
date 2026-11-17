
import numpy as np
import matplotlib.pyplot as plt

from fonctions.A_pm import A_pm


def tracer_A_pm(B, eta_cm, eta_ph_cm, valeurs_q, omega_cm_values,
                param_phys, para_swm, para_num, unites_phy, x_lim=None, y_lim=None, pas_x=10):

    # Calcul de la fonction spectrale
    A_pm_values = A_pm(
        B=B,
        eta_cm=eta_cm,
        eta_ph_cm=eta_ph_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        param_phys=param_phys,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    # Tracé de la fonction spectrale
    plt.figure(figsize=(9, 5))

    plt.plot(
        omega_cm_values,
        A_pm_values,
        color="blue",
        label=r"$A_{+-}(\omega)$"
    )

    plt.axvline(
        x=param_phys.omega_0_cm,
        color="red",
        linestyle="-",
        linewidth=1.5,
        label=rf"$\omega_0 = {param_phys.omega_0_cm}\ \mathrm{{cm}}^{{-1}}$"
    )

    plt.xlabel(r"$\omega$ (cm$^{-1}$)")
    plt.ylabel(r"$A_{+-}(\omega)$")

    if x_lim is not None:
        plt.xlim(x_lim)

    if y_lim is not None:
        plt.ylim(y_lim)

    if pas_x is not None:
        plt.xticks(np.arange(x_lim[0], x_lim[1] + pas_x, pas_x))

    plt.grid(True, alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.show()
   
    return A_pm_values
