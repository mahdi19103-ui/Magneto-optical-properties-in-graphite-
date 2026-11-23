import numpy as np
import matplotlib.pyplot as plt
from fonctions.Pi_pm import Pi_pm


def A_mp(B, eta_cm, eta_ph_cm, q_values, omega_cm_values, para_swm, para_num,
         para_phy, unites_phy,interacting=True):

    # Conversion des frequences et des elargissements en eV
    omega_ev = np.asarray(omega_cm_values) * unites_phy.cm_to_eV

    omega_0 = para_phy.omega_0_cm * unites_phy.cm_to_eV

    eta_ph = eta_ph_cm * unites_phy.cm_to_eV

    # Calcul de Pi_-+(omega)
    Pi_pm_values = Pi_pm(
        B=B,
        eta_cm=eta_cm,
        q_values=q_values,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        para_phy=para_phy,
        unites_phy=unites_phy
    )

    # Fonction de Green du phonon libre D0_+-
    D0_mp_values = (
        2 * omega_0
        / ((omega_ev + 1j * eta_ph)**2 - omega_0**2)
    )

    # Fonction de Green du phonon avec interaction electronique
    if interacting:
        D_mp_values = (
        D0_mp_values
        / (1 - D0_mp_values * Pi_pm_values)
    )
    else:
        D_mp_values= D0_mp_values

    # Fonction spectrale A_+-
    A_mp_values = -1 / np.pi * np.imag(D_mp_values)

    return A_mp_values



def plot_A_mp(B, eta_cm, eta_ph_cm, q_values, omega_cm_values, para_swm, para_num,
                para_phy, unites_phy, x_lim=None, y_lim=None, pas_x=10):

    # Calcul de la fonction spectrale
    A_mp_values = A_mp(
        B=B,
        eta_cm=eta_cm,
        eta_ph_cm=eta_ph_cm,
        q_values=q_values,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        para_phy=para_phy,
        unites_phy=unites_phy
    )

    # Tracé de la fonction spectrale
    plt.figure(figsize=(9, 5))

    plt.plot(
        omega_cm_values,
        A_mp_values,
        color="blue",
        label=r"$A_{-+}(\omega)$"
    )

    plt.axvline(
        x=para_phy.omega_0_cm,
        color="red",
        linestyle="-",
        linewidth=1.5,
        label=rf"$\omega_0 = {para_phy.omega_0_cm}\ \mathrm{{cm}}^{{-1}}$"
    )

    plt.xlabel(r"$\omega$ (cm$^{-1}$)")
    plt.ylabel(r"$A_{-+}(\omega)$")

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

    return A_mp_values