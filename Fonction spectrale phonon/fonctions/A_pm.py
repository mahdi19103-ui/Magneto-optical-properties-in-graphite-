
import numpy as np
import matplotlib.pyplot as plt
from fonctions.Pi_mp import Pi_mp


def A_pm(B, eta_cm, eta_ph_cm, q_values, omega_cm_values,
        para_swm, para_num,para_phy, unites_phy, interacting=True):




    """
Calculate the spectral function A_{+-} from the phonon Green's function.

The spectral function is obtained from the free phonon Green's function
and the electronic polarization Pi_{-+}. The electronic interaction
can be included or neglected through the interacting parameter.

Returns:
    A_pm_values: Spectral function as a function of omega.
    """

    # Conversion des frequences et des elargissements en eV
    omega_ev = np.asarray(omega_cm_values) * unites_phy.cm_to_eV

    omega_0 = para_phy.omega_0_cm * unites_phy.cm_to_eV

    eta_ph = eta_ph_cm * unites_phy.cm_to_eV

    # Calcul de Pi_-+(omega)
    Pi_mp_values = Pi_mp(
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
    D0_pm_values = (
        2 * omega_0
        / ((omega_ev + 1j * eta_ph)**2 - omega_0**2)
    )

    # Fonction de Green du phonon avec interaction electronique
    if interacting:
        D_pm_values = (
        D0_pm_values
        / (1 - D0_pm_values * Pi_mp_values)
                    )
    else:
        D_pm_values= D0_pm_values

    # Fonction spectrale A_+-
    A_pm_values = -1 / np.pi * np.imag(D_pm_values)

    return A_pm_values



def plot_A_pm(B, eta_cm, eta_ph_cm, q_values, omega_cm_values,
                 para_swm, para_num,para_phy, unites_phy, x_lim=None, y_lim=None, pas_x=10):

    """
Plot the spectral function A_{+-} as a function of omega
for a fixed magnetic field B.

The plot includes the phonon frequency omega_0 as a reference
and allows the axis limits to be adjusted.

Returns:
    A_pm_values: Spectral function as a function of omega.
    """

    # Calcul de la fonction spectrale
    A_pm_values = A_pm(
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
        A_pm_values,
        color="blue",
        label=r"$A_{+-}(\omega)$"
    )

    plt.axvline(
        x=para_phy.omega_0_cm,
        color="red",
        linestyle="-",
        linewidth=1.5,
        label=rf"$\omega_0 = {para_phy.omega_0_cm}\ \mathrm{{cm}}^{{-1}}$"
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
