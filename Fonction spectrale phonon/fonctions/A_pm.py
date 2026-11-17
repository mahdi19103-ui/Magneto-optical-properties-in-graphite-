
import numpy as np
from fonctions.Pi_mp import Pi_mp


def A_pm(B, eta_cm, eta_ph_cm, valeurs_q, omega_cm_values,
         param_phys, para_swm, para_num, unites_phy, interacting=True):

    # Conversion des frequences et des elargissements en eV
    omega_ev = np.asarray(omega_cm_values) * unites_phy.cm1_to_eV

    omega_0 = param_phys.omega_0_cm * unites_phy.cm1_to_eV

    eta_ph = eta_ph_cm * unites_phy.cm1_to_eV

    # Calcul de Pi_-+(omega)
    Pi_mp_values = Pi_mp(
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
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
