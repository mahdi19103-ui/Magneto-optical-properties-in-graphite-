
import numpy as np
from fonctions.A_pm import A_pm


def A_pm_B(
    eta_cm,
    eta_ph_cm,
    valeurs_q,
    omega_cm_values,
    param_phys,
    para_swm,
    para_num,
    unites_phy
):
    # Valeurs du champ magnetique entre 0 et 5 T
    B_values = np.linspace(0, 5, para_num.N_B)

    # Tableau pour stocker les fonctions spectrales
    A_values = np.zeros(
        (para_num.N_B, len(omega_cm_values)),
        dtype=float
    )

    # Calcul de A_pm pour chaque champ magnetique
    for i, B in enumerate(B_values):

        A_values[i, :] = A_pm(
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

    return B_values, A_values
