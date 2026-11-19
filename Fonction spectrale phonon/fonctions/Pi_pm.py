
import numpy as np
from fonctions.p_n_pm import P_n_pm


def Pi_pm(
    B, eta_cm, valeurs_q, omega_cm_values,
    para_swm, para_num, para_phy, unites_phy
):

    n_max = para_num.n_max

    # Initialiser la polarisation totale
    Pi_pm_values = np.zeros(
        len(omega_cm_values),
        dtype=complex
    )

    # Sommer les contributions de n=2 à n_max
    for n in range(2, n_max + 1):

        resultat = P_n_pm(
            n=n,
            B=B,
            eta_cm=eta_cm,
            valeurs_q=valeurs_q,
            omega_cm_values=omega_cm_values,
            para_swm=para_swm,
            para_num=para_num,
            para_phy=para_phy,
            unites_phy=unites_phy
        )

        P_n_pm_values = resultat[0]

        Pi_pm_values += P_n_pm_values

    return Pi_pm_values
