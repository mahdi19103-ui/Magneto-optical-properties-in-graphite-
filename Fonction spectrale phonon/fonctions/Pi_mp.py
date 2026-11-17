
import numpy as np
from fonctions.P_n_mp import P_n_mp


def Pi_mp(B, eta_cm, valeurs_q, omega_cm_values,
          para_swm, para_num, unites_phy):

    n_max = para_num.n_max

    Pi_mp_values = np.zeros(len(omega_cm_values), dtype=complex)

    for n in range(0, n_max - 1):

        resultat = P_n_mp(
            n=n,
            B=B,
            eta_cm=eta_cm,
            valeurs_q=valeurs_q,
            omega_cm_values=omega_cm_values,
            para_swm=para_swm,
            para_num=para_num,
            unites_phy=unites_phy
        )

        P_n_mp_values = resultat[0]

        Pi_mp_values += P_n_mp_values

    return Pi_mp_values
