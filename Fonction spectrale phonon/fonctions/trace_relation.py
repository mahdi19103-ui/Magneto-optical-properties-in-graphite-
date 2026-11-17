
import numpy as np
import matplotlib.pyplot as plt

from fonctions.Pi_pm import Pi_pm
from fonctions.Pi_mp import Pi_mp


def trace_Pi(
    B, eta_cm, valeurs_q, omega_cm_values,
    para_swm, para_num, unites_phy,
    pas_x=10, pas_y=None,
    x_lim=(10, 200), y_lim=None
):

    # 1. Calculer -Pi_+-(omega)
    Pi_pm_values = -Pi_pm(
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    # 2. Calculer Pi_-+(-omega)
    Pi_mp_values = Pi_mp(
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=-omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    # 3. Extraire les parties imaginaires
    Im_pm = np.imag(Pi_pm_values)
    Im_mp = np.imag(Pi_mp_values)

    # 4. Calculer la différence
    difference = Im_pm - Im_mp

    diff_max = np.max(np.abs(difference))

    print("Différence maximale =", diff_max)

    # 5. Première figure : les deux fonctions
    plt.figure(figsize=(10, 5))

    plt.plot(
        omega_cm_values,
        Im_pm,
        color="red",
        linewidth=2,
        linestyle="-",
        label=r"$-\mathrm{Im}\,\Pi_{+-}(\omega)$"
    )

    plt.plot(
        omega_cm_values,
        Im_mp,
        color="blue",
        linewidth=1.8,
        linestyle="--",
        label=r"$\mathrm{Im}\,\Pi_{-+}(-\omega)$"
    )

    plt.xlabel(r"$\omega$ (cm$^{-1}$)", fontsize=20)

    plt.title(
        rf"$B={B}$ T"
    )

    plt.xlim(x_lim)
    plt.xticks(
        np.arange(x_lim[0], x_lim[1] + pas_x, pas_x)
    )

    if y_lim is not None:
        plt.ylim(y_lim)

    if pas_y is not None:
        y_min, y_max = plt.ylim()
        plt.yticks(
            np.arange(y_min, y_max + pas_y * 0.5, pas_y)
        )

    plt.legend(fontsize=13)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
