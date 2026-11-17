import numpy as np
import matplotlib.pyplot as plt

from fonctions.Pi_mp import Pi_mp


def trace_Pi_mp(
    B, eta_cm, valeurs_q, omega_cm_values,
    para_swm, para_num, unites_phy,
    pas_x=10, pas_y=None,
    x_lim=(10, 200), y_lim=None
):

    # 1. Calculer la polarisation totale
    Pi_values = Pi_mp(
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    # 2. Créer la figure
    plt.figure(figsize=(12, 6))

    # 3. Tracer la partie imaginaire
    plt.plot(
        omega_cm_values,
        np.imag(Pi_values),
        color="red",
        linewidth=1.8
    )

    # 4. Mise en forme
    plt.xlabel(r"$\omega$ (cm$^{-1}$)")
    plt.ylabel(r"$\mathrm{Im}\,\Pi_{-+}(\omega)$")

    plt.title(
        rf"$B={B}$ T, $\eta={eta_cm}$ cm$^{{-1}}$"
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

    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    return Pi_values