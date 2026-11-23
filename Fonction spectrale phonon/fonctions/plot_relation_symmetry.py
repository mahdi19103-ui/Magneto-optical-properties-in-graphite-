
import numpy as np
import matplotlib.pyplot as plt

from fonctions.Pi_pm import Pi_pm
from fonctions.Pi_mp import Pi_mp


def plot_Pi_symmetry(
    B, q_values, omega_cm_values,eta_cm,unites_phy,
    para_swm, para_num, para_phy,
    pas_x=10, pas_y=None,
    x_lim=(10, 200), y_lim=None
):

    """
    Plot the imaginary parts of the two polarizations Pi_{+-}
    and Pi_{-+} for a fixed magnetic field B.

    The plot compares -Im(Pi_{+-}(omega)) and
    Im(Pi_{-+}(-omega)) and displays their maximum difference.

    Returns:
        Im_pm: Imaginary part of -Pi_{+-}.
        Im_mp: Imaginary part of Pi_{-+}.
        difference: Difference between the two curves.
    """

    # 1. Calculate -Pi_{+-}(omega)
    Pi_pm_values = -Pi_pm(
        B=B,
        q_values=q_values,eta_cm=eta_cm,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        para_phy=para_phy,unites_phy=unites_phy
    )

    # 2. Calculate Pi_{-+}(-omega)
    Pi_mp_values = Pi_mp(
        B=B,
        q_values=q_values,
        omega_cm_values=-omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        para_phy=para_phy,
        eta_cm=eta_cm,unites_phy=unites_phy
    )

    # 3. Extract imaginary parts
    Im_pm = np.imag(Pi_pm_values)
    Im_mp = np.imag(Pi_mp_values)

    # 4. Calculate the difference
    difference = Im_pm - Im_mp

    diff_max = np.max(np.abs(difference))

    print("Maximum difference =", diff_max)

    # 5. Plot both functions
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

    plt.title(rf"$B={B}$ T")

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

    return Im_pm, Im_mp, difference
