
import numpy as np
import matplotlib.pyplot as plt
from fonctions.p_n_pm import p_n_pm, plot_p_n_pm


def Pi_pm(
    B, q_values, omega_cm_values,unites_phy,eta_cm,
    para_swm, para_num, para_phy
):

    
    """
Calculate the total polarization Pi_{+-} by summing the contributions
P_n(omega) over Landau levels from n=2 to n_max.

The total polarization is obtained by adding the individual contributions
for each Landau level at a fixed magnetic field B.

Returns:
    Pi_pm_values: Total polarization as a function of omega.
"""
    n_max = para_num.n_max

    # Initialiser la polarisation totale
    Pi_pm_values = np.zeros(
        len(omega_cm_values),
        dtype=complex
    )

    # Sommer les contributions de n=2 à n_max
    for n in range(2, n_max + 1):

        resultat = p_n_pm(
            n=n,
            B=B,
            q_values=q_values,
            omega_cm_values=omega_cm_values,
            para_swm=para_swm,
            para_num=para_num,
            para_phy=para_phy,
            unites_phy=unites_phy,
            eta_cm=eta_cm
        )

        p_n_pm_values = resultat[0]

        Pi_pm_values += p_n_pm_values

    return Pi_pm_values


def plot_Pi_pm(
    B,  q_values, omega_cm_values,
    para_swm, para_num, para_phy,unites_phy,eta_cm,
    pas_x=10, pas_y=None,
    x_lim=(10, 200), y_lim=None
):

    
    """
    Plot the imaginary part of the total polarization Pi_{+-}
    as a function of omega for a fixed magnetic field B.

    The total polarization is calculated by summing the contributions
    from all Landau levels before plotting.

    The plot settings, including axis limits and tick spacing,
     can be adjusted through the function arguments.

    Returns:
    Pi_values: Total polarization as a function of omega.
    """

    # 1. Calculate the total polarization
    Pi_values = Pi_pm(
        B=B,
        q_values=q_values,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        para_phy=para_phy, unites_phy=unites_phy,eta_cm=eta_cm
    )

    # 2. Create the figure
    plt.figure(figsize=(10, 5))

    # 3. Plot the imaginary part
    plt.plot(
        omega_cm_values,
        np.imag(Pi_values),
        color="red",
        linewidth=1.8
    )

    # 4. Formatting
    plt.xlabel(r"$\omega$ (cm$^{-1}$)", fontsize=20)
    plt.ylabel(
        r"$\mathrm{Im}\,\Pi_{+-}(\omega)$",
        fontsize=19
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