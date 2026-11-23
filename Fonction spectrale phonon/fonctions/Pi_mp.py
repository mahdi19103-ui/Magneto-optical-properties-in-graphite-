
import numpy as np
import matplotlib.pyplot as plt
from fonctions.p_n_mp import p_n_mp


def Pi_mp(B, eta_cm, q_values, omega_cm_values,
          para_swm, para_num,para_phy, unites_phy):

    """
    Calculate the total polarization Pi_{-+} by summing the contributions
    P_n(omega) over Landau levels from n=0 to n_max.
    
    The total polarization is obtained by adding the individual contributions
    for each Landau level at a fixed magnetic field B.
    
    Returns:
        Pi_pm_values: Total polarization as a function of omega.
    """


    n_max = para_num.n_max

    Pi_mp_values = np.zeros(len(omega_cm_values), dtype=complex)

    for n in range(0, n_max - 1):

        resultat = p_n_mp(
            n=n,
            B=B,
            eta_cm=eta_cm,
            q_values=q_values,
            omega_cm_values=omega_cm_values,
            para_swm=para_swm,
            para_num=para_num,
            para_phy=para_phy,
            unites_phy=unites_phy
        )

        p_n_mp_values = resultat[0]

        Pi_mp_values += p_n_mp_values

    return Pi_mp_values


def plot_Pi_mp(
    B, eta_cm, q_values, omega_cm_values,
    para_swm, para_num,para_phy, unites_phy,
    pas_x=10, pas_y=None,
    x_lim=(10, 200), y_lim=None
):


    
    """
        Plot the imaginary part of the total polarization Pi_{-+}
        as a function of omega for a fixed magnetic field B.
    
        The total polarization is calculated by summing the contributions
        from all Landau levels before plotting.
    
        The plot settings, including axis limits and tick spacing,
         can be adjusted through the function arguments.
    
        Returns:
        Pi_values: Total polarization as a function of omega.
    """

    # 1. Calculer la polarisation totale
    Pi_values = Pi_mp(
        B=B,
        eta_cm=eta_cm,
        q_values=q_values,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        para_phy=para_phy,
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