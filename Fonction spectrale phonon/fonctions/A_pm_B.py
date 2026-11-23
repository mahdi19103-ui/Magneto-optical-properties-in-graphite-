

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from fonctions.A_pm import A_pm


def A_pm_B(
    eta_cm,
    eta_ph_cm,
    q_values,
    omega_cm_values,
    para_swm,
    para_num,
    para_phy,
    unites_phy
):
    """
    Calculate the spectral function A_{+-} for a range of magnetic fields.

    The magnetic field values are generated between 0 and 5 T.
    The spectral function is calculated for each magnetic field
    and stored as a two-dimensional array.

    Returns:
        B_values: Magnetic field values.
        A_pm_values: Spectral function for each magnetic field
                      as a function of omega.
    """

    # Generate magnetic field values between 0 and 5 T
    B_values = np.linspace(0, 5, para_num.N_B)

    # Initialize the spectral function array
    A_pm_values = np.zeros(
        (para_num.N_B, len(omega_cm_values)),
        dtype=float
    )

    # Calculate A_pm for each magnetic field
    for i, B in enumerate(B_values):

        A_pm_values[i, :] = A_pm(
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

    return B_values, A_pm_values


def plot_A_pm_B(
    B_values,
    A_pm_values,
    omega_cm_values,
    inverse_axes=False,
    B_lim=None,
    omega_lim=None,
    vmax=None,
    pas_omega=10,
    cmap="RdBu_r"
):
    """
    Plot the spectral function A_{+-} as a color map
    depending on the magnetic field B and frequency omega.

    The axes can be inverted to display either B or omega
    along the horizontal axis. The color scale is symmetric
    around zero and can be adjusted through vmax.

    Returns:
        fig: Figure object.
        ax: Axes object.
    """

    B_values = np.asarray(B_values)
    A_pm_values = np.asarray(A_pm_values)
    omega_cm_values = np.asarray(omega_cm_values)

    # Color scale
    if vmax is None:
        vmax = np.max(np.abs(A_pm_values))

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Invert or keep the axes
    if inverse_axes:

        image = ax.imshow(
            A_pm_values.T,
            origin="lower",
            aspect="auto",
            extent=[
                B_values[0], B_values[-1],
                omega_cm_values[0], omega_cm_values[-1]
            ],
            cmap=cmap,
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest"
        )

        ax.set_xlabel(r"$B$ (T)")
        ax.set_ylabel(r"$\omega$ (cm$^{-1}$)")

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(pas_omega))

        if B_lim is not None:
            ax.set_xlim(B_lim)

        if omega_lim is not None:
            ax.set_ylim(omega_lim)

    else:

        image = ax.imshow(
            A_pm_values,
            origin="lower",
            aspect="auto",
            extent=[
                omega_cm_values[0], omega_cm_values[-1],
                B_values[0], B_values[-1]
            ],
            cmap=cmap,
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest"
        )

        ax.set_xlabel(r"$\omega$ (cm$^{-1}$)")
        ax.set_ylabel(r"$B$ (T)")

        ax.xaxis.set_major_locator(MultipleLocator(pas_omega))

        if omega_lim is not None:
            ax.set_xlim(omega_lim)

        if B_lim is not None:
            ax.set_ylim(B_lim)

    # Color bar
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(r"$A_{+-}(\omega,B)$")

    ax.set_title(r"Spectral function $A_{+-}(\omega,B)$")

    plt.tight_layout()
    plt.show()

    return fig, ax

