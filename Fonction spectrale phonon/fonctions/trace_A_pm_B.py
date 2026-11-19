import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def tracer_A_pm_B(
    B_values,
    A_values,
    omega_cm_values,
    inverse_axes=False,
    B_lim=None,
    omega_lim=None,
    vmax=None,
    pas_omega=10,
    cmap="RdBu_r"
):

    B_values = np.asarray(B_values)
    A_values = np.asarray(A_values)
    omega_cm_values = np.asarray(omega_cm_values)

    # Echelle de couleur
    if vmax is None:
        vmax = np.max(np.abs(A_values))

    # Figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Inversion ou non des axes
    if inverse_axes:

        image = ax.imshow(
            A_values.T,
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

        ax.set_xlabel(r"$B\ (\mathrm{T})$")
        ax.set_ylabel(r"$\omega\ (\mathrm{cm}^{-1})$")

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(pas_omega))

        if B_lim is not None:
            ax.set_xlim(B_lim)

        if omega_lim is not None:
            ax.set_ylim(omega_lim)

    else:

        image = ax.imshow(
            A_values,
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

        ax.set_xlabel(r"$\omega\ (\mathrm{cm}^{-1})$")
        ax.set_ylabel(r"$B\ (\mathrm{T})$")

        if omega_lim is not None:
            ax.set_xlim(omega_lim)

        if B_lim is not None:
            ax.set_ylim(B_lim)

    # Barre de couleurs
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(r"$A_{+-}(\omega,B)$")

    ax.set_title(r"Fonction spectrale $A_{+-}(\omega,B)$")

    plt.tight_layout()
    plt.show()

    return fig, ax


