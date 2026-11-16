

import numpy as np
import matplotlib.pyplot as plt


def tracer_A_pm_B(
    B_values,
    A_values,
    omega_cm_values,
    x_lim=None,
    vmax=None,
    cmap="RdBu_r"
):

    B_values = np.asarray(B_values)
    A_values = np.asarray(A_values)
    omega_cm_values = np.asarray(omega_cm_values)

    # Figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # Echelle de couleur
    if vmax is None:
        vmax = np.max(np.abs(A_values))

    # Carte
    image = ax.imshow(
        A_values,
        origin="lower",
        aspect="auto",
        extent=[
            omega_cm_values[0],
            omega_cm_values[-1],
            B_values[0],
            B_values[-1]
        ],
        cmap=cmap,
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest"
    )

    # Barre de couleur
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(r"$A_{+-}(\omega,B)$")

    # Axes
    ax.set_xlabel(r"$\omega\ (\mathrm{cm}^{-1})$")
    ax.set_ylabel(r"$B\ (\mathrm{T})$")
    ax.set_title(r"Fonction spectrale $A_{+-}(\omega,B)$")

    if x_lim is not None:
        ax.set_xlim(x_lim)

    ax.set_ylim(B_values[0], B_values[-1])

    plt.tight_layout()
    plt.show()

    return fig, ax

