
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def tracer_A_pm_B_inverse(
    B_values,
    A_values,
    omega_cm_values,
    B_lim=None,
    omega_lim=None,
    vmax=None,pas_omega=10,
    cmap="RdBu_r"
):

    B_values = np.asarray(B_values)
    A_values = np.asarray(A_values)
    omega_cm_values = np.asarray(omega_cm_values)

    # Transposer A_values pour inverser les axes
    A_values_inverse = A_values.T

    # Définir l'échelle des couleurs
    if vmax is None:
        vmax = np.max(np.abs(A_values))

    # Créer la figure
    fig, ax = plt.subplots(figsize=(10, 6))

    image = ax.imshow(
    A_values_inverse,
    origin="lower",
    aspect="auto",
    extent=[
        B_values[0],
        B_values[-1],
        omega_cm_values[0],
        omega_cm_values[-1]
    ],
    cmap=cmap,
    vmin=-vmax,
    vmax=vmax,
    interpolation="nearest"
)

    # Barre de couleurs
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(r"$A_{+-}(\omega,B)$")

    # Axes
    ax.set_xlabel(r"$B\ (\mathrm{T})$")
    ax.set_ylabel(r"$\omega\ (\mathrm{cm}^{-1})$")
    ax.yaxis.set_major_locator(MultipleLocator(pas_omega))
    ax.set_title(r"Fonction spectrale $A_{+-}(\omega,B)$")

    # Limites des axes
    if B_lim is not None:
        ax.set_xlim(B_lim)

    if omega_lim is not None:
        ax.set_ylim(omega_lim)

    plt.tight_layout()
    plt.show()

    return fig, ax
