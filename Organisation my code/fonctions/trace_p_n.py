

import numpy as np
import matplotlib.pyplot as plt
from fonctions.p_n_pm import P_n_pm


def tracer_P_n_pm(n, B, eta_cm,valeurs_q, omega_cm_values,  para_swm, para_num, unites_phy, pas_x=10, pas_y=None,x_lim=(10, 200), 
    y_lim=None
):
    # 1. Calcul de la polarisation
    (
        P_n_pm_values,
        M_n,
        DeltaE_n,
        occupation_n,
        transitions_n
    ) = P_n_pm(
        n=n,
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    # 2. Étiquettes des branches
    def etiquettes_bloc(n):
        if n == 0:
            return {0: "0"}
        elif n == 1:
            return {0: "--", 1: "0", 2: "++"}
        else:
            return {0: "--", 1: "-", 2: "+", 3: "++"}

    labels_n = etiquettes_bloc(n)
    labels_n_moins_2 = etiquettes_bloc(n - 2)

    # 3. Regrouper les transitions par branches
    contributions = {}

    for i, t in enumerate(transitions_n): # i tres important donne quel trasnition pour q donnne
        bloc_initial, lambda_initial, bloc_final, lambda_final, q = t

        if bloc_initial == n - 2 and bloc_final == n:
            branche_initiale = labels_n_moins_2.get(
                lambda_initial, str(lambda_initial)
            )
            branche_finale = labels_n.get(
                lambda_final, str(lambda_final)
            )

        elif bloc_initial == n and bloc_final == n - 2:
            branche_initiale = labels_n.get(
                lambda_initial, str(lambda_initial)
            )
            branche_finale = labels_n_moins_2.get(
                lambda_final, str(lambda_final)
            )

        else:
            continue

        cle = (
            bloc_initial,
            lambda_initial,
            bloc_final,
            lambda_final,
            branche_initiale,
            branche_finale
        )

        contributions.setdefault(cle, []).append(i) # il donne la contriution pour i c'est a dire pour une type de trasnitions

    # 4. Calculer chaque contribution
    omega_ev = (
        np.asarray(omega_cm_values)
        * unites_phy.cm1_to_eV
    )
    eta_ev = eta_cm * unites_phy.cm1_to_eV

    def calculer_contribution(indices):
        indices = np.asarray(indices, dtype=int)

        D = 1.0 / (
            omega_ev[:, None]
            + DeltaE_n[indices][None, :]
            + 1j * eta_ev
        )

        return (2e-7 * B) * np.dot(
            D,
            M_n[indices] * occupation_n[indices]
        )

    # 5. Couleurs automatiques pour chaque combinaison
    palette = [
        "red", "blue", "green", "purple",
        "orange", "brown", "magenta", "cyan",
        "olive", "pink", "gray", "teal"
    ]

    couleurs = {}

    # 6. Tracer les contributions
    plt.figure(figsize=(12, 6))

    for cle, indices in contributions.items(): # boucle dans contributions pour chaqwu type des trasnitions

        (
            bloc_initial,
            lambda_initial,
            bloc_final,
            lambda_final,
            branche_initiale,
            branche_finale
        ) = cle

        Pi_partielle = calculer_contribution(indices) # on recoit une type donc on fait las omme sur touts les q


        cle_couleur = (
            bloc_initial, branche_initiale,
            bloc_final, branche_finale
            )

        if cle_couleur not in couleurs:
            couleurs[cle_couleur] = palette[len(couleurs) % len(palette)]

        couleur = couleurs[cle_couleur]

        label = rf"${bloc_initial}^{{{branche_initiale}}}\rightarrow{bloc_final}^{{{branche_finale}}}$"

        # Garder toutes les courbes, mais cacher -- et ++
        if "--" in (branche_initiale, branche_finale) or \
                "++" in (branche_initiale, branche_finale):
            label = "_nolegend_"

        plt.plot(
            omega_cm_values,
            np.imag(Pi_partielle),
            color=couleur,
            linewidth=1.5,
            alpha=0.8,
            label=label
        )

    plt.xlabel(r"$\omega$ (cm$^{-1}$)")
    plt.ylabel(r"$\mathrm{Im}\,\Pi_{+-}^{(n)}$")
    plt.title(
        rf"$n={n}$, $B={B}$ T, $\eta={eta_cm}$ cm$^{{-1}}$"
    )

    plt.xlim(x_lim)
    plt.xticks(np.arange(x_lim[0], x_lim[1] + pas_x, pas_x))

    if y_lim is not None:
        plt.ylim(y_lim)

    if pas_y is not None:
        y_min, y_max = plt.ylim()
        plt.yticks(np.arange(y_min, y_max + pas_y * 0.5, pas_y))

    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.show()

    return P_n_pm_values
