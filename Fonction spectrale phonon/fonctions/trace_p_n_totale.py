
import numpy as np
import matplotlib.pyplot as plt

from fonctions.p_n_pm import P_n_pm
from fonctions.P_n_mp import P_n_mp

# n commence par 2
def tracer_P_n_totale(
    n, B, eta_cm, valeurs_q, omega_cm_values,
    para_swm, para_num, unites_phy,
    pas_x=10, pas_y=None, x_lim=(10, 200),
    y_lim=None
):

    # --------------------------------------------------
    # 1. Calcul des deux polarisations
    # --------------------------------------------------

    resultat_pm = P_n_pm(
        n=n,
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    resultat_mp = P_n_mp(
        n=n-2,
        B=B,
        eta_cm=eta_cm,
        valeurs_q=valeurs_q,
        omega_cm_values=omega_cm_values,
        para_swm=para_swm,
        para_num=para_num,
        unites_phy=unites_phy
    )

    # Récupérer les résultats
    P_pm, M_pm, DeltaE_pm, occupation_pm, transitions_pm = resultat_pm
    P_mp, M_mp, DeltaE_mp, occupation_mp, transitions_mp = resultat_mp


    # --------------------------------------------------
    # 2. Étiquettes des branches
    # --------------------------------------------------

    def etiquettes_bloc(indice):
        if indice == 0:
            return {0: "0"}
        elif indice == 1:
            return {0: "--", 1: "0", 2: "++"}
        else:
            return {0: "--", 1: "-", 2: "+", 3: "++"}


    # --------------------------------------------------
    # 3. Fonction commune pour tracer les contributions
    # --------------------------------------------------

    def tracer_contributions(
        M, DeltaE, occupation, transitions,
        indice_bas, indice_haut,
        omega_ev, eta_ev,
        ax, couleur_ligne, polarisation
    ):

        labels_bas = etiquettes_bloc(indice_bas)
        labels_haut = etiquettes_bloc(indice_haut)

        # Regrouper les transitions par branches
        contributions = {}

        for i, t in enumerate(transitions):

            bloc_initial, lambda_initial, bloc_final, lambda_final, q = t

            if bloc_initial == indice_bas and bloc_final == indice_haut:

                branche_initiale = labels_bas.get(
                    lambda_initial, str(lambda_initial)
                )
                branche_finale = labels_haut.get(
                    lambda_final, str(lambda_final)
                )

            elif bloc_initial == indice_haut and bloc_final == indice_bas:

                branche_initiale = labels_haut.get(
                    lambda_initial, str(lambda_initial)
                )
                branche_finale = labels_bas.get(
                    lambda_final, str(lambda_final)
                )

            else:
                continue

            cle = (
                bloc_initial, lambda_initial,
                bloc_final, lambda_final,
                branche_initiale, branche_finale
            )

            contributions.setdefault(cle, []).append(i)


        # Calculer et tracer chaque contribution
        palette = [
            "red", "blue", "green", "purple",
            "orange", "brown", "magenta", "cyan",
            "olive", "pink", "gray", "teal"
        ]

        couleurs = {}

        for cle, indices in contributions.items():

            (
                bloc_initial, lambda_initial,
                bloc_final, lambda_final,
                branche_initiale, branche_finale
            ) = cle

            indices = np.asarray(indices, dtype=int)

            D = 1.0 / (
                omega_ev[:, None]
                + DeltaE[indices][None, :]
                + 1j * eta_ev
            )

            Pi_partielle = C_B * np.dot(
                D,
                M[indices] * occupation[indices]
            )

            # Couleur associée à la combinaison de branches
            cle_couleur = (
                bloc_initial, branche_initiale,
                bloc_final, branche_finale
            )

            if cle_couleur not in couleurs:
                couleurs[cle_couleur] = palette[
                    len(couleurs) % len(palette)
                ]

            couleur = couleurs[cle_couleur]

            # Étiquette de la courbe
            label = (
                rf"${bloc_initial}^{{{branche_initiale}}}"
                rf"\rightarrow {bloc_final}^{{{branche_finale}}}$"
            )

            if "--" in (branche_initiale, branche_finale) or \
                    "++" in (branche_initiale, branche_finale):
                label = "_nolegend_"

            ax.plot(
                omega_cm_values,
                np.imag(Pi_partielle),
                color=couleur,
                linestyle=couleur_ligne,
                linewidth=1.5,
                alpha=0.8,
                label=label
            )

        return len(transitions), len(contributions)


    # --------------------------------------------------
    # 4. Préparer les unités
    # --------------------------------------------------

    omega_ev = (
        np.asarray(omega_cm_values)
        * unites_phy.cm1_to_eV
    )

    eta_ev = eta_cm * unites_phy.cm1_to_eV


    # --------------------------------------------------
    # 5. Créer une seule figure
    # --------------------------------------------------

    fig, ax = plt.subplots(figsize=(12, 6))


    # PM : transitions entre n et n-2
    nombre_pm = tracer_contributions(
        M_pm, DeltaE_pm, occupation_pm, transitions_pm,
        n-2, n,
        omega_ev, eta_ev,
        ax, "-", r"\Pi_{+-}"
    )

    # MP : transitions entre n-2 et n
    nombre_mp = tracer_contributions(
        M_mp, DeltaE_mp, occupation_mp, transitions_mp,
        n-2, n,
        omega_ev, eta_ev,
        ax, "--", r"\Pi_{-+}"
    )


    # --------------------------------------------------
    # 6. Mise en forme de la figure
    # --------------------------------------------------


    ax.set_xlabel(r"$\omega$ (cm$^{-1}$)")
    ax.set_ylabel(
        r"$\mathrm{Im}\,\Pi_{+-}^{(n)}$ et "
        r"$\mathrm{Im}\,\Pi_{-+}^{(n-2)}$"
    )

    ax.set_title(
        rf"$n={n}$, $B={B}$ T, $\eta={eta_cm}$ cm$^{{-1}}$"
    )

    ax.set_xlim(x_lim)
    ax.set_xticks(
        np.arange(x_lim[0], x_lim[1] + pas_x, pas_x)
    )

    if y_lim is not None:
        ax.set_ylim(y_lim)

    if pas_y is not None:
        y_min, y_max = ax.get_ylim()
        ax.set_yticks(
            np.arange(y_min, y_max + pas_y * 0.5, pas_y)
        )

    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")

    plt.tight_layout()
    plt.show()

    return P_pm, P_mp
