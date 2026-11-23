import numpy as np
from fonctions.hamiltonien import diagonalisation
from fonctions.fermi_energy import energie_fermi
import matplotlib.pyplot as plt

def p_n_mp( n, B,eta_cm,  q_values, omega_cm_values, para_swm, para_num,para_phy,unites_phy):


    """
    Calculate the contribution P_n(omega) to the polarization Pi_{-+}
    associated with transitions between Landau levels n and n+2.

    The calculation takes into account the electronic energy differences,
    occupation factors at zero temperature, and matrix elements between
    the corresponding eigenstates.

    The final polarization is expressed in eV and includes the renormalization
    prefactor proportional to g3^2 * B, where g3 is the electron-phonon
    coupling constant.

    Returns:
        P_n_mp_values: Polarization contribution as a function of omega.
        M_tableau_n: Matrix elements for the allowed transitions.
        DeltaE_tableau_n: Energy differences between the states.
        occupation_tableau_n: Occupation factors.
        transitions_liste_n: List identifying the corresponding transitions.
    """

    if n < 0:
        raise ValueError("n doit être supérieur ou égal à 0.")

    # Énergie de Fermi
    EF = energie_fermi(
        B=B,
        para_swm=para_swm,
        para_num=para_num,para_phy=para_phy
    )

    # Listes pour stocker les contributions
    M_liste_n = []
    DeltaE_liste_n = []
    occupation_liste_n = []
    transitions_liste_n = []

    # Parcours des valeurs de q
    for q in q_values:

        Gamma = para_swm.Gamma(q)

        # Diagonalisation des blocs n et n+2
        energies_n, vecteurs_n = diagonalisation(
            n, B, q, para_swm,para_phy
        )

        energies_n_plus_2, vecteurs_n_plus_2 = diagonalisation(
            n + 2, B, q, para_swm,para_phy
        )

        # Parcours des bandes lambda et lambda'
        for indice_lambda in range(len(energies_n)):

            for indice_lambda_prime in range(
                len(energies_n_plus_2)
            ):

                E_n = energies_n[indice_lambda]

                E_n_plus_2 = energies_n_plus_2[
                    indice_lambda_prime
                ]

                # Occupations de Fermi à température nulle
                f_n = 1.0 if E_n < EF else 0.0

                f_n_plus_2 = (
                    1.0 if E_n_plus_2 < EF else 0.0
                )

                # Différence d'occupation
                facteur_occupation = f_n_plus_2 - f_n

                # Même occupation : aucune contribution
                if facteur_occupation == 0:
                    continue

                # Composantes des vecteurs propres
                x1 = vecteurs_n[
                    0, indice_lambda
                ]

                x4 = vecteurs_n_plus_2[
                    3, indice_lambda_prime
                ]

                # Élément de matrice : formule conservée
                M_n = Gamma**2 * np.abs(
                    np.conj(x4) * x1
                )**2

                # Différence d'énergie : formule conservée
                DeltaE = E_n_plus_2 - E_n

                # Stockage pour le calcul de la polarisation
                M_liste_n.append(M_n)
                DeltaE_liste_n.append(DeltaE)
                occupation_liste_n.append(facteur_occupation)

                # Enregistrement du sens physique de la transition
                if facteur_occupation == 1:
                    # Bloc n+2 occupé -> bloc n vide
                    transition = ( n + 2, indice_lambda_prime, n, indice_lambda,  q )
                else:
                    # Bloc n occupé -> bloc n+2 vide
                    transition = ( n, indice_lambda, n + 2,  indice_lambda_prime, q )
                transitions_liste_n.append(transition)

    # Conversion des listes en tableaux NumPy
    M_tableau_n = np.array(M_liste_n, dtype=float)
    DeltaE_tableau_n = np.array(DeltaE_liste_n, dtype=float)
    occupation_tableau_n = np.array(
        occupation_liste_n, dtype=float
    )

    # Conversion des fréquences et de eta de cm^-1 vers eV
    omega_ev = np.asarray(omega_cm_values) * unites_phy.cm_to_eV
    eta_ev = eta_cm * unites_phy.cm_to_eV

    # Dénominateur 
    D = 1.0 / (
        omega_ev[:, None]
        + DeltaE_tableau_n[None, :]
        + 1j * eta_ev
    )

    # Polarisation : formule conservée
    p_n_mp_values = (2.032e-7*(para_phy.g3)**2 * B) * np.dot(
        D,
        M_tableau_n * occupation_tableau_n
    ) # calcule final M fois occupation deux vecteur on multippe terme par tereme 

    

    return (
    p_n_mp_values,
    M_tableau_n,
    DeltaE_tableau_n,
    occupation_tableau_n,
    transitions_liste_n
)




def plot_p_n_mp(
    n, B, eta_cm, q_values, omega_cm_values,
    para_swm, para_num, para_phy, unites_phy,
    pas_x=10, pas_y=None, x_lim=(10, 200), y_lim=None
):
    """
    Plot the imaginary part of the polarization contribution P_n(omega)
    for transitions between Landau levels n and n-2 and n-2 to n.

    Each transition between energy branches is plotted separately.
    The total polarization P_n(omega) is also calculated and returned.
        """

    
    # 1. Calculate the polarization
    (
        p_n_mp_values,
        M_n,
        DeltaE_n,
        occupation_n,
        transitions_n
    ) = p_n_mp(
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

    # 2. Branch labels
    def etiquettes_bloc(n):
        if n == 0:
            return {0: "0"}
        elif n == 1:
            return {0: "--", 1: "0", 2: "++"}
        else:
            return {0: "--", 1: "-", 2: "+", 3: "++"}

    labels_n = etiquettes_bloc(n)
    labels_n_plus_2 = etiquettes_bloc(n + 2)

    # 3. Group transitions by branches
    contributions = {}

    for i, t in enumerate(transitions_n):
        bloc_initial, lambda_initial, bloc_final, lambda_final, q = t

        if bloc_initial == n + 2 and bloc_final == n:
            branche_initiale = labels_n_plus_2.get(
                lambda_initial, str(lambda_initial)
            )
            branche_finale = labels_n.get(
                lambda_final, str(lambda_final)
            )

        elif bloc_initial == n and bloc_final == n + 2:
            branche_initiale = labels_n.get(
                lambda_initial, str(lambda_initial)
            )
            branche_finale = labels_n_plus_2.get(
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

        contributions.setdefault(cle, []).append(i)

    # 4. Calculate each contribution
    omega_ev = (
        np.asarray(omega_cm_values)
        * unites_phy.cm_to_eV
    )
    eta_ev = eta_cm * unites_phy.cm_to_eV

    def calculer_contribution(indices):
        indices = np.asarray(indices, dtype=int)

        D = 1.0 / (
            omega_ev[:, None]
            + DeltaE_n[indices][None, :]
            + 1j * eta_ev
        )

        return (2.032e-7 * para_phy.g3**2 * B) * np.dot(
            D,
            M_n[indices] * occupation_n[indices]
        )

    # 5. Automatic colors for each transition
    palette = [
        "red", "blue", "green", "purple",
        "orange", "brown", "magenta", "cyan",
        "olive", "pink", "gray", "teal"
    ]

    couleurs = {}

    # 6. Plot contributions
    plt.figure(figsize=(12, 6))

    for cle, indices in contributions.items():

        (
            bloc_initial,
            lambda_initial,
            bloc_final,
            lambda_final,
            branche_initiale,
            branche_finale
        ) = cle

        p_partielle = calculer_contribution(indices)

        cle_couleur = (
            bloc_initial, branche_initiale,
            bloc_final, branche_finale
        )

        if cle_couleur not in couleurs:
            couleurs[cle_couleur] = palette[
                len(couleurs) % len(palette)
            ]

        couleur = couleurs[cle_couleur]

        label = (
            rf"${bloc_initial}^{{{branche_initiale}}}"
            rf"\rightarrow{bloc_final}^{{{branche_finale}}}$"
        )

        # Keep all curves but hide -- and ++ from the legend
        if "--" in (branche_initiale, branche_finale) or \
           "++" in (branche_initiale, branche_finale):
            label = "_nolegend_"

        plt.plot(
            omega_cm_values,
            np.imag(p_partielle),
            color=couleur,
            linewidth=1.5,
            alpha=0.8,
            label=label
        )

    # 7. Plot formatting
    plt.xlabel(r"$\omega$ (cm$^{-1}$)")
    plt.ylabel(r"$\mathrm{Im}\,\Pi_{-+}^{(n)}$")

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

    return p_n_mp_values
