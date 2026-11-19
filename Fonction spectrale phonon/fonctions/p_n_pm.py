

import numpy as np

from fonctions.hamiltonien import diagonalisation
from fonctions.fermi_energy import energie_fermi


def P_n_pm( n, B,eta_cm,  valeurs_q, omega_cm_values, para_swm, para_num,para_phy,unites_phy):
    """
    Calcule la contribution P_n à la polarisation Pi_{+-}.

    Les blocs étudiés sont n et n-2.

    Formule conservée :
        DeltaE = E_(n-2) - E_n
        facteur_occupation = f_(n-2) - f_n

        Pi_n(omega) = C_B * somme [
            M_n * (f_(n-2) - f_n)
            / (omega + DeltaE + i*eta)
        ]
    """

    if n < 2:
        raise ValueError("n doit être supérieur ou égal à 2.")

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
    for q in valeurs_q:

        Gamma = para_swm.Gamma(q)

        # Diagonalisation des blocs n et n-2
        energies_n, vecteurs_n = diagonalisation(
            n, B, q, para_swm,para_phy
        )

        energies_n_moins_2, vecteurs_n_moins_2 = diagonalisation(
            n - 2, B, q, para_swm,para_phy
        )

        # Parcours des bandes lambda et lambda'
        for indice_lambda in range(len(energies_n)):

            for indice_lambda_prime in range(
                len(energies_n_moins_2)
            ):

                E_n = energies_n[indice_lambda]

                E_n_moins_2 = energies_n_moins_2[
                    indice_lambda_prime
                ]

                # Occupations de Fermi à température nulle
                f_n = 1.0 if E_n < EF else 0.0

                f_n_moins_2 = (
                    1.0 if E_n_moins_2 < EF else 0.0
                )

                # Différence d'occupation
                facteur_occupation = f_n_moins_2 - f_n

                # Même occupation : aucune contribution
                if facteur_occupation == 0:
                    continue

                # Composantes des vecteurs propres
                x1 = vecteurs_n_moins_2[
                    0, indice_lambda_prime
                ]

                x4 = vecteurs_n[
                    3, indice_lambda
                ]

                # Élément de matrice : formule conservée
                M_n = Gamma**2 * np.abs(
                    np.conj(x1) * x4
                )**2

                # Différence d'énergie : formule conservée
                DeltaE = E_n_moins_2 - E_n

                # Stockage pour le calcul de la polarisation
                M_liste_n.append(M_n)
                DeltaE_liste_n.append(DeltaE)
                occupation_liste_n.append(facteur_occupation)

                # Enregistrement du sens physique de la transition
                if facteur_occupation == 1:
                    # Bloc n-2 occupé -> bloc n vide
                    transition = ( n - 2, indice_lambda_prime, n, indice_lambda,  q )
                else:
                    # Bloc n occupé -> bloc n-2 vide
                    transition = ( n, indice_lambda, n - 2,  indice_lambda_prime, q )
                transitions_liste_n.append(transition)

    # Conversion des listes en tableaux NumPy
    M_tableau_n = np.array(M_liste_n, dtype=float)
    DeltaE_tableau_n = np.array(DeltaE_liste_n, dtype=float)
    occupation_tableau_n = np.array(
        occupation_liste_n, dtype=float
    )

    # Conversion des fréquences et de eta de cm^-1 vers eV
    omega_ev = np.asarray(omega_cm_values) * unites_phy.cm1_to_eV
    eta_ev = eta_cm * unites_phy.cm1_to_eV

    # Dénominateur 
    D = 1.0 / (
        omega_ev[:, None]
        + DeltaE_tableau_n[None, :]
        + 1j * eta_ev
    )

    # Polarisation : formule conservée
    P_n_pm_values = (2e-10 * B) * np.dot(
        D,
        M_tableau_n * occupation_tableau_n
    ) # calcule final M fois occupation deux vecteur on multippe terme par tereme 


    return (
    P_n_pm_values,
    M_tableau_n,
    DeltaE_tableau_n,
    occupation_tableau_n,
    transitions_liste_n
)




