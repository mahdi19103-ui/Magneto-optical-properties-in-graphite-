"""
Objectives
---------
Definition of the graphite Hamiltonian in the presence of a magnetic field B,
and calculation of its eigenvalues and eigenvectors.

Remark1:
The effect of gamma3 is neglected in the electronic structure calculation.
Including gamma3 (trigonal warping) would require a more complex,
potentially infinite-dimensional Hamiltonian representation.
"""




import numpy as np


""""
n: Landau level index (n = 0, 1, 2, ...).

para_swm and para_phy: instances of the classes defined in parameters.py.
"""


def H_n(n, B, q, para_swm,para_phy):

    
    eps = para_phy.epsilon(B)

    Gamma = para_swm.Gamma(q)
    Gamma1 = para_swm.Gamma1(q)
    Gamma2 = para_swm.Gamma2(q)
    Gamma5 = para_swm.Gamma5(q)

    # on prend gamma4 , c'est le terme qui contientgamma4
    terme = eps * (para_swm.gamma4 / para_swm.gamma0) * Gamma

    if n == 0:
        return np.array([
            [Gamma2]
        ])

    elif n == 1:
        return np.array([
            [Gamma2, eps, terme],
            [eps, Gamma5, Gamma1],
            [terme, Gamma1, Gamma5]
        ])

    elif n >= 2:
        return np.array([
            [Gamma2, eps*np.sqrt(n), terme*np.sqrt(n), 0],
            [eps*np.sqrt(n), Gamma5, Gamma1, terme*np.sqrt(n-1)],
            [terme*np.sqrt(n), Gamma1, Gamma5, eps*np.sqrt(n-1)],
            [0, terme*np.sqrt(n-1), eps*np.sqrt(n-1), Gamma2]
        ])


# Diagonalsiation
def diagonalisation(n, B, q, para_swm, para_phy):
    H = H_n(n, B, q, para_swm,para_phy)

    energies, etats_propres = np.linalg.eigh(H)

    return  energies, etats_propres

## si j'ai besoin juste des valeurs propres comme dans lorsque je trace structure de bande
def energies_seules(n, B, q, para_swm,para_phy):
    matrice = H_n(n, B, q, para_swm,para_phy)
    energies = np.linalg.eigvalsh(matrice)
    return energies
