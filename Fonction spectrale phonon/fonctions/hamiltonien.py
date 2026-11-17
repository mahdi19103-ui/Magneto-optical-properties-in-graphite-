
import numpy as np  
def H_n(n, B, q, parametres):

    eps = parametres.epsilon(B)

    Gamma = parametres.Gamma(q)
    Gamma1 = parametres.Gamma1(q)
    Gamma2 = parametres.Gamma2(q)
    Gamma5 = parametres.Gamma5(q)

    # on prend gamma4 , c'est le terme qui contientgamma4
    terme = eps * (parametres.gamma4 / parametres.gamma0) * Gamma

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
def diagonalisation(n, B, q, parametres):
    H = H_n(n, B, q, parametres)

    energies, etats_propres = np.linalg.eigh(H)

    return  energies, etats_propres

## si j'ai besoin juste des valeurs propres comme dans lorsque je trace structure de bande
def energies_seules(n, B, q, parametres):
    matrice = H_n(n, B, q, parametres)
    energies = np.linalg.eigvalsh(matrice)
    return energies
