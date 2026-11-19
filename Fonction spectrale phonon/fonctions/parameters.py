import numpy as np

class Paraswm:
    def __init__(self,gamma0=3.37,gamma1=0.363,gamma2=-0.02,gamma3=0.31,gamma4=0.07,gamma5=0.05):
        self.gamma0 = gamma0 
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3 
        self.gamma4 = gamma4
        self.gamma5 = gamma5


    def Gamma(self,q): # avec q entre -0.5 et 0.5
        return 2*np.cos(np.pi*q)

    # Fonctions dépendant de q
    def Gamma1(self, q):
        return  self.gamma1*self.Gamma(q)

    def Gamma2(self, q):
        return  self.gamma2 *(self.Gamma(q)**2-2)

    def Gamma5(self, q):
        return  self.gamma5*(self.Gamma(q)**2-2)

class Paranum:
    def __init__( self, Nq=401, N_B=80,  n_max=30,  B=5,   N_omega=1200, omega_cm_min=-200.0,
                  omega_cm_max=200.0):
        self.Nq = Nq
        self.N_B = N_B
        self.n_max = n_max
        self.B = B
        self.N_omega = N_omega
        self.omega_cm_min = omega_cm_min
        self.omega_cm_max = omega_cm_max    

class Paraphy:
    def __init__(self,omega_0_cm=44.0,eta_cm=0.1,g3=0.01 ):
        # Fréquence propre du phonon, en cm⁻¹
        self.omega_0_cm = omega_0_cm

        # Élargissement, en cm⁻¹
        self.eta_cm = eta_cm

        # Couplage électron-phonon en eV/A
        self.g3 = g3
        
   # Echelle energetique qui depend du champ ce terme apparait dan sl'hamiltonien     
    def epsilon(self, B):
        return 0.036 * np.sqrt(B)  # eV et B en Tesla 

class Unitesphy():
    def __init__(self):
        self.e = 1.602176634e-19     # C
        self.eV_to_cm1 = 8065.544
        self.cm1_to_eV = 1 / 8065.544
        self.hbar = 1.054571817e-34  # J.s
        self.k_B = 1.380649e-23      # J/K
        self.J_to_eV = 1 / 1.602176634e-19    