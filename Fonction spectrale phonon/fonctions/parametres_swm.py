import numpy as np

class ParametresSWM:
    def __init__(self,gamma0=3.37,gamma1=0.363,gamma2=-0.02,gamma3=0.31,gamma4=0.07,gamma5=0.05):
        self.gamma0 = gamma0 #eV
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3 # eV
        self.gamma4 = gamma4
        self.gamma5 = gamma5



        # on met epsilon qui depend de B dans cette classe aussi
    def epsilon(self, B):
        return 0.036 * np.sqrt(B) # epsilon en eV et B en tesla

    def Gamma(self,q): # avec q entre -0.5 et 0.5
        return 2*np.cos(np.pi*q)

    # Fonctions dépendant de q
    def Gamma1(self, q):
        return  self.gamma1*self.Gamma(q)

    def Gamma2(self, q):
        return  self.gamma2 *(self.Gamma(q)**2-2)

    def Gamma5(self, q):
        return  self.gamma5*(self.Gamma(q)**2-2)


    