class ParametresPhysiques:
    def __init__(self,omega_0_cm=44.0,eta_cm=0.1,g3=1 ):
        # Fréquence propre du phonon, en cm⁻¹
        self.omega_0_cm = omega_0_cm

        # Élargissement, en cm⁻¹
        self.eta_cm = eta_cm

        # Couplage électron-phonon
        self.g3 = g3