class Paramnum:
    def __init__(
        self,
        Nq=401,
        N_B=80,
        n_max=30,
        B=5,
        N_omega=1200,
        omega_cm_min=-200.0, #cm
        omega_cm_max=200.0
    ):
        self.Nq = Nq
        self.N_B = N_B
        self.n_max = n_max
        self.B = B
        self.N_omega = N_omega
        self.omega_cm_min = omega_cm_min
        self.omega_cm_max = omega_cm_max