from astropy import units as u


class DensityModel:

    def __init__(self, name=None):
        self.name = name

    def evaluate(self, r):
        raise NotImplementedError


class SphericalBetaModel(DensityModel):
    def __init__(self, n0=0.0121 * u.cm ** -3, r_c=25 * u.kpc, beta=0.655):
        super().__init__("Spherical Beta Model")
        self.n0 = n0
        self.r_c = r_c
        self.beta = beta

    def evaluate(self, r):
        return self.n0 * (1 + (r / self.r_c) ** 2) ** (-3 * self.beta / 2)


class MB2015ModifiedBeta(DensityModel):
    def __init__(self, n0=0.0121 * u.cm ** -3, r_c=25 * u.kpc, beta=0.655):
        super().__init__("Miller and Bregman (2015) Modified Beta Model")
        self.n0 = n0
        self.r_c = r_c
        self.beta = beta

    def evaluate(self, r):
        return self.n0 * self.r_c ** (3 * self.beta) / r ** (3 * self.beta)
