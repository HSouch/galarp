from astropy import units as u
import numpy as np


__all__ = ['SigmaGas']


class SigmaGas:

    def __init__(self, sigma, nparticles = 1, units=u.Msun/u.pc**2):
        if isinstance(sigma, SigmaGas):
            self.sigma = sigma.sigma
            self.nparticles = sigma.nparticles
        else:
            self.sigma = sigma
            self.nparticles = nparticles
        self.units = units

        if not isinstance(self.sigma, u.Quantity):
            self.sigma *= self.units
        if self.sigma.isscalar:
            self.sigma = [self.sigma.to(self.units).value for i in range(nparticles)] * self.units

        self.sigma = self.sigma.to(self.units)
        self.nparticles = len(self.sigma)

    @staticmethod
    def load_spherical(mass, radius, n_particles=1):
        return SigmaGas(sigma = mass / (np.pi * radius**2), nparticles = n_particles)
