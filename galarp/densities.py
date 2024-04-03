from astropy import units as u
import numpy as np


__all__ = [
    "Density",
    "ExponentialDensity",
    "InterpolatedDensity",
]

class Density:
    """Parent class to represent time-variable densities."""

    def __init__(self, rho):
        self.rho = rho
        assert type(rho) == u.Quantity, "Density must be a Quantity"

    def evaluate(self, t):
        return self.rho.to(u.g / u.cm**3)
    
    def evaluate_arr(self, ts):
        return np.array([self.evaluate(t).value for t in ts])


class ExponentialDensity(Density):
    def __init__(self, rho, t0, width):
        super().__init__(rho)
        self.t0 = t0
        self.width = width

        assert type(t0) == u.Quantity, "t0 must be a Quantity"
        assert type(width) == u.Quantity, "width must be a Quantity"

    def evaluate(self, t):
        return self.rho * np.exp(-(((t - self.t0) / self.width) ** 2))


class InterpolatedDensity(Density):
    """Density class that returns an interpolated value for the density at a given time t.
    Helpful for if you are determining the density through a GALA orbit and want it to be as realistic
    as possible.
    """

    def __init__(self, interp):
        super().__init__(
            0
        )  # The "rho" can just be 0 here, as it will be overwritten by the interpolation

        self.interp = interp

    def evaluate(self, t):
        return self.interp(t)