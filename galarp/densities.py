from astropy import units as u
import numpy as np

from astropy.table import Table
from scipy.interpolate import interp1d


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
    """Class to represent an exponential density profile."""
    def __init__(self, rho, t0, width):
        super().__init__(rho)
        self.t0 = t0
        self.width = width

        assert type(t0) == u.Quantity, "t0 must be a Quantity"
        assert type(width) == u.Quantity, "width must be a Quantity"

    def evaluate(self, t):
        return self.rho * np.exp(-(((t - self.t0) / self.width) ** 2))


class InterpolatedDensity(Density):
    """Class to represent an interpolated density profile."""
    def __init__(self, interp=None, units=u.g/u.cm**3, **kwargs):
        super().__init__(rho=0 * u.g/u.cm**3, **kwargs)
        self.interp = interp

    def evaluate(self, t):
        """ Return the density at time t. """
        return self.interp(t) * u.g /u.cm **3

    @staticmethod
    def from_table(fn, t_key, rho_key,
                    t_units=u.s, rho_units=u.g/u.cm**3,
                    format="ascii", **kwargs):
        """
        Create an InterpolatedDensity object from a table.

        Parameters:
        fn (str): The filename of the table.
        t_key (str): The column name for the time values.
        rho_key (str): The column name for the density values.
        t_units (astropy.units.Unit, optional): The units for the time values. 
            Default is seconds (u.s).
        rho_units (astropy.units.Unit, optional): The units for the density values. 
            Default is grams per cubic centimeter (u.g/u.cm**3).
        format (str, optional): The format of the table. Default is "ascii".
        **kwargs: Additional keyword arguments to be passed to the InterpolatedDensity constructor.

        Returns:
        InterpolatedDensity: An InterpolatedDensity object created from the table.
        """
        t = Table.read(fn, format=format)

        ts = t[t_key] * t_units.to(u.Myr)

        rho = t[rho_key] * rho_units

        interp = interp1d(ts, rho, bounds_error=False, fill_value="extrapolate")

        return InterpolatedDensity(interp=interp, **kwargs)