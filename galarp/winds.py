# Third Party
import os
import pickle
import astropy.units as u
from astropy.units import Quantity
import numpy as np

from tqdm import tqdm

from matplotlib import pyplot as plt

# This module
from . import shadows, utils

import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic


#import integrate as gi            # If debugging
import gala.integrate as gi


__all__ = ['RPWind', 'LorentzianWind', 'StepFunctionWind', 'Density', 'ExponentialDensity', 'InterpolatedDensity']



class RPWind:
    """ Class to represent a ram pressure wind

        Usage:
            ```
            wind = gn.RPWind([50, 0, 50] * (u.km / u.s), units=gn.galactic)  # Initializes a wind at 45 degrees in the x-z plane.

            wind = gn.RPWind(units=galactic)                                 # Achieves the same thing.
            wind.init_from_inc(np.deg2rad(45), 100 * u.km  / u.s)          
            ```
    """

    def __init__(self, vector=None, units=None):
        """
        Args:
            vector (astropy Quantity, optional): Wind vector. Can also be initialized using self.init_from_inc(). Defaults to None.
            units (_type_, optional): _description_. Defaults to None.
        """
        if vector is None:
            self.vector = Quantity([0, 0, 0]) * (u.km / u.s)
        # self.vector = Quantity(vector)

        self.units = units
        if self.units is not None:
            if type(self.vector) == Quantity:
                self.vector = self.vector.to(self.units["length"] / self.units["time"])
            else:
                self.vector *= (self.units["length"] / self.units["time"])

    def evaluate(self, t):
        """ Return the wind vector at time t. For the default wind vector, this is JUST the vector converted to kpc/Myr """
        return self.vector.to(u.kpc / u.Myr).value

    def initialize_vector(self):
        if self.units is not None:
            self.vector = self.vector.to(self.units["length"] / self.units["time"])

    def init_from_inc(self, inclination, strength):
        """ Initialize the wind vector from an inclination and strength """
        x = strength * np.cos(inclination)
        z = strength * np.sin(inclination)
        self.vector = Quantity([x, 0 * x.unit, z]).to(self.units["length"] / self.units["time"])

    def wind_strength(self):
        """ Return the length (strength) of the wind vector """
        return np.sqrt(sum(self.vector ** 2))

    def normalized(self):
        """ Return the normalized wind vector """
        return self.vector / self.wind_strength()
    
    def vector_to_units(self, units):
        return self.vector.to(units)
    
    def vector_as_value(self):
        if type(self.vector) == Quantity:
            return self.vector.value
        else:
            return self.vector

    def inclination(self):
        x,y,z = self.vector_as_value()
        return np.arctan2(z, np.sqrt(x**2 + y**2))
    
    def __repr__(self):
        return f"<RP Wind Vector={self.vector}  Inclination={np.round(np.rad2deg(self.inclination()), 2):.2f}  >"


class LorentzianWind(RPWind):
    """ Ram pressure wind that is damped by a Lorentzian profile. 
        Reaches max (which is just the unadjusted RP wind value) at t0.

    Args:
        gn (_type_): _description_
    """
    def __init__(self, t0=0 * u.Myr, width=200 * u.Myr, **kwargs):
        super().__init__(**kwargs)
        self.t0 = t0.to(u.Myr).value                 # Units are in Myr
        self.width = width.to(u.Myr).value              # Units are in Myr
        
    def evaluate(self, t):
        """ Return the wind vector damped by a Lorentzian profile """
        return super().evaluate(t) * 1 / ((2 * (t - self.t0) / self.width)**2 + 1)


class StepFunctionWind(RPWind):
    """ A wind represented by a step function (no wind immediately into full wind at time t0)
        Good for testing that time variability is working.

    """
    def __init__(self, t0=0 * u.Myr, **kwargs):
        super().__init__(**kwargs)
        self.t0 = t0.to(u.Myr).value
    
    def evaluate(self, t):
        factor = t > self.t0
        return super().evaluate(t) * factor

    
class Density:
    """ Parent class to represent time-variable densities. 
    """
    def __init__(self, rho):
        self.rho = rho
        assert type(rho) == u.Quantity, "Density must be a Quantity"

    def evaluate(self, t):
        return self.rho.to(u.g / u.cm**3)
    

class ExponentialDensity(Density):
    def __init__(self, rho, t0, width):
        super().__init__(rho)
        self.t0 = t0
        self.width = width

        assert type(t0) == u.Quantity, "t0 must be a Quantity"
        assert type(width) == u.Quantity, "width must be a Quantity"

    def evaluate(self, t):
        return self.rho * np.exp(-((t - self.t0) / self.width) ** 2)
    

class InterpolatedDensity(Density):
    """ Density class that returns an interpolated value for the density at a given time t.
        Helpful for if you are determining the density through a GALA orbit and want it to be as realistic
        as possible.
    """
    def __init__(self, rho_interp):
        super().__init__(0)         # The "rho" can just be 0 here, as it will be overwritten by the interpolation

        self.rho_interp = rho_interp
        

    def evaluate(self, t):
        return np.interp(t, self.t, self.values)
    
