
from astropy.table import Table
import astropy.units as u
from astropy.units import Quantity
import numpy as np

from gala.units import galactic

from scipy.interpolate import interp1d

import warnings

__all__ = ["RPWind", "ConstantWind", "LorentzianWind", "StepFunctionWind", "InterpolatedStrengthWind"]



class RPWind:
    """ Base class for ram pressure winds.

    This class represents a base class for ram pressure winds. It provides methods to initialize the wind properties,
    evaluate the wind at a given time, update the wind properties, and provide a string representation of the wind.

    Parameters:
        inclination (float): The inclination angle of the wind in radians.
        strength (Quantity or float): The strength of the wind. If a Quantity object is provided, it should have units of length/time.
            If a float is provided, it is assumed to be in units of length/time.
        units (dict, optional): A dictionary containing the units for length and time. If not provided, default galactic units are used.
        **kwargs: Additional keyword arguments.

    Attributes:
        units (dict): A dictionary containing the units for length and time.
        inclination (float): The inclination angle of the wind in radians.
        strength (Quantity): The strength of the wind in units of length/time.
        unit_vector (ndarray): A numpy array representing the unit vector of the wind direction.

    Methods:
        evaluate(t): Evaluates the wind at a given time.
        evaluate_arr(ts): Evaluates the wind at multiple times.
        update(inclination, strength): Updates the wind properties.
        __repr__(): Returns a string representation of the wind.

    """

    def __init__(self, inclination=0, strength=800 * u.km / u.s, units=None, **kwargs):
        self.units = galactic if units is None else units
        self.inclination = inclination

        if isinstance(strength, (u.Quantity)):
            self.strength = strength.to(self.units["length"] / self.units["time"])
        else:
            self.strength = strength * self.units["length"] / self.units["time"]

        if self.strength > 5 * u.kpc / u.Myr:
            warnings.warn(f"The wind strength ({self.strength}) is quite high. Are your units correct?")

        self.unit_vector = np.array([np.sin(self.inclination), 0, np.cos(self.inclination)])

    def evaluate(self, t):
        return NotImplementedError("This method must be implemented by the subclass")
    
    def evaluate_arr(self, ts):
        try:
            return np.array([self.evaluate(t) for t in ts])
        except Exception:
            return np.array([self.evaluate(t) for t in ts.value])
        
    def update(self, inclination=None, strength=None):
        if inclination is not None:
            self.inclination = inclination
            self.unit_vector = np.array([np.sin(self.inclination), 0, np.cos(self.inclination)])
        if strength is not None:
            self.strength = strength

    def __repr__(self) -> str:
        return f'<RP Wind Vector: Inc={np.rad2deg(self.inclination)}>'
    


class ConstantWind(RPWind):

    def evaluate(self, t):
        return self.unit_vector * self.strength.value


class LorentzianWind(ConstantWind):
    """ A ram pressure wind shaped by the Lorentzian profile, which is a popular profile for winds for a 
        galaxy orbiting through a cluster.

    """
    def __init__(self, t0=200 * u.Myr, width=200 * u.Myr, **kwargs):
        """
            Parameters:
            - t0 (Quantity): The time at peak wind strength. Default is 200 Myr.
            - width (Quantity): The width of the Lorentzian profile. Default is 200 Myr.
        """
        super().__init__(**kwargs)
        
        self.t0 = t0.to(self.units["time"])
        self.width = width.to(self.units["time"])
    
    def evaluate(self, t):

        return super().evaluate(t) * 1 / ((2 * (t - self.t0.value) / self.width.value) ** 2 + 1)
    

class StepFunctionWind(ConstantWind):
    """ Representation of a "step function wind" where the wind turns on at the time t0.
        Useful for testing that time variability is working, or for allowing the particles to evolve in
        isolation before a constant wind hits.

    """
    def __init__(self, t0= 100 * u.Myr, **kwargs):
        """
            Parameters:
            - t0 (Quantity): The time at which the wind turns on. Default is 100 Myr.
        """
        super().__init__(**kwargs)
        self.t0 = t0.to(self.units["time"])
    
    def evaluate(self, t):
        return super().evaluate(t) * (t > self.t0)


class InterpolatedStrengthWind(RPWind):
    """
    A class representing an interpolated wind with strength.

    This class extends the RPWind class and provides methods to evaluate the wind strength at a given time.

    Attributes:
    - interp (callable): The interpolation function used to evaluate the wind strength.

    Methods:
    - __init__(self, interp=None, **kwargs): Initializes an InterpolatedStrengthWind object.
    - evaluate(self, t): Evaluates the wind strength at a given time.
    - from_table(fn, t_key, vel_keys, format="ascii", t_units=u.s, v_units=u.cm/u.s, verbose=False, **kwargs): 
      Creates an InterpolatedStrengthWind object from an input data table.

    Usage:
    >>> wind = InterpolatedStrengthWind.from_table("wind_data.txt", "time", 
        ["velocity_x", "velocity_y", "velocity_z"], units=galactic, inc=np.deg2rad(45))
    >>> print(wind.evaluate(10 * u.Myr))
    """

    def __init__(self, interp=None, **kwargs):
        super().__init__(**kwargs)
        self.interp = interp

    def evaluate(self, t):
        """
        Evaluates the wind strength at a given time.

        Parameters:
        - t: The time at which to evaluate the wind strength.

        Returns:
        - The wind strength at the given time.
        """
        return self.unit_vector * self.interp(t)
    
    @staticmethod
    def from_table(fn, t_key, vel_keys, format="ascii",
                   t_units=u.s, v_units=u.cm/u.s,
                   verbose=False, **kwargs):
        """
        Create an InterpolatedStrengthWind object from a table.

        Parameters:
        - fn (str): The filename of the table.
        - t_key (str): The key for the time column in the table.
        - vel_keys (str or list): The key(s) for the velocity column(s) in the table.
        - format (str, optional): The format of the table. Default is "ascii".
        - t_units (astropy.units.Unit, optional): The units for the time column. Default is u.s.
        - v_units (astropy.units.Unit, optional): The units for the velocity column(s). Default is u.cm / u.s.
        - verbose (bool, optional): Whether to print verbose output. Default is False.
        - **kwargs: Additional keyword arguments to be passed to the InterpolatedWind constructor.

        Returns:
        - InterpolatedWind: An InterpolatedWind object created from the table.
        Usage:
                >>> wind = InterpolatedWind.from_table("wind_data.txt", "time", 
                    ["velocity_x", "velocity_y", "velocity_z"], units=galactic, inc=np.deg2rad(45))
                >>> print(wind.evaluate(10 * u.Myr))
        """
        t = Table.read(fn, format=format)

        units = kwargs.get("units", galactic)

        if verbose:
            print(f'Loaded table with {len(t)} rows and keys: {t.keys()}')
        
        ts = t[t_key] * t_units.to(units["time"])

        if not isinstance(vel_keys, list):
            vel_keys = [vel_keys]
        
        vels = np.array([t[key] for key in vel_keys])

        vel_tot = np.sqrt(np.sum(vels ** 2, axis=0)) * v_units.to(units["length"] / units["time"])

        interp = interp1d(ts, vel_tot, bounds_error=False, fill_value="extrapolate")

        return InterpolatedStrengthWind(interp=interp, **kwargs)