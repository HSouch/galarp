
from astropy.table import Table
import astropy.units as u
from astropy.units import Quantity
import numpy as np

from gala.units import galactic

from scipy.interpolate import interp1d



__all__ = [
    "RPWind",
    "LorentzianWind",
    "StepFunctionWind",
    "InterpolatedWind",
]


class RPWind:
    """Class to represent a ram pressure wind

    Usage:
        ```
        wind = gn.RPWind([50, 0, 50] * (u.km / u.s), units=gn.galactic)  # Initializes a wind at 45 degrees in the x-z plane.

        wind = gn.RPWind(units=galactic)                                 # Achieves the same thing.
        wind.init_from_inc(np.deg2rad(45), 100 * u.km  / u.s)
        ```
    """

    def __init__(self, vector=None, units=None, **kwargs):
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
                self.vector *= self.units["length"] / self.units["time"]

    def evaluate(self, t):
        """Return the wind vector at time t. For the default wind vector, this is JUST the vector converted to kpc/Myr"""
        return self.vector.to(u.kpc / u.Myr).value
    
    def evaluate_arr(self, ts):
        return np.array([self.evaluate(t) for t in ts])

    def initialize_vector(self):
        if self.units is not None:
            self.vector = self.vector.to(self.units["length"] / self.units["time"])

    def init_from_inc(self, inclination, strength):
        """Initialize the wind vector from an inclination and strength"""

        x = strength * np.cos(inclination)
        z = strength * np.sin(inclination)
        self.vector = Quantity([x, 0 * x.unit, z]).to(
            self.units["length"] / self.units["time"]
        )

    def wind_strength(self):
        """Return the length (strength) of the wind vector"""
        return np.sqrt(sum(self.vector**2))

    def normalized(self):
        """Return the normalized wind vector"""
        return self.vector / self.wind_strength()

    def vector_to_units(self, units):
        return self.vector.to(units)

    def vector_as_value(self):
        if type(self.vector) == Quantity:
            return self.vector.value
        else:
            return self.vector

    def inclination(self):
        x, y, z = self.vector_as_value()
        return np.arctan2(z, np.sqrt(x**2 + y**2))

    def __repr__(self):
        return f"<RP Wind Vector={self.vector}  Inclination={np.round(np.rad2deg(self.inclination()), 2):.2f}  >"


class LorentzianWind(RPWind):
    """Ram pressure wind that is damped by a Lorentzian profile.
        Reaches max (which is just the unadjusted RP wind value) at t0.

    Args:
        gn (_type_): _description_
    """

    def __init__(self, t0=0 * u.Myr, width=200 * u.Myr, **kwargs):
        super().__init__(**kwargs)
        self.t0 = t0.to(u.Myr).value  # Units are in Myr
        self.width = width.to(u.Myr).value  # Units are in Myr

    def evaluate(self, t):
        """Return the wind vector damped by a Lorentzian profile"""
        return super().evaluate(t) * 1 / ((2 * (t - self.t0) / self.width) ** 2 + 1)


class StepFunctionWind(RPWind):
    """A wind represented by a step function (no wind immediately into full wind at time t0)
    Good for testing that time variability is working.

    """

    def __init__(self, t0=0 * u.Myr, **kwargs):
        super().__init__(**kwargs)
        self.t0 = t0.to(u.Myr).value

    def evaluate(self, t):
        factor = t > self.t0
        return super().evaluate(t) * factor
    

class InterpolatedWind(RPWind):
    """Wind that is interpolated from a table of values. """

    def __init__(self, interp=None, units=galactic, **kwargs):
        """
            Args:
                interp (callable): Interpolation function that takes a time value as input and returns the wind velocity at that time.
                units (dict, optional): Dictionary specifying the units of the wind vector. Defaults to galactic units.
                inc (float, optional): Inclination angle of the wind vector in radians. Defaults to 0.
            """
        super().__init__(**kwargs)
        self.interp = interp
        self.inc = kwargs.get("inc", 0)

        self.unit_vector = [np.cos(self.inc), 0, np.sin(self.inc)]
    
    def evaluate(self, t):
        """ Evaluate the wind vector at a given time by multiplying the unit vector by the interpolated
            wind strength at that time.
        """
        return self.unit_vector * self.interp(t)

    @staticmethod
    def from_table(fn, t_key, vel_keys, format="ascii",
                   t_units = u.s, v_units = u.cm / u.s, 
                   verbose=False, **kwargs):
        """
        Create an InterpolatedWind object from a table.

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

        if verbose:
            print(f'Loaded table with {len(t)} rows and keys: {t.keys()}')

        ts = t[t_key] * t_units.to(u.Myr)

        if not isinstance(vel_keys, list): vel_keys = [vel_keys]

        vels = np.array([t[key] for key in vel_keys])
        v_tot = np.sqrt(np.sum(vels**2, axis=0)) * v_units.to(u.kpc/u.Myr)

        interp = interp1d(ts, v_tot, bounds_error=False, fill_value="extrapolate")

        return InterpolatedWind(interp=interp, **kwargs)
