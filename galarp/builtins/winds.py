
import numpy as np
from scipy.interpolate import interp1d
from gala.units import galactic
from astropy import units as u

from ..winds import InterpolatedStrengthWind



def RB2006_Wind(inc = 45, peak = 800):
    
    peak = peak * (u.km/u.s).to(u.kpc / u.Myr)

    i_x, i_y = np.linspace(0, 50, 50), np.linspace(0, peak, 50)

    i_x = np.concatenate([i_x, np.linspace(50, 1000, 950)])
    i_y = np.concatenate([i_y, peak * np.ones(950)])

    interp = interp1d(i_x, i_y, bounds_error=False, fill_value="extrapolate")

    wind = InterpolatedStrengthWind(interp=interp, inclination=np.deg2rad(inc), units=galactic)

    return wind
