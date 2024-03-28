# Plotting routines and convenience functions.

# 3rd Party
from astropy import constants as c
from astropy import units as u
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.patches as patches

from scipy.interpolate import interp1d

import pickle


__all__ = ['v_circ', 'velocity', 'gen_mass_profile', 'm0', 'get_orbit_data', 'pickle_obj']


#############################################
######## Convenience Functions ##############
#############################################


def v_circ(M, R):
    return np.sqrt(c.G * M / R).to(u.km / u.s)


def velocity(mass_profile, R, theta):
    vx = -v_circ(mass_profile(R) * u.M_sun, R) * np.sin(theta)
    vy = v_circ(mass_profile(R) * u.M_sun, R) * np.cos(theta)
    return [vx.value, vy.value, 0] * (u.km / u.s)


def gen_mass_profile(potential, lims=[2e-2, 80]):
    pos = np.zeros((3, 100)) * u.kpc
    pos[0] = np.linspace(lims[0], lims[1], 100) * u.kpc
    m_profile = potential.mass_enclosed(pos)
    return interp1d(pos[0], m_profile, bounds_error=False, fill_value=0)


def m0(rho_d0):
    """Return the normalization constant for the Burkert Dark Matterpotential

    Args:
        rho_d0 (float): Central density of the Burkert potential (in g/cm^3)
    """
    return 1e9 * (rho_d0 / 1.46e-24) ** (-7 / 2)


def get_orbit_data(o):
    pos, vel = o.pos, o.vel

    x, y, z = pos.xyz.value
    x, y, z = x.T, y.T, z.T

    vx, vy, vz = vel.d_xyz.to(u.km / u.s).value
    vx, vy, vz = vx.T, vy.T, vz.T

    return x, y, z, vx, vy, vz


def pickle_obj(obj, name="obj.out"):
    with open(name, "wb") as f:
        pickle.dump(obj, f)
