# Plotting routines and convenience functions.

# 3rd Party
from astropy.cosmology import WMAP9 as cosmo
from astropy import constants as c
from astropy import units as u
import numpy as np


from scipy.interpolate import interp1d

import pickle


__all__ = ['v_circ', 'velocity', 'gen_mass_profile', 'm0', 'get_orbit_data', 'pickle_obj', 'Rvir', 'R100', 'R500',
           'rotate', 'rotate_yaxis']


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


def get_orbit_data(o, transposed=True):
    pos, vel = o.pos, o.vel

    x, y, z = pos.xyz.value
    if transposed:
        x, y, z = x.T, y.T, z.T

    vx, vy, vz = vel.d_xyz.to(u.km / u.s).value
    if transposed:
        vx, vy, vz = vx.T, vy.T, vz.T

    return x, y, z, vx, vy, vz


def pickle_obj(obj, name="obj.out"):
    with open(name, "wb") as f:
        pickle.dump(obj, f)


def Rvir(mass):
    """ Returns the virial radius of the potential """
    
    mass = mass.to(u.Msun)

    rho_c = 3 * cosmo.H(0.)**2 / (8*np.pi*c.G)
    rho_c = rho_c.to(u.Msun/u.kpc**3)

    return np.cbrt(mass / (200*rho_c.value) / (4/3*np.pi))

def R100(mass):
    """ Returns the radius at which the density is 100 times the critical density """
    
    mass = mass.to(u.Msun)

    rho_c = 3 * cosmo.H(0.)**2 / (8*np.pi*c.G)
    rho_c = rho_c.to(u.Msun/u.kpc**3)

    return np.cbrt(mass / (100*rho_c.value) / (4/3*np.pi))

def R500(mass):
    """ Returns the radius at which the density is 500 times the critical density """
    
    mass = mass.to(u.Msun)

    rho_c = 3 * cosmo.H(0.)**2 / (8*np.pi*c.G)
    rho_c = rho_c.to(u.Msun/u.kpc**3)

    return np.cbrt(mass / (500*rho_c.value) / (4/3*np.pi))


def rotate(vec, alpha=np.deg2rad(0), beta=np.deg2rad(0), gamma=np.deg2rad(0)):
    """ Rotate a set of positions by the given angles

    Args:
        vec (array_like): Array of positions to rotate 
        alpha (float): Rotation angle about the x-axis
        beta (float): Rotation angle about the y-axis
        gamma (float): Rotation angle about the z-axis 

    Returns:
        _type_: _description_
    """

    rot_matrix_1 = [np.cos(beta)*np.cos(gamma), np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma), 
                    np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma)]
    rot_matrix_2 = [np.cos(beta)*np.sin(gamma), np.sin(alpha)*np.sin(beta)*np.sin(gamma) + np.cos(alpha)*np.cos(gamma),
                     np.cos(alpha)*np.sin(beta)*np.sin(gamma) - np.sin(alpha)*np.cos(gamma)]
    rot_matrix_3 = [-np.sin(beta), np.sin(alpha)*np.cos(beta), np.cos(alpha)*np.cos(beta)]
    rot_matrix = np.array([rot_matrix_1, rot_matrix_2, rot_matrix_3])
    return np.dot(rot_matrix, vec)


def rotate_yaxis(vec, beta=np.deg2rad(0)):
    """ Rotate a set of positions by the given angles

    Args:
        vec (array_like): Array of positions to rotate 
        beta (float): Rotation angle about the y-axis

    Returns:
        _type_: _description_
    """

    rot_matrix_1 = [np.cos(beta), 0, np.sin(beta)]
    rot_matrix_2 = [0, 1, 0]
    rot_matrix_3 = [-np.sin(beta), 0, np.cos(beta)]
    rot_matrix = np.array([rot_matrix_1, rot_matrix_2, rot_matrix_3])
    return np.dot(rot_matrix, vec)


def handle_quantity(a, unit):
    if not isinstance(a, u.Quantity):
        return a * unit
    else:
        return a.to(unit)
    