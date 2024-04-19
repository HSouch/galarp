from astropy import units as u
from gala.units import galactic
import numpy as np

from ..utils import gen_mass_profile

from .. import builtins

from .. import winds, particle_grids, shadows, rampressure


__all__ = ['ExampleUniformEvent', 'ExampleLorentzianEvent']


def ExampleUniformEvent(inc=45., strength=600.):
    """ Returns a set of orbits for a simple uniform wind event.

    Args:
        inc (float, optional): Inclination of the wind vector (in deg). Defaults to 45.
        strength (float, optional): Strength of the wind vector (in km/s). Defaults to 600.

    Returns:
        galarp OrbitContainer: The orbit output.
    """
    pot = builtins.JZ2023_Satellite()
    mass_profile = gen_mass_profile(pot)
    
    wind = winds.RPWind(units=galactic)
    wind.init_from_inc(inclination=np.deg2rad(inc), strength= strength *u.km/u.s)
    
    shadow = shadows.UniformShadow()
    shadow.init_from_wind(wind)

    particles = particle_grids.UniformGrid(Rmax= 10, n_particles= 20, z_start= 0.5)
    particles.generate(mass_profile)

    sim = rampressure.RPSim(wind=wind, potential=pot, shadow=shadow, rho_icm=1e-26 * u.g/ u.cm**3)
    orbits = sim.run(particles, rho_icm=1e-26 * u.g/ u.cm**3, integration_time=1000 * u.Myr, dt=2*u.Myr)
    return orbits


def ExampleLorentzianEvent(inc=45, strength=600):
    """ Returns a set of orbits for a simple Lorentzian wind event.

    Args:
        inc (float, optional): Inclination of the wind vector (in deg). Defaults to 45.
        strength (float, optional): Strength of the wind vector (in km/s). Defaults to 600.

    Returns:
        galarp OrbitContainer: The orbit output.
    """

    pot = builtins.JZ2023_Satellite()
    mass_profile = gen_mass_profile(pot)
    wind = winds.LorentzianWind(t0=500 * u.Myr, width = 500 * u.Myr, units=galactic)
    wind.init_from_inc(inclination=np.deg2rad(inc), strength= strength *u.km/u.s)
    
    shadow = shadows.UniformShadow()
    shadow.init_from_wind(wind)

    particles = particle_grids.UniformGrid(Rmax= 10, n_particles= 20, z_start= 0.5)
    particles.generate(mass_profile)

    sim = rampressure.RPSim(wind=wind, potential=pot, shadow=shadow, rho_icm=1e-26 * u.g/ u.cm**3)
    orbits = sim.run(particles, rho_icm=1e-26 * u.g/ u.cm**3, integration_time=1000 * u.Myr, dt=2*u.Myr)
    return orbits