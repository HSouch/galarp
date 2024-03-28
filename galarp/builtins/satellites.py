from astropy import units as u
from gala import potential as gp

from gala.units import galactic


def JZ2023_Satellite():

    dm = gp.BurkertPotential(rho=5.93e-25 * u.g / u.cm**3, r0=11.87 * u.kpc, units=galactic)

    stars = gp.MiyamotoNagaiPotential(m=10**9.7 * u.Msun, a=2.5 * u.kpc, b=0.5 * u.kpc, units=galactic)
    gas = gp.MiyamotoNagaiPotential(m=10**9.7 * u.Msun, a=3.75 * u.kpc, b=0.75 * u.kpc, units=galactic)

    return dm + stars + gas

