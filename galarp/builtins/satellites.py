from astropy import units as u
from gala import potential as gp

from gala.units import galactic


def JZ2023_Satellite():

    dm = gp.BurkertPotential(rho=5.93e-25 * u.g / u.cm**3, r0=11.87 * u.kpc, units=galactic)

    stars = gp.MiyamotoNagaiPotential(m=10**9.7 * u.Msun, a=2.5 * u.kpc, b=0.5 * u.kpc, units=galactic)
    gas = gp.MiyamotoNagaiPotential(m=10**9.7 * u.Msun, a=3.75 * u.kpc, b=0.75 * u.kpc, units=galactic)

    return dm + stars + gas


def NA2023_Satellite():



    dm = gp.BurkertPotential(rho=5.93e-25 * u.g / u.cm**3, r0 = 17.36 * u.kpc, units=galactic)
    stars = gp.MiyamotoNagaiPotential(1e11 * u.Msun, a= 5.94 * u.kpc, b = 0.58 * u.kpc, units=galactic)
    gas = gp.MiyamotoNagaiPotential(1e10 * u.Msun, a= 10.1 * u.kpc, b = 0.87 * u.kpc, units=galactic)
