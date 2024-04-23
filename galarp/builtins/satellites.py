from astropy import units as u
from gala import potential as gp

from gala.units import galactic


def JZ2023_Satellite():
    """ Satellite potential for Zhu_+2023
        See https://iopscience.iop.org/article/10.3847/1538-4357/acfe6f
    """

    dm = gp.BurkertPotential(rho=5.93e-25 * u.g / u.cm**3, r0=11.87 * u.kpc, units=galactic)

    stars = gp.MiyamotoNagaiPotential(m=10**9.7 * u.Msun, a=2.5 * u.kpc, b=0.5 * u.kpc, units=galactic)
    gas = gp.MiyamotoNagaiPotential(m=10**9.7 * u.Msun, a=3.75 * u.kpc, b=0.75 * u.kpc, units=galactic)

    return gp.CompositePotential(dm=dm, stars=stars, gas=gas)


def NA2023_Satellite():
    """ Satellite potential for Akerman+2023 and Akerman+2024
        See https://iopscience.iop.org/article/10.3847/1538-4357/acbf4d
    """
    dm = gp.BurkertPotential.from_r0(r0=17.36 * u.kpc, units=galactic)
    stars = gp.MiyamotoNagaiPotential(1e11 * u.Msun, a= 5.94 * u.kpc, b = 0.58 * u.kpc, units=galactic)
    gas = gp.MiyamotoNagaiPotential(1e10 * u.Msun, a= 10.1 * u.kpc, b = 0.87 * u.kpc, units=galactic)

    return gp.CompositePotential(dm=dm, stars=stars, gas=gas)


def RB2006_Satellite():
    """ Satellite potential for Roediger & Brüggen 2006
        See https://academic.oup.com/mnras/article/369/2/567/3101476
    """

    dm = gp.BurkertPotential.from_r0(23 * u.kpc, units=galactic)

    stars = gp.MiyamotoNagaiPotential(m=1e11 * u.Msun, a=4 * u.kpc, b=0.25 * u.kpc, units=galactic)
    bulge = gp.HernquistPotential(m=1e10 * u.Msun, c=0.4 * u.kpc, units=galactic)
    
    gas = gp.MiyamotoNagaiPotential(m=10**10 * u.Msun, a=7 * u.kpc, b=0.4 * u.kpc, units=galactic)

    return gp.CompositePotential(dm=dm, stars=stars, bulge=bulge, gas=gas)
