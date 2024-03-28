import numpy as np

from astropy import units as u
from astropy import constants as const

from gala import potential as gp

from gala.units import galactic

from astropy.cosmology import WMAP9 as cosmo


def JZ2023_1e12():
    """ 10^12 Msun NFW host halo from Zhu+2023
        (https://iopscience.iop.org/article/10.3847/1538-4357/acfe6f)

    Returns:
        gala potential: Gala NFW potential
    """
    c_NFW = 10**0.945
    M_host = 1e12 * u.Msun
 
    return  gp.NFWPotential.from_M200_c(M200=M_host, 
                                        c=c_NFW, 
                                        units=galactic)


def JZ2023_1e13():
    """ 10^13 Msun NFW host halo from Zhu+2023
        (https://iopscience.iop.org/article/10.3847/1538-4357/acfe6f)

    Returns:
        gala potential: Gala NFW potential
    """
    c_NFW = 10**0.85
    M_host = 1e13 * u.Msun
 
    return  gp.NFWPotential.from_M200_c(M200=M_host, 
                                        c=c_NFW, 
                                        units=galactic)

def JZ2023_1e14():
    """ 10^14 Msun NFW host halo from Zhu+2023
        (https://iopscience.iop.org/article/10.3847/1538-4357/acfe6f)
    
    Returns:
        gala potential: Gala NFW potential
    """
    c_NFW = 10**0.75
    M_host = 1e14 * u.Msun
 
    return  gp.NFWPotential.from_M200_c(M200=M_host, 
                                        c=c_NFW, 
                                        units=galactic)
