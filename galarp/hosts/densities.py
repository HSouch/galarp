from astropy import units as u


def SphericalBetaModel(r, n0 = 0.0121 * u.cm ** -3, 
                       r_c=25 * u.kpc, 
                       beta = 0.655):
    """ The spherical beta model for the Milky Way halo.

    Args:
        r (astropy.quantity): The input radius
        n0 (Astropy quantity, optional): _description_. Defaults to 0.0121 .
        r_c (astropy.quantity, optional): _description_. Defaults to 25*u.kpc.
        beta (float, optional): _description_. Defaults to 0.655.

    Returns:
        astropy.quantity : The number density at the input radius
    """
    return n0 * (1 + (r / r_c) ** 2) ** (-3 * beta / 2)


def MB2015ModifiedBeta(r, n0 = 0.0121 * u.cm ** -3, 
                       r_c=25 * u.kpc,
                       beta =0.655):
    """ The Miller and Bregman (2015) modified beta model for the Milky Way halo.

    Args:
        r (astropy.quantity): _description_
        n0 (Astropy quantity, optional): _description_. Defaults to 0.0121 .
        r_c (_type_, optional): _description_. Defaults to 25*u.kpc.
        beta (float, optional): _description_. Defaults to 0.655.

    Returns:
        _type_: _description_
    """
    return n0 * r_c ** (3 * beta)  / r ** (3 * beta)