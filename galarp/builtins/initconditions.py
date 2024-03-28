from astropy import units as u
from ..hosts.initconditions import InitConditions


def JZ2023_1e12_IC():
    V200 = 143 
    V_phi = 0.655 * V200
    V_r = 0.832 * V200

    pos = [211, 0, 0] * u.kpc
    vel = [-V_r, V_phi, 0] * u.km / u.s

    return InitConditions(pos=pos, vel=vel)


def JZ2023_1e13_IC():
    V200 = 308 
    V_phi = 0.603 * V200
    V_r = 0.786 * V200

    pos = [455, 0, 0] * u.kpc
    vel = [-V_r, V_phi, 0] * u.km / u.s

    return InitConditions(pos=pos, vel=vel)


def JZ2023_1e14_IC():
    V200 = 663 
    V_phi = 0.53 * V200
    V_r = 0.782 * V200

    pos = [949, 0, 0] * u.kpc
    vel = [-V_r, V_phi, 0] * u.km / u.s

    return InitConditions(pos=pos, vel=vel)