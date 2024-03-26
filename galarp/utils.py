# Plotting routines and convenience functions.

# 3rd Party
from astropy import constants as c
from astropy import units as u
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.patches as patches

from scipy.interpolate import interp1d

import pickle


__all__ = [
    "ellipse_coords",
    "plot_axis_symbol",
    "plot_wind_vector",
    "get_orbit_data",
    "plot_rotation_vector",
    "plot_disk",
    "v_circ",
    "velocity",
    "gen_mass_profile",
    "m0",
    "diagnostic_3ax",
]


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


#############################################
######## Plotting routines ##################
#############################################


def ellipse_coords(x, y, a, b, theta, num_points=100, b_is_ellipticity=False):
    # Generate angles for sampling
    angles = np.linspace(0, 2 * np.pi, num_points)

    # If b is supplied as an ellipticity (and not as the semiminor axis directly, calculate semiminor axis)
    if b_is_ellipticity:
        b = a * (1 - b)

    # Parametric equation for the ellipse
    x_coords = (
        x + a * np.cos(angles) * np.cos(theta) - b * np.sin(angles) * np.sin(theta)
    )
    y_coords = (
        y + a * np.cos(angles) * np.sin(theta) + b * np.sin(angles) * np.cos(theta)
    )
    z_coords = np.zeros(num_points)
    # Return the sampled coordinates as a NumPy array
    return x_coords, y_coords, z_coords


def normalize_vector(vec):
    return vec / np.sqrt(np.sum(vec**2))


def plot_axis_symbol(ax, val, x0, y0, color, length):
    if abs(val) > 1e-1:
        if val > 0:
            ax.text(
                x0,
                y0,
                r"$\odot$",
                ha="center",
                va="center",
                size=length * 3,
                color=color,
            )
        else:
            ax.text(
                x0,
                y0,
                r"$\otimes$",
                ha="center",
                va="center",
                size=length * 3,
                color=color,
            )


def _limits_and_span(ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xspan, yspan = xmax - xmin, ymax - ymin
    return xmin, xmax, xspan, ymin, ymax, yspan


def plot_wind_vector(vector, ax, length=5, loc=(0, 0, 0), color="black"):
    vector = normalize_vector(vector) * length
    length = 5

    x0, y0, z0 = loc

    dx, dy, dz = vector * length

    xmin, xmax, xspan, ymin, ymax, yspan = _limits_and_span(ax[0])

    if abs(dx) > 1e-1 or abs(dy) > 1e-1:
        ax[0].arrow(x0, y0, dx, dy, head_width=1, color=color, zorder=20)
    else:
        plot_axis_symbol(ax[0], dy, x0, y0, color=color, length=length)

    xmin, xmax, xspan, ymin, ymax, yspan = _limits_and_span(ax[1])

    if abs(dx) > 1e-1 or abs(dz) > 1e-1:
        ax[1].arrow(x0, z0, dx, dz, head_width=1, color=color, zorder=20)
    else:
        plot_axis_symbol(ax[1], dy, x0, z0, color=color, length=length)

    xmin, xmax, xspan, ymin, ymax, yspan = _limits_and_span(ax[2])
    if abs(dy) > 1e-1 or abs(dz) > 1e-1:
        ax[2].arrow(y0, z0, dy, dz, head_width=1, color=color, zorder=20)
    else:
        plot_axis_symbol(ax[2], dx, y0, z0, color=color, length=length)


def plot_rotation_vector(ax, color="red", clockwise=False):
    xmin, xmax, xspan, ymin, ymax, yspan = _limits_and_span(ax)
    x0, y0 = xmax - xspan / 4, ymin + yspan / 5

    if clockwise:
        dx, dy = -xspan / 5, -yspan / 5
    else:
        dx, dy = xspan / 5, xspan / 5

    style = "Simple, tail_width=0.8, head_width=8, head_length=8"
    kw = dict(arrowstyle=style, color=color)
    a3 = patches.FancyArrowPatch(
        (x0, y0), (x0 + dx, y0 + dy), zorder=20, connectionstyle="arc3,rad=.25", **kw
    )
    ax.add_patch(a3)


def plot_disk(ax, R, lw=1, color="black"):
    """Plot a disk of a given radius in a 3-panel plot"""
    ell_xs, ell_ys, ell_zs = ellipse_coords(0, 0, R, R, 0)
    ax[0].plot(ell_xs, ell_ys, lw=lw, color=color)
    ax[1].plot(ell_xs, ell_zs, lw=lw, color=color)
    ax[2].plot(ell_ys, ell_zs, lw=lw, color=color)


def diagnostic_3ax(rpsim, length=3, wind_vec_loc=(0, 0, 0), outname=None):
    """Diagnostic plot for 3D wind vector and its projections onto the xy, xz, and yz planes.

    Args:
        rpsim (_type_): _description_
        length (int, optional): _description_. Defaults to 3.
        wind_vec_loc (tuple, optional): _description_. Defaults to (0,0,0).
    """

    fig, ax = plt.subplots(1, 3, facecolor="white", figsize=(12, 4))

    ell_xs, ell_ys, ell_zs = ellipse_coords(0, 0, 10, 10, 0, 100)

    ax[0].plot(ell_xs, ell_ys, color="black", lw=1)
    ax[1].plot(ell_xs, ell_zs, color="black", lw=1)
    ax[2].plot(ell_ys, ell_zs, color="black", lw=1)

    vector = rpsim.wind.normalized().value

    x0 = y0 = z0 = 0
    x1, y1, z1 = vector * length

    plot_wind_vector(vector, ax, length=length, loc=wind_vec_loc)

    xlabels, ylabels = ["x", "x", "y"], ["y", "z", "z"]

    for i, axis in enumerate(ax.flatten()):
        axis.set_xlim(-12, 12)
        axis.set_ylim(-12, 12)
        axis.set_xlabel(xlabels[i])
        axis.set_ylabel(ylabels[i])

    plt.tight_layout()
    if outname is not None:
        plt.savefig(outname, dpi=200)
    else:
        plt.show()
