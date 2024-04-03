"""Classes for implementing orbits under ram pressure acceleration."""

# Third Party
import os
import pickle
import astropy.units as u
import numpy as np

from .postprocessing import utils

# This module
from . import shadows, winds, densities

import gala.dynamics as gd
from gala.units import galactic


# import integrate as gi            # If debugging
import gala.integrate as gi


__all__ = ["F_RPS", "RPSim", "OrbitContainer"]


def F_RPS(
    t, w, potential, shadow, wind, rho, r_cloud, m_cloud, wind_on=True, debug=False
):
    # position units are in kpc
    # velocity units here are in kpc/Myr
    x, y, z, vx, vy, vz = w
    q = np.stack((x, y, z), axis=1)
    p = np.stack((vx, vy, vz), axis=1)

    # compute acceleration from potential:
    _t = np.array([0.0])
    acc = -potential._gradient(q, _t).T

    # Compute acceleration from ram pressure
    # a_ram = pi * rho * v_perp^2 * r_cloud^2 / m_cloud
    # TODO remove units altogether
    if wind_on:
        v_perp = p - wind.evaluate(t)
        a_ram = (np.pi * rho.evaluate(t) * r_cloud**2 / m_cloud).to(1 / u.kpc).value * (
            v_perp**2
        )

        a_ram = a_ram.T

        # If wind is on and shadow exists, apply shadow to appropriate particles
        if shadow is not None:
            shadow = shadow.evaluate(q, t).T
            a_ram *= shadow

        acc += a_ram

    return np.vstack((p.T, acc))


class RPSim:
    def __init__(
        self,
        wind,
        potential,
        rho_icm=None,
        shadow=None,
        potential_name="",
        method=F_RPS,
    ):
        self.method = method

        self.wind = wind
        self.potential = potential

        self.shadow = shadow

        self.potential_name = potential_name

        self.rho_icm = rho_icm

        self.inclination = wind.inclination()
        if (
            np.rad2deg(self.inclination) < 10
            and type(self.shadow) != shadows.EdgeOnShadow
        ):
            print(
                "Warning: Inclination is less than 10 degrees. You should be using an EdgeOnShadow."
            )

        self.sim_results = []

    def run(
        self,
        particles,
        rho_icm=2e-27 * (u.g / u.cm**3),
        t0=100 * u.Myr,
        r_cloud=50 * u.pc,
        m_cloud=1e5 * u.Msun,
        wind_on=True,
        integration_time=500 * u.Myr,
        dt=5 * u.Myr,
        printout=True,
        wind=None,
        outdir=None,
        debug=False,
    ):

        if type(rho_icm) is u.Quantity:
            self.rho_icm = densities.Density(rho_icm)
        else:
            self.rho_icm = rho_icm

        # Allow for user to switch out wind in the run method
        if wind is not None:
            self.wind = wind

        if printout:
            printout_width = 80
            print("".center(printout_width, "-"))
            print(f" Running GALA sim with  {self.wind.wind_strength():.2e}  wind at  {self.wind.inclination():.2f}  degrees ".center(
                    printout_width, "-"))
            print(f" Running for {integration_time}  at a timestep of  {dt}  ({integration_time / dt:.1f} steps) ".center(
                    printout_width, "-"))
            print("".center(printout_width, "-"))
            print(f" Container contains {len(particles.container)} particles ".center(
                    printout_width, "-"))
            print(f" Each cloud has a mass of  {m_cloud.value:.1e} solar masses  and a radius of  {r_cloud:.2f}   ".center(
                    printout_width, "-"))
            print(f" This corresponds to a density of {(m_cloud / r_cloud **3).to(u.g/u.cm**3):.2e}   ".center(
                    printout_width, "-"))
            print("".center(printout_width, "-"))

        # potential, shadow, wind_vec, rho, r_cloud, m_cloud, damping=1, wind_on=True)

        integrator = gi.RK5Integrator(
            self.method,
            func_args=(
                self.potential,
                self.shadow,
                self.wind,
                self.rho_icm,
                r_cloud,
                m_cloud,
                wind_on,
            ),
            func_units=galactic,
            progress=not debug,
        )

        orbits = integrator.run(
            gd.combine(particles.container), dt=dt, t1=0, t2=integration_time
        )
        self.sim_results.append(orbits)

        metadata = {
            "WIND": self.wind,
            "POTENTIAL": self.potential,
            "POTENTIAL_NAME": self.potential_name,
            "SHADOW": self.shadow,
            "RHO_ICM": self.rho_icm,
            "PARTICLE_GRID": particles,
            "R_CLOUD": r_cloud,
            "M_CLOUD": m_cloud,
            "WIND_ON": wind_on,
            "INTEGRATION_TIME": integration_time,
            "DT": dt,
        }

        out = OrbitContainer(orbits, units=galactic, metadata=metadata)

        if outdir is not None:
            os.makedirs(outdir, exist_ok=True)
            save_orbits(out, outdir=outdir)

        return out

    def __repr__(self) -> str:
        return f"Ram Pressure Sim \n Wind: \t\t{self.wind.__repr__()}\n Shadow: \t{self.shadow.__repr__()}\n Potential: \t{self.potential}"


def autogen_fn(md, suffix="orbits"):
    fn = ""

    if md["POTENTIAL_NAME"] != "":
        fn += f'{md["POTENTIAL_NAME"]}_'
    if md["WIND"] is not None:
        wind = md["WIND"]
        fn += f"WIND_{wind.wind_strength().to(u.km/u.s).value:.0f}_{np.rad2deg(wind.inclination()):.0f}_"
    if md["SHADOW"] is not None:
        shadow = md["SHADOW"]
        fn += f"_{shadow.shadow_name}_{shadow.damping}DAMP_"
    fn += f'{md["M_CLOUD"].value:.0e}Msun_{md["R_CLOUD"].value:.0f}pc'

    fn += f".{suffix}"

    return fn


def save_orbits(orbits, name="auto", outdir=""):
    metadata = orbits.metadata
    if name == "auto":
        name = autogen_fn(metadata)

    with open(outdir + name, "wb") as f:
        pickle.dump(orbits, f)


class OrbitContainer:
    """ Container for storing orbits and metadata from a completed GalaRP run.
    """
    def __init__(self, data, units=None, metadata={}):
        self.data = data
        self.units = units

        self.metadata = metadata

    def plot(self, plot_3d=False, plot_orbits=False):
        if plot_3d:
            utils.k3d_plot([self])

        if plot_orbits:
            utils.plot_orbits(
                self.data, wind=self.metadata["WIND"], shadow=self.metadata["SHADOW"]
            )
