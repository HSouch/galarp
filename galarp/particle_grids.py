import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt

from .postprocessing.plotting import general_plots
import gala.dynamics as gd
from . import utils


from tqdm import tqdm


__all__ = [
    "ParticleSet",
    "UniformGrid",
    "ExponentialGrid",
    "generate_exponential_positions",
]


class ParticleSet:
    def __init__(self, name="PG"):
        self.container = []  # List of PhaseSpacePosition objects
        self.name = name

        self.masses = None
        self.sigmas = None
        self.radii = None

    def generate(self):
        raise NotImplementedError(
            "This is an abstract class. Please use a \
            subclass."
        )

    def plot_phase_space(self, ax=None, color="black", outname=None):
        if ax is None:
            fig, ax = plt.subplots(1, 2, facecolor="white", figsize=(8, 4))
        for particle in self.container:
            ax[0].scatter(particle.x, particle.y, color=color)
            ax[1].scatter(particle.v_x, particle.v_y, color=color)

        plt.tight_layout()
        if outname is not None:
            plt.savefig(outname, dpi=200)
        else:
            plt.show()
    
    def get_xyz(self):
        x, y, z = [], [], []
        for particle in self.container:
            x.append(particle.x.value)
            y.append(particle.y.value)
            z.append(particle.z.value)
        return x, y, z


    def plot_grid(self, **kwargs):
        """ Plot the grid of particles in the xy and xz planes. """
        particle_size = kwargs.get("particle_size", 1)

        fig, ax = plt.subplots(1, 2, facecolor="white", figsize=(8, 4))
        xs, ys, zs = self.get_xyz()

        ax[0].scatter(xs, ys, color="black", s=particle_size)
        ax[1].scatter(xs, zs, color="black", s=particle_size)
        
        ax[0].set_xlabel("x [kpc]")
        ax[0].set_ylabel("y [kpc]")
        ax[1].set_xlabel("x [kpc]")
        ax[1].set_ylabel("z [kpc]")
        
        ax[0].set_xlim(-self.Rmax.value, self.Rmax.value)
        ax[0].set_ylim(-self.Rmax.value, self.Rmax.value)
        ax[1].set_xlim(-self.Rmax.value, self.Rmax.value)
        ax[1].set_ylim(-self.Rmax.value, self.Rmax.value)

        plt.tight_layout()
        plt.show()
    

    def plot_density(self, gridsize=40, outname=None, **kwargs):
        """ Generate hexbin plots of the x-y and x-z particle densities.

        Args:
            gridsize (int): Number of bins in the hexbin plots.
            outname (str): Filename for saving the output plot. Will just run plt.show() if None.
        """
        general_plots.plot_density(self.get_xyz(), gridsize=gridsize, outname=outname,  **kwargs)


    def save(self, outname):
        utils.pickle_obj(self.container, outname)


class UniformGrid(ParticleSet):
    """Class to initialize a grid of particles in the xy plane with circular
    velocities based on a given gravitational potential.

    Usage:
        ```
        mass_profile = gn.gen_mass_profile(potential)
        grid = gn.ParticleGrid()
        grid.generate(mass_profile)
        ```
    """

    def __init__(self, Rmax=10, n_particles=50, z_start=0.0, veldisp=10.0, 
                 masses=None, sigmas=None, radii=None):
        """Generate a grid of particles in the xy plane with circular
            velocities based on a given gravitational potential.
            Can also include a Gaussian velocity dispersion (helpful for
            including a scale height in the disk).

        Args:
            mass_profile (interp1d): The mass profile of the potential in scipy
            interp1d format. Defaults to None.
            Rmax (_type_): The maximum radius of the grid particles.
            n_particles (int, optional): Number of particles between -Rmax and
                Rmax. Defaults to 50.
            z_start (float, optional): The z-level to initialize the grid in.
                Defaults to 0.
            veldisp (float, optional): _description_. Defaults to 10.
            velocities (bool, optional): If false, particles have no velocities.
                Defaults to True.
        """

        super().__init__(name="UG")

        self.Rmax = Rmax
        self.n_particles = n_particles
        self.z_start = z_start
        self.veldisp = veldisp

        self.masses = masses
        self.sigmas = sigmas
        self.radii = radii

    def generate(self, mass_profile, velocities=True):
        particle_range = np.linspace(-self.Rmax, self.Rmax, self.n_particles) * u.kpc

        for x in particle_range:
            for y in particle_range:
                p_x0 = u.Quantity([x, y, self.z_start * u.kpc])
                R = np.sqrt(x**2 + y**2)
                if R > self.Rmax * u.kpc:
                    continue

                if velocities:
                    theta = np.arctan2(y, x)

                    # Initialize with circular velocity
                    p_v0 = utils.velocity(mass_profile, R, theta)
                    if self.veldisp is not None:
                        # Add random velocity dispersion
                        p_v0 += np.random.normal(scale=self.veldisp, size=3) * (
                            u.km / u.s
                        )
                else:
                    p_v0 = [0.0, 0.0, 0.0] * (u.km / u.s)

                w0 = gd.PhaseSpacePosition(pos=p_x0, vel=p_v0.to(u.kpc / u.Myr))
                self.container.append(w0)


class ExponentialGrid(ParticleSet):
    def __init__(self, h_R=4 * u.kpc, h_z=0.5 * u.kpc, n_particles=500, veldisp=10.0 * u.km / u.s,
        Rmax=None, zmax=None, masses=None, sigmas=None, radii=None):
        super().__init__(name="EG")

        self.h_R = h_R
        self.h_z = h_z

        self.Rmax = Rmax if Rmax is not None else self.h_R * 4
        self.zmax = zmax if zmax is not None else self.h_z * 4

        self.n_particles = n_particles
        self.veldisp = veldisp

        self.masses = masses
        self.sigmas = sigmas
        self.radii = radii

    def generate(self, mass_profile, velocities=True, positions=None):
        if positions is None:  # Generate positions if nothing given
            xs, ys, zs = generate_exponential_positions(
                h_R=self.h_R,
                h_z=self.h_z,
                n_particles=self.n_particles,
                Rmax=self.Rmax,
                zmax=self.zmax,
                outname=None,
            )
        elif isinstance(positions, str):    # If string assume this is a filename
            xs, ys, zs = self.load_positions(positions)
        else:
            xs, ys, zs = positions

        # Create PhaseSpacePositions for all points, and initialize velocities
        for i in range(self.n_particles):
            x, y, z = xs[i], ys[i], zs[i]

            p_x0 = u.Quantity([x, y, z] * u.kpc)

            R = np.sqrt(x**2 + y**2 + z**2) * u.kpc

            if velocities:
                theta = np.arctan2(y, x)
                vx, vy = self.velocity(mass_profile, R, theta)  # Initialize with circular velocity

                inc = np.arctan2(z, np.sqrt(x**2 + y**2))
                vtot = np.sqrt(vx**2 + vy**2)
                vz = (vtot * np.sin(inc))  # TODO improve this (better than nothing but not ideal. Leads to "sloshing")

                p_v0 = u.Quantity([vx, vy, vz])

                if self.veldisp is not None:
                    p_v0 += np.random.normal(scale=self.veldisp.value, size=3) * (u.km / u.s)  # Add random vel disp.
            else:
                p_v0 = [0.0, 0.0, 0.0] * (u.km / u.s)

            w0 = gd.PhaseSpacePosition(pos=p_x0, vel=p_v0.to(u.kpc / u.Myr))
            self.container.append(w0)

    def load_positions(self, filename):
        """Create the grid of particles from a file of positions (to save on time not needed to rerun the MC sim)
        Args:
            filename (numpy savetxt file): Text file of positions
        """
        positions = np.load(filename)
        return positions

    def velocity(self, mass_profile, R, theta):
        """ Get the velocity of a particle in the xy plane given an angle and the mass profile
        
        Args:
            mass_profile (interp1d): The mass profile of the potential in scipy interp1d format.
            R (float): The radius of the particle.
            theta (float): The angle of the particle.
        """
        vx = -utils.v_circ(mass_profile(R) * u.M_sun, R) * np.sin(theta)
        vy = utils.v_circ(mass_profile(R) * u.M_sun, R) * np.cos(theta)
        return vx, vy


def generate_exponential_positions(
    n_particles,
    h_R=4 * u.kpc,
    h_z=0.5 * u.kpc,
    Rmax=15 * u.kpc,
    Rmin=0.1 * u.kpc,
    zmax=2 * u.kpc,
    outname=None,
):
    def remap(value, a, b):
        return (value * (b - a)) + a

    n_0 = (
        1 / (4 * np.pi * h_R**2 * h_z).value
    )  # This normalizes the probability distribution

    h_R, h_z = h_R.to(u.kpc).value, h_z.to(u.kpc).value
    rmax = Rmax.to(u.kpc).value
    rmin = Rmin.to(u.kpc).value
    zmax = zmax.to(u.kpc).value

    Rs = []
    zs = []

    for n in tqdm(range(n_particles)):
        good = False
        while not good:  # TODO maybe add a max_tries here with a warning message
            rand_R, rand_z = np.random.random(2)
            R, z = remap(rand_R, rmin, rmax), remap(rand_z, -zmax, zmax)

            v = n_0 * np.exp(-R / h_R) * np.exp(-np.abs(z) / h_z)

            if np.random.random() <= v:
                Rs.append(R)
                zs.append(z)

                good = True

    thetas = np.random.random(len(Rs)) * (2 * np.pi)
    xs, ys = Rs * np.cos(thetas), Rs * np.sin(thetas)

    if outname is not None:
        np.save(outname, np.array([xs, ys, zs]))

    return xs, ys, zs
