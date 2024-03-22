
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt

import gala.dynamics as gd
from . import utils


__all__ = ["ParticleGrid"]


class ParticleGrid:
    """ Class to initialize a grid of particles in the xy plane with circular velocities based on a given gravitational potential.

        Usage:
            ```
            mass_profile = gn.gen_mass_profile(potential)
            grid = gn.ParticleGrid()
            grid.generate(mass_profile)
            ```
    """

    def __init__(self, Rmax=10, n_particles=50, z_start=0., veldisp=10.):
        """ Generate a grid of particles in the xy plane with circular velocities based on a given gravitational potential.
            Can also include a Gaussian velocity dispersion (helpful for including a scale height in the disk).

        Args:
            mass_profile (interp1d): The mass profile of the potential in scipy interp1d format. Defaults to None.
            Rmax (_type_): The maximum radius of the grid particles.
            n_particles (int, optional): Number of particles between -Rmax and Rmax. Defaults to 50.
            z_start (float, optional): The z-level to initialize the grid in. 0 is . Defaults to 0.
            veldisp (float, optional): _description_. Defaults to 10.
            velocities (bool, optional): If false, particles have no velocities. Defaults to True.
        """
        self.container = []             # List of PhaseSpacePosition objects
        self.Rmax = Rmax
        self.n_particles = n_particles
        self.z_start = z_start
        self.veldisp = veldisp

    
    def generate(self, mass_profile, velocities=True):
        particle_range = np.linspace(-self.Rmax, self.Rmax, self.n_particles) * u.kpc

        for x in particle_range:
            for y in particle_range:
                p_x0 =  u.Quantity([x, y, self.z_start * u.kpc])
                R = np.sqrt(x ** 2 + y **2)
                if R > self.Rmax * u.kpc:
                    continue
                
                if velocities:
                    theta = np.arctan2(y, x)
                    
                    p_v0 = utils.velocity(mass_profile, R, theta)          # Initialize with circular velocity
                    if self.veldisp is not None:
                        p_v0 += np.random.normal(scale=self.veldisp, size=3) * (u.km / u.s)    # Add random velocity dispersion
                else:
                    p_v0 = [0., 0., 0.] * (u.km / u.s)
                
                w0 = gd.PhaseSpacePosition(pos=p_x0, vel=p_v0.to(u.kpc / u.Myr))
                self.container.append(w0)
        

    def plot_phase_space(self, ax, color="black", outname=None):
        fig, ax = plt.subplots(1, 2, facecolor="white", figsize=(8, 4))
        for particle in self.container:
            ax[0].scatter(particle.x, particle.y, color=color)
            ax[1].scatter(particle.v_x, particle.v_y, color=color)
        
        plt.tight_layout()
        if outname is not None:
            plt.savefig(outname, dpi=200)
        else:
            plt.show()
