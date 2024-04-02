import numpy as np
from astropy import units as u
from astropy import constants as const
from gala.units import UnitSystem

from gala import dynamics as gd
from gala import potential as gp

from matplotlib import pyplot as plt



class HostOrbit:
    """ Class for a host orbit to determine a ram pressure profile.
    """
    def __init__(self, potential, init_conditions=None, 
                 density=None, temp=2.51e6 * u.K):
        self.potential = potential

        self.init_conditions = init_conditions
        self.density = density
        self.temp = temp

        self.orbit = None
        self.velocities = None

    

    def integrate(self, n_steps=10000, dt=1):
        assert self.init_conditions is not None

        w0 = gd.PhaseSpacePosition(pos=self.init_conditions.pos, vel=self.init_conditions.vel)

        self.orbit = gp.Hamiltonian(self.potential).integrate_orbit(w0, dt=dt, n_steps=n_steps)
        

    def velocities(self, units = u.km/u.s):
        dx, dy, dz = self.orbit.vel
        v = np.sqrt(dx**2 + dy**2 + dz**2).to(units)
        self.velocities = v

    def gen_ENZO_filename(self, units = UnitSystem(u.cm, u.s, u.radian, u.g)
):
        
        dx, dy, dz = self.orbit.vel.d_xyz


    def plot_infall_orbit(self, outname=None):
        fig, ax = plt.subplots(2, 3, figsize=(13,8))
        x,y,z = self.orbit.pos.xyz
        r = np.sqrt(x**2 + y**2 + z**2)
        
        t = self.orbit.t

        dx, dy, dz = self.orbit.vel.d_xyz

        vel = np.sqrt(dx**2 + dy**2 + dz**2).to(u.km/u.s)
        color, lw = "black", 2

        ax[0][0].plot(x,y, color=color, lw=lw)
        ax[0][0].scatter(0, 0, marker="+", color="black", s=400)
        ax[0][0].set_ylabel('y [kpc]')
        ax[0][0].set_xticks([])

        ax[1][0].plot(x,z, color=color, lw=lw)
        ax[1][0].scatter(0, 0, marker="+", color="black", s=400)
        ax[1][0].set_xlabel('x [kpc]')
        ax[1][0].set_ylabel('z [kpc]')

        ax[0][1].plot(t, r, color=color, lw=lw)
        ax[0][1].set_ylabel('r [kpc]')
        ax[0][1].set_xticks([])

        ax[0][2].plot(t, vel, color=color, lw=lw)
        ax[0][2].set_ylabel('v [km/s]')
        ax[0][2].set_xticks([])

        if self.density is not None:
            dens = self.density.evaluate(r) * const.u.cgs
            ax[1][1].plot(t, dens, color=color, lw=lw)
            ax[1][1].set_ylabel('Density [g/cm^3]')

            rp_profile = (dens * vel**2).to(u.g / u.cm / u.s**2)
            

            ax[1][2].plot(t, rp_profile, color=color, lw=lw)
            ax[1][2].set_ylabel('Ram Pressure [g/cm/s^2]')

        ax[1][1].set_xlabel('t')
        ax[1][1].set_ylabel('Density [1e-3 cm^-3]')

        labels = ["XY", "Radius", "Total Velocity", "XZ", "Density", "RP Profile"]

        for i, axis in enumerate(ax.flatten()):
            xmin, xmax = axis.get_xlim()
            ymin, ymax = axis.get_ylim()
            dx, dy = xmax - xmin, ymax - ymin

            axis.text(xmax  - 0.3 * dx, ymin + 0.9 * dy, 
                    labels[i], fontsize=12, color='black', weight='bold')


        plt.tight_layout()
        if outname is not None:
            plt.savefig(outname)
        else:
            plt.show()