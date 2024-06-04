import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt

from scipy.interpolate import interp1d

from . import utils

from .postprocessing import analysis

__all__ = ["DynamicShadow",  "UniformShadow", "ExponentialShadow", "EdgeOnShadow", 
           "UniformLinearZVariableShadow", "UniformExponentialZVariableShadow"]


def _shadow_tangent(z, phi):
    return z / np.tan(phi)


class ShadowBase:
    """
    Base class for creating shadow objects.

    Parameters:
    - damping (float): The damping factor for the shadow.
    - R_disk (Quantity or float): The radius of the shadow disk.
    - shadow_name (str): The name of the shadow.

    Methods:
    - evaluate(xyz): Evaluates the shadow at the given coordinates.
    - init_from_wind(wind): Initializes the shadow from a wind object.
    - plot_shadow(ax, wind, color, outname, x0, y0, z0): Plots the shadow on the given axes.
    - __repr__(): Returns a string representation of the shadow object.
    """

    def __init__(self, damping=0.5, R_disk=10, shadow_name="ShadowBase", **kwargs):
        self.damping = damping
        if isinstance(R_disk, u.Quantity):
            R_disk = R_disk.to(u.kpc).value
        self.R_disk = R_disk
        self.shadow_name = shadow_name

        self.dynamic_shadow = kwargs.get("dynamic", False)

    def evaluate(self, q, p, t):
        """
        Evaluates the shadow at the given coordinates.

        Parameters:
        - xyz (array-like): The coordinates at which to evaluate the shadow.

        Returns:
        - result (array-like): The evaluated shadow values.
        """
        raise NotImplementedError

    def init_from_wind(self, wind):
        """
        Initializes an angled shadow from a wind object.

        Parameters:
        - wind (Wind): The wind object from which to initialize the shadow.
        """
        self.phi = wind.inclination

    def plot_shadow(
        self, ax=None, wind=None, color="black", outname=None, x0=0, y0=0, z0=0
    ):
        """
        Plots the shadow on the given axes.

        Parameters:
        - ax (array-like, optional): The axes on which to plot the shadow. If not provided, a new figure will be created.
        - wind (Wind, optional): The wind object to plot along with the shadow.
        - color (str, optional): The color of the shadow.
        - outname (str, optional): The name of the output file to save the plot.
        - x0, y0, z0 (float, optional): The coordinates of the origin of the shadow.

        Returns:
        - None
        """
        R_disk = self.R_disk.value if type(self.R_disk) is u.Quantity else self.R_disk
        R_plot = R_disk * 1.5

        # Create figure if ax is not provided
        if ax is None:
            fig, ax = plt.subplots(1, 3, facecolor="white", figsize=(12, 4))

        for axis in ax.flatten():
            axis.set_xlim(-R_plot, R_plot)
            axis.set_ylim(-R_plot, R_plot)

        x = np.linspace(-R_disk, R_disk, 100)
        y = np.linspace(-R_disk, R_disk, 100)
        X, Y = np.meshgrid(x, y)
        XY = self.evaluate(np.array([X, Y, np.zeros(X.shape) + z0]).T, 0)
        XZ = self.evaluate(np.array([X, np.zeros(X.shape) + y0, Y]).T, 0)
        YZ = self.evaluate(np.array([np.zeros(X.shape) + x0, X, Y]).T, 0)

        im1 = ax[0].imshow(
            1 - XY,
            origin="lower",
            cmap="Greys",
            alpha=0.5,
            extent=(-R_disk, R_disk, -R_disk, R_disk),
        )
        im2 = ax[1].imshow(
            1 - XZ,
            origin="lower",
            cmap="Greys",
            alpha=0.5,
            extent=(-R_disk, R_disk, -R_disk, R_disk),
        )
        im3 = ax[2].imshow(
            1 - YZ,
            origin="lower",
            cmap="Greys",
            alpha=0.5,
            extent=(-R_disk, R_disk, -R_disk, R_disk),
        )

        plt.colorbar(mappable=im1, ax=ax[0], location="top")
        plt.colorbar(mappable=im2, ax=ax[1], location="top")
        plt.colorbar(mappable=im3, ax=ax[2], location="top")

        if wind is not None:
            utils.plot_wind_vector(
                wind.normalized().value, ax, length=0.5, loc=(-R_disk, -R_disk, -R_disk)
            )

        if ax is None:
            plt.tight_layout()
            if outname is not None:
                plt.savefig(outname, dpi=200)
            else:
                plt.show()
    
    def plot_shadow_xz(self, ax=None, **kwargs):
        """ Plot the shadow in the x-z plane for easy debugging"""
        xrange = kwargs.get('xrange', (-20, 40))
        zrange = kwargs.get('zrange', (-10, 70))

        t = kwargs.get('t', 0)

        zs, xs = np.mgrid[zrange[0]:zrange[1], xrange[0]:xrange[1]]
        ys = np.zeros_like(xs)

        evaluation = self.evaluate(np.array([xs,ys,zs]).T, t)
        
        plt.figure(figsize=kwargs.get('figsize', (6, 6)))
        plt.imshow(evaluation, origin="lower", cmap=kwargs.get("cmap", "Greys_r"), vmin=0, vmax=1, 
                extent=[xrange[0], xrange[1], zrange[0], zrange[1]])
        
        plt.xlabel("X [kpc]")
        plt.ylabel("Z [kpc]")
        plt.colorbar(label="Shadowing    [0 = Complete Reduction]")

        plt.tight_layout()
        
        outname = kwargs.get('outname', None)
        if outname is not None:
            plt.savefig(outname)
        else:
            plt.show()

    def __repr__(self):
        """
        Returns a string representation of the shadow object.

        Returns:
        - repr (str): The string representation of the shadow object.
        """
        return f"<{self.shadow_name} RP Shadow:  Phi={np.rad2deg(self.phi):.2f}  Damping={self.damping}  R_Disk={self.R_disk}  zmin={self.zmin} >"


class UniformShadow(ShadowBase):
    """A class representing a uniform angled shadow.

    This class inherits from the ShadowBase class.

    Args:
        damping (float, optional): The damping factor. Defaults to 0.5.
        R_disk (Quantity, optional): The radius of the disk. Defaults to 10 kpc.
        zmin (Quantity, optional): The minimum value of z above the disk. Defaults to 0 kpc.
        phi (float, optional): The angle in radians. Defaults to 20 degrees.
    """

    def __init__(self, damping=0.5, R_disk=10, zmin=0.5, phi=np.deg2rad(20), **kwargs ):
        super().__init__(damping=damping, R_disk=R_disk, shadow_name="Uniform", **kwargs)
        if isinstance(zmin, u.Quantity):
            zmin = zmin.to(u.kpc).value
        self.zmin = zmin
        self.phi = phi

        self.frac = kwargs.get("frac", 0.9)
        self.Rmax = kwargs.get("Rmax", 20)
        self.zmax = kwargs.get("zmax", 2)
        self.debug = kwargs.get("debug", False)

        self.Rdisks = []

    def evaluate(self, q, p, t):
        x, y, z = q.T

        if self.dynamic_shadow:
            self.R_disk = analysis.calculate_rstrip(q.T, frac=self.frac, rmax=self.Rmax, zmax=self.zmax)
        if self.debug:
            self.Rdisks.append(self.R_disk)

        cent = _shadow_tangent(z, self.phi)
        dist = np.sqrt((x - cent) ** 2 + y**2)

        out = np.ones(dist.shape)
        in_disk = np.logical_and((z > self.zmin), (dist < self.R_disk))
        out[in_disk] = self.damping
        return out

    def plot_shadow(
        self, ax=None, wind=None, color="black", outname=None, x0=0, y0=0, z0=None
    ):
        z0 = self.zmin + 0.5 if z0 is None else z0
        super().plot_shadow(
            ax=ax, wind=wind, color=color, outname=outname, x0=z0, y0=y0, z0=z0
        )


class ExponentialShadow(ShadowBase):
    """A class representing an exponential drop-off angled shadow.

    This class inherits from the ShadowBase class.

    Args:
        damping (float, optional): The damping factor. Defaults to 0.5.
        R_disk (Quantity, optional): The radius of the disk. Defaults to 10 kpc.
        zmin (Quantity, optional): The minimum value of z above the disk. Defaults to 0 kpc.
        phi (float, optional): The angle in radians. Defaults to 20 degrees.
    """

    def __init__(
        self, damping=0.5, R_disk=10 * u.kpc, zmin=0 * u.kpc, phi=np.deg2rad(20)
    ):
        super().__init__(damping=damping, R_disk=R_disk, shadow_name="Exponential")

        if isinstance(zmin, u.Quantity):
            zmin = zmin.to(u.kpc).value
        self.zmin = zmin
        self.phi = phi

    def evaluate(self, q, p, t):
        x, y, z = q.T
        cent = _shadow_tangent(z, self.phi)
        dist = np.sqrt((x - cent) ** 2 + y**2)
        out = np.exp(-dist / self.R_disk)
        in_disk = np.logical_and((z > self.zmin), (dist < self.R_disk))
        out[in_disk] *= self.damping
        return out

    def plot_shadow(
        self, ax=None, wind=None, color="black", outname=None, x0=0, y0=0, z0=None
    ):
        z0 = self.zmin + 0.5 if z0 is None else z0
        super().plot_shadow(
            ax=ax, wind=wind, color=color, outname=outname, x0=z0, y0=y0, z0=z0
        )


class EdgeOnShadow(ShadowBase):
    def __init__(self, damping=0.5, R_disk=10 * u.kpc, Z_disk=2 * u.kpc, x0=0 * u.kpc, **kwargs):
        super().__init__(damping=damping, R_disk=R_disk, shadow_name="EdgeOn", **kwargs)

        # Always assume kpc for xyz inputs (galactic coordinate system)
        self.x0 = x0.to(u.kpc).value
        self.R_disk = R_disk.to(u.kpc).value
        self.Z_disk = Z_disk.to(u.kpc).value

        self.frac = kwargs.get("frac", 0.9)
        self.Rmax = kwargs.get("Rmax", 20)
        self.zmax = kwargs.get("zmax", 2)

    def evaluate(self, q, p, t):
        x, y, z = q.T

        if self.dynamic_shadow:
            self.R_disk = analysis.calculate_rstrip(q.T, frac=self.frac, rmax=self.Rmax, zmax=self.zmax)

        in_ellipsoid = (y / self.R_disk) ** 2 + (z / self.Z_disk) ** 2 < 1

        out = np.ones(x.shape)

        out[np.logical_and((x >= 0), in_ellipsoid)] = self.damping

        return out

    def plot_shadow(
        self, ax=None, wind=None, color="black", outname=None, x0=0, y0=0, z0=0
    ):
        super().plot_shadow(
            ax=ax, wind=wind, color=color, outname=outname, x0=z0, y0=y0, z0=z0
        )


class UniformLinearZVariableShadow(ShadowBase):
    def __init__(self, damping=0.5, R_disk=10, zmin=0.5, phi=np.deg2rad(20), z_dropoff=10, **kwargs ):
        super().__init__(damping=damping, R_disk=R_disk, shadow_name="Uniform", **kwargs)
        if isinstance(zmin, u.Quantity):
            zmin = zmin.to(u.kpc).value
        self.zmin = zmin
        self.phi = phi
        self.z_dropoff = z_dropoff

        self.frac = kwargs.get("frac", 0.9)
        self.Rmax = kwargs.get("Rmax", 20)
        self.zmax = kwargs.get("zmax", 2)
        self.debug = kwargs.get("debug", False)

        self.Rdisks = []

    def evaluate(self, q, p, t):
        x, y, z = q.T

        if self.dynamic_shadow:
            self.R_disk = analysis.calculate_rstrip(q.T, frac=self.frac, rmax=self.Rmax, zmax=self.zmax)
        if self.debug:
            self.Rdisks.append(self.R_disk)

        cent = _shadow_tangent(z, self.phi)
        dist = np.sqrt((x - cent) ** 2 + y**2)

        out = np.ones(dist.shape)
        in_disk = np.logical_and((z > self.zmin), (dist < self.R_disk))
        out[in_disk] = self.damping + z[in_disk] / self.z_dropoff

        out[out > 1] = 1

        return out


class UniformExponentialZVariableShadow(ShadowBase):
    def __init__(self, damping=0.5, R_disk=10, zmin=0.5, phi=np.deg2rad(20), z_dropoff=10, **kwargs ):
        super().__init__(damping=damping, R_disk=R_disk, shadow_name="Uniform", **kwargs)
        if isinstance(zmin, u.Quantity):
            zmin = zmin.to(u.kpc).value
        self.zmin = zmin
        self.phi = phi
        self.z_dropoff = z_dropoff

        self.frac = kwargs.get("frac", 0.9)
        self.Rmax = kwargs.get("Rmax", 20)
        self.zmax = kwargs.get("zmax", 2)
        self.debug = kwargs.get("debug", False)

        self.dynamic_ratio = kwargs.get("dynamic_ratio", 1.5)

        self.Rdisks = []

    def evaluate(self, q, p, t):
        x, y, z = q.T

        if self.dynamic_shadow:
            self.R_disk = analysis.calculate_rstrip(q.T, frac=self.frac, rmax=self.Rmax, zmax=self.zmax)
        if self.debug:
            self.Rdisks.append(self.R_disk)

        cent = _shadow_tangent(z, self.phi)
        dist = np.sqrt((x - cent) ** 2 + y**2)

        out = np.ones(dist.shape)
        in_disk = np.logical_and((z > self.zmin), (dist < self.R_disk))
        out[in_disk] = self.damping +  (1 - np.exp(-z[in_disk] / self.z_dropoff))
        
        out[out > 1] = 1

        return out


class DynamicShadow:
    
    def __init__(self, wind, a_disk=20, b_disk=2,
                  rho=grp.Density(1e-26 * u.g / u.cm ** 3), 
                 tau_wind=1 * u.Myr, masses=None, radii=None, **kwargs):
        self.wind = wind
        self.rho = rho

        self.masses = masses
        self.radii = radii * u.pc

        self.shadow_name = "Dynamic"

        self.a_disk = a_disk
        self.y_range = kwargs.get("y_range", (-self.a_disk, self.a_disk))
        
        self.b_disk = observed_minor_axis(a_disk, b_disk, wind.inclination)
        self.z_range = kwargs.get("z_range", (-self.b_disk, self.b_disk))
        
        self.n_bins = kwargs.get("n_bins", 20)

        self.L_bin = (self.y_range[1] - self.y_range[0]) / self.n_bins * u.kpc
        self.W_bin = (self.z_range[1] - self.z_range[0]) / self.n_bins * u.kpc

        self.A_bin = self.L_bin * self.W_bin
        self.tau_wind = tau_wind

        self.x_range = kwargs.get("x_range", (-20, 50))
        self.n_bins_wind_direction = kwargs.get("n_bins_wind_direction", 101)

        self.debug = kwargs.get("debug", False)
        self.interps = []


    def evaluate(self, q, p, t):
        """
        Evaluate the shadowing effect on particles dynamically, by rotating particles to the wind-frame of reference
        and then computing line-of sight depths.

        NOTE: In the wind frame, the Y-Z plane is what the wind "sees", and the X-axis is the distance along the
        wind direction.

        Parameters:
        - q: numpy array, shape (3, N)
            Array of particle positions.
        - t: float
            Time parameter (Redundant for this shadow, but might be useful if we implement a changing wind inclination).

        Returns:
        - shadowing: numpy array, shape (N,)
            Array of shadowing values for each particle.
        """

        q, p = q.T, p.T
    
        # STEP 1: Rotate the particles to the wind frame    

        xyz_rotated = grp.rotate_yaxis(q, beta=np.pi/2 - self.wind.inclination)
        xyz_rotated_T = xyz_rotated.T
        v_xyz_rotated = grp.rotate_yaxis(p, beta=np.pi/2 - self.wind.inclination)

        v_wind = np.sqrt(np.sum(wind.evaluate(t) ** 2)) * u.kpc / u.Myr

        shadowing = np.ones(len(q[0]))

        # STEP 2: Bin the particles along the y-z axis to get distribution along wind's line of sight
        # Create bins along the y and z axis
        y_bins = np.linspace(self.y_range[0], self.y_range[1], self.n_bins + 1)
        z_bins = np.linspace(self.z_range[0], self.z_range[1], self.n_bins + 1)

        L_slab = self.tau_wind * v_wind        
        A_clouds = np.pi * self.radii ** 2

        # This is the momentum deposited by the slab into each cloud
        # (p_dep / p_slab) = 2(1 + R_cl/L_slab)(1 - v_x/v_wind)(A_cloud/A_slab)
        normalized_p = 2 * (1 + self.radii / L_slab) 
        normalized_p *= (A_clouds / self.A_bin).cgs * (1 - v_xyz_rotated[0] * u.km / u.s / v_wind)
        

        # Bin particles along y-z axis
        y_bin_indices = np.digitize(xyz_rotated[1], y_bins)
        z_bin_indices = np.digitize(xyz_rotated[2], z_bins)

        # Get particles in the bin, and take only the z component
        # We take a cumulative sum of the histogram to get the integrated number of particles at each distance
        for i in range(self.n_bins):
            for j in range(self.n_bins):
                
                # Determine which particles are in this current bin              
                in_bin = np.bitwise_and(y_bin_indices == i + 1, z_bin_indices == j + 1)

                particles_in_z_bin = xyz_rotated_T[in_bin].T[0]     # Take x component of the rotated particles in bin
                sorted_indices = np.argsort(particles_in_z_bin)     # Sort the particles by x component

                if len(sorted_indices) == 0:
                    continue
                
                cum_p = np.cumsum(normalized_p[sorted_indices])
                cum_p = np.clip(cum_p, 0, 1)

                interp = interp1d(particles_in_z_bin[sorted_indices], 1 - cum_p, kind="linear", 
                                  fill_value=1, bounds_error=False, assume_sorted=True)

                # Get a cumulative histogram of the particle x directions and make it into an interp object
                # NOTE: This might be improved in the future to speed things up, but good for now

                # Apply shadowing to the particles in this bin
                shadowing[in_bin] = interp(xyz_rotated_T[in_bin].T[0])
                
                if self.debug:
                    self.interps.append({"i": i, "j": j, "interp": interp})
                    # fig, ax = plt.subplots(1, 2)
                    # ax[0].hist(normalized_p, bins=30)
                    # ax[1].hist(shadowing)
                    # plt.show()

        return shadowing


    def debug_evaluate(self, q, p, t, **kwargs):
        shadowing = self.evaluate(q, p, t)
        xyz_rotated = grp.rotate_yaxis(q.T, beta=np.pi/2 - self.wind.inclination)
        v_xyz_rotated = grp.rotate_yaxis(p.T, beta=np.pi/2 - self.wind.inclination)


        size = kwargs.get("size", 1)

        outname = kwargs.get("outname", None)
        cmap = kwargs.get("cmap", "magma")

        fig, ax = plt.subplots(1, 3, figsize=(12, 5))

        m1 = ax[0].scatter(xyz_rotated[0], xyz_rotated[1], c=v_xyz_rotated[0], cmap="rainbow", s=size)
        ax[0].set(xlabel="X [kpc]", ylabel="Y [kpc]", xlim=(-20, 20), ylim=(-20, 20))

        m2 = ax[1].scatter(xyz_rotated[0], xyz_rotated[2], c=shadowing, cmap=cmap, s=size)
        ax[1].set(xlabel="X [kpc]", ylabel="Z [kpc]", xlim=(-20, 20), ylim=(-20, 20))

        m3 = ax[2].scatter(xyz_rotated[1], xyz_rotated[2], c=shadowing, cmap=cmap, s=size)
        ax[2].set(xlabel="Y [kpc]", ylabel="Z [kpc]", xlim=(-20, 20), ylim=(-20, 20))

        plot_grid(ax[2], 
                  np.linspace(self.y_range[0], self.y_range[1], self.n_bins + 1), 
                  np.linspace(self.z_range[0], self.z_range[1], self.n_bins + 1))

        plt.colorbar(mappable=m1, ax=ax[0], label="X-Velocity [km s$^{-1}$]", orientation="horizontal", location="top")
        plt.colorbar(mappable=m2, ax=ax[1], label="Wind Strength", orientation="horizontal", location="top")
        plt.colorbar(mappable=m3, ax=ax[2], label="Wind Strength", orientation="horizontal", location="top")

        plt.tight_layout()
        if outname is not None:
            plt.savefig(outname, dpi=kwargs.get("dpi", 200))
            plt.close()