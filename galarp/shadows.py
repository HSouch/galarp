import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt

from . import utils

from .postprocessing import analysis

__all__ = ["UniformShadow", "ExponentialShadow", "EdgeOnShadow", 
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

    def evaluate(self, xyz, t):
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
        self.phi = wind.inclination()

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
        XY = self.evaluate(np.array([X, Y, np.zeros(X.shape) + z0]).T)
        XZ = self.evaluate(np.array([X, np.zeros(X.shape) + y0, Y]).T)
        YZ = self.evaluate(np.array([np.zeros(X.shape) + x0, X, Y]).T)

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

        utils.plot_disk(ax, R_disk, lw=1, color="black")

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

    def evaluate(self, q, t):
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
        z0 = self.zmin + 0.5 * u.kpc if z0 is None else z0
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

    def evaluate(self, q, t):
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
        z0 = self.zmin + 0.5 * u.kpc if z0 is None else z0
        super().plot_shadow(
            ax=ax, wind=wind, color=color, outname=outname, x0=z0, y0=y0, z0=z0
        )


class EdgeOnShadow(ShadowBase):
    def __init__(self, damping=0.5, R_disk=10 * u.kpc, Z_disk=2 * u.kpc, x0=0 * u.kpc):
        super().__init__(damping=damping, R_disk=R_disk, shadow_name="EdgeOn")

        # Always assume kpc for xyz inputs (galactic coordinate system)
        self.x0 = x0.to(u.kpc).value
        self.R_disk = R_disk.to(u.kpc).value
        self.Z_disk = Z_disk.to(u.kpc).value

    def evaluate(self, q, t):
        x, y, z = q.T

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

    def evaluate(self, q, t):
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

        self.Rdisks = []

    def evaluate(self, q, t):
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
