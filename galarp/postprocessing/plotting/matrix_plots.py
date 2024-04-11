
from ..analysis import calculate_rstrip

from ..utils import ellipse_coords
from ...utils import get_orbit_data


import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


class MatrixPlot:
    """ Base class to create matrix plots. """
    def __init__(self, **kwargs):
        self.cmap = kwargs.get("cmap", "viridis")
        self.figsize = kwargs.get("figsize", (10, 10))
        self.nrows = kwargs.get("nrows", 4)
        self.ncols = kwargs.get("ncols", 4)

        self.sharex = kwargs.get("sharex", True)
        self.sharey = kwargs.get("sharey", True)

    def create(self):
        self.fig, self.ax = plt.subplots(nrows=self.nrows, ncols=self.ncols, figsize=self.figsize, 
                               sharex=self.sharex, sharey=self.sharey)
        return self.fig, self.ax
    
    def populate(self, method, *args):
        index = 0
        for i in range(self.nrows):
            for j in range(self.ncols):
                method(i, j, index, *args)
                index += 1
    
    def add_xy_labels(self, xlabel, ylabel):
        for i in range(self.nrows):
            self.ax[i, 0].set_ylabel(ylabel)
        for j in range(self.ncols):
            self.ax[self.nrows - 1, j].set_xlabel(xlabel)
    
    def add_colorbar(self,i,j, mappable, loc="lower left"):
        ia = inset_axes(self.ax[i,j], width="50%", height="5%", loc=loc)
        plt.colorbar(mappable=mappable, cax=ia, orientation="horizontal")
        ia.xaxis.set_ticks_position("top")
    


def density_matrix(orbits, x_ind=0, y_ind=1, nrows=4, ncols=4, outname=None, **kwargs):
    """ Create a matrix plot of the density of the orbits for a given principle axis.

    Args:
        orbits (OrbitContainer): The orbits to plot.
        x_ind (int, optional): The x-axis index. Defaults to 0.
        y_ind (int, optional): The y-axis index. Defaults to 1.
        nrows (int, optional): The number of rows in the matrix. Defaults to 4.
        ncols (int, optional): The number of columns in the matrix. Defaults to 4.
        outname (str, optional): The name of the output file. Defaults to None.
    
    Usage:
        ```
        density_matrix(orbits, x_ind=0, y_ind=1, nrows=4, ncols=4, outname="test_out.pdf", **kwargs)
        ```
    """
    
    times = kwargs.get("times", np.linspace(0, len(orbits.data.t) - 1, nrows * ncols).astype(int))
    xextent = kwargs.get("xextent", (-40., 40.))
    xextent = (-xextent, xextent) if isinstance(xextent, (int, float)) else xextent
    yextent = kwargs.get("yextent", (-40., 40.))
    yextent = (-yextent, yextent) if isinstance(yextent, (int, float)) else yextent

    hlines, vlines = kwargs.get("hlines", []), kwargs.get("vlines", [])

    cbar_loc = kwargs.get("cbar_loc", "lower left")

    vmin, vmax = kwargs.get("vmin", 1), kwargs.get("vmax", 500)

    gridsize = kwargs.get("gridsize", 30)

    xlabel=kwargs.get("xlabel", "X (kpc)")
    ylabel=kwargs.get("ylabel", "Y (kpc)")
    
    mplot = MatrixPlot(**kwargs)
    fig, ax = mplot.create()

    data = get_orbit_data(orbits.data)

    x, y = data[x_ind], data[y_ind]

    mappables = []
    index = 0
    def method(i, j, index, x, y):
        this_x, this_y = x.T[times[index]], y.T[times[index]]

        hb = ax[i, j].hexbin(this_x, this_y, bins="log", cmap=mplot.cmap, gridsize=(gridsize, gridsize),
                        extent=[xextent[0], xextent[1], yextent[0], yextent[1]], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        mappables.append(hb)
        xmin, xmax = ax[i, j].get_xlim()
        ymin, ymax = ax[i, j].get_ylim()
        dx, dy = xmax - xmin, ymax - ymin
        
        ax[i, j].text(xmin + dx/20, ymax - dy/10, f'{orbits.data.t[times[index]]}')
        for hline in hlines:
            ax[i, j].axhline(hline, color="Grey", lw=1, zorder=1)
        for vline in vlines:
            ax[i, j].axvline(vline, color="Grey", lw=1, zorder=1)

    mplot.populate(method, x, y)
    mplot.add_xy_labels(xlabel, ylabel)

    mplot.add_colorbar(i=nrows-1, j=ncols-1, mappable=mappables[-1], loc=cbar_loc)

    plt.subplots_adjust(hspace=0, wspace=0, left=0.07, right=0.99, top=0.99, bottom=0.05)

    if outname is not None:
        plt.savefig(outname)
    else:
        plt.show()


def rstrip_check(orbits, x_ind=0, y_ind=1, nrows=4, ncols=4, outname=None, **kwargs):
    """
    Check the computation of the stripping radius 

    Args:
        orbits: An object containing orbit data.
        x_ind: The index of the x-coordinate in the orbit data. Default is 0.
        y_ind: The index of the y-coordinate in the orbit data. Default is 1.
        nrows: The number of rows in the plot grid. Default is 4.
        ncols: The number of columns in the plot grid. Default is 4.
        outname: The name of the output file to save the plot. Default is None.
    Returns:
    - None

    """
    
    times = kwargs.get("times", np.linspace(0, len(orbits.data.t) - 1, nrows * ncols).astype(int))
    xextent = kwargs.get("xextent", (-40., 40.))
    xextent = (-xextent, xextent) if isinstance(xextent, (int, float)) else xextent
    yextent = kwargs.get("yextent", (-40., 40.))
    yextent = (-yextent, yextent) if isinstance(yextent, (int, float)) else yextent

    hlines, vlines = kwargs.get("hlines", []), kwargs.get("vlines", [])

    cbar_loc = kwargs.get("cbar_loc", "lower left")

    vmin, vmax = kwargs.get("vmin", 1), kwargs.get("vmax", 500)

    gridsize = kwargs.get("gridsize", 30)

    xlabel=kwargs.get("xlabel", "X (kpc)")
    ylabel=kwargs.get("ylabel", "Z (kpc)")
    
    mplot = MatrixPlot(**kwargs)
    fig, ax = mplot.create()

    x,y,z, *_ = get_orbit_data(orbits.data)

    mappables = []
    index = 0
    def method(i, j, index, x, y, z):
        this_x, this_y, this_z = x.T[times[index]], y.T[times[index]], z.T[times[index]]

        rstrip = calculate_rstrip([this_x, this_y, this_z], 
                                  rmax=kwargs.get("rmax", 20),
                                  zmax=kwargs.get("zmax", 2),
                                  frac=kwargs.get("frac", 0.9))
        
        ellipse = ellipse_coords(0, 0, rstrip, rstrip / 5, 0)
        ellipse_shifted = ellipse_coords(np.median(this_x), np.median(this_y), rstrip, rstrip / 5, 0)

        hb = ax[i, j].hexbin(this_x, this_z, bins="log", cmap=mplot.cmap, gridsize=(gridsize, gridsize),
                        extent=[xextent[0], xextent[1], yextent[0], yextent[1]], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        mappables.append(hb)
        xmin, xmax = ax[i, j].get_xlim()
        ymin, ymax = ax[i, j].get_ylim()
        dx, dy = xmax - xmin, ymax - ymin
        
        ax[i, j].text(xmin + dx/20, ymax - dy/10, f'{orbits.data.t[times[index]]}')
        for hline in hlines:
            ax[i, j].axhline(hline, color="Grey", lw=1, zorder=1)
        for vline in vlines:
            ax[i, j].axvline(vline, color="Grey", lw=1, zorder=1)
        
        ax[i, j].plot(ellipse[0], ellipse[1], color="red", lw=1, zorder=10)
        ax[i, j].plot(ellipse_shifted[0], ellipse_shifted[1], color="blue", lw=1, zorder=10)

    mplot.populate(method, x, y, z)
    mplot.add_xy_labels(xlabel, ylabel)

    mplot.add_colorbar(i=nrows-1, j=ncols-1, mappable=mappables[-1], loc=cbar_loc)

    plt.subplots_adjust(hspace=0, wspace=0, left=0.07, right=0.99, top=0.99, bottom=0.05)

    if outname is not None:
        plt.savefig(outname)
    else:
        plt.show()
