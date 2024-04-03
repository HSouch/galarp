
from ...utils import get_orbit_data

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


class MatrixPlot:
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
    


def density_xy_matrix(orbits, x_ind=0, y_ind=1, nrows=4, ncols=4, outname=None, **kwargs):
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
