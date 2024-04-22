import numpy as np

from astropy import units as u

from matplotlib import pyplot as plt

from scipy import stats

from .. import analysis

from ...utils import get_orbit_data
from ...rampressure import OrbitContainer


def rstrip_plot(orbits, zmax=2 * u.kpc, rstrip_frac = 0.8, rmax=20 * u.kpc, 
                close_plot=False, title=None, **kwargs):
    """
    Plot the stripping radius for an interaction as a function of time.

    Parameters:
    - orbits: str or OrbitContainer
        The orbits to be analyzed. If a string is provided, it is assumed to be the path to a file containing the orbits.
    - zmax: Quantity, optional
        The maximum absolute value of the z-coordinate to consider. Default is 2 kpc.
    - rstrip_frac: float, optional
        The fraction of particles to consider when calculating the r-strip radius. Default is 0.8.
    - rmax: Quantity, optional
        The maximum radius to consider. Default is 20 kpc.
    - close_plot: bool, optional
        Whether to close the plot after saving or showing. Default is False.
    - title: str, optional
        The title of the plot. Default is None.
    - **kwargs: keyword arguments
        Additional keyword arguments to be passed to the function.

    Returns:
    None
    """
    
    if isinstance(orbits, str):
        orbits = OrbitContainer.load(orbits)

    outname = kwargs.get("outname", None)
    
    x,y,z, *_ = get_orbit_data(orbits.data)
    x,y,z = x.T, y.T, z.T
    
    r = np.sqrt(x**2 + y**2 + z**2)
    times = orbits.data.t

    
    fig, ax = plt.subplots(1, 2, figsize = (10, 5), facecolor="white")
    
    for i in np.linspace(0, len(r) - 1, 6).astype(int):
        this_r, this_z = r[i], z[i]
        this_r_cut = this_r[np.abs(this_z) < zmax.value]
        this_r_cut = this_r_cut[this_r_cut < rmax.value]
        
        cdf = ax[0].ecdf(this_r_cut, lw=2)
        xdata, ydata = cdf.get_data()
        
        closest = np.argmin(np.abs(ydata - rstrip_frac))
        rstrip = xdata[closest]
        
        ax[0].axvline(xdata[closest], color=cdf.get_color(), alpha=.3)
        
        cdf.set(label = f't = {times[i]} Myr, ' + r'$R_{strip} = $' + f'{rstrip:.2f} kpc',)
        
        ax[1].scatter(times[i], rstrip, color=cdf.get_color(), zorder=3, marker="s", s=50, alpha=0.8)
    
    strip_times, rstrips = [], []
    for i in range(0, len(r), 5):
        this_r, this_z = r[i], z[i]
        this_r_cut = this_r[np.abs(this_z) < zmax.value]
        this_r_cut = this_r_cut[this_r_cut < rmax.value]
        
        cdf = stats.ecdf(this_r_cut)
        cdf_xs, cdf_vals = cdf.cdf.quantiles, cdf.cdf.probabilities
        
        closest = np.argmin(np.abs(cdf_vals - rstrip_frac)) 
        strip_times.append(times[i].value)
        rstrips.append(cdf_xs[closest])
        
    ax[1].scatter(strip_times, rstrips, color="black", s=10, marker="s")

    ax[0].axhline(rstrip_frac, color="black", alpha=0.3)
    ax[0].set_ylim(0, 1.05)
    ax[0].set_ylabel(r'$F(< R)$')
    ax[1].set_ylabel(f'R({rstrip_frac}' + r'$ N_{tot} < R)$')

    ax[0].legend(loc="lower right", fontsize=8)

    ax[0].set_xlabel(r'$R$    [kpc]')
    ax[1].set_xlabel(r'$t$    [Myr]')

    xmin, xmax = ax[1].get_xlim()
    ymin, ymax = ax[1].get_ylim()
    dx, dy = xmax - xmin, ymax - ymin
    
    ax[1].text(xmax - dx/2, ymax - 1 * dy/10, r'$f_{strip}$ = ' + f'{rstrip_frac}', fontsize=15)
    ax[1].text(xmax - dx/2, ymax - 1.6 * dy/10, r'$z_{max}$ = ' + f'{zmax}', fontsize=15)
    
    if title is not None:
        plt.suptitle(title)
        
    plt.tight_layout()
    
    if outname is not None:
        plt.savefig(outname)
    else:
        plt.show()
    if close_plot:
        plt.close()


def stripped_plot(orbits, **kwargs):
    stripped_frac = analysis.stripped(orbits)
    outname = kwargs.get("outname", None)

    plt.figure(figsize=kwargs.get("figsize", (6, 5)))

    plt.plot(orbits.data.t, stripped_frac)
    plt.ylim(-0.1, 1.1)
    plt.xlabel("Time (Myr)")
    plt.ylabel("Fraction of stripped particles")

    plt.tight_layout()
    if outname is not None:
        plt.savefig(outname, dpi = kwargs.get("dpi", 200))
    else:
        plt.show()
