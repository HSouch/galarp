import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

from astropy import units as u

from .. import analysis 
from .. import utils 


def animated_hexbin_plot(orbits, x_ind=1, y_ind=2, n_frames=100, **kwargs):
    """
    Create an animated hexbin plot of the given orbits.

    Parameters:
    - orbits: The orbits data.
    - x_ind: The index of the x-coordinate in the orbit data. Default is 1.
    - y_ind: The index of the y-coordinate in the orbit data. Default is 2.
    - n_frames: The number of frames in the animation. Default is 100.
    - **kwargs: Additional keyword arguments for customizing the plot.

    Returns:
    - None
    """
    
    outname = kwargs.get("outname", "animated_hexbin.gif")
    close = kwargs.get("close", True)
    cmap = kwargs.get("cmap", "viridis")

    figsize = kwargs.get("figsize", (10, 10))

    gridsize = kwargs.get("gridsize", 30)
    xextent = kwargs.get("xextent", (-40., 40.))
    xextent = (-xextent, xextent) if isinstance(xextent, (int, float)) else xextent
    yextent = kwargs.get("yextent", (-40., 40.))
    yextent = (-yextent, yextent) if isinstance(yextent, (int, float)) else yextent
    vextent = kwargs.get("vextent", (-100, 300))
    vextent = (-vextent, vextent) if isinstance(vextent, (int, float)) else vextent

    vmin, vmax = kwargs.get("vmin", 1), kwargs.get("vmax", 500)

    def add_labels(ax):
        ax[1][0].set_xlabel("X [kpc]")
        ax[0][0].set_ylabel("Y [kpc]")
        ax[1][0].set_ylabel("Z [kpc]")
        ax[0][1].set_xlabel("Z [kpc]")
        ax[0][1].set_ylabel(r"$V_z$ [km/ s$^{-1}$]")
        ax[1][1].set_xlabel("Time [Myr]")
        ax[1][1].set_ylabel(r"RP Profile [g/ (cm s$^2$)]")

    
    fig, ax = plt.subplots(2, 2, facecolor="white", figsize=figsize)

    frames = np.linspace(0, len(orbits.data.t) - 1, n_frames).astype(int)

    x,y,z,vx,vy,vz = orbits.get_orbit_data(transposed=False)   

    wind_profile = np.sqrt(np.sum(orbits.metadata["WIND"].evaluate_arr(orbits.data.t) ** 2, axis=1)) * u.kpc/u.Myr
    wind_profile = wind_profile.to(u.cm/u.s)

    dens_profile = orbits.metadata["RHO_ICM"].evaluate_arr(orbits.data.t) * (u.g / u.cm**3)
    rp_profile = wind_profile ** 2 * dens_profile
    
    
    ax[0][0].hexbin(x[0], y[0], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                    extent=[xextent[0], xextent[1], xextent[0], xextent[1]], 
                    vmin=vmin, vmax=vmax, zorder = 5)
    ax[1][0].hexbin(x[0], z[0], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                    extent=[xextent[0], xextent[1], yextent[0], yextent[1]], 
                    vmin=vmin, vmax=vmax, zorder = 5)
    ax[0][1].hexbin(z[0], vz[0], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                    extent=[vextent[0], vextent[1], -400, 400], 
                    vmin=vmin, vmax=vmax, zorder = 5)
    
    ax[1][1].plot(orbits.data.t, rp_profile, color="black", lw=2)
    
    add_labels(ax)

    plt.tight_layout()

    
    def animate(i):
        this_x, this_y, this_z = x[i], y[i], z[i]
        rstrip, strip_med = analysis.rstrip_and_medians([this_x, this_y, this_z], frac=0.9)


        ell = utils.ellipse_coords(strip_med[0], strip_med[1], rstrip, rstrip, 0, 100)

        for axis in ax.flatten()[:]:
            axis.cla()

        ax[0][0].hexbin(x[i], y[i], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                        extent=[xextent[0], xextent[1], xextent[0], xextent[1]], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        
        
        ax[1][0].hexbin(x[i], z[i], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                        extent=[xextent[0], xextent[1], yextent[0], yextent[1]], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        ax[0][1].hexbin(z[i], vz[i], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                        extent=[vextent[0], vextent[1], -400, 400], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        ylims = ax[0][0].get_ylim()

        ax[0][0].plot(ell[0], ell[1], color="black", zorder=6)
        ax[1][0].plot(ell[0], ell[2], color="black", zorder=6)

        ax[1][1].plot(orbits.data.t, rp_profile, color="black", lw=2)

        current_time = orbits.data.t[i].value
        ylims = ax[1][1].get_ylim()
        ax[1][1].plot([current_time, current_time], [ylims[0], ylims[1]], color="red", lw=1, ls="dashed")

        add_labels(ax)

    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=100)
    ani.save(outname, writer='pillow', fps=24)

    if close:
        plt.close(fig)


def r_vr(orbits, **kwargs):
    n_frames = kwargs.get("n_frames", 200)
    outname = kwargs.get("outname", "r_vr.gif")

    x,y,z, vx, vy, vz = orbits.get_orbit_data(transposed=False)

    r = np.sqrt(x**2 + y**2 + z**2)
    vr = (x*vx + y*vy + z*vz) / r

    gridsize = 30
    xextent = (-10, 100)
    yextent = (-50, 800)
    vmin, vmax = 1, 200

    fig, ax = plt.subplots(1, 1)

    ax.hexbin(r[0], vr[0], bins="log", cmap=kwargs.get("cmap", "viridis"), gridsize=(gridsize, gridsize),
              extent=[xextent[0], xextent[1], yextent[0], yextent[1]], 
              vmin=vmin, vmax=vmax, zorder = 5)
    rs = np.linspace(*xextent, 500)

    pot = orbits.metadata["POTENTIAL"]
    v_escs = np.sqrt(2 * np.abs(pot.energy(np.array([np.zeros(len(rs)), np.zeros(len(rs)), rs])))).to(u.km/u.s)
    
    ax.plot(rs, v_escs,
              color="k", lw=1, zorder=10)
    
    def add_labels(ax):
        ax.set_xlabel("R [kpc]")
        ax.set_ylabel(r"$v_R$ [km s$^{-1}$]")

    add_labels(ax)

    def animate(i):
        this_r, this_vr = r[i], vr[i]

        stripped = this_vr > np.sqrt(2 * np.abs(pot.energy(np.array([x[i], y[i], z[i]])))).to(u.km/u.s).value
        ax.cla()

        ax.plot(rs, v_escs, color="k", lw=1, zorder=10)
        
        ax.hexbin(this_r, this_vr, bins="log", cmap=kwargs.get("cmap", "viridis"), gridsize=(gridsize, gridsize),
                  extent=[xextent[0], xextent[1], yextent[0], yextent[1]], 
                  vmin=vmin, vmax=vmax, zorder = 5)

        xlims, ylims = ax.get_xlim(), ax.get_ylim()
        dx, dy = xlims[1] - xlims[0], ylims[1] - ylims[0]

        ax.text(xlims[1] - 0.2 * dx, ylims[1] - 0.05 * dy, f"t = {orbits.data.t[i]:.1f}", fontsize=8)
        ax.text(xlims[1] - 0.2 * dx, ylims[1] - 0.1 * dy, f"Stripped: {stripped.sum() / len(stripped):.2f}", fontsize=8)

        add_labels(ax)
    
    ani = animation.FuncAnimation(fig, animate, frames=n_frames, interval=100)
    ani.save(outname, writer='pillow', fps=24)