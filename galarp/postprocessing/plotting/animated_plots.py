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

    imin = kwargs.get("imin", 0)
    imax = kwargs.get("imax", len(orbits.data.t) - 1)
    frames = np.linspace(imin, imax, n_frames).astype(int)

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


def dynamic_shadow_animated_plot(orbits, shadow, **kwargs):
    """
    Create an animated plot showing the effects of a given dynamic shadow on a GalaRP orbit container.

    Parameters:
    - orbits: An object representing the orbits.
    - shadow: An object representing the shadow.
    - **kwargs: Additional keyword arguments for customization.

    Keyword Arguments:
    - figsize: Tuple specifying the figure size. Default is (10, 10).
    - outname: String specifying the output filename. Default is "animated_hexbin.gif".
    - close: Boolean indicating whether to close the figure after saving the animation. Default is True.
    - gridsize: Integer specifying the grid size for hexbin plots. Default is 50.
    - xextent: Tuple or float specifying the x-axis extent. Default is (-40., 40.).
    - yextent: Tuple or float specifying the y-axis extent. Default is (-40., 40.).
    - imax: Integer specifying the maximum index of the orbits data. Default is len(orbits.data.t) - 1.
    - vmin: Float specifying the minimum value (of particles) for color mapping. Default is 1.
    - vmax: Float specifying the maximum value (of particles) for color mapping. Default is 500.
    - n_frames: Integer specifying the number of frames for the animation. Default is 100.

    Returns:
    - None
    """
    x,y,z, *_ = orbits.get_orbit_data(orbits.data, transposed=False)

    fig, ax = plt.subplots(2, 2, figsize=kwargs.get("figsize", (10, 10)))

    outname = kwargs.get("outname", "animated_hexbin.gif")
    close = kwargs.get("close", True)

    gridsize = kwargs.get("gridsize", 50)
    xextent = kwargs.get("xextent", (-40., 40.))
    xextent = (-xextent, xextent) if isinstance(xextent, (int, float)) else xextent
    yextent = kwargs.get("yextent", (-40., 40.))
    yextent = (-yextent, yextent) if isinstance(yextent, (int, float)) else yextent

    imax = kwargs.get("imax", len(orbits.data.t) - 1)

    vmin, vmax = kwargs.get("vmin", 1), kwargs.get("vmax", 500)

    n_frames = min(kwargs.get("n_frames", 100), imax)

    frames = np.linspace(0, imax, n_frames).astype(int)


    def add_labels(ax):
        ax[0][0].set(ylabel="Y [kpc]", xlim=xextent, ylim=xextent)
        ax[1][0].set(xlabel="Y [kpc]", ylabel="Z [kpc]", xlim=xextent, ylim=yextent)


    cmaps = ["Greys", "Blues", "Greens", "Reds"]
    colors = ["Grey", "Blue", "Green", "Red"]
    damped = [1.0, 0.9, 0.5, 0.1]

    bins = np.linspace(-0.1, 1.1, 23)

    def animate(i):
        this_x, this_y, this_z = x[i], y[i], z[i]

        this_damping = shadow.evaluate(np.array([this_x, this_y, this_z]).T, 0)
        hist, bin_edges = np.histogram(this_damping, bins=bins)

        for axis in ax.flatten()[:]:
            axis.cla()

        ax[1][1].stairs(hist, bin_edges, color="black", lw=2)

        xlims = ax[1][1].get_xlim()
        ylims = ax[1][1].get_ylim()
        dx, dy = xlims[1] - xlims[0], ylims[1] - ylims[0]
        current_time = orbits.data.t[i].value

        ax[1][1].text(xlims[1] - dx /2, ylims[1] - dy / 10,  f"Time: {current_time:.0f} Myr", ha="left")

        for i in range(len(damped)):
            this_x_damped = this_x[this_damping < damped[i]]
            this_y_damped = this_y[this_damping < damped[i]]
            this_z_damped = this_z[this_damping < damped[i]]

            ax[0][0].hexbin(this_x_damped, this_y_damped, bins="log", cmap=cmaps[i], gridsize=(gridsize, gridsize),
                            extent=[*xextent, *xextent], 
                            vmin=vmin, vmax=vmax, zorder = i + 5)
            
            ax[0][1].hexbin(this_x_damped, this_z_damped, bins="log", cmap=cmaps[i], gridsize=(gridsize, gridsize),
                            extent=[*xextent, *yextent], 
                            vmin=vmin, vmax=vmax, zorder = i + 5)
            ax[1][0].hexbin(this_y_damped, this_z_damped, bins="log", cmap=cmaps[i], gridsize=(gridsize, gridsize),
                            extent=[*xextent, *yextent], 
                            vmin=vmin, vmax=vmax, zorder = i + 5)


            ax[1][1].text(xlims[1] - dx/2,
                          ylims[1] - (i + 2) * dy / 10, 
                          f"d < {damped[i]}: {len(this_x_damped)}", 
                          ha="left", color=colors[i])

        add_labels(ax)

    animate(0)
    plt.tight_layout()
    
    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=100)
    ani.save(outname, writer='pillow', fps=24)

    if close:
        plt.close(fig)


def animated_surface_density_plot(orbits, masses, x_ind=1, y_ind=2, n_frames=100, **kwargs):
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
        ax[0][1].set_xlabel("R [kpc]")
        ax[0][1].set_ylabel(r"$\Sigma(r)$")
        ax[1][1].set_xlabel("Time [Myr]")
        ax[1][1].set_ylabel(r"RP Profile [g/ (cm s$^2$)]")

    
    fig, ax = plt.subplots(2, 2, facecolor="white", figsize=figsize)

    frames = np.linspace(0, len(orbits.data.t) - 1, n_frames).astype(int)

    x,y,z,vx,vy,vz = orbits.get_orbit_data(transposed=False)   

    wind_profile = np.sqrt(np.sum(orbits.metadata["WIND"].evaluate_arr(orbits.data.t) ** 2, axis=1)) * u.kpc/u.Myr
    wind_profile = wind_profile.to(u.cm/u.s)

    dens_profile = orbits.metadata["RHO_ICM"].evaluate_arr(orbits.data.t) * (u.g / u.cm**3)
    rp_profile = wind_profile ** 2 * dens_profile
    
    rs, surface_density_profiles = analysis.calc_surface_density_profile(orbits, masses)

    mdisk_evol = analysis.Mdisk(orbits, masses=masses)[1]
    print(len(mdisk_evol))
    
    ax[0][0].hexbin(x[0], y[0], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                    extent=[xextent[0], xextent[1], xextent[0], xextent[1]], 
                    vmin=vmin, vmax=vmax, zorder = 5)
    
    ax[1][0].hexbin(x[0], z[0], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                    extent=[xextent[0], xextent[1], xextent[0], xextent[1]], 
                    vmin=vmin, vmax=vmax, zorder = 5)
    
    ax[0][1].plot(rs, surface_density_profiles[0], color="black", lw=2, label=f"M_disk={mdisk_evol[0]:.2e} Msun")
    ax[0][1].set(ylim=kwargs.get("sd_ylim", (1e-6, 5e3)), yscale="log")
    ax[0][1].legend(loc="upper right")  

    ax[1][1].plot(orbits.data.t, rp_profile, color="black", lw=2)
    
    add_labels(ax)

    plt.tight_layout()

    
    def animate(i):
        this_x, this_y, this_z = x[i], y[i], z[i]

        for axis in ax.flatten()[:]:
            axis.cla()

        ax[0][0].hexbin(x[i], y[i], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                        extent=[xextent[0], xextent[1], xextent[0], xextent[1]], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        
        ax[1][0].hexbin(x[i], z[i], bins="log", cmap=cmap, gridsize=(gridsize, gridsize),
                        extent=[xextent[0], xextent[1], xextent[0], xextent[1]], 
                        vmin=vmin, vmax=vmax, zorder = 5)
        
       
        ax[0][1].plot(rs, surface_density_profiles[i], color="black", lw=2, label=f"M_disk={mdisk_evol[i]:.2e} Msun")
        ax[0][1].set(ylim=kwargs.get("sd_ylim", (1e-6, 5e3)), yscale="log")
        ax[0][1].legend(loc="upper right")
        

        ylims = ax[0][0].get_ylim()
        ax[1][1].plot(orbits.data.t, rp_profile, color="black", lw=2)

        current_time = orbits.data.t[i].value
        ylims = ax[1][1].get_ylim()
        ax[1][1].plot([current_time, current_time], [ylims[0], ylims[1]], color="red", lw=1, ls="dashed")

        add_labels(ax)

    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=100)
    ani.save(outname, writer='pillow', fps=24)

    if close:
        plt.close(fig)