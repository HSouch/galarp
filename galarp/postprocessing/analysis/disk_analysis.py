from ...utils import get_orbit_data
from ...rampressure import OrbitContainer


import numpy as np
from scipy import stats

from astropy import units as u
from astropy.modeling import models, fitting


def calculate_rstrip(xyz, rmax=20, zmax=2, frac=0.8):
    x,y,z = xyz
    r = np.sqrt(x**2 + y**2 + z**2)
    this_r_cut = r[np.abs(z) < zmax]
    this_r_cut = this_r_cut[this_r_cut < rmax]
    
    try:
        cdf = stats.ecdf(this_r_cut)
        cdf_xs, cdf_vals = cdf.cdf.quantiles, cdf.cdf.probabilities
        return cdf_xs[np.argmin(np.abs(cdf_vals - frac))]
    except Exception:
        return 0
    

def rstrip_and_medians(xyz, rmax=20, zmax=2, frac=0.9):
    x,y,z = xyz
    r = np.sqrt(x**2 + y**2 + z**2)

    rcut = (np.abs(z) < zmax) & (r < rmax)
    
    try:
        cdf = stats.ecdf(r[rcut])
        cdf_xs, cdf_vals = cdf.cdf.quantiles, cdf.cdf.probabilities
        rstrip =  cdf_xs[np.argmin(np.abs(cdf_vals - frac))]
    except Exception:
        rstrip = 0
    
    median_x, median_y, median_z = np.median([x[rcut], y[rcut], z[rcut]], axis=1)
    
    return rstrip, [median_x, median_y, median_z]


def rstrip(orbits, zmax=2 * u.kpc, r_strip_frac=0.8, rmax=20 * u.kpc):
    """
    Calculate the evolution of the stripping radius for a given OrbitContainer.

    Parameters:
    - orbits (str or OrbitContainer): The orbits to analyze. If a string is provided, it is assumed to be the path to 
                                      a file containing the orbits.
    - zmax (Quantity, optional): The maximum absolute value of z-coordinate to consider. Defaults to 2 kpc.
    - r_strip_frac (float, optional): The fraction of particles in the cumultative distribution function
                                      that sets the stripping radius. Defaults to 0.8 (80%).
    - rmax (Quantity, optional): The maximum radius to consider. This is a global radius to ensure that edge-on events
                                 (which don't incude strong vertical orbits) are still computed reasonably using this
                                 method. Defaults to 20 kpc.

    Returns:
    - strip_times (list): The times at which the stripping radius was evaluated.
    - rstrips (list): The stripping radii at each measured time.
    """

    if isinstance(orbits, str): 
        orbits = OrbitContainer.load(orbits)

    x,y,z, *_ = get_orbit_data(orbits.data)
    
    x,y,z = x.T, y.T, z.T
    
    r = np.sqrt(x**2 + y**2 + z**2)
    
    times = orbits.data.t
    
    strip_times, rstrips = [], []
    for i in range(0, len(r), 5):
        strip_times.append(times[i].value)
        rstrips.append(calculate_rstrip((x[i], y[i], z[i]), rmax=rmax.value, zmax=zmax.value, frac=r_strip_frac))
    return strip_times, rstrips



def final_rstrip(orbits, zmax=2 * u.kpc, r_strip_frac=0.8, rmax=20 * u.kpc):
    """
    For a given OrbitContainer, calculate the final stripping radius by fitting a Gaussian to the overall 
    time series.

    Parameters:
    - orbits (str or OrbitContainer): The orbits to analyze. If a string is provided, it is assumed to be the path to 
                                      a file containing the orbits.
    - zmax (Quantity, optional): The maximum absolute value of z-coordinate to consider. Defaults to 2 kpc.
    - r_strip_frac (float, optional): The fraction of particles in the cumultative distribution function
                                      that sets the stripping radius. Defaults to 0.8 (80%).
    - rmax (Quantity, optional): The maximum radius to consider. This is a global radius to ensure that edge-on events
                                 (which don't incude strong vertical orbits) are still computed reasonably using this
                                 method. Defaults to 20 kpc.

    Returns:
    - The final value for the stripping radius.
    """

    xs, ys = rstrip(orbits, zmax=zmax, r_strip_frac=r_strip_frac, rmax=rmax)

    @models.custom_model
    def PedastalGauss(x, amplitude=1, mean=0, stddev=1, pedastal=0):
        return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2)) + pedastal
    
    mod = PedastalGauss(amplitude = np.nanmax(ys) - np.nanmin(ys), mean = 0, stddev = np.nanmax(xs) / 5, 
                        pedastal = np.nanmin(ys))

    fit = fitting.LevMarLSQFitter()

    mod_fit = fit(mod, xs, ys)

    rstrip_final = mod_fit.pedastal.value       # The final stripping radius is the pedastal of the Gaussian model

    return rstrip_final


def stripped(orbits, method="vrad"):
    """
    Calculate the number of stripped particles in the given orbits.

    Parameters:
    - orbits (OrbitCollection): The collection of orbits.

    Returns:
    - stripped_particles (numpy.ndarray): An array containing the number of stripped particles for each time step.

    """
    x, y, z, vx, vy, vz = get_orbit_data(orbits.data, transposed=False)

    r = np.sqrt(x**2 + y**2 + z**2)

    pot = orbits.metadata["POTENTIAL"]
    energy = pot.energy([x, y, z])

    v_escape = np.sqrt(2 * np.abs(energy)).to(u.km/u.s).value

    if method == "vrad":
        v_particle = (x * vx + y * vy + z * vz) / r
    elif method == "vtot":
        v_particle = np.sqrt(vx**2 + vy**2 + vz**2)

    stripped = v_particle > v_escape

    return np.sum(stripped, axis=1) / len(x[0])


def Mdisk(orbits, masses=None, rmax = 27 * u.kpc, thickness=1 * u.kpc, **kwargs):

    x,y,z, *_ = orbits.get_orbit_data(transposed=False)

    r = np.sqrt(x**2 + y**2)

    in_disk = np.logical_and(r < rmax.value, np.abs(z) < thickness.value)

    evol = []
    for i in range(len(orbits.data.t)):
        if masses is None:
            evol.append(np.sum(in_disk[i]) * orbits.metadata["M_CLOUD"].to(u.Msun).value)
        else:
            evol.append(np.sum(masses[in_disk[i]]))
    
    return orbits.data.t, np.array(evol)


def mass_profile(orbits, masses, z_width = 5 * u.kpc, rmin=0, rmax=40, nbins=21, **kwargs):
    x,y,z, *_ = orbits.get_orbit_data(transposed=False)

    rs = np.linspace(rmin, rmax, nbins)
    diff = np.mean(np.diff(rs))

    mass_profiles = []

    for i in range(0, len(orbits.data.t), 10):
        this_x, this_y, this_z = x[i], y[i], z[i]
        z_cut = np.abs(this_z) < z_width.value

        this_x, this_y = this_x[z_cut], this_y[z_cut]
        this_mass = masses[z_cut]

        r = np.sqrt(this_x**2 + this_y**2)

        binned_rs = np.digitize(r, rs)

        m_r = []
        for i in range(1, len(rs)):
            in_bin = binned_rs == i
            m_r.append(np.sum(this_mass[in_bin]))
        mass_profiles.append(m_r)

    
    return np.array(mass_profiles)
        

def calc_surface_density_profile(orbits, masses, z_width = 5 * u.kpc, rmin=0, rmax=40, nbins=51, **kwargs):
    x, y, z, *_ = orbits.get_orbit_data(transposed=False)

    rs = np.linspace(rmin, rmax, nbins)
    diff = np.mean(np.diff(rs))

    surface_density_profiles = []

    for i in range(0, len(orbits.data.t)):
        this_x, this_y, this_z = x[i], y[i], z[i]
        z_cut = np.abs(this_z) < z_width.value

        this_x, this_y = this_x[z_cut], this_y[z_cut]
        
        this_mass = masses[z_cut]
       # this_mass = np.copy(masses)

        r = np.sqrt(this_x**2 + this_y**2)

        binned_rs = np.digitize(r, rs)

        sd_r = []
        for i in range(1, len(rs)):
            in_bin = binned_rs == i
            sd_r.append(np.sum(this_mass[in_bin]) / (np.pi * (rs[i] ** 2 - rs[i-1] ** 2)))
        surface_density_profiles.append((sd_r * u.Msun / u.kpc**2).to(u.Msun / u.pc**2).value)
    
    return rs[:-1] + diff/2, np.array(surface_density_profiles)
