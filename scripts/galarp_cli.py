import os

import galarp as grp

from gala import potential as gp
from gala.units import galactic

from astropy import units as u

import argparse

import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Run GalaRP in a CLI')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    return parser.parse_args()


def create_potential():
    potential_list = ["Burkert", "NFW", "Miyamoto-Nagi", "Hernquist", "Kepler", None, None, None, None, None, None]
    

    for i, potential in enumerate(potential_list):
        if potential is not None:
            print(f"    {potential} ({i})")
    print(f"\n    None (9)")

    def add_individual(name="dark matter", default=1):

        selection = int(input(f"Enter the  potential to add for the {name} (default={potential_list[default]}, {default}):") or default)

        if selection == 0:          # Handle Burkert potential
            r0 = float(input("  Enter the scale radius r0 in kpc (8):") or 8.)
            return gp.BurkertPotential.from_r0(r0=r0 * u.kpc, units=galactic)
        elif selection == 1:        # Handle NFW potential
            m200 = float(input("  Enter the virial mass M200 in Msun (1.2e15):")) or 1.2e15
            c = float(input("  Enter the concentration c (9):")) or 9
            print(f"Creating NFW potential with M200={m200:.2e} Msun and c={c}")
            return gp.NFWPotential.from_M200_c(M200=m200 * u.Msun, c=c, units=galactic)
        elif selection == 2:        # Handle Miyamoto-Nagai potential
            mass = float(input("  Enter the mass (Msun) (1e11):") or 1e11)
            a = float(input("  Enter the scale radius a (kpc) (5):") or 5)
            b = float(input("  Enter the scale height b (kpc) (0.5):") or 0.5)
            print(f"Creating Miyamoto-Nagai potential with mass={mass:.2e} Msun, a={a} kpc, and b={b} kpc")
            return gp.MiyamotoNagaiPotential(m=mass * u.Msun, a=a * u.kpc, b=b * u.kpc, units=galactic)
        elif selection == 3:        # Handle Hernquist potential
            mass = float(input("  Enter the mass (Msun) (1e10):") or 1e10)
            c = float(input("  Enter the scale radius c (kpc) (0.4):") or 0.4)
            print(f"Creating Hernquist potential with mass={mass:.2e} Msun and c={c} kpc")
            return gp.HernquistPotential(m=mass * u.Msun, c=c * u.kpc, units=galactic)
        elif selection == 4:        # Handle Kepler potential
            mass = float(input("  Enter the mass (Msun) (1e9):") or 1e9)
            print(f"Creating Kepler potential with mass={mass:.2e} Msun")
            return gp.KeplerPotential(m=mass * u.Msun, units=galactic)
        elif selection == 9:
            return gp.NullPotential(units=galactic)
        else:
            print("Invalid selection, try again.")
            return add_individual(name, default)
    

    dm = add_individual("dark matter", default=0)
    gas = add_individual("gas", default=2)
    stars = add_individual("stars", default=2)
    bulge = add_individual("stellar bulge", default=9)

    combined = gp.CompositePotential(dm=dm, gas=gas, stars=stars, bulge=bulge)
    return combined


def create_wind():    
    wind_list = ["Uniform", "Lorentzian", "Interpolated", "None", None, None, None, None, None, None, None]


    for i, wind in enumerate(wind_list):
        if wind is not None:
            print(f"    {wind} ({i})")
    print(f"\n    None (9)")

    selection = input("Enter the number of the wind model to use (default=0):") or 0

    if selection in [0, 1]:
        strength = float(input("  Enter the wind strength in km/s (800)") or 800.)
        inclination = float(input("  Enter the wind inclination in degrees (0: face-on)") or 90.)
        
        if selection == 0:
            return grp.ConstantWind(inclination=np.deg2rad(inclination), strength=strength * u.km/u.s, units=galactic)
        elif selection == 1:
            t0 = float(input("  Enter the time at peak wind strength in Myr (200)") or 500.)
            width = float(input("  Enter the width of the Lorentzian profile in Myr (200)") or 200.)
            return grp.LorentzianWind(t0=t0 * u.Myr, width=width * u.Myr, 
                                      inclination=np.deg2rad(inclination), strength=strength * u.km/u.s, units=galactic)
        elif selection == 2:
            fn = str(input("  Enter the path to the interpolated wind file:"))
            if fn is None:
                print("No file provided, defaulting to None.")
                return None
            
            return grp.InterpolatedStrengthWind.from_table(fn)


def create_particles(potential):
    mass_profile = grp.gen_mass_profile(potential)

    position_list = ["UniformGrid", "ExponentialGrid", "None", None, None, None, None, None, None, "From File", None]

    for i, position in enumerate(position_list):
        if position is not None:
            print(f"    {position} ({i})")
    print(f"\n    None (9)")

    selection = int(input("Enter the format of the particle distribution to use (0):") or 0)

    if selection == 0:
        Rmax = float(input("    Enter the maximum radius in kpc (10):") or 10.)
        spacing = float(input("    Enter the spacing in kpc (0.5):") or 0.5)

        positions = grp.PlaneDistribution(Rmax=Rmax, spacing=spacing)
    elif selection == 1:
        h_R = float(input("    Enter the radial scale length in kpc (4):") or 4.)
        h_z = float(input("    Enter the vertical scale length in kpc (0.5):") or 0.5)
        n_particles = int(input("    Enter the number of particles to generate (100):") or 100)

        positions = grp.ExponentialDistribution(h_R=h_R * u.kpc, h_z=h_z * u.kpc, 
                                                n_particles=n_particles)
    elif selection == 8:
        fn = str(input("    Enter the path to the particle file:"))
        if fn is None:
            print(    "No file provided, defaulting to None.")
            return None
        return grp.ParticleDistribution.from_file(fn)
    
    
    



    pset = grp.ParticleSet(particles=positions,  units=galactic)
    pset.generate(mass_profile)

    return pset
    


###########################################
##      Main function      ################
###########################################

def run_galarp():
    args = parse_args()
    
    print("Welcome to the GalaRP CLI, let's get started.\n")
    
    config = str(input("Enter the name of the configuration file (None):") or None)
    
    if config is None:
        print("No configuration file provided, defining components manually.")
    
    # # Add potential
    # print("\nDefine the gravitational potential:")
    # potential = create_potential()

    # # Define the wind
    # print("\nDefine the wind model:")
    # wind = create_wind()

    # Define the particles
    print("\nDefine the constituent particles:")
    particles = create_particles(grp.builtins.JZ2023_Satellite())

    # Define the shadow
    print("\nDefine the shadow model:")

    # Extra stuff (output dirs and whatnot)
    outdir = str(input("Enter the output directory (output):") or "output")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    


    # Run the damn thing
    # sim = grp.RPSim(potential=potential, wind=wind)
    # orbits = sim.run(particles)

    
    if args.debug:
        print("Config:", config)



if __name__ == "__main__":
    run_galarp()