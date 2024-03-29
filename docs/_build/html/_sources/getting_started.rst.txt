Getting Started with GalaRP
===========================

Getting started with GalaRP is pretty simple. The following steps are needed for a complete simulation:

1. Create a satellite galaxy potential.
2. Determine the ram pressure wind model.
3. Run the simulation.

and that's it. The following sections will guide you through these steps.

Step 1: Create a Satellite Galaxy Potential
-------------------------------------------


To create a satellite galaxy potential, you can follow this example. Here, we are initializing 
a potential that is similar to the dwarf galaxy used in
`Zhu+2023 <https://iopscience.iop.org/article/10.3847/1538-4357/acfe6f/pdf/>`_.

>>> import galarp as grp

GalaRP uses Gala for all of its gravitational potentials, though there are some built-in potentials available as 
defaults, based on various papers.

>>> from gala import potential as gp
>>> from gala.units import galactic

In this case the Zhu+2023 potential is already available as a built-in, but it is just a combination of an NFW DM halo,
a Miyamoto-Nagai stellar disk, and an exponential gas disk. You define these in the exact same way as one would in
`Gala <https://gala-astro.readthedocs.io/en/latest/potential/compositepotential.html>`_.

>>> pot = grp.builtins.JZ2023_Satellite()


Step 2: Define the RP Interaction Parameters
--------------------------------------------

A ram pressure interaction in GalaRP is primarily done by assuming spherical gas blobs, and the Gunn and Gott (1972)
stripping criterion.


.. math::
   P_{\text{ram}} = \rho_{\text{ICM}} v^2

where :math:`\rho_{\text{ICM}}` is the ICM density, and :math:`v` is the relative velocity between the satellite and 
the medium it is travelling through.

For this example, we will define a simple time-variable *Lorentzian wind* that follows a Lorentzian profile, similar to
the computed ram pressure profiles for a galaxy's orbit through a cluster.

.. math::
   v_{i}(t) = \frac{v_{i, \text{max}}}{(2(t - t_0) / w)^2 + 1}

where :math:`t_0` is the time of peak ram pressure, which we set as 400 Myr into the simulation,
and :math:`w` is the characteristic width of the Lorentzian profile, which we set as a Gyr. The angle of 
the wind is 45 degrees, and the peak velocity is 700 km/s. We first initialize the object, and define the wind
vector using :code:`init_from_inc`.


>>> phi = 45.0
>>> strength = 700
>>> 
>>> wind = gn.LorentzianWind(units=galactic, t0=400 * u.Myr, width=1000 * u.Myr)
>>> wind.init_from_inc(np.deg2rad(phi), strength * u.km  / u.s)


We also will specify an RP shadow, which models the reduction in ram pressure strength because the wind runs into the 
disk. We will specify that the wind strength is reduced to 40% of the total when inside the shadow, which is
effectively a "tilted cylinder" along the wind direction. 

>>> shadow = gn.UniformShadow(damping=0.2)
>>> shadow.init_from_wind(wind)


Step 3: Define the Parameter Grid
---------------------------------

Lastly, we need to define the particles that we are simulating. In this case, we will use a simple uniform grid of 
evenly-spaced particles, which are helpful for showing `where` in the disk particles will get perturbed or stripped.
To initialize the particles we need a mass profile from the potential, which we can get using GalaRP.

>>> mass_profile = grp.gen_mass_profile(pot)
>>> uniform = gn.UniformGrid(n_particles=60)
>>> uniform.generate(mass_profile=mass_profile)


Step 4: Run the Simulation
--------------------------

With all the components in hand we call the RPSim object to run the ram pressure interaction. We will run the 
simulation for 2 Gyr at 2 Myr spacing.


>>> sim = grp.RPSim(wind=wind, potential=pot, shadow=shadow)
>>> orbits = sim.run(test, wind_on=True, rho_icm=1e-26 * u.g / u.cm**3,
                     integration_time=2000 * u.Myr, dt=2 * u.Myr, 
                     m_cloud=1e6 * u.Msun, printout=True)

The RP simulation spits out an OrbitContainer object which contains both the integrated Gala orbits, as well as a
metadata dictionary object that contains the parameters used in the simulation. These can be plotted simply, or 
passed to GalaRP's postpcrocessing tools for more detailed analysis.

>>> orbits.data.plot()

