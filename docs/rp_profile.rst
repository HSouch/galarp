Determining the Ram Pressure Profile
====================================

The key factor for a ram pressure interaction is the *ram pressure profile*. An infalling galaxy will experience a 
time-variable wind strength. The ram pressure goes as :math:`P_\text{ram} = \rho v^2`. As a galaxy approaches 
pericenter, both the local medium density and the velocity increases to a peak at pericenter, where it then steadily
decreases. Therefore, implementing a time-variable strength of ram pressure on a galaxy is considerably more physically
realistic than a more simple *constant wind*. 

Getting a Ram Pressure Profile from an Infalling Orbit
------------------------------------------------------

**GalaRP** has built-in classes to make this proceess seamless and easy for users, whether they want to initialize a 
ram pressure profile from a Gala-integrated infalling orbit, or manually defining the profile 

A typical density profile for a galaxy cluster is given by the Beta profile:

.. math::
   n(r) = n_0\biggl[1 + \biggl(\frac{r}{r_c}\biggr)^2\biggr]^{-3\beta / 2}

where :math:`n_0` is the central density and :math:`r_c` is the core radius. This can be converted into a mass density
by assuming a gas mass, commonly by multiplying by an atomic mass unit (in grams) to assume ionized atomic Hydrogen.

To make a ram pressure profile, we can use the **GalaRP** HostOrbit class. To make complete use of this class we have
to specify the host potential, the initial position and velocity conditions, and the host density model. 
Like with some of our other examples, we 
will use initial conditions given in `Zhu+2023 <https://iopscience.iop.org/article/10.3847/1538-4357/acfe6f/pdf/>`_.


.. code-block:: python

    import galarp as grp

    test_host_potential = grp.builtins.JZ2023_1e14()
    init_conditions = grp.builtins.JZ2023_1e14_IC()
    
    host_1e14 = grp.HostOrbit(potential=test_host_potential, 
                            init_conditions=init_conditions,
                            density=grp.SphericalBetaModel())


To generate our orbit, simply call :code:`host.integrate()` on your HostOrbit object. You can also generate a simple
plot showing the orbit, along with the computed density and ram pressure profiles.

.. code-block:: python

    host_1e14.integrate(n_steps=10000)          # Integrate orbits
    host_1e14.plot_infall_orbit()               # Plot orbits

.. raw:: html

   <img src="_static/plots/infall_orbit.png" width="80%%"  style="margin-bottom: 32px;"/>


