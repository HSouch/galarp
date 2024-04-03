Matrix Plots ``plotting.matrixplots`` 
=====================================

Density Plots
-------------

After :doc:`obtaining the results of an RPSim <getting_started>` , it is often helpful to explore the distribution of particles, both along
principle axes as well as in phase-space. To do this, we can generate a matrix plot, which shows this distribution and
how it evolves at different "snapshots" of a ``GalaRP`` simulation. This is *especially* useful when comparing to 
snapshots of a hydrodynamical simulation.

.. code-block:: python

    import galarp as grp

    orbits = sim.run(...)

    grp.density_matrix(orbits, outname="matrix_xy.png")


This will plot the x-y plane by default. You can adjust the x-axis and y-axis plot by changing ``x_ind`` and ``y_ind``, 
where the index of what is plotted corresponds to the output of ``grp.get_orbit_data(orbits)``, where they are
organized as :math:`(x,y,z,v_x,v_y,v_z)`. Thus, to plot the :math:`x-z` plane, you would set ``y_ind=2``, and to
plot the :math:`z-v_z` phase space, you would set ``x_ind=2`` and ``y_ind=5``.

.. raw:: html

   <img src="_static/plots/density_xy_matrix.png" width="45%%"  style="margin-bottom: 32px;"/>
   <img src="_static/plots/density_zvz_matrix.png" width="45%%"  style="margin-bottom: 32px;"/>

