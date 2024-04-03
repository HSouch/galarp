Plotting Results
================

``GalaRP`` has a wealth of plotting functions to plot all stages of the RP interaction, 
from the initial infall conditions, to the RP simulation setup, to the particle grids,
to the final orbits. 

``GalaRP`` has a default pyplot style which you can load as follows:

.. code-block:: python

    import galarp as grp
    grp.pyplot_style()

You can also load in the colormap used for some of these figures with the following:

.. code-block:: python

    cmap = grp.lavender_cmap(step_1=50)

which can be tuned to your liking by adjusting the ``step_1`` parameter.

.. toctree::
   :maxdepth: 1

   genplots
   matrixplots
   aestheticplots
