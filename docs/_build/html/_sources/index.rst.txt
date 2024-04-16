.. GalaRP documentation master file, created by
   sphinx-quickstart on Tue Mar 26 15:11:59 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. raw:: html

   <img src="_static/Galarp-Logo.png" width="50%"  style="margin-bottom: 32px;"/>

.. module:: gala

******
GalaRP
******

Ram pressure is a physical process that a galaxy experiences when it travels through a 
gaseous medium, 
such as the intergalactic medium or the intracluster medium. The gas is stripped from the galaxy, 
which can lead to the formation of a tail of gas and stars. This process is known as ram pressure stripping. 
This process leads to the creation of "RPS Galaxies": some of the most visually striking objects in galaxy evolution.

``GalaRP`` leverages the `Gala python package <http://gala.adrian.pw/en/latest/>`_, combined with nonconservative orbit 
integration,
to simulate the larger-scale effects of ram pressure absent hydrodynamical instabilities. Its primary purpose is to 
be an efficient tool for finding optimal setup conditions for hydrodynamical wind-tunnel simulations of ram pressure
stripping.

The central purpose of ``GalaRP`` is to simulate (to first order) the expected large-scale conditions of an RP 
interaction from a hydrodynamical simulation. 

Development is active on `Github <https://github.com/HSouch/galarp/>`_, and the Authors encourage
updates and suggestions to improve and enhance its functionality. Furthermore, some videos showing what ``GalaRP`` can 
do are available (and can be used in presentations with credit given) on 
`Youtube <https://www.youtube.com/@GalaRPS>`_.

.. toctree::
   :maxdepth: 1

   installation
   getting_started


.. toctree::
   :maxdepth: 2

   running



.. toctree::
   :maxdepth: 2

   postprocess



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
