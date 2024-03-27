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

