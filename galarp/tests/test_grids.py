from .helpers import ParticleGridTestBase
from ..particle_grids import *

from ..builtins.satellites import JZ2023_Satellite
from ..utils import gen_mass_profile


class TestUniformGrid(ParticleGridTestBase):
    name = "UniformGrid"
    grid = UniformGrid(Rmax=10, n_particles=50, z_start=0.0, veldisp=10.0)
    

    example_potential = JZ2023_Satellite()
    mass_profile = gen_mass_profile(example_potential)

    grid.generate(mass_profile, velocities=True)