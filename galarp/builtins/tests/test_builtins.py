
from .. import hosts
from .. import satellites

from ...utils import Rvir, R100, R500

from astropy import units as u
import gala.dynamics as gd
import gala.potential as gp


class BuiltInPotentialTestBase:

    name = None
    potential = None
    show_plots = False

    pos = [1, 0, 0]
    vel = [10, 10, 0]

    def test_creation(self):
        assert self.potential is not None
    
    def test_potential_integration(self, nsteps=1000):
        w0 = gd.PhaseSpacePosition(pos=self.pos, vel=self.vel)
        orbit = gp.Hamiltonian(potential=self.potential).integrate_orbit(w0, dt=1, n_steps=nsteps)
        
        assert len(orbit.pos) == nsteps + 1
        assert type(orbit) == gd.Orbit
    
    def test_rs(self):
        mass = self.potential.mass_enclosed([500, 0., 0.] * u.Mpc)
        


class TestJZ2023_1e12(BuiltInPotentialTestBase):

    name = "JZ2023_1e12"
    potential = hosts.JZ2023_1e12()


class TestJZ2023_1e13(BuiltInPotentialTestBase):
    
        name = "JZ2023_1e13"
        potential = hosts.JZ2023_1e13()


class TestJZ2023_1e14(BuiltInPotentialTestBase):
         
        name = "JZ2023_1e14"
        potential = hosts.JZ2023_1e14()


class TestJZ2023_Satellite(BuiltInPotentialTestBase):
    
    name = "JZ2023_Satellite"
    potential = satellites.JZ2023_Satellite()
