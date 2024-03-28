from gala import dynamics as gd
from gala import integrate as gi 
from gala import potential as gp


class HostOrbit:
    """ Class for a host orbit to determine a ram pressure profile.
    """
    def __init__(self, potential, init_conditions=None):
        self.potential = potential

        self.init_conditions = init_conditions

        self.orbit = None
    

    def integrate(self, n_steps=10000):
        assert self.init_conditions is not None

        w0 = gd.PhaseSpacePosition(pos=self.init_conditions.pos, vel=self.init_conditions.vel)

        self.orbit = gp.Hamiltonian(self.potential).integrate_orbit(w0, dt=1, n_steps=n_steps)
        
