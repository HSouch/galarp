import galarp as grp

from gala.units import galactic

from astropy import units as u
from numpy import deg2rad


class TestRPSim:
    potential = grp.builtins.JZ2023_Satellite()
    wind = grp.LorentzianWind(units=galactic)
    wind.init_from_inc(inclination=deg2rad(45), strength = 800 * u.km / u.s)

    particles = grp.UniformGrid(n_particles=12)
    particles.generate(grp.gen_mass_profile(potential))

    density = grp.Density(rho=1e-25 * u.g / u.cm ** 3)
    def test_sim(self):
        sim = grp.RPSim(wind=self.wind, potential=self.potential)

        sim.run(particles=self.particles)

    def test_sim_with_density(self):
        sim = grp.RPSim(wind=self.wind, potential=self.potential)
        sim.run(particles=self.particles, rho_icm=self.density)
        
        