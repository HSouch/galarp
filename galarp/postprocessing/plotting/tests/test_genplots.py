from .. import general_plots as genplots

import galarp as grp

class TestDensityPlots:
    pot = grp.builtins.JZ2023_Satellite()
    mass_profile = grp.gen_mass_profile(pot)
    particles = grp.ExponentialGrid(n_particles=200)
    particles.generate(mass_profile=mass_profile)

    def test_2ax_plot(self, tmp_path):
        d = tmp_path / "sub"
        d.mkdir()
        genplots.plot_density(self.particles.get_xyz(), gridsize=10, outname=f'{d}_2ax_test.png')
    
    def test_3ax_plot(self, tmp_path):
        d = tmp_path / "sub"
        d.mkdir()
        genplots.plot_density_3ax(self.particles.get_xyz(), gridsize=10, outname=f'{d}_3ax_test.png')
