from ..general_plots import *

import galarp as grp

import pytest

# class TestK3DPlot():

#     #@pytest.mark.filterwarnings("ignore")
#     def test_K3D_plot(self, tmp_path):

#         d = tmp_path / "sub"
#         d.mkdir()

#         orbits = grp.builtins.ExampleUniformEvent()

#         k3d_plot([orbits], outname=f'{d}.k3D_test.html')


class TestDensityPlots:
    pot = grp.builtins.JZ2023_Satellite()
    mass_profile = grp.gen_mass_profile(pot)
    particles = grp.ExponentialGrid(n_particles=200)
    particles.generate(mass_profile=mass_profile)

    def test_2ax_plot(self, tmp_path):
        d = tmp_path / "sub"
        d.mkdir()
        plot_density(self.particles.get_xyz(), gridsize=10, outname=f'{d}_2ax_test.png')
    
    def test_3ax_plot(self, tmp_path):
        d = tmp_path / "sub"
        d.mkdir()
        plot_density_3ax(self.particles.get_xyz(), gridsize=10, outname=f'{d}_3ax_test.png')