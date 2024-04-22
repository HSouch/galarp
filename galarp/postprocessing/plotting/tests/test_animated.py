
import galarp as grp

class TestAnimatedPlots:
    orbits = grp.builtins.ExampleLorentzianEvent()

    def test_animated_hexbin(self, tmp_path):
        grp.plotting.animated_hexbin_plot(self.orbits, n_frames=25, outname=f'{tmp_path}animated_hexbin.gif')

    def test_r_vr(self, tmp_path):
        grp.plotting.r_vr(self.orbits, n_frames=25, outname=f'{tmp_path}r_vr.gif')