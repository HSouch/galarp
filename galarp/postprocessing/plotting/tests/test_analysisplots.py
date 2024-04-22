import galarp as grp


class TestAnalysisPlots:
    orbits = grp.builtins.ExampleLorentzianEvent()

    def test_rstrip(self, tmp_path):
        grp.plotting.rstrip_plot(self.orbits, outname=f'{tmp_path}rstrip.png')
    
    def test_stripped_plot(self, tmp_path):
        grp.plotting.stripped_plot(self.orbits, outname=f'{tmp_path}stripped.png')