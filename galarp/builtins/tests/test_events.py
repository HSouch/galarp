from .. import events

import galarp as grp


class TestBuiltinUniformEvent:
    
    def test_orbit(self):
        orbits = events.ExampleUniformEvent()

        assert orbits is not None
        assert isinstance(orbits.metadata["WIND"], grp.RPWind)


class TestBuiltinLorentzianEvent:
    
    def test_orbit(self):
        orbits = events.ExampleLorentzianEvent()

        assert orbits is not None
        assert isinstance(orbits.metadata["WIND"], grp.LorentzianWind)