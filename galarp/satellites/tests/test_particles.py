from astropy import units as u

from .. import particles


class SigmaGasTestBase:
    sg = None

    def test_units(self):
        assert self.sg.units == u.solMass / u.pc**2
        assert self.sg.sigma.unit == self.sg.units

    def test_sigma(self):
        assert self.sg.sigma.isscalar == False

    def test_length(self):   
        assert len(self.sg.sigma) == self.sg.nparticles


class TestSigmaGasNoUnits(SigmaGasTestBase):
    # Test class with just a single scalar value input (with no units)
    sg = particles.SigmaGas(sigma=20., nparticles=10)


class TestSigmaGasWithUnits(SigmaGasTestBase):
    # Test class with just a single scalar value input (with units)
    sg = particles.SigmaGas(sigma=4e-3 * u.g / u.cm **2, nparticles=10)


class TestSigmaGasList(SigmaGasTestBase):
    # Test class with a list of scalar values input (with no units)
    sg = particles.SigmaGas(sigma=[20, 50, 80])


class TestSigmaGasListWithUnits(SigmaGasTestBase):
    # Test class with a list of scalar values input (with units)
    sg = particles.SigmaGas(sigma=[1e-3, 1e-4, 1e-5] * u.g / u.cm **2)
    

class TestSigmaGasClassInput(SigmaGasTestBase):
    # Test class with another SigmaGas object as input
    sg = particles.SigmaGas(sigma=particles.SigmaGas(sigma=20, nparticles=1), nparticles=10)