
from .. import hosts
from .. import satellites

import gala.potential as gp

class BuiltInPotentialTestBase:

    name = None
    potential = None
    show_plots = False

    def test_creation(self):
        assert self.potential is not None


class TestJZ2023_10e12(BuiltInPotentialTestBase):

    name = "JZ2023_10e12"
    potential = hosts.JZ2023_10e12()


class TestJZ2023_10e13(BuiltInPotentialTestBase):
    
        name = "JZ2023_10e13"
        potential = hosts.JZ2023_10e13()

class TestJZ2023_10e14(BuiltInPotentialTestBase):
         
        name = "JZ2023_10e14"
        potential = hosts.JZ2023_10e14()
