import galarp as grp

from . import helpers

class TestInitConditions(helpers.TestHostBase):

    init_conditions = grp.builtins.JZ2023_1e14_IC()

    