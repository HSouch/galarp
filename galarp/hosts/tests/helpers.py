class TestHostBase:
    init_conditions = None

    def test_init_conditions(self):
        assert self.init_conditions is not None