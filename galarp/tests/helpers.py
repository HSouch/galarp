import pytest
from ..shadows import ShadowBase


class ShadowTestBase:
    name = None
    shadow = None
    show_plots = False

    xyz = None

    def test_evaluation(self):
        eval = self.shadow.evaluate(self.xyz)

        assert eval is not None
        assert eval.shape == self.xyz.shape
