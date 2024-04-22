import numpy as np

from .helpers import ShadowTestBase
from .. import shadows 

import pytest


class TestShadowBase(ShadowTestBase):
    shadow = shadows.ShadowBase()
    
    @pytest.mark.skip(reason="Abstract class")
    def test_evaluation(self):
        return super().test_evaluation()
    
    @pytest.mark.skip(reason="Abstract class")
    def test_plotting(self):
        return super().test_plotting()


class TestUniformShadow(ShadowTestBase):
    shadow = shadows.UniformShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0


class TestExponentialShadow(ShadowTestBase):
    shadow = shadows.ExponentialShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0


class TestEdgeOnShadow(ShadowTestBase):
    shadow = shadows.EdgeOnShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0
    

class TestUniformLinearZVariableShadow(ShadowTestBase):
    shadow = shadows.UniformLinearZVariableShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0


class TestUniformExponentialZVariableShadow(ShadowTestBase):
    shadow = shadows.UniformExponentialZVariableShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0