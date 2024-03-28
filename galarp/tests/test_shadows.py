import numpy as np

from .helpers import ShadowTestBase
from ..shadows import *

class TestUniformShadow(ShadowTestBase):
    shadow = UniformShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0

class TestExponentialShadow(ShadowTestBase):
    shadow = ExponentialShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0

class TestEdgeOnShadow(ShadowTestBase):
    shadow = EdgeOnShadow()
    xyz = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    t = 0
    