
import numpy as np
from .helpers import RPWindTestBase
from ..winds import RPWind, LorentzianWind, StepFunctionWind, InterpolatedWind


class TestRPWind(RPWindTestBase):
    wind = RPWind()
    t = 0
