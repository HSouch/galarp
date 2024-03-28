
import numpy as np
from astropy import units as u

from scipy.interpolate import interp1d

from .helpers import RPWindTestBase
from ..winds import RPWind, LorentzianWind, StepFunctionWind, InterpolatedWind


class TestRPWind(RPWindTestBase):
    wind = RPWind()
    t = 0


class TestLorentzianWind(RPWindTestBase):
    wind = LorentzianWind(t0=500 * u.Myr)
    t = 0


class TestStepFunctionWind(RPWindTestBase):
    wind = StepFunctionWind(t0=500 * u.Myr)
    

class TestInterpolatedWind(RPWindTestBase):
    xs = np.linspace(0, 1000, 100) * u.Myr
    lor = LorentzianWind(t0=500 * u.Myr)
    ys = []
    for x in xs:
        ys.append(np.sqrt(np.sum(lor.evaluate(x.value) ** 2, axis=0)))
    ys = np.array(ys)

    wind = InterpolatedWind()
    wind.from_xy(xs, ys)
