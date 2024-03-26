# Testing utilities


# Third Party
import os
import numpy as np


# This module
import galarp as grp
from . import utils

from gala.units import galactic
import astropy.units as u


__all__ = ["test_uniform_shadow"]


def test_uniform_shadow(plot_dir="plots/shadow_tests/"):
    """Test that the angled shadows being created are correct"""

    os.makedirs(plot_dir, exist_ok=True)

    for angle in range(10, 90, 10):
        wind = grp.RPWind(units=galactic)
        wind.init_from_inc(np.deg2rad(angle), 300 * u.km / u.s)

        shadow = grp.UniformShadow()
        shadow.init_from_wind(wind)
        shadow.plot_shadow(wind=wind, outname=f"{plot_dir}shadow_{angle}.pdf")
