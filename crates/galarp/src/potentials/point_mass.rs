//! Point mass (Kepler) potential, Phi = -G M / r. Mostly for tests.

use super::{norm, radial, Param, Potential};
use crate::units::G;

#[derive(Clone, Debug)]
pub struct PointMass {
    pub mass: Param,
}

impl Potential for PointMass {
    fn potential(&self, x: &[f64; 3], t: f64) -> f64 {
        -G * self.mass.at(t) / norm(x)
    }

    fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3] {
        let r = norm(x);
        radial(x, r, G * self.mass.at(t) / (r * r * r))
    }
}

