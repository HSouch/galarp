//! Miyamoto-Nagai (1975) disk, symmetric about the z axis:
//! Phi = -G M / sqrt(R^2 + (a + sqrt(z^2 + b^2))^2).

use super::{Param, Potential};
use crate::units::G;

#[derive(Clone, Debug)]
pub struct MiyamotoNagai {
    pub mass: Param,
    /// Radial scale length.
    pub a: Param,
    /// Vertical scale height.
    pub b: Param,
}

impl Potential for MiyamotoNagai {
    fn potential(&self, x: &[f64; 3], t: f64) -> f64 {
        let (a, b) = (self.a.at(t), self.b.at(t));
        let s = a + (x[2] * x[2] + b * b).sqrt();
        -G * self.mass.at(t) / (x[0] * x[0] + x[1] * x[1] + s * s).sqrt()
    }

    fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3] {
        let (a, b) = (self.a.at(t), self.b.at(t));
        let zb = (x[2] * x[2] + b * b).sqrt();
        let s = a + zb;
        let d2 = x[0] * x[0] + x[1] * x[1] + s * s;
        let gm_d3 = G * self.mass.at(t) / (d2 * d2.sqrt());
        [-gm_d3 * x[0], -gm_d3 * x[1], -gm_d3 * x[2] * s / zb]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::test_util::{check_gradient, sample_points};

    #[test]
    fn acceleration_is_minus_gradient() {
        let p = MiyamotoNagai { mass: 5.0.into(), a: 3.0.into(), b: 0.3.into() };
        check_gradient(&p, &sample_points(3.0), 0.0);
    }

    #[test]
    fn far_field_is_keplerian() {
        let p = MiyamotoNagai { mass: 5.0.into(), a: 3.0.into(), b: 0.3.into() };
        let r = 1e5;
        let phi = p.potential(&[r, 0.0, 0.0], 0.0);
        assert!((phi / (-G * 5.0 / r) - 1.0).abs() < 1e-4);
    }
}
