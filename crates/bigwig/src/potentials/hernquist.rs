//! Hernquist (1990) bulge: Phi = -G M / (r + a), M(<r) = M r^2 / (r + a)^2.

use super::{norm, radial, Param, Potential};
use crate::units::G;

#[derive(Clone, Debug)]
pub struct Hernquist {
    /// Total mass.
    pub mass: Param,
    /// Scale radius.
    pub a: Param,
}

impl Potential for Hernquist {
    fn potential(&self, x: &[f64; 3], t: f64) -> f64 {
        -G * self.mass.at(t) / (norm(x) + self.a.at(t))
    }

    fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3] {
        let r = norm(x);
        let s = r + self.a.at(t);
        radial(x, r, G * self.mass.at(t) / (s * s * r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::test_util::{check_gradient, sample_points};

    #[test]
    fn acceleration_is_minus_gradient_and_mass_profile_matches() {
        let p = Hernquist { mass: 2.0.into(), a: 1.5.into() };
        check_gradient(&p, &sample_points(1.5), 0.0);
        for r in [0.5, 1.5, 10.0] {
            let a = p.acceleration(&[r, 0.0, 0.0], 0.0)[0];
            let m_enc = -a * r * r / G;
            assert!((m_enc - 2.0 * r * r / ((r + 1.5) * (r + 1.5))).abs() < 1e-12);
        }
    }
}
