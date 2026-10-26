//! Burkert (1995) cored halo: rho = rho_0 r_0^3 / ((r + r_0)(r^2 + r_0^2)).
//!
//! With x = r / r_0:
//! M(<r)  = pi rho_0 r_0^3 [ln((1 + x^2)(1 + x)^2) - 2 arctan x]
//! Phi(r) = -G M(<r) / r - 4 pi G rho_0 r_0^2 [pi/4 + ln(1 + x)/2 - ln(1 + x^2)/4 - arctan(x)/2]
//! where the second term is the contribution of shells outside r.

use std::f64::consts::PI;

use super::{norm, radial, Param, Potential};
use crate::units::G;

#[derive(Clone, Debug)]
pub struct Burkert {
    /// Central density.
    pub rho_0: Param,
    /// Core radius.
    pub r_0: Param,
}

impl Burkert {
    pub fn enclosed_mass(&self, r: f64, t: f64) -> f64 {
        let r0 = self.r_0.at(t);
        let x = r / r0;
        PI * self.rho_0.at(t) * r0 * r0 * r0 * (((1.0 + x * x) * (1.0 + x) * (1.0 + x)).ln() - 2.0 * x.atan())
    }
}

impl Potential for Burkert {
    fn potential(&self, x: &[f64; 3], t: f64) -> f64 {
        let r = norm(x);
        let (rho0, r0) = (self.rho_0.at(t), self.r_0.at(t));
        let s = r / r0;
        let outer = PI / 4.0 + 0.5 * (1.0 + s).ln() - 0.25 * (1.0 + s * s).ln() - 0.5 * s.atan();
        let inner = if r == 0.0 { 0.0 } else { G * self.enclosed_mass(r, t) / r };
        -inner - 4.0 * PI * G * rho0 * r0 * r0 * outer
    }

    fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3] {
        let r = norm(x);
        radial(x, r, G * self.enclosed_mass(r, t) / (r * r * r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::test_util::{check_gradient, sample_points};

    #[test]
    fn acceleration_is_minus_gradient() {
        let p = Burkert { rho_0: 0.02.into(), r_0: 5.0.into() };
        check_gradient(&p, &sample_points(5.0), 0.0);
    }

    #[test]
    fn enclosed_mass_matches_density_integral() {
        let p = Burkert { rho_0: 0.02.into(), r_0: 5.0.into() };
        let r_max = 12.0;
        let n = 100_000;
        let dr = r_max / n as f64;
        let m: f64 = (0..n)
            .map(|k| {
                let r = (k as f64 + 0.5) * dr;
                4.0 * PI * r * r * 0.02 * 125.0 / ((r + 5.0) * (r * r + 25.0)) * dr
            })
            .sum();
        assert!((m / p.enclosed_mass(r_max, 0.0) - 1.0).abs() < 1e-8);
    }
}
