//! NFW (1997) halo, in the same parameterization as Gala's `NFWPotential`:
//! Phi = -G m ln(1 + r/r_s) / r, with m = 4 pi rho_0 r_s^3 and
//! M(<r) = m [ln(1 + x) - x / (1 + x)], x = r / r_s.

use super::{norm, radial, Param, Potential};
use crate::units::G;

#[derive(Clone, Debug)]
pub struct Nfw {
    /// Scale mass m = 4 pi rho_0 r_s^3 (not the virial mass).
    pub m: Param,
    /// Scale radius.
    pub r_s: Param,
}

impl Nfw {
    /// From virial mass and concentration, where `m_vir` is the mass inside
    /// `r_vir = c r_s`.
    pub fn from_virial(m_vir: f64, r_vir: f64, concentration: f64) -> Self {
        let c = concentration;
        let m = m_vir / ((1.0 + c).ln() - c / (1.0 + c));
        Nfw { m: m.into(), r_s: (r_vir / c).into() }
    }

    /// Mass enclosed within radius `r`.
    pub fn enclosed_mass(&self, r: f64, t: f64) -> f64 {
        let x = r / self.r_s.at(t);
        self.m.at(t) * ((1.0 + x).ln() - x / (1.0 + x))
    }
}

impl Potential for Nfw {
    fn potential(&self, x: &[f64; 3], t: f64) -> f64 {
        let r = norm(x);
        let rs = self.r_s.at(t);
        if r == 0.0 {
            return -G * self.m.at(t) / rs;
        }
        -G * self.m.at(t) * (1.0 + r / rs).ln() / r
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
        let p = Nfw { m: 50.0.into(), r_s: 20.0.into() };
        check_gradient(&p, &sample_points(20.0), 0.0);
    }

    #[test]
    fn virial_constructor_encloses_virial_mass() {
        let p = Nfw::from_virial(100.0, 200.0, 8.0);
        assert!((p.enclosed_mass(200.0, 0.0) / 100.0 - 1.0).abs() < 1e-12);
        assert!((p.r_s.at(0.0) - 25.0).abs() < 1e-12);
    }
}
